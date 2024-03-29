version 1.0

workflow slalom {
    input {
        String zones
        String docker
        File phenolistfile
        String sumstats_pattern
        File annotation_script
        Boolean grch38
        Boolean align_alleles
    }

    Array[Array[String]] all_phenos = read_tsv(phenolistfile)

    scatter (pheno in all_phenos) {
        call preprocess {
            input: zones=zones, docker=docker, pheno=pheno[0], grch38=grch38, sumstats_pattern=sumstats_pattern
        }

        if (preprocess.had_results) {
            File zfile = preprocess.zfiles[0]
            scatter (zfile in preprocess.zfiles) {
                call annotate {
                    input: zones=zones, annotation_script=annotation_script,
                        snp=zfile, grch38=grch38, align_alleles=align_alleles,
                        n_samples_arg=pheno[1]
                }
            }
            call combine {
                input: zones=zones, docker=docker, pheno=pheno[0],
                    abf_snp=annotate.ld_snp, summary=annotate.summary
            }
        }
    }
}

task preprocess {
    input {
        String pheno
        String sumstats_pattern
        File sumstats = sub(sumstats_pattern,"\\{PHENO\\}",pheno)    
        String zones
        String docker
        String rsid_col
        String chromosome_col
        String position_col
        String allele1_col
        String allele2_col
        String freq_col
        String beta_col
        String se_col
        String p_col
        String delimiter
        Array[String]? extra_cols
        Float p_threshold
        Boolean? grch38
        Boolean x_chromosome
        Int window
    }

    command <<<

        make_finemap_inputs.py \
            --sumstats ~{sumstats} \
            --rsid-col "~{rsid_col}" \
            --chromosome-col "~{chromosome_col}" \
            --position-col "~{position_col}" \
            --allele1-col "~{allele1_col}" \
            --allele2-col "~{allele2_col}" \
            --freq-col "~{freq_col}" \
            --beta-col "~{beta_col}" \
            --se-col "~{se_col}" \
            --p-col "~{p_col}" \
            ~{true="--extra-cols " false="" defined(extra_cols)}~{sep=" " extra_cols} \
            --delimiter "~{delimiter}" \
            ~{true="--grch38 " false=" " grch38} \
            --exclude-MHC \
            --no-upload \
            --prefix ~{pheno} \
            --out ~{pheno} \
            --window ~{window} \
            --wdl \
            ~{true="--x-chromosome" false=" " x_chromosome} \
            --p-threshold ~{p_threshold}

        res=`cat ~{pheno}_had_results`

        if [ "$res" == "False" ]
        then
            touch ~{pheno}".z"
            touch ~{pheno}".lead_snps.txt"
            touch ~{pheno}".bed"
        fi

    >>>

    output {

        String out_pheno = pheno
        Array[File] zfiles = glob("*.z")
        File leadsnps = pheno + ".lead_snps.txt"
        File bed = pheno + ".bed"
        File log = pheno + ".log"
        Boolean had_results = read_boolean("~{pheno}_had_results")
    }

    runtime {
        docker: docker
        cpu: "1"
        memory: "104 GB"
        disks: "local-disk 10 HDD"
        zones: zones
        preemptible: 2
        noAddress: false
    }
}

task annotate {
    input {
        File snp
        File annotation_script
        String prefix = basename(snp, ".z")
        String n_samples_arg
        File? lead_variant
        Boolean grch38
        Boolean align_alleles
        String zones
        String docker
    }

    command <<<

        mem=$(grep MemTotal /proc/meminfo | awk '{printf("%dg\n", $2 / 1024 / 1024)}')

        cat << "__EOF__" > /usr/local/lib/python3.8/dist-packages/pyspark/conf/spark-defaults.conf
        spark.hadoop.google.cloud.auth.service.account.enable true
        spark.hadoop.fs.gs.requester.pays.mode AUTO
        spark.hadoop.fs.gs.requester.pays.project.id finngen-xavier
        __EOF__

        if ~{defined(lead_variant)}
        then
            awk '
            FNR == NR {
                lead_variants[$1] = $1
                next
            }
            FNR < NR && FNR == 1 {
                for (i = 1; i <= NF; i++) {
                    col[$i] = i
                }
            }
            FNR < NR && FNR > 1 {
                v = sprintf( \
                    "%s:%s:%s:%s", \
                    $col["chromosome"], \
                    $col["position"], \
                    $col["allele1"], \
                    $col["allele2"] \
                )
                if (v in lead_variants) {
                    print v > "lead_variant.txt"
                    exit
                }
            }
            ' ~{lead_variant} ~{snp}
        fi

        pip3 install cython && pip3 install gnomad
        PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=$mem pyspark-shell" \
        python3 ~{annotation_script} \
        --snp ~{snp} \
        --out ~{prefix}.ld.snp \
        --out-summary ~{prefix}.summary.txt \
        --annotate-cups \
        --annotate-consequence \
        --annotate-gnomad-freq \
        --export-r \
        --lead-variant-choice "prob" \
        --weighted-average-r ~{n_samples_arg} \
        --dentist \
        --abf \
        --summary \
        ~{true='--lead-variant $(cat lead_variant.txt)' false='' defined(lead_variant)} \
        ~{true='--reference-genome GRCh38 ' false='--reference-genome GRCh37 ' grch38} \
        ~{true='--align-alleles ' false='' align_alleles}

    >>>

    output {
        File ld_snp = prefix + ".ld.snp"
        File summary = prefix + ".summary.txt"
    }

    runtime {

        docker: docker
        cpu: "2"
        memory: "32 GB"
        disks: "local-disk 50 HDD"
        zones: zones
        preemptible: 2
        noAddress: false
    }
}

task combine {
    input {
        String pheno
        Array[File] abf_snp
        Array[File] summary
        String zones
        String docker
    }

    command <<<
        cpu=`grep -c ^processor /proc/cpuinfo`

        cat << "__EOF__" > combine_snp.awk
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
            gsub(" ", "\t")
            print "trait", "region", "variant", $0
        }
        FNR == 1 {
            match(FILENAME, /(chr[0-9X]+)\.([0-9]+-[0-9]+)\./, a)
            region = a[1]":"a[2]
        }
        FNR > 1 {
            chrom = $col["chromosome"]
            # sub(/^chr/, "", chrom)
            # sub(/^0/, "", chrom)
            v = sprintf( \
                "%s:%s:%s:%s", \
                chrom, \
                $col["position"], \
                $col["allele1"], \
                $col["allele2"] \
            )
            gsub(" ", "\t")
            print pheno, region, v, $0 | "sort -V -k2,3"
        }
        __EOF__

        # Combine abf .snp files
        awk -f combine_snp.awk -v pheno=~{pheno} ~{sep=" " abf_snp} | bgzip -c -@ $cpu > ~{pheno}.ABF.snp.bgz

        # Combine summary files
        awk -v pheno=~{pheno} '
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            print "trait", "region", $0
        }
        FNR == 1 {
            match(FILENAME, /(chr[0-9X]+)\.([0-9]+-[0-9]+)\./, a)
            region = a[1]":"a[2]
        }
        FNR > 1 {
            print pheno, region, $0 | "sort -V -k2,2"
        }
        ' ~{sep=" " summary} | bgzip -c -@ $cpu > ~{pheno}.SLALOM.summary.bgz

        tabix -s 5 -b 6 -e 6 -S 1 ~{pheno}.ABF.snp.bgz

    >>>

    output {

        File out_abf_snp = pheno + ".ABF.snp.bgz"
        File out_abf_snp_tbi = pheno + ".ABF.snp.bgz.tbi"
        File out_summary = pheno + ".SLALOM.summary.bgz"

    }

    runtime {

        docker: docker
        cpu: "1"
        memory: "7 GB"
        disks: "local-disk 30 HDD"
        zones: zones
        preemptible: 2
        noAddress: false
    }
}
