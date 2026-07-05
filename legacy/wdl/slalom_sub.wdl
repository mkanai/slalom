task annotate {
    File snp
    File annotation_script
    String prefix = basename(snp, ".z")
    String n_samples_arg
    File? lead_variant
    Boolean grch38
    String zones
    String docker

    command <<<

        mem=$(grep MemTotal /proc/meminfo | awk '{printf("%dg\n", $2 / 1024 / 1024)}')

        cat << "__EOF__" > /usr/local/lib/python3.6/dist-packages/pyspark/conf/spark-defaults.conf
        spark.hadoop.google.cloud.auth.service.account.enable true
        spark.hadoop.fs.gs.requester.pays.mode AUTO
        spark.hadoop.fs.gs.requester.pays.project.id encode-uk-biobank-restrict
        __EOF__

        if ${defined(lead_variant)}
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
            ' ${lead_variant} ${snp}
        fi

        pip3 install cython && pip3 install gnomad
        PYSPARK_SUBMIT_ARGS="--conf spark.driver.memory=$mem pyspark-shell" \
        python3 ${annotation_script} \
        --snp ${snp} \
        --out ${prefix}.ld.snp \
        --out-summary ${prefix}.summary.txt \
        --annotate-consequence \
        --annotate-gnomad-freq \
        --export-r \
        --lead-variant-choice "prob" \
        --weighted-average-r ${n_samples_arg} \
        --dentist \
        --abf \
        --summary \
        ${true='--lead-variant $(cat lead_variant.txt)' false='' defined(lead_variant)} \
        ${true='--reference-genome GRCh38 ' false='--reference-genome GRCh37 ' grch38}

    >>>

    output {
        File ld_snp = prefix + ".ld.snp"
        File summary = prefix + ".summary.txt"
    }

    runtime {

        docker: "${docker}"
        cpu: "2"
        memory: "104 GB"
        disks: "local-disk 50 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

task combine {
    String pheno
    Array[File] abf_snp
    Array[File] summary
    String zones
    String docker

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
        awk -f combine_snp.awk -v pheno=${pheno} ${sep=" " abf_snp} | bgzip -c -@ $cpu > ${pheno}.ABF.snp.bgz

        # Combine summary files
        awk -v pheno=${pheno} '
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
        ' ${sep=" " summary} | bgzip -c -@ $cpu > ${pheno}.SLALOM.summary.bgz

        tabix -s 5 -b 6 -e 6 -S 1 ${pheno}.ABF.snp.bgz

    >>>

    output {

        File out_abf_snp = pheno + ".ABF.snp.bgz"
        File out_abf_snp_tbi = pheno + ".ABF.snp.bgz.tbi"
        File out_summary = pheno + ".SLALOM.summary.bgz"

    }

    runtime {

        docker: "${docker}"
        cpu: "1"
        memory: "7 GB"
        disks: "local-disk 30 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

workflow finemap {

    String zones
    String docker
    String pheno
    Array[File] zfiles
    Boolean grch38
    String n_samples_arg

    File annotation_script
    File? lead_variant

    scatter (zfile in zfiles) {

        call annotate {
            input: zones=zones, annotation_script=annotation_script, lead_variant=lead_variant,
                snp=zfile, grch38=grch38, n_samples_arg=n_samples_arg
        }

    }

    call combine {
        input: zones=zones, docker=docker, pheno=pheno,
            abf_snp=annotate.ld_snp, summary=annotate.summary
    }
}
