import "slalom_sub.wdl" as sub

task preprocess {
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
    Float p_threshold
    Boolean? grch38
    Boolean x_chromosome
    Int window

    command <<<

        make_finemap_inputs.py \
            --sumstats ${sumstats} \
            --rsid-col "${rsid_col}" \
            --chromosome-col "${chromosome_col}" \
            --position-col "${position_col}" \
            --allele1-col "${allele1_col}" \
            --allele2-col "${allele2_col}" \
            --freq-col "${freq_col}" \
            --beta-col "${beta_col}" \
            --se-col "${se_col}" \
            --p-col "${p_col}" \
            --extra-cols "all_meta_AF" "p_het" "n_cases" "n_controls" "n_samples" "n_datasets" "n_biobanks" "is_strand_flip" "is_diff_AF_gnomAD" "n_afr" "n_amr" "n_eas" "n_fin" "n_nfe" \
            --delimiter "${delimiter}" \
            ${true='--grch38 ' false=' ' grch38} \
            --exclude-MHC \
            --no-upload \
            --prefix ${pheno} \
            --out ${pheno} \
            --window ${window} \
            --wdl \
            ${true='--x-chromosome' false=' ' x_chromosome} \
            --p-threshold ${p_threshold}

        res=`cat ${pheno}_had_results`

        if [ "$res" == "False" ]
        then
            touch ${pheno}".z"
            touch ${pheno}".lead_snps.txt"
            touch ${pheno}".bed"
        fi

    >>>

    output {

        String out_pheno = pheno
        Array[File] zfiles = glob("*.z")
        File leadsnps = pheno + ".lead_snps.txt"
        File bed = pheno + ".bed"
        File log = pheno + ".log"
        Boolean had_results = read_boolean("${pheno}_had_results")
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "104 GB"
        disks: "local-disk 10 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

workflow slalom {

    String zones
    String docker
    File phenolistfile
    String sumstats_pattern
    Boolean grch38

    Array[Array[String]] all_phenos = read_tsv(phenolistfile)

    scatter (pheno in all_phenos) {

        call preprocess {
            input: zones=zones, docker=docker, pheno=pheno[0], grch38=grch38, sumstats_pattern=sumstats_pattern
        }

        if (preprocess.had_results) {
            call sub.finemap {
                input: zones=zones, docker=docker, pheno=preprocess.out_pheno, zfiles=preprocess.zfiles, grch38=grch38,
                n_samples_arg=pheno[1]
            }
        }

    }
}
