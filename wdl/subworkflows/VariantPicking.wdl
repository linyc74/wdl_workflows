version 1.0


task VariantPicking {
    input {
        File refFa
        String sampleName
        File inFileVcfLofreq
        File inFileVcfMutect2
        File inFileVcfMuse
    }
 
    command <<<
        set -e -o pipefail
        python /usr/local/seqslab/omic variant-picking \
        --ref-fa ~{refFa} \
        --lofreq ~{inFileVcfLofreq} \
        --mutect2 ~{inFileVcfMutect2} \
        --muse ~{inFileVcfMuse} \
        --somatic-sniper None \
        --vardict None \
        --varscan None \
        --output-vcf ~{sampleName}_picked.vcf \
        --min-snv-caller 1 \
        --min-indel-callers 1
        bgzip \
        --stdout ~{sampleName}_picked.vcf > ~{sampleName}_picked.vcf.gz
        tabix \
        --preset vcf \
        ~{sampleName}_picked.vcf.gz
    >>>
 
    output {
        File outFileVcfGz = "~{sampleName}_picked.vcf.gz"
        File outFileVcfIndex = "~{sampleName}_picked.vcf.gz.tbi"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

