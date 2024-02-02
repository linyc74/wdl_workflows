version 1.0


task VariantPicking {
    input {
        File inFileVcfSS
        File inFileVcfMU
        File inFileVcfM2
        File inFileVcfLF
        File inFileVcfVD
        File infileVcfVS
        File refFa
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        python /usr/local/seqslab/omic variant-picking \
        --ref-fa ~{refFa} \
        --somatic-sniper None \
        --muse ~{inFileVcfMU} \
        --mutect2 ~{inFileVcfM2} \
        --lofreq ~{inFileVcfLF} \
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

