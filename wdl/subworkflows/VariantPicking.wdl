version 1.0


task VariantPicking {
    input {
        File inFileVcfLofreq
        File inFileVcfMutect2
        File inFileVcfMuse
        File inFileVcfSomaticsniper
        File inFileVcfVardict
        File inFileVcfVarscan
        File refFa
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        python /usr/local/seqslab/omic variant-picking \
        --ref-fa ~{refFa} \
        --somatic-sniper ~{inFileVcfSomaticsniper} \
        --muse ~{inFileVcfMuse} \
        --mutect2 ~{inFileVcfMutect2} \
        --lofreq ~{inFileVcfLofreq} \
        --vardict ~{inFileVcfVardict} \
        --varscan None ~{inFileVcfVarscan} \
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

