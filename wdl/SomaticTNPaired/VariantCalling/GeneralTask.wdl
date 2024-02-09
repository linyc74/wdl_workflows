version 1.0


task VariantFiltering {
    input {
        File inFileVcfGz
        String sampleName
        String callerName
        String flaggingCriteria = "\"low_depth: DP <= 10\""
        String removalFlags = "fragment,haplotype,normal_artifact,panel_or_normals,position,slippage,weak_evidence,map_qual,panel_of_normals,base_qual,strand_bias,multiallelic,orientation,contamination,clustered_events,MQ40,QD2,SOR3,MQRankSum-12.5,FS60,FS200,ReadPosRankSum-8,ReadPosRankSum-20,Bias,q22.5,Q10,Cluster0bp"
    }

    command <<<
        set -e -o pipefail
        zcat ~{inFileVcfGz} > in.vcf
        python /usr/local/seqslab/omic variant-filtering \
        --input-vcf in.vcf \
        --output-vcf ~{sampleName}_filtered.vcf \
        --variant-flagging-criteria ~{flaggingCriteria}  \
        --variant-removal-flags ~{removalFlags} \
        --only-pass
        bgzip --stdout ~{sampleName}_filtered.vcf > ~{sampleName}_~{callerName}_filtered.vcf.gz
        tabix --preset vcf ~{sampleName}_~{callerName}_filtered.vcf.gz
    >>>

    output {
        File outFileVcfGz = "~{sampleName}_~{callerName}_filtered.vcf.gz"
        File outFileVcfIndex = "~{sampleName}_~{callerName}_filtered.vcf.gz.tbi"
    }

    runtime {
        docker: 'nycu:latest'
    }
}


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


task ConcatSnvIndelVcfs {
    input {
        File inFileSnvVcf
        File inFileSnvVcfIndex
        File inFileIndelVcf
        File infileIndelVcfIndex
        String sampleName
        String callerName
    }

    command <<<
        set -e -o pipefail
        bcftools concat \
        --allow-overlaps \
        --output-type v \
        --output ~{sampleName}.vcf \
        ~{inFileSnvVcf} \
        ~{inFileIndelVcf}
        bgzip --stdout ~{sampleName}.vcf > ~{sampleName}_~{callerName}.vcf.gz
        tabix --preset vcf ~{sampleName}_~{callerName}.vcf.gz
    >>>

    output {
        File outFileVcfGz = "~{sampleName}_~{callerName}.vcf.gz"
        File outFileVcfIndex = "~{sampleName}_~{callerName}.vcf.gz.tbi"
    }

    runtime {
        docker: 'nycu:latest'
    }
}