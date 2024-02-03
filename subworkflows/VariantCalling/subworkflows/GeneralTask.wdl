version 1.0


task VariantFiltering {
    input {
        File inFileVcfGz
        String flaggingCriteria = "\"low_depth: DP <= 10\""
        String removalFlags = "fragment,haplotype,normal_artifact,panel_or_normals,position,slippage,weak_evidence,map_qual,panel_of_normals,base_qual,strand_bias,multiallelic,orientation,contamination,clustered_events,MQ40,QD2,SOR3,MQRankSum-12.5,FS60,FS200,ReadPosRankSum-8,ReadPosRankSum-20,Bias,q22.5,Q10,Cluster0bp"
        String sampleName
        String callerName
    }

    command <<<
        set -e -o pipefail
        zcat ~{inFileVcfGz} > in.vcf
        python /usr/local/seqslab/omic variant-filtering \
        --input-vcf in.vcf \
        --output-vcf ~{sampleName}_filtered.vcf \
        --variant-flagging-criteria ~{flaggingCriteria}  \
        --variant-removal-flags ~{removalFlags}
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