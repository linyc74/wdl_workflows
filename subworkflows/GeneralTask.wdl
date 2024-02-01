version 1.0


task BamStats {
    input {
        File inFileBam
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        samtools stats ~{inFileBam} > ~{sampleName}_stats.txt
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileBamStats = "~{sampleName}_stats.txt"
    }
}


task Concat {
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


task FastQC {
    input {
        File inFileFastqR1
        File inFileFastqR2
    }
 
    command <<<
        set -e -o pipefail
        mkdir -p fastqc_report
        fastqc ~{inFileFastqR1} ~{inFileFastqR2} -o fastqc_report
        tar cvf fastqc_report.tar fastqc_report
    >>>
 
    output {
        File outFileFastqcReportTar = "fastqc_report.tar"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}


task Mpileup {
    input {
        File inFileBam
        File inFileIntervalBed
        File refFa
        File refFai
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        samtools mpileup \
        --fasta-ref ~{refFa} \
        --positions ~{inFileIntervalBed} \
        --output ~{sampleName}_pileup.txt \
        ~{inFileBam}
    >>>
 
    output {
        File outFilePileup = "~{sampleName}_pileup.txt"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}


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
        --output-vcf ~{sampleName}_Pyfiltered.vcf \
        --variant-flagging-criteria ~{flaggingCriteria}  \
        --variant-removal-flags ~{removalFlags}
        bgzip --stdout ~{sampleName}_Pyfiltered.vcf > ~{sampleName}_~{callerName}_Pyfiltered.vcf.gz
        tabix --preset vcf ~{sampleName}_~{callerName}_Pyfiltered.vcf.gz
    >>>
 
    output {
        File outFileVcfGz = "~{sampleName}_~{callerName}_Pyfiltered.vcf.gz"
        File outFileVcfIndex = "~{sampleName}_~{callerName}_Pyfiltered.vcf.gz.tbi"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}