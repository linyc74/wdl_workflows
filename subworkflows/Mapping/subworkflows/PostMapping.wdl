version 1.0


workflow PostMapping {
    input {
        File inFileUnSortRawBam
        File inFileDbsnpVcf
        File inFileDbsnpVcfIndex
        File refFa
        File refFai
        File refDict
        String sampleName
    }
    
    call Sort { 
        input:
            inFileBam = inFileUnSortRawBam,
            sampleName = sampleName
    }

    call MarkDuplicates {
        input:
            inFileBam = Sort.outFileBam,
            sampleName = sampleName
    }

    call BaseRecalibrator {
        input:
            inFileBam = MarkDuplicates.outFileBam,
            sampleName = sampleName,
            inFileDbsnpVcf = inFileDbsnpVcf,
            inFileDbsnpVcfIndex = inFileDbsnpVcfIndex,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict
    }

    call ApplyBqsr {
        input:
            inFileBam = MarkDuplicates.outFileBam,
            inFileRecalibrationTable = BaseRecalibrator.outFileRecalibrationTable,
            sampleName = sampleName,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict
    }

    output {
        File outFileBam = ApplyBqsr.outFileBam
        File outFileBamIndex = ApplyBqsr.outFileBamIndex
        File outFileSortedRawBam = Sort.outFileBam
    }
}


task Sort {
    input {
        File inFileBam
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        samtools sort ~{inFileBam} > ~{sampleName}.bam
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileBam = "~{sampleName}.bam"
    }
}


task MarkDuplicates {
    input {
        File inFileBam
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        gatk MarkDuplicates \
        --INPUT ~{inFileBam} \
        --OUTPUT ~{sampleName}.bam \
        --METRICS_FILE ~{sampleName}_metrics.txt \
        --REMOVE_DUPLICATES false
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileBam = "~{sampleName}.bam"
        File outFileMetrics = "~{sampleName}_metrics.txt"
    }
}


task BaseRecalibrator {
    input {
        File inFileDbsnpVcf
        File inFileDbsnpVcfIndex
        File inFileBam
        File refFa
        File refFai
        File refDict
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        gatk BaseRecalibrator \
        --input ~{inFileBam} \
        --reference ~{refFa} \
        --known-sites ~{inFileDbsnpVcf} \
        --output ~{sampleName}.recalibration.table
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileRecalibrationTable = '~{sampleName}.recalibration.table'
    }
}


task ApplyBqsr {
    input {
        File inFileBam
        File inFileRecalibrationTable
        File refFa
        File refFai
        File refDict
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        gatk ApplyBQSR \
        --input ~{inFileBam} \
        --reference ~{refFa} \
        --bqsr-recal-file ~{inFileRecalibrationTable} \
        --output ~{sampleName}.bam
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileBam = "~{sampleName}.bam"
        File outFileBamIndex = "~{sampleName}.bai"
    }
}