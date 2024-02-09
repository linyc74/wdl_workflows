version 1.0


workflow QCStats {
    input {
        File inFileTumorFastqR1
        File inFileTumorFastqR2
        File inFileTumorBam
        String tumorSampleName

        File inFileNormalFastqR1
        File inFileNormalFastqR2
        File inFileNormalBam
        String normalSampleName
    }

    call FastQC as tumorFastqc {
        input:
            inFileFastqR1 = inFileTumorFastqR1,
            inFileFastqR2 = inFileTumorFastqR2,
            sampleName = tumorSampleName
    }

    call FastQC as normalFastqc {
        input:
            inFileFastqR1 = inFileNormalFastqR1,
            inFileFastqR2 = inFileNormalFastqR2,
            sampleName = normalSampleName
    }

    call BamStats as tumorBamStats {
        input:
            inFileBam = inFileTumorBam,
            sampleName = tumorSampleName
    }

    call BamStats as normalBamStats {
        input:
            inFileBam = inFileNormalBam,
            sampleName = normalSampleName
    }

    output {
        File outFileTumorFastqcReportTar = tumorFastqc.outFileFastqcReportTar
        File outFileNormalFastqcReportTar = normalFastqc.outFileFastqcReportTar
        File outFileTumorBamStats = tumorBamStats.outFileBamStats
        File outFileNormalBamStats = normalBamStats.outFileBamStats
    }
}


task FastQC {
    input {
        File inFileFastqR1
        File inFileFastqR2
        String sampleName
    }

    command <<<
        set -e -o pipefail
        mkdir -p fastqc_report
        fastqc ~{inFileFastqR1} ~{inFileFastqR2} -o fastqc_report
        tar cvf fastqc_report.tar fastqc_report
    >>>

    output {
        File outFileFastqcReportTar = "~{sampleName}_fastqc_report.tar"
    }

    runtime {
        docker: 'nycu:latest'
    }
}


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
