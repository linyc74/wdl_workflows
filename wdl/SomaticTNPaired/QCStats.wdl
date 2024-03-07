version 1.0


workflow QCStats {
    input {
        Array[File] inFileTumorFastqPair
        File inFileTumorBam
        String tumorSampleName

        Array[File] inFileNormalFastqPair
        File inFileNormalBam
        String normalSampleName
    }

    call FastQC as tumorFastqc {
        input:
            inFileFastqR1 = inFileTumorFastqPair[0],
            inFileFastqR2 = inFileTumorFastqPair[1],
            sampleName = tumorSampleName
    }

    call FastQC as normalFastqc {
        input:
            inFileFastqR1 = inFileNormalFastqPair[0],
            inFileFastqR2 = inFileNormalFastqPair[1],
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
        mkdir -p ~{sampleName}_fastqc_report
        fastqc ~{inFileFastqR1} ~{inFileFastqR2} -o ~{sampleName}_fastqc_report
        tar -cvf ~{sampleName}_fastqc_report.tar ~{sampleName}_fastqc_report
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
