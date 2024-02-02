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