version 1.0

import "GeneralTask.wdl" as general


workflow Varscan {
    input {
        File inFileTumorBam
        File inFileNormalBam
        File inFileIntervalBed
        File refFa
        File refFai
        String tumorSampleName
        String normalSampleName
        String sampleName
    }
 
    call Mpileup as tumorMpileup {
        input:
            inFileBam = inFileTumorBam,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            sampleName = tumorSampleName
    }
 
    call Mpileup as normalMpileup {
        input:
            inFileBam = inFileNormalBam,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            sampleName = normalSampleName
    }

    call VarscanSomatic {
        input:
            inFileTumorPileup = tumorMpileup.outFilePileup,
            inFileNormalPileup = normalMpileup.outFilePileup,
            sampleName = sampleName
    }

    call general.ConcatSnvIndelVcfs as concat {
        input:
            inFileSnvVcf = VarscanSomatic.outFileSnpVcfGz,
            inFileSnvVcfIndex = VarscanSomatic.outFileSnpVcfIndex,
            inFileIndelVcf = VarscanSomatic.outFileIndelVcfGz,
            infileIndelVcfIndex = VarscanSomatic.outFileIndelVcfIndex,
            sampleName = sampleName,
            callerName = "varscan"
    }

    call general.VariantFiltering as filter {
        input:
            inFileVcfGz = concat.outFileVcfGz,
            sampleName = sampleName,
            callerName = "varscan"
    }

    output {
        File outFileVcfGz = concat.outFileVcfGz
        File outFileVcfIndex = concat.outFileVcfIndex
        File outFileFilteredVcfGz = filter.outFileVcfGz
        File outFileFilteredVcfIndex = filter.outFileVcfIndex
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


task VarscanSomatic {
    input {
        File inFileTumorPileup
        File inFileNormalPileup
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        java -jar /usr/local/seqslab/VarScan.v2.3.7.jar somatic \
        ~{inFileNormalPileup} \
        ~{inFileTumorPileup} \
        --output-snp ~{sampleName}_snp.vcf \
        --output-indel ~{sampleName}_indel.vcf \
        --strand-filter 1 \
        --output-vcf 1
        bgzip --stdout ~{sampleName}_snp.vcf > ~{sampleName}_snp.vcf.gz
        bgzip --stdout ~{sampleName}_indel.vcf > ~{sampleName}_indel.vcf.gz        
        tabix --preset vcf ~{sampleName}_snp.vcf.gz
        tabix --preset vcf ~{sampleName}_indel.vcf.gz
    >>>
 
    output {
        File outFileSnpVcfGz = "~{sampleName}_snp.vcf.gz"
        File outFileIndelVcfGz = "~{sampleName}_indel.vcf.gz"
        File outFileSnpVcfIndex = "~{sampleName}_snp.vcf.gz.tbi"
        File outFileIndelVcfIndex = "~{sampleName}_indel.vcf.gz.tbi"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}