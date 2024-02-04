version 1.0

import "GeneralTask.wdl" as general


workflow Lofreq {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File inFileNormalBam
        File inFileNormalBamIndex
        File inFileIntervalBed
        File refFa
        File refFai
        String sampleName
    }
 
    call LofreqSomatic {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            sampleName = sampleName
    }

    call general.ConcatSnvIndelVcfs as concat {
        input:
            inFileSnvVcf = LofreqSomatic.outFileSnvVcf,
            inFileSnvVcfIndex = LofreqSomatic.outFileSnvVcfIndex,
            inFileIndelVcf = LofreqSomatic.outFileIndelVcf,
            infileIndelVcfIndex = LofreqSomatic.outFileIndelVcfIndex,
            sampleName = sampleName,
            callerName = "lofreq"
    }

    call general.VariantFiltering as filter {
        input:
            inFileVcfGz = concat.outFileVcfGz,
            sampleName = sampleName,
            callerName = "lofreq"
    }
   
    output {
        File outFileVcfGz = concat.outFileVcfGz
        File outFileVcfIndex = concat.outFileVcfIndex
        File outFileFilteredVcfGz = filter.outFileVcfGz
        File outFileFilteredVcfIndex = filter.outFileVcfIndex
    }
}


task LofreqSomatic {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File inFileNormalBam
        File inFileNormalBamIndex
        File inFileIntervalBed
        File refFa
        File refFai
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        lofreq somatic \
        --normal ~{inFileNormalBam} \
        --tumor ~{inFileTumorBam} \
        --ref ~{refFa} \
        --threads 1 \
        --call-indels \
        -l ~{inFileIntervalBed} \
        -o ~{sampleName}
    >>>
 
    output {
        File outFileSnvVcf = "~{sampleName}somatic_final.snvs.vcf.gz"
        File outFileSnvVcfIndex = "~{sampleName}somatic_final.snvs.vcf.gz.tbi"
        File outFileIndelVcf = "~{sampleName}somatic_final.indels.vcf.gz"
        File outFileIndelVcfIndex = "~{sampleName}somatic_final.indels.vcf.gz.tbi"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}