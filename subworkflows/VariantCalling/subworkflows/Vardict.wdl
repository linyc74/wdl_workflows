version 1.0

import "GeneralTask.wdl" as general


workflow Vardict {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File inFileNormalBam
        File inFileNormalBamIndex
        File inFileIntervalBed
        File refFa
        File refFai
        Float minimumAF
        String tumorSampleName
        String normalSampleName
        String sampleName
    }
 
    call VardictPaired {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            minimumAF = minimumAF,
            tumorSampleName = tumorSampleName
    }

    call TestsomaticR {
        input:
            inFileVardictStep1 = VardictPaired.outFileResult
    }

    call Var2Vcf {
        input:
            inFileVardictStep2 = TestsomaticR.outFileResult,
            minimumAF = minimumAF,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = sampleName
    }

    call general.VariantFiltering as filter {
        input:
            inFileVcfGz = Var2Vcf.outFileVcfGz,
            sampleName = sampleName,
            callerName = "vardict"
    }

    output {
        File outFileVcfGz = Var2Vcf.outFileVcfGz
        File outFileVcfIndex = Var2Vcf.outFileVcfIndex
        File outFileFilteredVcfGz = filter.outFileVcfGz
        File outFileFilteredVcfIndex = filter.outFileVcfIndex
    }
}


task VardictPaired {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File inFileNormalBam
        File inFileNormalBamIndex
        File inFileIntervalBed
        File refFa
        File refFai
        Float minimumAF
        String tumorSampleName
    }
 
    command <<<
        set -e -o pipefail
        VarDict \
        -G ~{refFa} \
        -f ~{minimumAF} \
        -N ~{tumorSampleName} \
        -b "~{inFileTumorBam}|~{inFileNormalBam}" \
        -c 1 \
        -S 2 \
        -E 3 \
        -g 4 \
        ~{inFileIntervalBed} \
        1> VardictStep1
    >>>
 
    output {
        File outFileResult = "VardictStep1"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}


task TestsomaticR {
    input {
        File inFileVardictStep1
    }
 
    command <<<
        set -e -o pipefail
        cat ~{inFileVardictStep1} \
        | \
        testsomatic.R \
        1> VardictStep2 
    >>>
 
    output {
        File outFileResult = "VardictStep2"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}


task Var2Vcf {
    input {
        File inFileVardictStep2
        Float minimumAF
        String tumorSampleName
        String normalSampleName
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        cat ~{inFileVardictStep2} \
        | \
        var2vcf_paired.pl \
        -N "~{tumorSampleName} | ~{normalSampleName}" \
        -f ~{minimumAF} \
        1> ~{sampleName}.vcf
        bgzip --stdout ~{sampleName}.vcf > ~{sampleName}_vardict.vcf.gz
        tabix --preset vcf ~{sampleName}_vardict.vcf.gz
    >>>
 
    output {
        File outFileVcfGz = "~{sampleName}_vardict.vcf.gz"
        File outFileVcfIndex = "~{sampleName}_vardict.vcf.gz.tbi"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}
