version 1.0

import "General.wdl" as general


workflow Mutect2 {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File? inFileNormalBam
        File? inFileNormalBamIndex
        File? refGermlineResourceVcfGz
        File? refGermlineResourceVcfIndex
        File? refPonVcfGz
        File? refPonVcfIndex
        File refIntervalBed
        File refFa
        File refFai
        File refDict
        String tumorSampleName
        String? normalSampleName
        String sampleName
        String? extraArgs
    }
 
    call Mutect2 {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            refGermlineResourceVcfGz = refGermlineResourceVcfGz,
            refGermlineResourceVcfIndex = refGermlineResourceVcfIndex,
            refPonVcfGz = refPonVcfGz,
            refPonVcfIndex = refPonVcfIndex,
            refIntervalBed = refIntervalBed,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = sampleName,
            extraArgs = extraArgs
    }

    call LearnReadOrientationModel {
        input:
            inFileF1R2 = Mutect2.outFileF1R2
    }

    call FilterMutectCalls {
        input:
            inFileArtifactPriors = LearnReadOrientationModel.outFileArtifactPriors,
            inFileVcfGz = Mutect2.outFileVcfGz,
            inFileVcfStats = Mutect2.outFileVcfStats,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,   
            sampleName = sampleName     
    }

    call general.VariantFiltering as filter {
        input:
            inFileVcfGz = FilterMutectCalls.outFileVcfGz,
            sampleName = sampleName,
            callerName = "mutect2"
    }

    output {
        File outFileVcfGz = FilterMutectCalls.outFileVcfGz
        File outFileVcfIndex = FilterMutectCalls.outFileVcfIndex
        File outFileFilteredVcfGz = filter.outFileVcfGz
        File outFileFilteredVcfIndex = filter.outFileVcfIndex
    }
}


task Mutect2 {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File? inFileNormalBam
        File? inFileNormalBamIndex
        File? refGermlineResourceVcfGz
        File? refGermlineResourceVcfIndex
        File? refPonVcfGz
        File? refPonVcfIndex
        File refIntervalBed
        File refFa
        File refFai
        File refDict
        String tumorSampleName
        String? normalSampleName
        String sampleName
        String? extraArgs
    }
 
    command <<<
        set -e -o pipefail
        gatk Mutect2 \
        --reference ~{refFa} \
        --intervals ~{refIntervalBed} \
        --input ~{inFileTumorBam} \
        ~{"--input " + inFileNormalBam} \
        --tumor-sample ~{tumorSampleName} \
        ~{"--normal-sample " + normalSampleName} \
        --output ~{sampleName}.vcf.gz \
        --f1r2-tar-gz ~{sampleName}_f1r2.tar.gz \
        --max-reads-per-alignment-start 0 \
        ~{"--germline-resource " + refGermlineResourceVcfGz} \
        ~{"--panel-of-normals " + refPonVcfGz} \
        ~{extraArgs}
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileVcfGz = "~{sampleName}.vcf.gz"
        File outFileF1R2 = "~{sampleName}_f1r2.tar.gz"
        File outFileVcfStats =  "~{sampleName}.vcf.gz.stats"
    }
}


task LearnReadOrientationModel {
    input {
        File inFileF1R2
    }
 
    command <<<
        set -e -o pipefail
        gatk LearnReadOrientationModel \
        --input ~{inFileF1R2} \
        --output artifact-prior-table.tar.gz
    >>>
 
    output {
        File outFileArtifactPriors = "artifact-prior-table.tar.gz"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}


task FilterMutectCalls {
    input {
        File inFileArtifactPriors
        File inFileVcfGz
        File inFileVcfStats
        File refFa
        File refFai
        File refDict
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        tabix -p vcf ~{inFileVcfGz}
        gatk FilterMutectCalls \
        --variant ~{inFileVcfGz} \
        --reference ~{refFa} \
        --output ~{sampleName}_filterMutectCalls.vcf.gz \
        --filtering-stats ~{inFileVcfStats} \
        --orientation-bias-artifact-priors ~{inFileArtifactPriors}
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileVcfGz = "~{sampleName}_filterMutectCalls.vcf.gz"
        File outFileVcfIndex = "~{sampleName}_filterMutectCalls.vcf.gz.tbi"
    }
}