version 1.0

import "TrimAndMap.wdl" as tm
import "PostMapping.wdl" as pm


workflow Preprocessing {
    input {
        Array[File] inFileTumorFastqPair
        Array[File] inFileNormalFastqPair
        File refDbsnpVcfGz
        File refDbsnpVcfIndex
        File refAmb
        File refAnn
        File refBwt
        File refPac
        File refSa
        File refFa
        File refFai
        File refDict
        String tumorSampleName
        String normalSampleName
    }

    call tm.TrimAndMap as trimAndMapTumor {
        input:
            inFileFastqPair = inFileTumorFastqPair,
            refAmb = refAmb,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refSa = refSa,
            refFa = refFa,
            refFai = refFai,
            sampleName = tumorSampleName
    }

    call tm.TrimAndMap as trimAndMapNormal {
        input:
            inFileFastqPair = inFileNormalFastqPair,
            refAmb = refAmb,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refSa = refSa,
            refFa = refFa,
            refFai = refFai,
            sampleName = normalSampleName,
    }

    call pm.PostMapping as postMappingTumor {
        input:
            inFileUnSortRawBam = trimAndMapTumor.outFileUnSortRawBam,
            refDbsnpVcfGz = refDbsnpVcfGz,
            refDbsnpVcfIndex = refDbsnpVcfIndex,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            sampleName = tumorSampleName
    }

    call pm.PostMapping as postMappingNormal {
        input:
            inFileUnSortRawBam = trimAndMapNormal.outFileUnSortRawBam,
            refDbsnpVcfGz = refDbsnpVcfGz,
            refDbsnpVcfIndex = refDbsnpVcfIndex,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            sampleName = normalSampleName   
    }
 
    output {
        Array[File] outFileTumorTrimmedFastqs = trimAndMapTumor.outFileTrimmedFastq
        Array[File] outFileNormalTrimmedFastqs = trimAndMapNormal.outFileTrimmedFastq

        File outFileTumorBam = postMappingTumor.outFileBam
        File outFileNormalBam = postMappingNormal.outFileBam

        File outFileTumorBamIndex = postMappingTumor.outFileBamIndex
        File outFileNormalBamIndex = postMappingNormal.outFileBamIndex

        File outFileTumorSortedRawBam = postMappingTumor.outFileSortedRawBam
        File outFileNormalSortedRawBam = postMappingNormal.outFileSortedRawBam
    }
}