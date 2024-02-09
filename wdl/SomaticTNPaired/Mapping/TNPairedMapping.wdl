version 1.0

import "TrimAndMapping.wdl" as trimAndMapping
import "PostMapping.wdl" as postMapping


workflow TNPairedMapping {
    input {
        Array[File] inFileTumorFastqs
        Array[File] inFileNormalFastqs
        File inFileDbsnpVcf
        File inFileDbsnpVcfIndex
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

    call trimAndMapping.TrimAndMapping as trimAndMapTumor {
        input:
            refAmb = refAmb,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refSa = refSa,
            refFa = refFa,
            refFai = refFai,
            sampleName = tumorSampleName,
            inFileFastqs = inFileTumorFastqs
    }

    call trimAndMapping.TrimAndMapping as trimAndMapNormal {
        input:
            refAmb = refAmb,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refSa = refSa,
            refFa = refFa,
            refFai = refFai,
            sampleName = normalSampleName,
            inFileFastqs = inFileNormalFastqs
    }

    call postMapping.PostMapping as postMappingTumor {
        input:
            inFileUnSortRawBam = trimAndMapTumor.outFileUnSortRawBam,
            inFileDbsnpVcf = inFileDbsnpVcf,
            inFileDbsnpVcfIndex = inFileDbsnpVcfIndex,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            sampleName = tumorSampleName
    }

    call postMapping.PostMapping as postMappingNormal {
        input:
            inFileUnSortRawBam = trimAndMapNormal.outFileUnSortRawBam,
            inFileDbsnpVcf = inFileDbsnpVcf,
            inFileDbsnpVcfIndex = inFileDbsnpVcfIndex,
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