version 1.0

import "subworkflows/GeneralTask.wdl" as general
import "subworkflows/Mapping/TNPairedMapping.wdl" as mapper
import "subworkflows/VariantCalling/TNPairedVariantCalling.wdl" as caller
import "subworkflows/VariantPicking.wdl" as pick
import "subworkflows/Annotate.wdl" as annot


workflow SomaticPipelineTumorNormalMode {
    input {
        Array[Array[File]] inFileTumorFastqs
        Array[Array[File]] inFileNormalFastqs
        File inFileDbsnpVcf
        File inFileDbsnpVcfIndex
        File inFileIntervalBed
        File inFileGermlineResource
        File inFileGermlineResourceIndex
        File inFilePON
        File inFilePONindex
        File inDirPCGRref
        File refAmb
        File refAnn
        File refBwt
        File refPac
        File refSa
        File refFa
        File refFai
        File refDict
        Array[String] tumorSampleNames
        Array[String] normalSampleNames
        Array[String] finalOutputNames
    }

    scatter (i in range(length(finalOutputNames))) {
        Array[File] iFTFs = inFileTumorFastqs[i]
        Array[File] iFNFs = inFileNormalFastqs[i]
        String tumorSampleName = tumorSampleNames[i]
        String normalSampleName = normalSampleNames[i]
        String finalOutputName = finalOutputNames[i]

        call general.FastQC as fastqcTumorFastq {
            input:
                inFileFastqR1 = iFTFs[0],
                inFileFastqR2 = iFTFs[1]
        }

        call general.FastQC as fastqcNormalFastq {
            input:
                inFileFastqR1 = iFNFs[0],
                inFileFastqR2 = iFNFs[1]
        }

        call mapper.TNPairedMapping as TNmapping {
            input:
                inFileTumorFastqs = iFTFs,
                inFileNormalFastqs = iFNFs,
                inFileDbsnpVcf = inFileDbsnpVcf,
                inFileDbsnpVcfIndex = inFileDbsnpVcfIndex,
                refAmb = refAmb,
                refAnn = refAnn,
                refBwt = refBwt,
                refPac = refPac,
                refSa = refSa,
                refFa = refFa,
                refFai = refFai,
                refDict = refDict,
                tumorSampleName = tumorSampleName,
                normalSampleName = normalSampleName
        }

        call general.BamStats as tumorBamStats {
            input:
                inFileBam = TNmapping.outFileTumorBam,
                sampleName = tumorSampleName
        }

        call general.BamStats as normalBamStats {
            input:
                inFileBam = TNmapping.outFileNormalBam,
                sampleName = normalSampleName
        }

        call caller.TNPairedVariantCalling as variantCalling {
            input:
                inFileTumorBam = TNmapping.outFileTumorBam,
                inFileTumorBamIndex = TNmapping.outFileTumorBamIndex,
                inFileNormalBam = TNmapping.outFileNormalBam,
                inFileNormalBamIndex = TNmapping.outFileNormalBamIndex,
                inFileIntervalBed = inFileIntervalBed,
                inFileGermlineResource = inFileGermlineResource,
                inFileGermlineResourceIndex = inFileGermlineResourceIndex,
                inFilePON = inFilePON,
                inFilePONindex = inFilePONindex,
                refFa = refFa,
                refFai = refFai,
                refDict = refDict,
                tumorSampleName = tumorSampleName,
                normalSampleName = normalSampleName,
                sampleName = finalOutputName,
                vardictMinimumAF = 0.01
        }

        call pick.VariantPicking as variantPicking {
            input:
                inFileVcfSS = variantCalling.outFileSomaticsniperFilteredVcfGz,
                inFileVcfMU = variantCalling.outFileMuseFilteredVcfGz,
                inFileVcfM2 = variantCalling.outFileMutect2FilteredVcfGz,
                inFileVcfLF = variantCalling.outFileLofreqFilteredVcfGz,
                inFileVcfVD = variantCalling.outFileVardictFilteredVcfGz,
                infileVcfVS = variantCalling.outFileVarscanFilteredVcfGz,
                refFa = refFa,
                sampleName = finalOutputName
        }

        call annot.Annotate as variantAnnotation {
            input:
                inFileVcfGz = variantPicking.outFileVcfGz,
                inFileVcfIndex = variantPicking.outFileVcfIndex,
                inDirPCGRref = inDirPCGRref,
                sampleName = finalOutputName
        }
    }
    
    output {
        Array[Array[File]] outFileTumorTrimmedFastqs = TNmapping.outFileTumorTrimmedFastqs
        Array[Array[File]] outFileNormalTrimmedFastqs = TNmapping.outFileNormalTrimmedFastqs

        Array[File] outFileTumorFastqcReportTar = fastqcTumorFastq.outFileFastqcReportTar
        Array[File] outFileNormalFastqcReportTar = fastqcNormalFastq.outFileFastqcReportTar

        Array[File] outFileTumorBam = TNmapping.outFileTumorBam
        Array[File] outFileNormalBam = TNmapping.outFileNormalBam
        Array[File] outFileTumorBamIndex = TNmapping.outFileTumorBamIndex
        Array[File] outFileNormalBamIndex = TNmapping.outFileNormalBamIndex

        Array[File] outFileTumorSortedRawBam = TNmapping.outFileTumorSortedRawBam
        Array[File] outFileNormalSortedRawBam = TNmapping.outFileNormalSortedRawBam

        Array[File] outFileTumorBamStats = tumorBamStats.outFileBamStats
        Array[File] outFileNormalBamStats = normalBamStats.outFileBamStats

        Array[File] outFileSomaticsniperFilteredVcfGz = variantCalling.outFileSomaticsniperFilteredVcfGz
        Array[File] outFileSomaticsniperFilteredVcfIndex = variantCalling.outFileSomaticsniperFilteredVcfIndex
        Array[File] outFileLofreqFilteredVcfGz = variantCalling.outFileLofreqFilteredVcfGz
        Array[File] outFileLofreqFilteredVcfIndex = variantCalling.outFileLofreqFilteredVcfIndex
        Array[File] outFileMuseFilteredVcfGz = variantCalling.outFileMuseFilteredVcfGz
        Array[File] outFileMuseFilteredVcfIndex = variantCalling.outFileMuseFilteredVcfIndex
        Array[File] outFileMutect2FilteredVcfGz = variantCalling.outFileMutect2FilteredVcfGz
        Array[File] outFileMutect2FilteredVcfIndex = variantCalling.outFileMutect2FilteredVcfIndex
        Array[File] outFileVardictFilteredVcfGz = variantCalling.outFileVardictFilteredVcfGz
        Array[File] outFileVardictFilteredVcfIndex = variantCalling.outFileVardictFilteredVcfIndex
        Array[File] outFileVarscanFilteredVcfGz = variantCalling.outFileVarscanFilteredVcfGz
        Array[File] outFileVarscanFilteredVcfIndex = variantCalling.outFileVarscanFilteredVcfIndex

        Array[File] outFilePCGRannotatedVcf = variantAnnotation.outFilePCGRannotatedVcf
        Array[File] outFilePCGRannotatedVcfIndex = variantAnnotation.outFilePCGRannotatedVcfIndex

        Array[File] outFileMaf = variantAnnotation.outFileMaf
        Array[File] outFileCsv = variantAnnotation.outFileCsv

        Array[File] outFilePCGRflexdbHtml = variantAnnotation.outFilePCGRflexdbHtml
        Array[File] outFilePCGRhtml = variantAnnotation.outFilePCGRhtml
    }
}