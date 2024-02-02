version 1.0

import "subworkflows/GeneralTask.wdl" as general
import "subworkflows/Mapping/TNPairedMapping.wdl" as mapper
import "subworkflows/VariantCalling/TNPairedVariantCalling.wdl" as caller
import "subworkflows/PickAndAnnotate.wdl" as annotate


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
        Array[String] tumorSampleName
        Array[String] normalSampleName
        Array[String] finalOutputName
    }

    scatter (i in range(length(finalOutputName))) {
        Array[File] iFTFs = inFileTumorFastqs[i]
        Array[File] iFNFs = inFileNormalFastqs[i]
        String tSN = tumorSampleName[i]
        String nSN = normalSampleName[i]
        String fON = finalOutputName[i]

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
                tumorSampleName = tSN,
                normalSampleName = nSN
        }

        call general.BamStats as tumorBamStats {
            input:
                inFileBam = TNmapping.outFileTumorBam,
                sampleName = tSN
        }

        call general.BamStats as normalBamStats {
            input:
                inFileBam = TNmapping.outFileNormalBam,
                sampleName = nSN
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
                tumorSampleName = tSN,
                normalSampleName = nSN,
                sampleName = fON,
                vardictMinimumAF = 0.01
        }

        call annotate.PickAndAnnotate as vcfAnnotate {
            input:
                inFileVcfSS = variantCalling.outFileSomaticsniperFilteredVcfGz,
                inFileVcfMU = variantCalling.outFileMuseFilteredVcfGz,
                inFileVcfM2 = variantCalling.outFileMutect2FilteredVcfGz,
                inFileVcfLF = variantCalling.outFileLofreqFilteredVcfGz,
                inFileVcfVD = variantCalling.outFileVardictFilteredVcfGz,
                infileVcfVS = variantCalling.outFileVarscanFilteredVcfGz,
                inDirPCGRref = inDirPCGRref,
                refFa = refFa,
                sampleName = fON
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

        Array[File] outFileStatsTumorBam = tumorBamStats.outFileBamStats
        Array[File] outFileStatsNormalBam = normalBamStats.outFileBamStats

        Array[File] outFileSomaticsniperFilteredVcfGz = variantCalling.outFileSomaticsniperFilteredVcfGz
        Array[File] outFileBamsomaticsniperFilteredVcfIndex = variantCalling.outFileBamsomaticsniperFilteredVcfIndex
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

        Array[File] outFileBamsomaticsniperVcfGz = variantCalling.outFileBamsomaticsniperVcfGz
        Array[File] outFileBamsomaticsniperVcfIndex = variantCalling.outFileBamsomaticsniperVcfIndex
        Array[File] outFileLofreqVcfGz = variantCalling.outFileLofreqVcfGz
        Array[File] outFileLofreqVcfIndex = variantCalling.outFileLofreqVcfIndex
        Array[File] outFileMuseVcfGz = variantCalling.outFileMuseVcfGz
        Array[File] outFileMuseVcfIndex = variantCalling.outFileMuseVcfIndex
        Array[File] outFileMutect2VcfGz = variantCalling.outFileMutect2VcfGz
        Array[File] outFileMutect2VcfIndex = variantCalling.outFileMutect2VcfIndex
        Array[File] outFileVardictVcfGz = variantCalling.outFileVardictVcfGz
        Array[File] outFileVardictVcfIndex = variantCalling.outFileVardictVcfIndex
        Array[File] outFileVarscanVcfGz = variantCalling.outFileVarscanVcfGz
        Array[File] outFileVarscanVcfIndex = variantCalling.outFileVarscanVcfIndex

        Array[File] outFilePCGRannotatedVcf = vcfAnnotate.outFilePCGRannotatedVcf
        Array[File] outFilePCGRannotatedVcfIndex = vcfAnnotate.outFilePCGRannotatedVcfIndex

        Array[File] outFileMaf = vcfAnnotate.outFileMaf
        Array[File] outFileCsv = vcfAnnotate.outFileCsv

        Array[File] outFilePCGRflexdbHtml = vcfAnnotate.outFilePCGRflexdbHtml
        Array[File] outFilePCGRhtml = vcfAnnotate.outFilePCGRhtml
    }
}