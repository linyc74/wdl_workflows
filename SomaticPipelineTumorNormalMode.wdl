version 1.0

import "subworkflows/GeneralTask.wdl" as general
import "subworkflows/Mapping/TNpairedMapping.wdl" as mapper
import "subworkflows/VariantsCalling/TNpairedVariantsCalling.wdl" as caller
import "subworkflows/PickAndAnnotate.wdl" as annotate

# WORKFLOW DEFINITION

# NYCU Dentistry somatic pipeline in Tumor-Normal paired mode
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

        call mapper.TNpairedMapping as TNmapping {
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

        call caller.TNpairedVariantsCalling as variantCalling {
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
                inFileVcfSS = variantCalling.outFileBamsomaticsniperPyVcfGz,
                inFileVcfMU = variantCalling.outFileMusePyVcfGz,
                inFileVcfM2 = variantCalling.outFileMutect2PyVcfGz,
                inFileVcfLF = variantCalling.outFileLofreqPyVcfGz,
                inFileVcfVD = variantCalling.outFileVardictPyVcfGz,
                infileVcfVS = variantCalling.outFileVarscanPyVcfGz,
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
        Array[File] outFileBamsomaticsniperPyVcfGz = variantCalling.outFileBamsomaticsniperPyVcfGz
        Array[File] outFileBamsomaticsniperPyVcfIndex = variantCalling.outFileBamsomaticsniperPyVcfIndex
        Array[File] outFileLofreqPyVcfGz = variantCalling.outFileLofreqPyVcfGz
        Array[File] outFileLofreqPyVcfIndex = variantCalling.outFileLofreqPyVcfIndex
        Array[File] outFileMusePyVcfGz = variantCalling.outFileMusePyVcfGz
        Array[File] outFileMusePyVcfIndex = variantCalling.outFileMusePyVcfIndex
        Array[File] outFileMutect2PyVcfGz = variantCalling.outFileMutect2PyVcfGz
        Array[File] outFileMutect2PyVcfIndex = variantCalling.outFileMutect2PyVcfIndex
        Array[File] outFileVardictPyVcfGz = variantCalling.outFileVardictPyVcfGz
        Array[File] outFileVardictPyVcfIndex = variantCalling.outFileVardictPyVcfIndex
        Array[File] outFileVarscanPyVcfGz = variantCalling.outFileVarscanPyVcfGz
        Array[File] outFileVarscanPyVcfIndex = variantCalling.outFileVarscanPyVcfIndex
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
        Array[File] outFilePCGRflexdbHtml = vcfAnnotate.outFilePCGRflexdbHtml
        Array[File] outFilePCGRhtml = vcfAnnotate.outFilePCGRhtml
    }
}