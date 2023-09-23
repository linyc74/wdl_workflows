version 1.0

import "subworkflows/GeneralTask.wdl" as general
import "subworkflows/Mapping/subworkflows/TrimAndMapping.wdl" as ponMapping
import "subworkflows/Mapping/subworkflows/PostMappingProcess.wdl" as ponPostMapping
import "subworkflows/Mapping/TNpairedMapping.wdl" as mapper
import "subworkflows/VariantsCalling/TNpairedVariantsCalling.wdl" as caller
import "subworkflows/PickAndAnnotate.wdl" as annotate
import "subworkflows/CreateMutect2PoN.wdl" as createPoN

# WORKFLOW DEFINITION

# 'Workflow description'
workflow ScatterInProcessPoNSomaticpipelineTumorNormalMode {
    input {
        Array[Array[File]] inFileTumorFastqs
        Array[Array[File]] inFileNormalFastqs
        Array[Array[File]] inFilePoNfastqs
        File inFileDbsnpVcf
        File inFileDbsnpVcfIndex
        File inFileIntervalBed
        File inFileGermlineResource
        File inFileGermlineResourceIndex
        File inDirPCGRref
        File refAmb
        File refAnn
        File refBwt
        File refPac
        File refSa
        File refFa
        File refFai
        File refDict
        String ponName
        Array[String] ponSampleName
        Array[String] tumorSampleName
        Array[String] normalSampleName
        Array[String] finalOutputName
    }
 
    scatter (i in range(length(ponSampleName))) {
        Array[File] iFPFs = inFilePoNfastqs[i]
        String pSN = ponSampleName[i]

        call general.FastQC as fastqcPoNfastq {
            input:
                inFileFastqR1 = iFPFs[0],
                inFileFastqR2 = iFPFs[1]
    }

        call ponMapping.TrimAndMapping {
            input:
                refAmb = refAmb,
                refAnn = refAnn,
                refBwt = refBwt,
                refPac = refPac,
                refSa = refSa,
                refFa = refFa,
                refFai = refFai,
                sampleName = pSN,
                inFileFastqs = iFPFs
        }

        call ponPostMapping.PostMappingProcess {
            input:
                inFileUnSortRawBam = TrimAndMapping.outFileUnSortRawBam,
                inFileDbsnpVcf = inFileDbsnpVcf,
                inFileDbsnpVcfIndex = inFileDbsnpVcfIndex,
                refFa = refFa,
                refFai = refFai,
                refDict = refDict,
                sampleName = pSN
        }

        call general.BamStats as ponBamStats {
            input:
                inFileBam = PostMappingProcess.outFileBam,
                sampleName = pSN
        }
    }

    call createPoN.CreateMutect2PoN {
        input:
            inFileNormalBams = PostMappingProcess.outFileBam,
            inFileNormalBamIndexs = PostMappingProcess.outFileBamIndex,
            inFileGermlineResource = inFileGermlineResource,
            inFileGermlineResourceIndex = inFileGermlineResourceIndex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            normalSampleName = ponSampleName,
            extraArgs = "--max-mnp-distance 0",
            ponName = ponName
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
                inFilePON = CreateMutect2PoN.outFilePoNvcf,
                inFilePONindex = CreateMutect2PoN.outFilePoNvcfIndex,
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
        Array[Array[File]] outFileTrimmedTumorFastqs = TNmapping.outFileTumorTrimmedFastqs
        Array[Array[File]] outFileTrimmedNormalFastqs = TNmapping.outFileNormalTrimmedFastqs
        Array[File] outFileTumorFastqcReportTar = fastqcTumorFastq.outFileFastqcReportTar
        Array[File] outFileNormalFastqcReportTar = fastqcNormalFastq.outFileFastqcReportTar
        Array[File] outFilePoNfastqcReportTar = fastqcPoNfastq.outFileFastqcReportTar
        Array[File] outFileTumorBam = TNmapping.outFileTumorBam
        Array[File] outFileNormalBam = TNmapping.outFileNormalBam
        Array[File] outFileTumorBamIndex = TNmapping.outFileTumorBamIndex
        Array[File] outFileNormalBamIndex = TNmapping.outFileNormalBamIndex
        Array[File] outFileTumorSortedRawBam = TNmapping.outFileTumorSortedRawBam
        Array[File] outFileNormalSortedRawBam = TNmapping.outFileNormalSortedRawBam
        Array[File] outFileStatsTumorBam = tumorBamStats.outFileBamStats
        Array[File] outFileStatsNormalBam = normalBamStats.outFileBamStats
        Array[File] outFileStatsPoNbam = ponBamStats.outFileBamStats        
        Array[File] outFileBamsomaticsniperPyVcfGz = variantCalling.outFileBamsomaticsniperPyVcfGz
        Array[File] outFileBamsomaticsniperPyVcfIndex = variantCalling.outFileBamsomaticsniperPyVcfIndex
        Array[File] outFileLofreqPyVcfGz = variantCalling.outFileLofreqPyVcfGz
        Array[File] outFileLofreqPyVcfIndex = variantCalling.outFileLofreqPyVcfIndex
        Array[File] outFileMusePyVcfGz = variantCalling.outFileMusePyVcfGz
        Array[File] outFileMusePyVcfIndex = variantCalling.outFileMusePyVcfIndex
        Array[File] outFileMutect2PyVcfGz = variantCalling.outFileMutect2PyVcfGz
        Array[File] outFileMutect2PyVcfIndex = variantCalling.outFileMutect2PyVcfIndex
        Array[File] outFileVardictPyVcfGz = variantCalling.outFileVardictPyVcfGz
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
        File outFilePON = CreateMutect2PoN.outFilePoNvcf
        File outFilePONindex = CreateMutect2PoN.outFilePoNvcfIndex
    }
}