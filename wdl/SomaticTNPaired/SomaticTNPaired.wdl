version 1.0

import "QCStats.wdl" as qc
import "Preprocessing/Preprocessing.wdl" as prep
import "VariantCalling/VariantCalling.wdl" as caller
import "Annotate.wdl" as annot


workflow SomaticTNPaired {
    input {
        Array[Array[File]] inFileTumorFastqs
        Array[Array[File]] inFileNormalFastqs
        Array[File] inFileIntervalBeds
        File inFileDbsnpVcf
        File inFileDbsnpVcfIndex
        File inFileGermlineResource
        File inFileGermlineResourceIndex
        File inFilePON
        File inFilePONindex
        File inFileVepRef
        File inDirPcgrRef
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
        File inFileIntervalBed = inFileIntervalBeds[i]
        String tumorSampleName = tumorSampleNames[i]
        String normalSampleName = normalSampleNames[i]
        String finalOutputName = finalOutputNames[i]

        call prep.Preprocessing as preprocessing {
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

        call qc.QCStats as qcstats {
            input:
                inFileTumorFastqR1 = iFTFs[0],
                inFileTumorFastqR2 = iFTFs[1],
                inFileTumorBam = preprocessing.outFileTumorBam,
                tumorSampleName = tumorSampleName,
                inFileNormalFastqR1 = iFNFs[0],
                inFileNormalFastqR2 = iFNFs[1],
                inFileNormalBam = preprocessing.outFileNormalBam,
                normalSampleName = normalSampleName
        }

        call caller.VariantCalling as variantCalling {
            input:
                inFileTumorBam = preprocessing.outFileTumorBam,
                inFileTumorBamIndex = preprocessing.outFileTumorBamIndex,
                inFileNormalBam = preprocessing.outFileNormalBam,
                inFileNormalBamIndex = preprocessing.outFileNormalBamIndex,
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
                sampleName = finalOutputName
        }

        call annot.Annotate as variantAnnotation {
            input:
                inFileVcfGz = variantCalling.outFilePickedVcfGz,
                inFileVcfIndex = variantCalling.outFilePickedVcfIndex,
                refFa = refFa,
                inFileVepRef = inFileVepRef,
                inDirPcgrRef = inDirPcgrRef,
                tumorSampleName = tumorSampleName,
                normalSampleName = normalSampleName,
        }
    }
    
    output {
        Array[Array[File]] outFileTumorTrimmedFastqs = preprocessing.outFileTumorTrimmedFastqs
        Array[Array[File]] outFileNormalTrimmedFastqs = preprocessing.outFileNormalTrimmedFastqs

        Array[File] outFileTumorBam = preprocessing.outFileTumorBam
        Array[File] outFileNormalBam = preprocessing.outFileNormalBam
        Array[File] outFileTumorBamIndex = preprocessing.outFileTumorBamIndex
        Array[File] outFileNormalBamIndex = preprocessing.outFileNormalBamIndex

        Array[File] outFileTumorSortedRawBam = preprocessing.outFileTumorSortedRawBam
        Array[File] outFileNormalSortedRawBam = preprocessing.outFileNormalSortedRawBam

        Array[File] outFileTumorFastqcReportTar = qcstats.outFileTumorFastqcReportTar
        Array[File] outFileNormalFastqcReportTar = qcstats.outFileTumorFastqcReportTar
        Array[File] outFileTumorBamStats = qcstats.outFileTumorBamStats
        Array[File] outFileNormalBamStats = qcstats.outFileNormalBamStats

        Array[File] outFileLofreqFilteredVcfGz = variantCalling.outFileLofreqFilteredVcfGz
        Array[File] outFileMutect2FilteredVcfGz = variantCalling.outFileMutect2FilteredVcfGz
        Array[File] outFileMuseFilteredVcfGz = variantCalling.outFileMuseFilteredVcfGz
        Array[File] outFilePickedVcfGz = variantCalling.outFilePickedVcfGz

        Array[File] outFileVepVcfGz = variantAnnotation.outFileVepVcfGz
        Array[File] outFileVepMaf = variantAnnotation.outFileVepMaf
        Array[File] outFileVepCsv = variantAnnotation.outFileVepCsv

        Array[File] outFilePcgrVcfGz = variantAnnotation.outFilePcgrVcfGz
        Array[File] outFilePcgrFlexdbHtml = variantAnnotation.outFilePcgrFlexdbHtml
        Array[File] outFilePcgrHtml = variantAnnotation.outFilePcgrHtml
    }
}