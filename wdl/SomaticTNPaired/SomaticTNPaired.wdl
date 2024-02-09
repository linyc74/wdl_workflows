version 1.0

import "QCStats.wdl" as qc
import "Preprocessing/Preprocessing.wdl" as prep
import "VariantCalling/VariantCalling.wdl" as caller
import "Annotate.wdl" as annot


workflow SomaticTNPaired {
    input {
        Array[Array[File]] inFileTumorFastqPairs
        Array[Array[File]] inFileNormalFastqPairs
        Array[File] refIntervalBeds
        File refDbsnpVcfGz
        File refDbsnpVcfIndex
        File refGermlineResourceVcfGz
        File refGermlineResourceVcfIndex
        File refPonVcfGz
        File refPonVcfIndex
        File refVepTarGz
        File refPcgrDir
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
    }

    scatter (i in range(length(tumorSampleNames))) {
        Array[File] inFileTumorFastqPair = inFileTumorFastqPairs[i]
        Array[File] inFileNormalFastqPair = inFileNormalFastqPairs[i]
        File refIntervalBed = refIntervalBeds[i]
        String tumorSampleName = tumorSampleNames[i]
        String normalSampleName = normalSampleNames[i]

        call prep.Preprocessing as preprocessing {
            input:
                inFileTumorFastqPair = inFileTumorFastqPair,
                inFileNormalFastqPair = inFileNormalFastqPair,
                refDbsnpVcfGz = refDbsnpVcfGz,
                refDbsnpVcfIndex = refDbsnpVcfIndex,
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
                inFileTumorFastqPair = inFileTumorFastqPair,
                inFileTumorBam = preprocessing.outFileTumorBam,
                tumorSampleName = tumorSampleName,
                inFileNormalFastqPair = inFileNormalFastqPair,
                inFileNormalBam = preprocessing.outFileNormalBam,
                normalSampleName = normalSampleName
        }

        call caller.VariantCalling as variantCalling {
            input:
                inFileTumorBam = preprocessing.outFileTumorBam,
                inFileTumorBamIndex = preprocessing.outFileTumorBamIndex,
                inFileNormalBam = preprocessing.outFileNormalBam,
                inFileNormalBamIndex = preprocessing.outFileNormalBamIndex,
                refIntervalBed = refIntervalBed,
                refGermlineResourceVcfGz = refGermlineResourceVcfGz,
                refGermlineResourceVcfIndex = refGermlineResourceVcfIndex,
                refPonVcfGz = refPonVcfGz,
                refPonVcfIndex = refPonVcfIndex,
                refFa = refFa,
                refFai = refFai,
                refDict = refDict,
                tumorSampleName = tumorSampleName,
                normalSampleName = normalSampleName
        }

        call annot.Annotate as annotate {
            input:
                inFileVcfGz = variantCalling.outFilePickedVcfGz,
                inFileVcfIndex = variantCalling.outFilePickedVcfIndex,
                refFa = refFa,
                refVepTarGz = refVepTarGz,
                refPcgrDir = refPcgrDir,
                tumorSampleName = tumorSampleName,
                normalSampleName = normalSampleName,
        }
    }
    
    output {
        Array[Array[File]] outFileTumorTrimmedFastqs = preprocessing.outFileTumorTrimmedFastqs
        Array[Array[File]] outFileNormalTrimmedFastqs = preprocessing.outFileNormalTrimmedFastqs
        Array[File] outFileTumorBam = preprocessing.outFileTumorBam
        Array[File] outFileNormalBam = preprocessing.outFileNormalBam
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

        Array[File] outFileVepVcfGz = annotate.outFileVepVcfGz
        Array[File] outFileVepMaf = annotate.outFileVepMaf
        Array[File] outFileVepCsv = annotate.outFileVepCsv
        Array[File] outFilePcgrVcfGz = annotate.outFilePcgrVcfGz
        Array[File] outFilePcgrFlexdbHtml = annotate.outFilePcgrFlexdbHtml
        Array[File] outFilePcgrHtml = annotate.outFilePcgrHtml
    }
}