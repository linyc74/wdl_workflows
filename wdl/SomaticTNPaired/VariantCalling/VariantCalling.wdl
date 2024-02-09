version 1.0

import "Lofreq.wdl" as lofreq
import "Mutect2.wdl" as mutect2
import "Muse.wdl" as muse
import "GeneralTask.wdl" as general


workflow VariantCalling {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File inFileNormalBam
        File inFileNormalBamIndex
        File inFileIntervalBed
        File inFileGermlineResource
        File inFileGermlineResourceIndex
        File? inFilePON
        File? inFilePONindex
        File refFa
        File refFai
        File refDict
        String tumorSampleName
        String normalSampleName
    }

    call lofreq.Lofreq as lofreq {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            sampleName = tumorSampleName
    }

    call mutect2.Mutect2 as mutect2 {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            inFileGermlineResource = inFileGermlineResource,
            inFileGermlineResourceIndex = inFileGermlineResourceIndex,
            inFilePON = inFilePON,
            inFilePONindex = inFilePONindex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = tumorSampleName
    }

    call muse.Muse as muse {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            refFa = refFa,
            refFai = refFai,
            sampleName = tumorSampleName
    }

    call general.VariantPicking as variantPicking {
        input:
            refFa = refFa,
            sampleName = tumorSampleName,
            inFileVcfLofreq = lofreq.outFileFilteredVcfGz,
            inFileVcfMutect2 = mutect2.outFileFilteredVcfGz,
            inFileVcfMuse = muse.outFileFilteredVcfGz
    }

    output {
        File outFileLofreqVcfGz = lofreq.outFileVcfGz
        File outFileLofreqVcfIndex = lofreq.outFileVcfIndex
        File outFileMuseVcfGz = muse.outFileVcfGz
        File outFileMuseVcfIndex = muse.outFileVcfIndex
        File outFileMutect2VcfGz = mutect2.outFileVcfGz
        File outFileMutect2VcfIndex = mutect2.outFileVcfIndex

        File outFileLofreqFilteredVcfGz = lofreq.outFileFilteredVcfGz
        File outFileLofreqFilteredVcfIndex = lofreq.outFileFilteredVcfIndex
        File outFileMuseFilteredVcfGz = muse.outFileFilteredVcfGz
        File outFileMuseFilteredVcfIndex = muse.outFileFilteredVcfIndex
        File outFileMutect2FilteredVcfGz = mutect2.outFileFilteredVcfGz
        File outFileMutect2FilteredVcfIndex = mutect2.outFileFilteredVcfIndex

        File outFilePickedVcfGz = variantPicking.outFileVcfGz
        File outFilePickedVcfIndex = variantPicking.outFileVcfIndex
    }
}