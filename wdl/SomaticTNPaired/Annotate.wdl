version 1.0


workflow Annotate {
    input {
        File inFileVcfGz
        File inFileVcfIndex
        File refFa
        File refVepTarGz
        File refPcgrDir
        String tumorSampleName
        String normalSampleName
    }

    call VEP {
        input:
            inFileVcfGz = inFileVcfGz,
            inFileVcfIndex = inFileVcfIndex,
            refFa = refFa,
            refVepTarGz = refVepTarGz,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName
    }

    call PCGR {
        input:
            inFileVcfGz = inFileVcfGz,
            inFileVcfIndex = inFileVcfIndex,
            refPcgrDir = refPcgrDir,
            sampleName = tumorSampleName
    }

    output {
        File outFileVepVcfGz = VEP.outFileVcfGz
        File outFileVepVcfIndex = VEP.outFileVcfIndex
        File outFileVepMaf = VEP.outFileMaf
        File outFileVepCsv = VEP.outFileCsv

        File outFilePcgrVcfGz = PCGR.outFileVcfGz
        File outFilePcgrVcfIndex = PCGR.outFileVcfIndex
        File outFilePcgrFlexdbHtml = PCGR.outFileFlexdbHtml
        File outFilePcgrHtml = PCGR.outFileHtml
    }
}


task VEP {
    input {
        File inFileVcfGz
        File inFileVcfIndex
        File refFa
        File refVepTarGz
        String tumorSampleName
        String normalSampleName
    }

    command <<<
        set -e -o pipefail
        mkdir ./vep-cache
        tar -xvzf ~{refVepTarGz} -C "./vep-cache"
        vep \
          --offline \
          --input_file ~{inFileVcfGz} \
          --fasta ~{refFa} \
          --dir_cache ./vep-cache \
          --merged \
          --everything \
          --fork 1 \
          --buffer_size 10000 \
          --vcf \
          --output_file ~{tumorSampleName}_vep.vcf
        vcf2maf.pl \
          --input-vcf ~{tumorSampleName}_vep.vcf \
          --ref-fasta ~{refFa} \
          --tmp-dir vcf2maf-tmp-dir \
          --output-maf ~{tumorSampleName}_vep.maf \
          --tumor-id ~{tumorSampleName} \
          --normal-id ~{normalSampleName} \
          --inhibit-vep \
          --species homo_sapiens \
          --ncbi-build GRCh38
        python /usr/local/seqslab/omic vcf2csv \
          --input-vcf ~{tumorSampleName}_vep.vcf \
          --output-csv ~{tumorSampleName}_vep.csv
        bgzip --stdout ~{tumorSampleName}_vep.vcf > ~{tumorSampleName}_vep.vcf.gz
        tabix --preset vcf ~{tumorSampleName}_vep.vcf.gz
    >>>

    output {
        File outFileVcfGz = "~{tumorSampleName}_vep.vcf.gz"
        File outFileVcfIndex = "~{tumorSampleName}_vep.vcf.gz.tbi"
        File outFileMaf = "~{tumorSampleName}_vep.maf"
        File outFileCsv = "~{tumorSampleName}_vep.csv"
    }

    runtime {
        docker: 'nycu:latest'
    }
}


task PCGR {
    input {
        File inFileVcfGz
        File inFileVcfIndex
        File refPcgrDir
        String sampleName
    }

    command <<<
        set -e -o pipefail
        pcgr \
        --input_vcf ~{inFileVcfGz} \
        --pcgr_dir ~{refPcgrDir} \
        --output_dir pcgr_output \
        --genome_assembly grch38 \
        --vep_buffer_size 30000 \
        --sample_id ~{sampleName}
    >>>

    output {
        File outFileVcfGz = "pcgr_output/~{sampleName}_pcgr.vcf.gz"
        File outFileVcfIndex = "pcgr_output/~{sampleName}_pcgr.vcf.gz.tbi"
        File outFileFlexdbHtml = "pcgr_output/~{sampleName}_pcgr.flexdb.html"
        File outFileHtml = "pcgr_output/~{sampleName}_pcgr.grch38.html"
    }

    runtime {
        docker: 'nycu:latest'
    }
}