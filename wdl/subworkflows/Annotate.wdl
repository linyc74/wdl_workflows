version 1.0


workflow Annotate {
    input {
        File inFileVcfGz
        File inDirPCGRref
        String sampleName
    }

    call PCGR {
        input:
            inFileVcfGz = inFileVcfGz,
            inDirPCGRref = inDirPCGRref,
            sampleName = sampleName
    }

    call Vcf2Csv {
        input:
            inFileVcfGz = PCGR.outFileVcfGz,
            sampleName = sampleName
    }

    output {
        File outFilePCGRannotatedVcf = PCGR.outFileVcfGz
        File outFilePCGRannotatedVcfIndex = PCGR.outFileVcfIndex
        File outFileMaf = PCGR.outFileMaf
        File outFilePCGRflexdbHtml = PCGR.outFileFlexdbHtml
        File outFilePCGRhtml = PCGR.outFileHtml
        File outFileCsv = Vcf2Csv.outFileCsv
    }
}


task PCGR {
    input {
        File inFileVcfGz
        File inDirPCGRref
        String sampleName
    }

    command <<<
        set -e -o pipefail
        pcgr \
        --input_vcf ~{inFileVcfGz} \
        --pcgr_dir ~{inDirPCGRref} \
        --output_dir pcgr_output \
        --genome_assembly grch38 \
        --vep_buffer_size 30000 \
        --sample_id ~{sampleName}
        zcat pcgr_output/~{sampleName}.pcgr_acmg.grch38.vcf.gz > file_for_vcf2maf.vcf
        vcf2maf.pl \
        --input-vcf file_for_vcf2maf.vcf \
        --output-maf ~{sampleName}.maf \
        --tumor-id ~{sampleName} \
        --inhibit-vep \
        --ncbi-build GRCh38 \
        --ref-fasta ~{inDirPCGRref}/data/grch38/.vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    >>>

    output {
        File outFileVcfGz = "pcgr_output/~{sampleName}.pcgr_acmg.grch38.vcf.gz"
        File outFileVcfIndex = "pcgr_output/~{sampleName}.pcgr_acmg.grch38.vcf.gz.tbi"
        File outFileMaf = "~{sampleName}.maf"
        File outFileFlexdbHtml = "pcgr_output/~{sampleName}.pcgr_acmg.grch38.flexdb.html"
        File outFileHtml = "pcgr_output/~{sampleName}.pcgr_acmg.grch38.html"
    }

    runtime {
        docker: 'nycu:latest'
    }
}



task Vcf2Csv {

    input {
        File inFileVcfGz
        String sampleName
    }

    command <<<
        set -e -o pipefail
        python /usr/local/seqslab/omic vcf2csv \
        --input-vcf ~{inFileVcfGz}
        --output-csv ~{sampleName}.csv
    >>>

    output {
        File outFileCsv = "~{sampleName}.csv"
    }

    runtime {
        docker: 'nycu:latest'
    }
}