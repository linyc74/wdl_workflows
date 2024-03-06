version 1.0


workflow TrimAndMap {
    input {
        Array[File] inFileFastqPair
        File refAmb
        File refAnn
        File refBwt
        File refPac
        File refSa
        File refFa
        File refFai
        String sampleName
    }

    call TrimGalore {
        input:
            sampleName = sampleName,
            inFileFastqR1 = inFileFastqPair[0],
            inFileFastqR2 = inFileFastqPair[1]
    }

    call BwaMemMapping as mapping {
        input:
            inFileFastqR1 = TrimGalore.outFileTrimmedFastqR1,
            inFileFastqR2 = TrimGalore.outFileTrimmedFastqR2,
            sampleName = sampleName,
            refAmb = refAmb,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refSa = refSa,
            refFa = refFa,
            refFai = refFai
    }
 
    output {
        Array[File] outFileTrimmedFastq = TrimGalore.outFileTrimmedFastqs
        File outFileUnSortRawBam = mapping.outFileUnSortRawBam
    }
}


task TrimGalore {
    input {
        File inFileFastqR1
        File inFileFastqR2
        String sampleName
    }

    command <<<
        set -e -o pipefail
        trim_galore \
        --paired \
        --quality 20 \
        --phred33 \
        --illumina \
        --length 20 \
        --max_n 0 \
        --trim-n \
        --gzip \
        --output_dir out \
        --basename ~{sampleName} \
        ~{inFileFastqR1} \
        ~{inFileFastqR2}
    >>>

    runtime {
        docker: 'nycu:latest'
    }

    output {
        File outFileTrimmedFastqR1 = "out/~{sampleName}_val_1.fq.gz"
        File outFileTrimmedFastqR2 = "out/~{sampleName}_val_2.fq.gz"
    }
}


task BwaMemMapping {
    input {
        File inFileFastqR1
        File inFileFastqR2
        String sampleName
        File refAmb
        File refAnn
        File refBwt
        File refPac
        File refSa
        File refFa
        File refFai
    }
 
    command <<<
        set -e -o pipefail
        bwa mem \
        -t 4 \
        -R "@RG\tID:~{sampleName}\tSM:~{sampleName}\tPL:ILLUMINA\tLB:~{sampleName}" \
        ~{refFa} \
        ~{inFileFastqR1} \
        ~{inFileFastqR2} \
        | \
        samtools view -b - > ~{sampleName}.bam  
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileUnSortRawBam = "~{sampleName}.bam"
    }
}
