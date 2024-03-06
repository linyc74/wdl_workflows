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


task TwistMapping {

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
        # (1) Convert Fastq to unaligned BAM
        picard FastqToSam \
            -O unaligned.bam \
            -F1 ~{inFileFastqR1} \
            -F2 ~{inFileFastqR2} \
            -SM ~{sampleName} \
            -LB Library1 \
            -PU Unit1 \
            -PL Illumina
        # (2) Extract UMI bases into unaligned BAM tag
        fgbio ExtractUmisFromBam \
            --input=unaligned.bam \
            --output=unaligned_umi_extracted.bam \
            --read-structure=5M2S+T 5M2S+T \
            --molecular-index-tags=ZA ZB \
            --single-tag=RX
        # (3) Convert unaligned BAM to Fastq
        picard SamToFastq \
            --INPUT unaligned_umi_extracted.bam \
            --FASTQ umi_extracted.fastq \
            --INTERLEAVE true
        # (4) Align Fastq
        bwa mem -p -t 1 \
            ~{refFa} \
            umi_extracted.fastq | \
        samtools sort -@ 1 \
            -o aligned_umi_extracted.bam
        # (5) Merge aligned BAM with unaligned BAM containing UMI tags
        picard MergeBamAlignment \
            --UNMAPPED_BAM unaligned_umi_extracted.bam \
            --ALIGNED_BAM aligned_umi_extracted.bam \
            --OUTPUT aligned_tag_umi.bam \
            --REFERENCE_SEQUENCE ~{refFa} \
            --CLIP_ADAPTERS false \
            --VALIDATION_STRINGENCY SILENT \
            --CREATE_INDEX true \
            --EXPECTED_ORIENTATIONS FR \
            --MAX_GAPS -1 \
            --SORT_ORDER coordinate \
            --ALIGNER_PROPER_PAIR_FLAGS false
        # (6) Optional: Collect Hs metrics prior to consensus calling
        # (7) Group reads by UMI
        fgbio GroupReadsByUmi \
            --strategy=Paired \
            --input=aligned_tag_umi.bam \
            --output=grouped_by_umi.bam \
            --raw-tag=RX \
            --min-map-q=10 \
            --edits=1
        # (8) Call consensus reads
        fgbio CallDuplexConsensusReads \
            --input=grouped_by_umi.bam \
            --output=unaligned_consensus.bam \
            --error-rate-pre-umi=45 \
            --error-rate-post-umi=30 \
            --min-input-base-quality=30 \
            --min-reads 2 1 1
        # (9) Align duplex consensus reads
        picard SamToFastq \
            --INPUT unaligned_consensus.bam \
            --FASTQ consensus.fastq \
            --INTERLEAVE true
        bwa mem -p -t 1 \
            ~{refFa} \
            -o aligned_consensus.bam \
            consensus.fastq
        # --------------------STUCK HERE--------------------
        # (10) Merge unaligned consensus BAM and aligned consensus BAM to retain UMI tag metadata
        picard MergeBamAlignment \
            --UNMAPPED_BAM unaligned_consensus.bam \
            --ALIGNED_BAM aligned_consensus.bam \
            --OUTPUT merged_no_read_group_consensus.bam \
            --REFERENCE_SEQUENCE ./chr9.fa \
            --CLIP_ADAPTERS false \
            --VALIDATION_STRINGENCY SILENT \
            --CREATE_INDEX true \
            --EXPECTED_ORIENTATIONS FR \
            --MAX_GAPS -1 \
            --SORT_ORDER coordinate \
            --ALIGNER_PROPER_PAIR_FLAGS false \
            --ATTRIBUTES_TO_RETAIN X0 \
            --ATTRIBUTES_TO_RETAIN ZS \
            --ATTRIBUTES_TO_RETAIN ZI \
            --ATTRIBUTES_TO_RETAIN ZM \
            --ATTRIBUTES_TO_RETAIN ZC \
            --ATTRIBUTES_TO_RETAIN ZN \
            --ATTRIBUTES_TO_RETAIN ad \
            --ATTRIBUTES_TO_RETAIN bd \
            --ATTRIBUTES_TO_RETAIN cd \
            --ATTRIBUTES_TO_RETAIN ae \
            --ATTRIBUTES_TO_RETAIN be \
            --ATTRIBUTES_TO_RETAIN ce
        picard AddOrReplaceReadGroups \
            I=merged_no_read_group_consensus.bam \
            O=~{sampleName}_twist.bam \
            RGID=~{sampleName} \
            RGLB=~{sampleName} \
            RGPL=Illumina \
            RGSM=~{sampleName} \
            RGPU=NA
    >>>

    runtime {
        docker: 'nycu:latest'
    }

    output {
        File outFile = "~{sampleName}_twist.bam"
    }
}
