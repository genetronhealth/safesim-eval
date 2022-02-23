#!/usr/bin/env bash
# consensus duplex sequencing data

outdir=$(pwd -P)

while getopts :i:j:l:s:o: opt; do
    case ${opt} in
    i)
        fq1="${OPTARG}"
        ;;
    j)
        fq2="${OPTARG}"
        ;;
    l)
        lumi="${OPTARG}"
        ;;
    o)
        outdir="${OPTARG}"
        ;;
    s)
        sample="${OPTARG}"
        ;;
    ?)
        echo "Usage: $0 -i R1.fq.gz -j R2.fq.gz -l umi_length -s sample [-o outdir] "
        exit 0
        ;;
    esac
done

test -d ${outdir} && echo "exist consensus results must be removed!"
mkdir -p ${outdir}
${java} -Xmx8G -jar $picard FastqToSam --FASTQ ${fq1} --FASTQ2 ${fq2} --OUTPUT ${outdir}/${sample}.unalign.bam --READ_GROUP_NAME ${sample} --SAMPLE_NAME ${sample} --LIBRARY_NAME ${sample} --PLATFORM_UNIT ${sample}.L1 --PLATFORM ILLUMINA --TMP_DIR $TMP
${java} -Xmx8G -jar $fgbio --tmp-dir=$TMP ExtractUmisFromBam -i ${outdir}/${sample}.unalign.bam -o ${outdir}/${sample}.umi.unalign.bam -r ${lumi}M+T ${lumi}M+T -t ZA ZB -s RX
samtools fastq ${outdir}/${sample}.umi.unalign.bam |bwa mem -t 8 -p ${REF} - |samtools view -@ 4 -b -o ${outdir}/${sample}.umi.align.bam

${java} -Xmx8G -jar $picard MergeBamAlignment -R ${REF} -ALIGNED ${outdir}/${sample}.umi.align.bam -UNMAPPED ${outdir}/${sample}.umi.unalign.bam -O ${outdir}/${sample}.merge.align.bam --MAX_GAPS -1 --ATTRIBUTES_TO_RETAIN XS,XA --ALIGNER_PROPER_PAIR_FLAGS true --TMP_DIR $TMP
${java} -Xmx8G -jar $fgbio --tmp-dir=$TMP GroupReadsByUmi -s Paired -i ${outdir}/${sample}.merge.align.bam -o ${outdir}/${sample}.group.bam -t RX
${java} -Xmx8G -jar $fgbio --tmp-dir=$TMP CallDuplexConsensusReads -i ${outdir}/${sample}.group.bam -o ${outdir}/${sample}.consensus.bam
samtools fastq ${outdir}/${sample}.consensus.bam |bwa mem -t 8 -p ${REF} -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPU:${sample}.L1\tPL:ILLUMINA" - |samtools sort -@ 4 -o ${outdir}/${sample}.consensus.align.bam
samtools index ${outdir}/${sample}.consensus.align.bam
