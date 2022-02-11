#!/usr/bin/env bash

outdir=$(pwd -P)

Usage(){
    echo -e "\nUsage: $0 [options]
    -n nbam          normal bam file, required
    -s snv           snv mutations file to simulate, required
    -i indel         indel mutations file to simulate, required
    -r region.bed    target region, required
    -l umi_length    UMI length (Duplex), run consensus if provided
    -S sample        sample ID, required
    -o outdir        result directory, default is current directory
    -h help          show this help message and exit"
    exit 0
}

while getopts :n:s:i:r:l:S:o: opt; do
    case ${opt} in
    n)
        nbam="${OPTARG}"
        ;;
    s)
        snv="${OPTARG}"
        ;;
    i)
        indel="${OPTARG}"
        ;;
    r)
        bed="${OPTARG}"
        ;;
    l)
        lumi="${OPTARG}"
        ;;
    o)
        outdir="${OPTARG}"
        ;;
    S)
        sample="${OPTARG}"
        ;;
    ?)
        Usage
        ;;
    esac
done

if [[ -z "${nbam}" ]] || [[ -z "${snv}" ]] || [[ -z "${indel}" ]] || [[ -z "${sample}" ]] || [[ -z "${bed}" ]]; then
echo "command is imcomplete!"
Usage
fi

mkdir -p ${outdir}

# run bamsurgeon
cd $outdir

${python3} ${addsnv} -v ${snv} -f ${nbam} -r ${ref} -o ${sample}.simsnv.bam -p 8 --aligner mem --picardjar ${picard} --mindepth 10 --maxdepth 1000000 --coverdiff 0.5 --tmpdir addsnv.tmp
${python3} ${addindel} -v ${indel} -f ${nbam} -r ${ref} -o ${sample}.simindel.bam -p 8 --aligner mem --picardjar ${picard} --mindepth 10 --maxdepth 1000000 --tmpdir addindel.tmp
samtools sort -@ 8 -o ${sample}.simsnv.sort.bam ${sample}.simsnv.bam
samtools sort -@ 8 -o ${sample}.simindel.sort.bam ${sample}.simindel.bam
samtools index -@ 8 ${sample}.simsnv.sort.bam
samtools index -@ 8 ${sample}.simindel.sort.bam

if [[ -z "${lumi}" ]];then
# bcftools mpileup
bcftools mpileup -R ${bed} --max-depth 99999 --max-idepth 99999 --min-BQ 0 -f ${ref} -a AD,ADF,ADR,DP,SP ${sample}.simsnv.sort.bam |bcftools norm -f ${ref} -m- -Oz -o ${sample}_sim.snv.mpileup.vcf.gz
bcftools index ${sample}_sim.snv.mpileup.vcf.gz
bcftools mpileup -R ${bed} --max-depth 99999 --max-idepth 99999 --min-BQ 0 -f ${ref} -a AD,ADF,ADR,DP,SP ${sample}.simindel.sort.bam |bcftools norm -f ${ref} -m- -Oz -o ${sample}_sim.indel.mpileup.vcf.gz
bcftools index ${sample}_sim.indel.mpileup.vcf.gz
bcftools concat -aD ${sample}_sim.snv.mpileup.vcf.gz ${sample}_sim.indel.mpileup.vcf.gz -Oz -o ${sample}_sim.mpileup.vcf.gz
bcftools index ${sample}_sim.mpileup.vcf.gz

else
# consensus & mpileup
$(dirname "$0")/consensus.sh -i ${sample}_sim_1.fq.gz -j ${sample}_sim_2.gz -l ${lumi} -s ${sample}_sim
samtools fastq ${sample}_sim.consensus.bam |bwa mem -t 8 -p ${ref} - |samtools sort -@ 4 -o ${sample}_sim.consensus.align.bam
samtools index ${sample}_sim.consensus.align.bam
bcftools mpileup -R ${bed} --max-depth 99999 --max-idepth 99999 --min-BQ 0 -f ${ref} -a AD,ADF,ADR,DP,SP ${sample}_sim.consensus.align.bam |bcftools norm -f ${ref} -m- -Oz -o ${sample}_sim.mpileup.vcf.gz
bcftools index ${sample}_sim.mpileup.vcf.gz
fi
