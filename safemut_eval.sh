#!/usr/bin/env bash

outdir=$(pwd -P)

Usage(){
    echo -e "\nUsage: $0 [options]
    -n nbam          normal bam file, required
    -v vcf           mutations vcf to simulate, required
    -r region.bed    target region, required
    -l umi_length    UMI length (Duplex), run consensus if provided
    -s sample        sample ID, required
    -o outdir        result directory, default is current directory
    -h help          show this help message and exit"
    exit 0
}

while getopts :n:v:r:s:l:o: opt; do
    case ${opt} in
    n)
        nbam="${OPTARG}"
        ;;
    v)
        vcf="${OPTARG}"
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
    s)
        sample="${OPTARG}"
        ;;
    ?)
        Usage
        ;;
    esac
done

if [[ -z "${nbam}" ]] || [[ -z "${vcf}" ]] || [[ -z "${sample}" ]] || [[ -z "${bed}" ]]; then
echo "command is imcomplete!"
Usage
fi

mkdir -p ${outdir}

# run safemut
safemut -b ${nbam} -v ${vcf} -1 ${outdir}/${sample}_sim_1.fq.gz -2 ${outdir}/${sample}_sim_2.fq.gz
gzip -d ${outdir}/${sample}_sim_1.fq.gz
gzip -d ${outdir}/${sample}_sim_2.fq.gz
fastq-sort -n -S 10G ${outdir}/${sample}_sim_1.fq > ${outdir}/${sample}_sim_1.sort.fq
fastq-sort -n -S 10G ${outdir}/${sample}_sim_2.fq > ${outdir}/${sample}_sim_2.sort.fq
gzip ${outdir}/${sample}_sim_1.sort.fq
gzip ${outdir}/${sample}_sim_2.sort.fq

if [[ -z "${lumi}" ]];then
# mapping & mpileup
bwa mem -t 8 -R "@RG\tID:${sample}_sim\tSM:${sample}_sim\tLB:${sample}_sim\tPU:${sample}_sim.L1\tPL:ILLUMINA" ${ref} ${outdir}/${sample}_sim_1.sort.fq.gz ${outdir}/${sample}_sim_2.sort.fq.gz |samtools sort -@ 4 -o ${outdir}/${sample}_sim.bam
samtools index -@ 4 ${outdir}/${sample}_sim.bam
bcftools mpileup -R ${bed} --max-depth 99999 --max-idepth 99999 --min-BQ 0 -f ${ref} -a AD,ADF,ADR,DP,SP ${outdir}/${sample}_sim.bam |bcftools norm -f ${ref} -m- -Oz -o ${outdir}/${sample}_sim.mpileup.vcf.gz
bcftools index ${outdir}/${sample}_sim.mpileup.vcf.gz

else
# consensus & mpileup
$(dirname "$0")/consensus.sh -i ${outdir}/${sample}_sim_1.sort.fq.gz -j ${outdir}/${sample}_sim_2.sort.gz -l ${lumi} -s ${sample}_sim -f ${ref} -o ${outdir}
samtools fastq ${outdir}/${sample}_sim.consensus.bam |bwa mem -t 8 -p ${ref} - |samtools sort -@ 4 -o ${outdir}/${sample}_sim.consensus.align.bam
samtools index ${outdir}/${sample}_sim.consensus.align.bam
bcftools mpileup -R ${bed} --max-depth 99999 --max-idepth 99999 --min-BQ 0 -f ${ref} -a AD,ADF,ADR,DP,SP ${outdir}/${sample}_sim.consensus.align.bam |bcftools norm -f ${ref} -m- -Oz -o ${outdir}/${sample}_sim.mpileup.vcf.gz
bcftools index ${outdir}/${sample}_sim.mpileup.vcf.gz
fi
