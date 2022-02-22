#!/usr/bin/env bash

outdir=$(pwd -P)

Usage(){
    echo -e "\nUsage: $0 [options]
    -n nbam          normal bam file, required
    -m mut           mutations file to simulate, required
    -r region.bed    target region, required
    -l umi_length    UMI length (Duplex), run consensus if provided
    -s sample        sample ID, required
    -o outdir        result directory, default is current directory
    -h help          show this help message and exit"
    exit 0
}

while getopts :n:m:r:s:l:o: opt; do
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

if [[ -z "${nbam}" ]] || [[ -z "${mut}" ]] || [[ -z "${sample}" ]] || [[ -z "${bed}" ]]; then
echo "command is imcomplete!"
Usage
fi

mkdir -p ${outdir}

# run varben
samtools view -@ 8 -b -F 0x800 -o ${outdir}/${sample}.rmhc.bam $nbam
samtools index -@ 8 ${outdir}/${sample}.rmhc.bam
${python2} ${muteditor} -m ${vcf} -b ${outdir}/${sample}.rmhc.bam -r ${REF} --mindepth 10 --minmutreads 1 -p 8 --aligner bwa --alignerIndex ${REF} --seqer illumina -o ${outdir}

if [[ -z "${lumi}" ]];then
# bcftools mpileup
bcftools mpileup -R ${bed} --max-depth 99999 --max-idepth 99999 --min-BQ 0 -f ${REF} -a AD,ADF,ADR,DP,SP ${outdir}/edit.sorted.bam |bcftools norm -f ${REF} -m- -Oz -o ${outdir}/${sample}_sim.mpileup.vcf.gz
bcftools index -f ${outdir}/${sample}_sim.mpileup.vcf.gz

else
# consensus & mpileup
#bam to fastq
$(dirname "$0")/consensus.sh -i ${outdir}/${sample}_sim_1.fq.gz -j ${outdir}/${sample}_sim_2.gz -l ${lumi} -s ${sample}_sim -o ${outdir}
bcftools mpileup -R ${bed} --max-depth 99999 --max-idepth 99999 --min-BQ 0 -f ${REF} -a AD,ADF,ADR,DP,SP ${outdir}/${sample}_sim.consensus.align.bam |bcftools norm -f ${REF} -m- -Oz -o ${outdir}/${sample}_sim.mpileup.vcf.gz
bcftools index -f ${outdir}/${sample}_sim.mpileup.vcf.gz
fi
