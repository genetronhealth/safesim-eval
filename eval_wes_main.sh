#!/usr/bin/env sh
scriptdir=$(dirname $(readlink -e "$0"))/scripts
simfiles=$(dirname $(readlink -e "$0"))/simfiles
init=$(dirname $(readlink -e "$0"))/init_eval.sh

inputdir=$(dirname $(readlink -e "$0"))/Rawdata/WES
outdir=$(dirname $(readlink -e "$0"))/Results/WES

while getopts :i:o: opt; do
    case ${opt} in
    i)
        inputdir="${OPTARG}"
        ;;
    o)
        outdir="${OPTARG}"
        ;;
    ?)
        echo "\nUsage: ./process_wes.sh -i inputdir -o outdir"
        ;;
    esac
done

mkdir -p "${outdir}/scriptlog"
mkdir -p "${outdir}/realdata"
mkdir -p "${outdir}/nbam"

cat $(dirname "$0")/tn_sample.list |grep ^W |while read c t n
do

# realdata
tfq1="${inputdir}/${t}_1.fastq.gz"
tfq2="${inputdir}/${t}_2.fastq.gz"

echo "source ${init} WES
bwa mem -t 8 -R \"@RG\\tID:${t}\\tSM:${t}\\tLB:${t}\\tPU:${t}.L1\\tPL:ILLUMINA\" \${REF} ${tfq1} ${tfq2} |samtools sort -@ 4 -o ${outdir}/realdata/${t}.bam
samtools index -@ 4 ${outdir}/realdata/${t}.bam
bcftools mpileup -R ${simfiles}/wes_simmut.bed --max-depth 99999 --max-idepth 99999 --min-BQ 0 -f \${REF} -a AD,ADF,ADR,DP,SP ${outdir}/realdata/${t}.bam |bcftools norm -f \${REF} -m- -Oz -o ${outdir}/realdata/${t}.mpileup.vcf.gz
bcftools index ${outdir}/realdata/${t}.mpileup.vcf.gz
\${python2} ${scriptdir}/extract_info.py ${outdir}/realdata/${t}.mpileup.vcf.gz ${c} |awk '{print \$0\"\\tSEQC2\"}' > ${outdir}/realdata/${t}_simmut.info.txt" > "${outdir}/scriptlog/${t}.sh"

# simulation
nfq1="${inputdir}/${n}_1.fastq.gz"
nfq2="${inputdir}/${n}_2.fastq.gz"

echo "source ${init} WES
bwa mem -t 8 -R \"@RG\\tID:${n}\\tSM:${n}\\tLB:${n}\\tPU:${n}.L1\\tPL:ILLUMINA\" \${REF} ${nfq1} ${nfq2} |samtools sort -@ 4 -o ${outdir}/nbam/${n}.bam
samtools index -@ 4 ${outdir}/nbam/${n}.bam" > "${outdir}/scriptlog/${n}_bwa.sh"

echo "source ${init} WES
bash ${scriptdir}/safemut_sim.sh -n ${outdir}/nbam/${n}.bam -v ${simfiles}/wes_simmut_safemut.vcf.gz -r ${simfiles}/wes_simmut.bed -s ${n} -o ${outdir}/safemut/${n}
\${python2} ${scriptdir}/extract_info.py ${outdir}/safemut/${n}/${n}_sim.mpileup.vcf.gz ${c} |awk '{print \$0\"\\tSimulation\"}' > ${outdir}/safemut/${n}/${n}_simmut.info.txt" > "${outdir}/scriptlog/${n}.safemut.sh"

echo "source ${init} WES
bash ${scriptdir}/varben_sim.sh -n ${outdir}/nbam/${n}.bam -v ${simfiles}/wes_simmut_varben.tsv -r ${simfiles}/wes_simmut.bed -s ${n} -o ${outdir}/varben/${n}
\${python2} ${scriptdir}/extract_info.py ${outdir}/varben/${n}/${n}_sim.mpileup.vcf.gz ${c} |awk '{print \$0\"\\tSimulation\"}' > ${outdir}/varben/${n}/${n}_simmut.info.txt" > "${outdir}/scriptlog/${n}.varben.sh"

echo "source ${init} WES
bash ${scriptdir}/bamsurgeon_sim.sh -n ${outdir}/nbam/${n}.bam -s ${simfiles}/wes_simmut_bamsurgeon_snv.tsv -i ${simfiles}/wes_simmut_bamsurgeon_indel.tsv -r ${simfiles}/wes_simmut.bed -S ${n} -o ${outdir}/bamsurgeon/${n}
\${python2} ${scriptdir}/extract_info.py ${outdir}/bamsurgeon/${n}/${n}_sim.mpileup.vcf.gz ${c} |awk '{print \$0\"\\tSimulation\"}' > ${outdir}/bamsurgeon/${n}/${n}_simmut.info.txt" >"${outdir}/scriptlog/${n}.bamsurgeon.sh"

done

# summary & plot
mkdir -p "${outdir}/plot"

echo "source ${init} WES
paste ${outdir}/realdata/*_simmut.info.txt |cut -f 1,2,6,10,14,18,22,26 > ${outdir}/realdata/realdata.txt

for tool in safemut varben bamsurgeon
do
paste ${outdir}/\${tool}/SRR*/*_simmut.info.txt |cut -f 1,2,6,10,14,18,22,26 > ${outdir}/\${tool}/simulation.txt
\${Rscript} ${scriptdir}/plot_wes.r ${simfiles}/wes_simmut_gs.txt ${outdir}/realdata ${outdir}/\${tool} ${outdir}/plot/\${tool}
done" > "${outdir}/scriptlog/Summary_plot.sh"
