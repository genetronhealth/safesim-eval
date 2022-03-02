#!/usr/bin/env sh

scriptdir=$(dirname $(readlink -e "$0"))/scripts
simfiles=$(dirname $(readlink -e "$0"))/simfiles
init=$(dirname $(readlink -e "$0"))/init_eval.sh

inputdir=$(dirname $(readlink -e "$0"))/Rawdata/ctDNA
outdir=$(dirname $(readlink -e "$0"))/Results/ctDNA

while getopts :i:o: opt; do
    case ${opt} in
    i)
        inputdir="${OPTARG}"
        ;;
    o)
        outdir="${OPTARG}"
        ;;
    ?)
        echo "\nUsage: ./process_ct.sh -i inputdir -o outdir"
        ;;
    esac
done

mkdir -p "${outdir}/scriptlog"
mkdir -p "${outdir}/realdata"
mkdir -p "${outdir}/rmdup-sim/consensus"
mkdir -p "${outdir}/sim-rmdup/nbam"

cat $(dirname "$0")/tn_sample.list |grep ^I |while read c t n
do
if [[ "${c}" == "IDT" ]];then 
    lumi=3
elif [[ "${c}" == "ILM" ]];then
    lumi=8
else
    echo "Errror: unknown umi length"; exit
fi

### realdata
tfq1="${inputdir}/${t}_1.fastq.gz"
tfq2="${inputdir}/${t}_2.fastq.gz"

echo "source ${init} ctDNA
bash ${scriptdir}/consensus.sh -i ${tfq1} -j ${tfq2} -l ${lumi} -s ${t} -o ${outdir}/realdata
bcftools mpileup -R ${simfiles}/${c}_simmut.bed --max-depth 99999 --max-idepth 99999 --min-BQ 0 -f \${REF} -a AD,ADF,ADR,DP,SP ${outdir}/realdata/${t}.consensus.align.bam |bcftools norm -f \${REF} -m- -Oz -o ${outdir}/realdata/${t}.mpileup.vcf.gz
bcftools index ${outdir}/realdata/${t}.mpileup.vcf.gz
\${python2} ${scriptdir}/extract_info.py ${outdir}/realdata/${t}.mpileup.vcf.gz ${c} |awk -v c=${c} '{print \$0\"\\t\"c\"\\tLbx_High\"}' > ${outdir}/realdata/${t}_simmut.info.txt" > "${outdir}/scriptlog/${t}.sh"

# rmdup-sim
mkdir -p "${outdir}/rmdup-sim"
nfq1="${inputdir}/${n}_1.fastq.gz"
nfq2="${inputdir}/${n}_2.fastq.gz"

echo "source ${init} ctDNA
bash ${scriptdir}/consensus.sh -i ${nfq1} -j ${nfq2} -l ${lumi} -s ${n} -o ${outdir}/rmdup-sim/consensus" > "${outdir}/scriptlog/${n}_consensus.sh"

echo "source ${init} ctDNA
bash ${scriptdir}/safemut_sim.sh -n ${outdir}/rmdup-sim/consensus/${n}.consensus.align.bam -v ${simfiles}/${c}_simmut_safemut.vcf.gz -r ${simfiles}/${c}_simmut.bed -s ${n} -o ${outdir}/rmdup-sim/safemut/${n}
\${python2} ${scriptdir}/extract_info.py ${outdir}/rmdup-sim/safemut/${n}/${n}_sim.mpileup.vcf.gz ${c} |awk -v c=${c} '{print \$0\"\\t\"c\"\\tSimulation\"}' > ${outdir}/rmdup-sim/safemut/${n}/${n}_simmut.info.txt" > "${outdir}/scriptlog/${n}_rmdup-sim.safemut.sh"

echo "source ${init} ctDNA
bash ${scriptdir}/varben_sim.sh -n ${outdir}/rmdup-sim/consensus/${n}.consensus.align.bam -v ${simfiles}/${c}_simmut_varben.tsv -r ${simfiles}/${c}_simmut.bed -s ${n} -o ${outdir}/rmdup-sim/varben/${n}
\${python2} ${scriptdir}/extract_info.py ${outdir}/rmdup-sim/varben/${n}/${n}_sim.mpileup.vcf.gz ${c} |awk -v c=${c} '{print \$0\"\\t\"c\"\\tSimulation\"}' > ${outdir}/rmdup-sim/varben/${n}/${n}_simmut.info.txt" > "${outdir}/scriptlog/${n}_rmdup-sim.varben.sh"

echo "source ${init} ctDNA
bash ${scriptdir}/bamsurgeon_sim.sh -n ${outdir}/rmdup-sim/consensus/${n}.consensus.align.bam -s ${simfiles}/${c}_simmut_bamsurgeon_snv.tsv -i ${simfiles}/${c}_simmut_bamsurgeon_indel.tsv -r ${simfiles}/${c}_simmut.bed -S ${n} -o ${outdir}/rmdup-sim/bamsurgeon/${n}
\${python2} ${scriptdir}/extract_info.py ${outdir}/rmdup-sim/bamsurgeon/${n}/${n}_sim.mpileup.vcf.gz ${c} |awk -v c=${c} '{print \$0\"\\t\"c\"\\tSimulation\"}' > ${outdir}/rmdup-sim/bamsurgeon/${n}/${n}_simmut.info.txt" > "${outdir}/scriptlog/${n}_rmdup-sim.bamsurgeon.sh"

# sim-rmdup
mkdir -p "${outdir}/sim-rmdup/nbam"

echo "source ${init} ctDNA
bwa mem -t 8 -R \"@RG\\tID:${n}\\tSM:${n}\\tLB:${n}\\tPU:${n}.L1\\tPL:ILLUMINA\" \${REF} ${nfq1} ${nfq2} |samtools sort -@ 4 -o ${outdir}/sim-rmdup/nbam/${n}.bam
samtools index -@ 4 ${outdir}/sim-rmdup/nbam/${n}.bam " > "${outdir}/scriptlog/${n}_bwa.sh"

echo "source ${init} ctDNA
bash ${scriptdir}/safemut_sim.sh -n ${outdir}/sim-rmdup/nbam/${n}.bam -v ${simfiles}/${c}_simmut_safemut.vcf.gz -r ${simfiles}/${c}_simmut.bed -s ${n} -o ${outdir}/sim-rmdup/safemut/${n} -l ${lumi}
\${python2} ${scriptdir}/extract_info.py ${outdir}/sim-rmdup/safemut/${n}/${n}_sim.mpileup.vcf.gz ${c} |awk -v c=${c} '{print \$0\"\\t\"c\"\\tSimulation\"}' > ${outdir}/sim-rmdup/safemut/${n}/${n}_simmut.info.txt" > "${outdir}/scriptlog/${n}_sim-rmdup.safemut.sh"

echo "source ${init} ctDNA
bash ${scriptdir}/varben_sim.sh -n ${outdir}/sim-rmdup/nbam/${n}.bam -v ${simfiles}/${c}_simmut_varben.tsv -r ${simfiles}/${c}_simmut.bed -s ${n} -o ${outdir}/sim-rmdup/varben/${n} -l ${lumi}
\${python2} ${scriptdir}/extract_info.py ${outdir}/sim-rmdup/varben/${n}/${n}_sim.mpileup.vcf.gz ${c} |awk -v c=${c} '{print \$0\"\\t\"c\"\\tSimulation\"}' > ${outdir}/sim-rmdup/varben/${n}/${n}_simmut.info.txt" > "${outdir}/scriptlog/${n}_sim-rmdup.varben.sh"

done

# summary & plot
mkdir -p "${outdir}/plot"

echo "source ${init} ctDNA

paste ${outdir}/realdata/*_simmut.info.txt |cut -f 1,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62,66,70,74,78,82,86,90,94 > ${outdir}/realdata/realdata.txt

for tool in safemut varben bamsurgeon
do
paste ${outdir}/rmdup-sim/\${tool}/SRR*/*_simmut.info.txt |cut -f 1,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62,66,70,74,78,82,86,90,94 > ${outdir}/rmdup-sim/\${tool}/simulation.txt
\${Rscript} ${scriptdir}/plot_ct.r ${simfiles}/ctDNA_simmut_gs.txt ${outdir}/realdata ${outdir}/rmdup-sim/\${tool} ${outdir}/plot/\${tool}_rmdup-sim
done

for tool in safemut varben
do
paste ${outdir}/sim-rmdup/\${tool}/SRR*/*_simmut.info.txt |cut -f 1,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62,66,70,74,78,82,86,90,94 > ${outdir}/sim-rmdup/\${tool}/simulation.txt
\${Rscript} ${scriptdir}/plot_ct.r ${simfiles}/ctDNA_simmut_gs.txt ${outdir}/realdata ${outdir}/sim-rmdup/\${tool} ${outdir}/plot/\${tool}_sim-rmdup
done" > "${outdir}/scriptlog/Summary_plot.sh"

