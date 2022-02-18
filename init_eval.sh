data=$1   # ctDNA / WES

declare -r rootdir=$(dirname $(readlink -e ${BASH_SOURCE[0]}))
extdir=${rootdir}/ext
ctDNA_ref=${rootdir}/genome/hs37d5.fa
wes_ref=${rootdir}/genome/GCA_000001405.15_GRCh38_full_analysis_set.fna

# global dependencies
export TMP=/sxyf_PR/user/guojingyu/tmp
export python2=python2
export python3=/bionfsdate/ctDNA/experiment/guojingyu/software/miniconda2/envs/facets/bin/python
export Rscript=Rscript
export PATH=${extdir}/bin:$PATH   # samtools, bcftools, bwa, wgsim

if [[ $data == "WES" ]];then

echo "Initiate Evaluation on WES Data!"
export REF=${wes_ref}
elif [[ $data == "ctDNA" ]];then
echo "Initiate Evaluation on ctDNA Data!"
export REF=${ctDNA_ref}
# run consensus dependencies
export java=java
export fgbio=${extdir}/fgbio-1.5.0.jar
else
echo "Please check input datatype!"
exit 0

fi

# safemut dependencies
export PATH=${extdir}/safesim-2022-0210:$PATH

# varben dependencies
export muteditor=${extdir}/VarBen/bin/muteditor.py

# bamsurgeon dependencies
export addsnv=${extdir}/bamsurgeon-1.3/bin/addsnv.py
export addindel=${extdir}/bamsurgeon-1.3/bin/addindel.py
export picard=${extdir}/picard.jar

