data=$1   # ctDNA / WES

export rootdir=$(dirname $(readlink -e ${BASH_SOURCE[0]}))
extdir=${rootdir}/ext
ctDNA_ref=${rootdir}/genome/hs37d5.fa
wes_ref=${rootdir}/genome/GCA_000001405.15_GRCh38_full_analysis_set.fna

# global dependencies
export TMP=${TMPDIR}
mkdir -p ${TMP} && chmod uga+w ${TMP}/
export python2=python2
export python3=python3
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
export debarcode=${extdir}/uvc-0.8.0/bin/debarcode
elif [[ $data == "reference" ]];then
export java=java
echo "Prepare reference genome"
else
echo "Please check input datatype!"
fi

# safemut dependencies
export PATH=${extdir}/safesim-2022-0210:$PATH

# varben dependencies
export muteditor=${extdir}/VarBen/bin/muteditor.py

# bamsurgeon dependencies
export addsnv=${extdir}/bamsurgeon-1.3/bin/addsnv.py
export addindel=${extdir}/bamsurgeon-1.3/bin/addindel.py
export picard=${extdir}/picard.jar

# bedtools and safemut
export bedtools=${extdir}/bedtools2/bin/bedtools

