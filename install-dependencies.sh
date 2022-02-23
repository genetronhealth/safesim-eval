#!/usr/bin/env sh

### NOTICE: please set python2 and python3 path
# Python 2.7.15+ is required for VarBen installation
# Python 3.8.10+ is required for Bamsurgeon installation 
python2=python2
python3=python3

wd=`pwd -P`
extdir=$(dirname "$0")/ext
test -d ${exdir} && rm -rf ${extdir}
mkdir -p "${extdir}/bin" && cd "${extdir}"

### install samtools-1.11
wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
tar jxvf samtools-1.11.tar.bz2
cd samtools-1.11
./configure --prefix=`pwd -P` --disable-bz2 --disable-lzma # which will make some CRAM files produced elsewhere unreadable
make && make install
cp samtools ../bin
cp misc/wgsim ../bin
cd ..

### install bcftools-1.11
wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2
tar jxvf bcftools-1.11.tar.bz2
cd bcftools-1.11
./configure --prefix=`pwd -P` --disable-bz2 --disable-lzma # which will make some CRAM files produced elsewhere unreadable
make && make install
cp bcftools ../bin
cd ..

### install bwa-0.7.17
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar jxvf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
cp bwa ../bin
cd ..

### install safemut
# TODO

### install varben commit 0f66e35dc85b80938df2beafa5919330c9356953
# NOTICE: please check python2 configuration at the start of current script !
${python2} -m pip install pysam
${python2} -m pip install numpy
git clone https://github.com/nccl-jmli/VarBen.git

### install bamsurgeon
# Step1: install external prerequisites: picard, exonerate, Velvet, pysam
wget https://github.com/broadinstitute/picard/releases/download/2.26.10/picard.jar

wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz
tar zxvf velvet_1.2.10.tgz
make -C velvet_1.2.10
cp velvet_1.2.10/velvetg ./bin
cp velvet_1.2.10/velveth ./bin

git clone https://github.com/adamewing/exonerate.git
cd exonerate
git checkout v2.4.0
autoreconf -i
./configure && make && make check && make install
cd ..

# NOTICE: please check python3 configuration at the start of current script !
${python3} -m pip install cython
${python3} -m pip install pysam

# Step2: install bamsurgeon
wget https://github.com/adamewing/bamsurgeon/archive/refs/tags/1.3.tar.gz
mv 1.3.tar.gz bamsurgeon-1.3.tar.gz
tar zxvf bamsurgeon-1.3.tar.gz
cd bamsurgeon-1.3
${python3} setup.py build
${python3} setup.py install
cd ..

### install fgbio
wget https://github.com/fulcrumgenomics/fgbio/releases/download/1.5.0/fgbio-1.5.0.jar

cd ${wd}
