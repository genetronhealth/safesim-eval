This README describes how to evaluate SafeMut with BamSurgeon and VarBen. 

# Install Dependencies
Java8, Python 2.7.15+, Python 3.8.10+ are presumed be installed, you may use conda to do python version management.
To install from scratch, please set python path at start of ./install-dependencies.sh and run: `./install-dependencies.sh`.
Bamsurgeon may be failed during exonerate installation. you may install exonerate using conda: `conda install -c bioconda exonerate`
R 4.0+ is required for plotting. The R packages `ggplot2`, `gridExtra` and `this.path` need to be installed. These packages can be installed from CRAN.

1. samtools, bcftools, bwa.

If these dependencies were already installed, please run following commands:

```
cp /path/to/samtools ./ext/bin
cp /path/to/bcftools ./ext/bin
cp /path/to/bwa      ./ext/bin
cp /path/to/wgism    ./ext/bin
```

2. safemut, VarBen, Bamsurgeon.

For manual installation, see following instructions.

Safemut ...
VarBen requires python2 and external prerequisites. Please see https://github.com/nccl-jmli/VarBen for installation.
Bamsurgeon requires python3 and external prerequisites. Please see https://github.com/adamewing/bamsurgeon for installation.
If these dependencies were installed manually, please assign to corresponding environment variables in ./init_eval.sh

3. fgbio

fgbio is used for calling consensus by UMIs in ctDNA data processing.
If fgbio was already installed, please run `cp /path/to/fgbio.jar ./ext`.

4. picard

picard is requiered for bamsurgeon and calling consensus.
If picard was already installed, please run `cp /path/to/picard.jar ./ext`.

---------------
# Download Data

1. Reference genome

please run `./reference.sh` to download and build index for hs37d5.fa(ctDNA_ref) and GCA_000001405.15_GRCh38_full_analysis_set.fna(wes_ref).
If reference genomes were already prepared, please configure `ctDNA_ref` and `wes_ref` environment variables in ./init_eval.sh 

2. WES dataset and ctDNA dataset

SRA Accession numbers of WES samples are listed in ./tn_sample.list.
Data will be dowoloaded in ./RawData default.

./RawData structure (the fastq.gz files should be pre-downloaded using your favorite tool such as ascp, wget, or curl):
```
./
|--WES
|  |--SRRXXXXX_1.fastq.gz
|  |--SRRXXXXX_2.fastq.gz
|--ctDNA
|  |--SRRXXXXX_1.fastq.gz
|  |--SRRXXXXX_2.fastq.gz
```
------------------------
# Run Evaluation Scripts

**NOTICE:**
Before evaluation, please check configuration in `./init_eval.sh` and make sure dependencies and dataset are configured correctly.

1. Evaluation on WES dataset

Evaluation on WES dataset using safemut, VarBen and Bamsurgeon.

Please run `bash ./eval_wes_main.sh` to generate scripts.
Default input directory is ./Rawdata/WES, and default output directory is Results/WES.
You can use customized path by set options `-i` and `-o`. Please note the input directory must have the same structure as the above example.

Run scripts:
Run the evaluation and plot of all samples: please run `ls ./Results/WES/scriptlog/*.sh |awk '{print "bash " $1}' |sh > ./Results/WES/scriptlog/runall.logo 2>runall.loge`.
 
For one sample (tumor-normal pair): please run `bash ./Results/WES/scriptlog/tumor.sh ; ls ./Results/WES/scriptlog/normal*.sh |awk '{print "bash " $1}' |sh >./Results/WES/scriptlog/tumor-normal.logo 2> ./Results/WES/scriptlog/tumor-normal.loge`. tumor is SRA accession number of tumor sample and normal is SRA accession number of matched normal sample. (this is only for code testing.)


Optionally, you can also use `parallel` to process data in multithreads and you can just process one sample each time in step 1,2,3.
```
# Step1: process tumor
ls ./Results/WES/scriptlog/SRR*[0-9].sh |awk '{print "bash " $1 " >" $1 ".logo 2>" $1 ".loge"}' |parallel -j $job_number

# Step2: prepare normal bam
ls ./Results/WES/scriptlog/SRR*_bwa.sh |awk '{print "bash " $1 " >" $1 ".logo 2>" $1 ".loge"}' |parallel -j $job_number

# Step3: mutation simulation
ls ./Results/WES/scriptlog/SRR*.*.sh |awk '{print "bash " $1 " >" $1 ".logo 2>" $1 ".loge"}' |parallel -j $job_number

# Step4: summary and plot
bash ./Results/WES/scriptlog/Summary_plot.sh
```

2. Evaluation on ctDNA dataset

Evaluation on ctDNA dataset consists of two methods: 
1.simulation before calling consensus (safemut, VarBen)
2.simulation after calling consensus  (safemut, VarBen, Bamsurgeon)

Please run `bash ./eval_ct_main.sh` to generate scripts.
Default input directory is ./Rawdata/ctDNA, and default output directory is Results/ctDNA.
Using customized path by set options `-i` and `-o`. Please note the input directory must have the same structure as example above.

Run scripts: 
Run all samples evaluation and plot: please run `ls ./Results/ctDNA/scriptlog/*.sh |awk '{print "bash " $1}' |sh > ./Results/ctDNA/scriptlog/runall.logo 2>./Results/ctDNA/scriptlog/runall.loge`.

For one sample (tumor-normal pair): please run `bash ./Results/ctDNA/scriptlog/tumor.sh ; ls ./Results/ctDNA/scriptlog/normal*.sh |awk '{print "bash " $1}' |sh > ./Result/ctDNA/scriptlog/tumor-normal.logo 2> ./Result/ctDNA/scriptlog/tumor-normal.loge`. tumor is SRA accession number of tumor sample and normal is SRA accession number of matched normal sample (this is only for code testing).
 
Optionaly, you can also use `parallel` to process data in multithreads and you can just process one sample each time in step 1,2,3.
```
# Step1: process tumor
ls ./Results/ctDNA/scriptlog/SRR*[0-9].sh |awk '{print "bash " $1 " >" $1 ".logo 2>" $1 ".loge"}' |parallel -j $job_number

# Step2: prepare consensus bam
ls ./Results/ctDNA/scriptlog/SRR*_consensus.sh |awk '{print "bash " $1 " >" $1 ".logo 2>" $1 ".loge"}' |parallel -j $job_number

# Step3: simulation on consensus bam
ls ./Results/ctDNA/scriptlog/SRR*_rmdup-sim.*.sh |awk '{print "bash " $1 " >" $1 ".logo 2>" $1 ".loge"}' |parallel -j $job_number

# Step4: prepare normal bam
ls ./Results/ctDNA/scriptlog/SRR*_bwa.sh |awk '{print "bash " $1 " >" $1 ".logo 2>" $1 ".loge"}' |parallel -j $job_number

# Step5: simulation on original normal bam
ls ./Results/ctDNA/scriptlog/SRR*_sim-rmdup.*.sh |awk '{print "bash " $1 " >" $1 ".logo 2>" $1 ".loge"}' |parallel -j $job_number

# Step6: summary and plot
bash ./Results/ctDNA/scriptlog/Summary_plot.sh
```

3. Output directory contents
3.1 Output directory structure

```
./Results
|--WES
|  |--realdata   # call variants in tumor data
|  |--nbam       # normal bam for mutaion simulation
|  |--safemut
|  |--varben
|  |--bamsurgeon
|  |--plot
|  |--scriptlog
|--ctDNA
|  |--realdata
|  |--rmdup-sim
|  |  |--consensus  # normal consensus for mutation simulation
|  |  |--safemut
|  |  |--varben
|  |  |--bamsurgeon
|  |--sim-rmdup
|  |  |--safemut
|  |  |--varben
|  |--plot
|  |--scriptlog
```

3.2 key results of each mutation-simulation software

safemut:    sample_sim_1.sort.fastq.gz sample_sim_2.sort.fastq.gz sample_sim.mpileup.vcf.gz
VarBen:     edited.sort.bam sample_sim.mpileup.vcf.gz
Bamsurgeon: sample_simsnv.sort.bam sample_simindel.sort.bam sample_sim.mpileup.vcf.gz

3.3 evaluation plot

The plot directory includes two classes of plots. 
Class 1: XXX.mut.pdf shows mutation Allele Fraction Or/And Depth of mutations both realdata and simulation data.
Class 2: XXX.qqplot.pdf, XXX.mean.pdf, XXX.var.plot shows comparison of statistic paramaters (Z-score, Mean, Variance) between realdata and simulation data.
