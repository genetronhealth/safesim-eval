mkdir -p genome
cd ./genome

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gzip -d hs37d5.fa.gz

bwa index hs37d5.fa
samtools faidx hs37d5.fa
${java} -jar ${picard} -R hs37d5.fa

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
gzip -d GCA_000001405.15_GRCh38_full_analysis_set.fna.gz

bwa index GCA_000001405.15_GRCh38_full_analysis_set.fna
samtools index GCA_000001405.15_GRCh38_full_analysis_set.fna
${java} -jar ${picard} -R GCA_000001405.15_GRCh38_full_analysis_set.fna
