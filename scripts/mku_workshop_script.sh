# MKU workshop Script
#!/bin/bash

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh

conda update --yes conda

# Install some conda channels
# A channel is where conda looks for packages
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create an environment
conda create -n ngs python=3

# activate the environment
conda activate ngs

#------------------------------------------------------------------------------
# create a directory you work in
mkdir analysis

# change into the directory
cd analysis

#create a directory for the data
mkdir data

# change into the data directory
cd data

#Copy the sequencing data files (compressed FASTQs) into the data directory
# cp data/*.fastq.gz .
cp -r /home/wanjau/Downloads/MKU_bioinfomatics_workshop_2022-20220412T083435Z-001/MKU_bioinfomatics_workshop_2022/data/FASTQ/FASTQ_files/*.fastq.gz . .

# change into the directory
#cd analysis
cd ..
# create env and install tools
conda create --yes -n qc fastp fastqc multiqc

# activate env
conda activate qc

# create directory for the placing the trimmed output
mkdir trimmed

# trimming the data
fastp --detect_adapter_for_pe \
      --overrepresentation_analysis \
      --correction \
      --cut_right \
      --thread 2 \
      --html trimmed/GSM1275862.fastp.html --json trimmed/GSM1275862.fastp.json \
      -i data/GSM1275862_R1.fastq.gz -I data/GSM1275862_R2.fastq.gz \
      -o trimmed/GSM1275862_R1.fastq.gz -O trimmed/GSM1275862_R2.fastq.gz

# create fastqc_report directory
mkdir fastqc_report

# check data quality
fastqc -o fastqc_report data/GSM1275862_R1.fastq.gz

# check quality of all samples
fastqc -o fastqc_report data/*.fastq.gz

# check multiqc report
mkdir multiqc_report
multiqc -o multiqc_report fastqc_report

# Create a mapping result directory
mkdir mappings

# Create the conda environment and Install HISAT2
conda create --yes -n mapping samtools hisat2 qualimap r-base
conda activate mapping

# Align the reads to the reference chromosome 2 genome
hisat2 -q -x genome_index/hg38_chr2_genome \
       -1 data/trimmed/GSM1275862_R1.fastq.gz \
       -2 data/trimmed/GSM1275862_R2.fastq.gz \
       -S mappings/GSM1275862.sam

# Sort the alignment
# first install samtools using conda
# conda install -y samtools
samtools sort -n -O sam mappings/GSM1275862.sam | samtools fixmate -m -O bam \
- mappings/GSM1275862.fixmate.bam

# remove the sam file
rm mappings/GSM1275862.sam

# Sorting
# Convert to bam file and sort
samtools sort -O bam -o mappings/GSM1275862.sorted.bam mappings/GSM1275862.fixmate.bam

# Once it successfully finished, delete the fixmate file to save space
rm mappings/GSM1275862.fixmate.bam

# Remove duplicates
samtools markdup -r -S mappings/GSM1275862.sorted.bam mappings/GSM1275862.sorted.dedup.bam

# If it worked, delete the original file
rm mappings/GSM1275862.sorted.bam

# Mapping statistics using SAMtools
samtools flagstat mappings/GSM1275862.sorted.dedup.bam

# For indepth statistics run
samtools index mappings/GSM1275862.sorted.dedup.bam
samtools idxstats mappings/GSM1275862.sorted.dedup.bam

samtools depth mappings/GSM1275862.sorted.dedup.bam | gzip > mappings/GSM1275862.depth.txt.gz

# Statistics with QualiMap
# Run QualiMap
qualimap bamqc -bam mappings/GSM1275862.sorted.dedup.bam

# Once finished open result page using your web browser.
# You may also install firefox in your conda env and use it to visualize the results
firefox mappings/GSM1275862.sorted.dedup_stats/qualimapReport.html

# Sub-selecting reads
# Selecting concordantly mapped reads
samtools view -h -b -f 3 mappings/GSM1275862.sorted.dedup.bam > mappings/GSM1275862.sorted.dedup.concordant.bam

#-----------------------------------------------------------------------------------------------------------------
# Unmapped reads
#Get unmapped reads
samtools view -b -f 4 mappings/GSM1275862.sorted.dedup.bam > mappings/GSM1275862.sorted.unmapped.bam

# Delete the original to save space, however, in reality you might want to save it to investigate later
rm mappings/GSM1275862.sorted.dedup.bam

# count the unmapped reads
samtools view -c mappings/GSM1275862.sorted.unmapped.bam

# Extract reads in R1 fastq and R2 fastq files
samtools fastq -1 mappings/GSM1275862.sorted.unmapped.R1.fastq.gz \
               -2 mappings/GSM1275862.sorted.unmapped.R2.fastq.gz \
               mappings/GSM1275862.sorted.unmapped.bam

# delete files not needed
rm mappings/GSM1275862.sorted.unmapped.bam

#----------------------------------------------------------------------------------------------------

# Expression quantification and normalization

cd analysis
conda create --yes -n htseq_count htseq samtools
conda activate htseq_count
mkdir expression

# index the concordant file
samtools index mappings/GSM1275862.sorted.dedup.concordant.bam

htseq-count \
     -f bam \
     -s no \
     -t exon \
     -i gene_id \
     --additional-attr=gene_name \
     mappings/GSM1275862.sorted.dedup.concordant.bam \
     genomes_and_annotations/hg38_chr2_annotations.gtf > expression/GSM1275862_gene_counts.txt

#-----------------------------------------------------------------------------------------------------
