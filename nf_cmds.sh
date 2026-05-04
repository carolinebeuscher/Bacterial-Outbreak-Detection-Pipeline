#!/bin/bash

##########################################
# Nextflow Commands #
# Caroline Beuscher # 
##########################################

### Download and convert SRA files ###

# Set up and activate conda environment
source ~/miniforge3/etc/profile.d/conda.sh

# conda create -n ex3 -c bioconda -c conda-forge entrez-direct sra-tools fastqc trimmomatic skesa spades pigz tree -y
conda activate ex3

mkdir -pv ./raw_fastq
cd ./raw_fastq


#use for loop to get all accessions from ncbi 
for accession in SRR2584863 SRR9094324 ERR019289 SRR1172848 SRR2093876; do
    prefetch "${accession}"
done

# convert to all sra files to fastq files using for loop 
for accession in SRR2584863 SRR9094324 ERR019289 SRR1172848 SRR2093876; do
  fasterq-dump \
   "${accession}" \
   --outdir . \
   --split-files \
   --skip-technical
done

# compress files: fastq -> fastq.gz
pigz -9 *.fastq

conda deactivate
cd ..


### Remove low quality reads ###

# Set up and activate conda environment
# conda create -n read_cleaning_env -c bioconda -c conda-forge fastp trimmomatic pigz tree -y
conda activate read_cleaning_env

mkdir ./clean_fastq
cd ./clean_fastq

# Clean with fastp
for read in ../raw_fastq/*_1.fastq.gz; do
  sample="$(basename ${read} _1.fastq.gz)"
   fastp \
   -i "${read}" \
   -I "${read%_1.fastq.gz}_2.fastq.gz" \
   -o "${sample}.R1.fq.gz" \
   -O "${sample}.R2.fq.gz" \
   --json "${sample}.json" \
   --html "${sample}.html"
done


conda deactivate 
cd ..

### Genome Assembly ###

# Set up and activate conda environment
conda activate ex3
mkdir ./Assemblies
cd ./Assemblies

# assembly with skesa
for read in ../clean_fastq/*.R1.fq.gz; do
  sample="$(basename ${read} .R1.fq.gz)"
  skesa \
   --reads "${read}","${read%R1.fq.gz}R2.fq.gz" \
   --cores 4 \
   --min_contig 1000 \
   --contigs_out "${sample}".fna
done


conda deactivate
cd ../clean_fastq

## get trimmed fastq metrics ##

#conda create -n seqkit_env -c bioconda seqkit -y
conda activate seqkit_env

#-a All Statistics (including N50, Q30)
#-b grabs basename

seqkit stats -a -b *.fq.gz > fastq_stats.tsv

conda deactivate
cd ..
