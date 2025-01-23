#!/bin/bash

################################################### Prep files & Environment ###########################################################
# Set up the conda environment and install the necessary tools
echo "Setting up the environment and preparing files..."

conda create -n ngs_analysis -y python=3.9
conda activate ngs_analysis
conda install -y -c bioconda fastqc multiqc trimmomatic bwa samtools gatk4 picard
conda install -y -c conda-forge pandas

# Define directories and file paths
project_dir=/media/my_drive/Bioinfo/
data_dir=$project_dir/rawdata/
ref=$project_dir/Genome/hg38.fa
known_sites=$project_dir/supporting_files/Homo_sapiens_assembly38.dbsnp138.vcf
aligned_reads=$project_dir/aligned_reads
reads=$project_dir/reads
results=$project_dir/results
mutect2_supporting_files=$project_dir/supporting_files/mutect2_supporting_files


# Download reference files
wget -P /media/my_drive/Bioinfo/Genome/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip /media/my_drive/Bioinfo/Genome/hg38.fa.gz

# Create reference index and dictionary
samtools faidx "$ref"
samtools dict "$ref" > "${ref%.fa}.dict"

# Download known sites files for BQSR from GATK resource bundle
wget -P /media/my_drive/Bioinfo/supporting_files/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P /media/my_drive/Bioinfo/supporting_files/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# Input raw data files
tumor_r1="$data_dir/PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq.gz"
tumor_r2="$data_dir/PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz"
normal_r1="$data_dir/PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz"
normal_r2="$data_dir/PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz"

# Create necessary directories
mkdir -p "$aligned_reads" "$reads" "$results"

###################################################### VARIANT CALLING STEPS ###########################################################

# -------------------
# STEP 1: QC - Run fastqc 
# -------------------

echo "STEP 1: QC - Running FastQC"

fastqc "$tumor_r1" -o "$reads/"
fastqc "$tumor_r2" -o "$reads/"
fastqc "$normal_r1" -o "$reads/"
fastqc "$normal_r2" -o "$reads/"

# Generate a combined QC report using MultiQC
echo "Generating MultiQC report..."
multiqc "$reads" -o "$results/"


# -----------------------------
# STEP 2: Trimming adapters 
# -----------------------------

#echo "STEP 2: Trimming adapters with Trimmomatic"

trimmed_tumor_r1="$reads/PA220KH-Tumor_trimmed_R1.fastq.gz"
trimmed_tumor_r2="$reads/PA220KH-Tumor_trimmed_R2.fastq.gz"
trimmed_normal_r1="$reads/PA221MH-Norm_trimmed_R1.fastq.gz"
trimmed_normal_r2="$reads/PA221MH-Norm_trimmed_R2.fastq.gz"

trimmomatic PE \
 -threads 4 \
 -phred33 \
 "$tumor_r1" "$tumor_r2" \
 "$trimmed_tumor_r1" "$reads/unpaired_tumor_r1.fastq.gz" \
 "$trimmed_tumor_r2" "$reads/unpaired_tumor_r2.fastq.gz" \
 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

trimmomatic PE \
 -threads 4 \
 -phred33 \
 "$normal_r1" "$normal_r2" \
 "$trimmed_normal_r1" "$reads/unpaired_normal_r1.fastq.gz" \
 "$trimmed_normal_r2" "$reads/unpaired_normal_r2.fastq.gz" \
 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# -----------------------------
# STEP 3: Alignment using BWA-MEM
# -----------------------------

echo "STEP 3: Aligning reads to reference using BWA-MEM"

# Index reference genome
bwa index "$ref"

# Align tumor samples
bwa mem "$ref" "$trimmed_tumor_r1" "$trimmed_tumor_r2" | samtools view -Sb - | samtools sort -o "$aligned_reads/tumor_sorted.bam"

# Align normal samples
bwa mem "$ref" "$trimmed_normal_r1" "$trimmed_normal_r2" | samtools view -Sb - | samtools sort -o "$aligned_reads/normal_sorted.bam"



# -----------------------------
# STEP 4: Add read group info
# -----------------------------

echo "STEP 4: Adding read group information"

picard AddOrReplaceReadGroups \
  I="$aligned_reads/normal_sorted.bam" \
  O="$aligned_reads/normal_sorted_rg.bam" \
  RGID=NormalSample RGLB=library1 RGPL=ILLUMINA RGPU=unit1 RGSM=NormalSample

picard AddOrReplaceReadGroups \
  I="$aligned_reads/tumor_sorted.bam" \
  O="$aligned_reads/tumor_sorted_rg.bam" \
  RGID=TumorSample RGLB=library1 RGPL=ILLUMINA RGPU=unit1 RGSM=TumorSample


# -----------------------------------------
# STEP 5: Mark Duplicates - Picard
# -----------------------------------------
echo "STEP 5: Mark Duplicates - Picard"

picard MarkDuplicates I="$aligned_reads/normal_sorted_rg.bam" O="$aligned_reads/normal_sorted_dedup.bam" M="$aligned_reads/normal_dedup_metrics.txt"
picard MarkDuplicates I="$aligned_reads/tumor_sorted_rg.bam" O="$aligned_reads/tumor_sorted_dedup.bam" M="$aligned_reads/tumor_dedup_metrics.txt"

# -----------------------------------------
# STEP 6: Base Quality Recalibration
# -----------------------------------------
echo "STEP 6: Base Quality Recalibration"

# Build the model
gatk BaseRecalibrator -I "$aligned_reads/normal_sorted_dedup.bam" -R "$ref" --known-sites "$known_sites" -O "$aligned_reads/normal_recal_data.table"
gatk BaseRecalibrator -I "$aligned_reads/tumor_sorted_dedup.bam" -R "$ref" --known-sites "$known_sites" -O "$aligned_reads/tumor_recal_data.table"

# Apply the model
gatk ApplyBQSR -I "$aligned_reads/normal_sorted_dedup.bam" -R "$ref" --bqsr-recal-file "$aligned_reads/normal_recal_data.table" -O "$aligned_reads/normal_sorted_dedup_bqsr.bam"
gatk ApplyBQSR -I "$aligned_reads/tumor_sorted_dedup.bam" -R "$ref" --bqsr-recal-file "$aligned_reads/tumor_recal_data.table" -O "$aligned_reads/tumor_sorted_dedup_bqsr.bam"


# -----------------------------------------------
# STEP 7: Collect Alignment & Insert Size Metrics
# -----------------------------------------------
echo "STEP 7: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics R="$ref" I="$aligned_reads/normal_sorted_dedup_bqsr.bam" O="$aligned_reads/normal_alignment_metrics.txt"
gatk CollectInsertSizeMetrics INPUT="$aligned_reads/normal_sorted_dedup_bqsr.bam" OUTPUT="$aligned_reads/normal_insert_size_metrics.txt" HISTOGRAM_FILE="$aligned_reads/normal_insert_size_histogram.pdf"



gatk CollectAlignmentSummaryMetrics R="$ref" I="$aligned_reads/tumor_sorted_dedup_bqsr.bam" O="$aligned_reads/tumor_alignment_metrics.txt"
gatk CollectInsertSizeMetrics INPUT="$aligned_reads/tumor_sorted_dedup_bqsr.bam" OUTPUT="$aligned_reads/tumor_insert_size_metrics.txt" HISTOGRAM_FILE="$aligned_reads/tumor_insert_size_histogram.pdf"
