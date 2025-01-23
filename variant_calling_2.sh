#!/bin/bash

# Define paths to directories and files
project_dir=/media/my_drive/Bioinfo
data_dir=$project_dir/rawdata
ref=$project_dir/Genome/hg38.fa
known_sites=$project_dir/supporting_files/Homo_sapiens_assembly38.dbsnp138.vcf
aligned_reads=$project_dir/aligned_reads
reads=$project_dir/reads
results=$project_dir/results
mutect2_supporting_files=$project_dir/supporting_files/mutect2_supporting_files

# Create necessary directories if they do not exist
mkdir -p "${results}" "${mutect2_supporting_files}" "${project_dir}/tmp"

###################################################
# Download Mutect2 Supporting Files (only once)
###################################################

echo "Downloading Mutect2 supporting files..."
# gnomAD
wget -P ${mutect2_supporting_files} https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz
wget -P ${mutect2_supporting_files} https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi

# PoN
wget -P ${mutect2_supporting_files} https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
wget -P ${mutect2_supporting_files} https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi

# Intervals
wget -P ${mutect2_supporting_files} https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list
wget -P ${mutect2_supporting_files} https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list

# Funcotator Data Sources
wget -P ${mutect2_supporting_files} https://storage.googleapis.com/broad-public-datasets/funcotator/funcotator_dataSources.v1.8.hg19.20230908g.tar.gz
tar -xvzf ${mutect2_supporting_files}/funcotator_dataSources.v1.8.hg19.20230908g.tar.gz -C ${mutect2_supporting_files}

gatk FuncotatorDataSourceDownloader --germline --validate-integrity --hg19 --extract-after-download
gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --hg38 --extract-after-download

###################################################
# Mutect2 Pipeline: Create Panel of Normals (PoN)
###################################################

# Run Mutect2 in tumor-only mode on normal BAM
echo "Creating Panel of Normals (PoN)..."
gatk Mutect2 \
    -R ${ref} \
    -I ${aligned_reads}/normal_sorted_dedup_bqsr.bam \
    -tumor TumorSample \
    -L ${mutect2_supporting_files}/wgs_calling_regions.hg38.interval_list \
    -O ${mutect2_supporting_files}/tumor_mutect2.vcf.gz

gatk CreateSomaticPanelOfNormals \
   --java-options "-Xmx50G" \
   -V ${mutect2_supporting_files}/tumor_mutect2.vcf.gz \
   -O ${mutect2_supporting_files}/PON.vcf.gz

###################################################
# Mutect2 Variant Calling
###################################################

# Step 8: Call Variants with Mutect2
echo "Calling somatic variants with Mutect2..."
gatk Mutect2 \
    -R ${ref} \
    -I ${aligned_reads}/tumor_sorted_dedup_bqsr.bam \
    -I ${aligned_reads}/normal_sorted_dedup_bqsr.bam \
    -tumor TumorSample \
    -normal NormalSample \
    --germline-resource ${mutect2_supporting_files}/af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals ${mutect2_supporting_files}/PON.vcf.gz \
    -O ${results}/somatic_variants_mutect2.vcf.gz \
    --f1r2-tar-gz ${results}/f1r2.tar.gz

# Step 9: Estimate Contamination
echo "Estimating cross-sample contamination..."
gatk GetPileupSummaries \
    -I ${aligned_reads}/tumor_sorted_dedup_bqsr.bam \
    -V ${mutect2_supporting_files}/af-only-gnomad.hg38.vcf.gz \
    -L ${mutect2_supporting_files}/wgs_calling_regions.hg38.interval_list \
    -O ${results}/tumor_getpileupsummaries.table

gatk GetPileupSummaries \
    -I ${aligned_reads}/normal_sorted_dedup_bqsr.bam \
    -V ${mutect2_supporting_files}/af-only-gnomad.hg38.vcf.gz \
    -L ${mutect2_supporting_files}/wgs_calling_regions.hg38.interval_list \
    -O ${results}/normal_getpileupsummaries.table

gatk CalculateContamination \
    -I ${results}/tumor_getpileupsummaries.table \
    --matched ${results}/normal_getpileupsummaries.table \
    -O ${results}/pair_calculatecontamination.table

# Step 10: Read Orientation Artifacts
echo "Estimating read orientation artifacts..."
gatk LearnReadOrientationModel \
    -I ${results}/f1r2.tar.gz \
    -O ${results}/read-orientation-model.tar.gz

# Step 11: Filter Variants
echo "Filtering Mutect2 calls..."
gatk FilterMutectCalls \
    -V ${results}/somatic_variants_mutect2.vcf.gz \
    -R ${ref} \
    --contamination-table ${results}/pair_calculatecontamination.table \
    --ob-priors ${results}/read-orientation-model.tar.gz \
    -O ${results}/somatic_variants_filtered_mutect2.vcf


###################################################
# Varscan SNP and Indel Calling
###################################################

# Create mpileup file
echo "Creating mpileup file..."
samtools mpileup -f ${ref} \
    ${aligned_reads}/tumor_sorted_dedup_bqsr.bam \
    ${aligned_reads}/normal_sorted_dedup_bqsr.bam \
    -o ${results}/tumor_normal.mpileup

# Run Varscan
echo "Running Varscan..."
java -jar /home/pari/miniconda3/envs/myenv/share/varscan-2.4.6-0/VarScan.jar somatic \
    ${results}/tumor_normal.mpileup \
    ${results}/varscan_output \
    --strand-filter 0 --min-coverage 8 --min-var-freq 0.1 --p-value 0.99 --mpileup 1 --output-vcf 1

###################################################
# Filter, Compress, Index, and Merge VCF Files
###################################################

# Filter Mutect2 and Varscan VCF Files
echo "Filtering Mutect2 and Varscan VCF files..."
grep -E '^#|PASS' ${results}/somatic_variants_filtered_mutect2.vcf > ${results}/mutect_filtered_somatic.vcf
grep -E '^#|SOMATIC' ${results}/varscan_output.snp.vcf > ${results}/varscan_filtered_somatic.vcf

# Compress and index VCF files
echo "Compressing and indexing VCF files..."
bgzip -c ${results}/mutect_filtered_somatic.vcf > ${results}/mutect_filtered_somatic.vcf.gz
tabix -p vcf ${results}/mutect_filtered_somatic.vcf.gz

bgzip -c ${results}/varscan_filtered_somatic.vcf > ${results}/varscan_filtered_somatic.vcf.gz
tabix -p vcf ${results}/varscan_filtered_somatic.vcf.gz

# Merge VCF files
echo "Merging VCF files..."
bcftools merge \
    ${results}/mutect_filtered_somatic.vcf.gz \
    ${results}/varscan_filtered_somatic.vcf.gz \
    -o ${results}/merged_variants.vcf

echo "Pipeline execution completed successfully."


###################################################
#  Annotation of the final filtered and merged vcf
###################################################

# Step 12: Annotate Variants with Funcotator
echo "Annotating variants with Funcotator..."
gatk Funcotator \
    --variant ${results}/merged_variants.vcf \
    --reference ${ref} \
    --ref-version hg38 \
    --data-sources-path ${mutect2_supporting_files}/funcotator_dataSources.v1.8.hg38.20230908s \
    --output ${results}/somatic_variants_funcotated.vcf \
    --output-file-format VCF
