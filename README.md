# Somatic_mutation
In this hands-on tutorial, we will call somatic short mutations, using the GATK4 Mutect2 and Varscan

This repository contains code for processing and analyzing raw sequencing data from NGS samples. The task is divided into two part:


### Part 1: Variant Calling 1

1. **Environment Setup**: 
   - A conda environment (`ngs_analysis`) is created and activated.
   - Tools like `fastqc`, `multiqc`, `trimmomatic`, `bwa`, `samtools`, `gatk4`, and `picard` are installed.

2. **Directory and File Setup**: 
   - Various directories are defined for raw data, results, reference genome, and supporting files.
   - Reference genome (`hg38.fa`) and known sites file (`dbsnp138.vcf`) are downloaded, uncompressed, and indexed.

3. **Preprocessing**: 
   - Raw fastq files for tumor and normal samples are defined.
   - Directories for aligned reads, raw reads, and results are created.

4. **Quality Control (QC)**:
   - FastQC is run on raw fastq files, followed by a MultiQC report for a combined overview of the quality.

5. **Trimming**:
   - Adapters are trimmed from the fastq files using `Trimmomatic`.

6. **Alignment**:
   - BWA is used for aligning the tumor and normal samples to the reference genome.

7. **Read Group Addition**:
   - `Picard` tools are used to add read group information to the aligned BAM files.

8. **Mark Duplicates**:
   - Duplicates are marked using `Picard`'s `MarkDuplicates` tool.

9. **Base Quality Recalibration**:
   - GATK `BaseRecalibrator` and `ApplyBQSR` are used to perform base quality score recalibration.

10. **Metrics Collection**:
    - GATK `CollectAlignmentSummaryMetrics` and `CollectInsertSizeMetrics` are used to collect various metrics for the aligned and deduplicated BAM files.


### Part 2: Variant Calling 2

1. **Paths Setup**: 
   - Paths to directories, data, and reference files are defined, similar to the first script.

2. **Download Mutect2 Supporting Files**:
   - Mutect2 supporting files such as the Panel of Normals (PoN), gnomAD data, and functional annotation files are downloaded.

3. **Mutect2 PoN Creation**:
   - A Panel of Normals (PoN) is created using Mutect2 in tumor-only mode on the normal BAM file.

4. **Variant Calling with Mutect2 and Varscan**:
   - Somatic variants are called using `Mutect2` in a tumor-only mode, followed by the creation of the PoN.
   - Further variant calling steps are defined with `Mutect2` and Varscan.

5. **Filtering, merging and annotation of the variants**:
   - Filter the Mutect2 VCF file to retain only those entries with the "PASS" filter using grep -E '^#|PASS'.
   - Similarly,filter the Varscan VCF file to retain entries with the "SOMATIC" tag using grep -E '^#|SOMATIC'.
   - The filtered Mutect2 and Varscan VCF files are then compressed using bgzip and indexed using tabix to make them compatible with downstream tools.
   - The compressed and indexed VCF files from Mutect2 and Varscan are merged using bcftools merge to create a combined VCF file.
   - Aannotate the merged VCF file using GATK's Funcotator tool to add functional annotations based on the provided reference genome (hg38) 
