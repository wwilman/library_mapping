# SGRNA Mapping and Gene Expression Analysis Pipeline

This pipeline performs the mapping of sgRNA sequences to the GRCh38 reference genome and analyzes gene expression based on the mapping results. The pipeline is implemented using Nextflow and consists of several steps outlined below.

## Step 1: Download Reference Genome

The reference genome (GRCh38) is downloaded using the wget command and extracted using gunzip. The downloaded reference genome file is saved as GCA_000001405.15_GRCh38_full_analysis_set.fna.

## Step 2: Build Index with Bowtie2

The reference genome is used to build an index for mapping using Bowtie2. The index files are generated with the prefix human_ref and saved as human_ref.1.bt2, human_ref.2.bt2, human_ref.3.bt2, human_ref.4.bt2, human_ref.rev.1.bt2, and human_ref.rev.2.bt2.

## Step 3: Map sgRNA Library to Reference Genome

The sgRNA library sequences (library.fa) are mapped to the reference genome using Bowtie2. The resulting alignments are saved in the SAM format as alignment.sam.

## Step 4: Sort Aligned Data with Samtools

The aligned data in SAM format is sorted using Samtools to produce a sorted BAM file (alignment.bam).

## Step 5: Extract Gene Information

From the sorted BAM file, the pipeline extracts chromosome names, start and end positions, and strands for each mapped sgRNA sequence. The extracted information is saved in a tab-separated file called gene_info.txt.

## Step 6: Extract Gene Names from Fasta Descriptions

The gene names provided within the description lines of the fasta file (library.fa) are extracted and saved in a file called gene_names.txt.

## Step 7: Map Gene Names Based on Mapping Results

The gene names are mapped based on the mapping results obtained from the aligned data. For each mapped sgRNA, the pipeline checks if the gene name is present in the provided gene names file. The mapping results along with the assigned gene names (if found) are saved in a file called mapped_genes.txt.

## Step 8: Check Differences in Annotation

This step compares the gene names obtained from the fasta file descriptions with the mapped gene names. It identifies any differences in annotation between the two and outputs the results to the console.

## Step 9: Retrieve Gene Expression Matrix for Sample TCGA-A7-A13D-01A-13R-A12P-07

The pipeline retrieves the gene expression matrix for the gene IDs obtained from the previous step for the sample "TCGA-A7-A13D-01A-13R-A12P-07" from the TCGA-BRCA dataset. The gene expression matrix is obtained from the provided GDC file (c9f641bc-d9fe-4cd6-8e73-bbf6935aee52.rna_seq.augmented_star_gene_counts.tsv) and saved as matched_lines_TCGA-A7-A13D-01A-13R-A12P-07.txt.

## Step 10: Retrieve Gene Expression Matrix for Sample TCGA-E9-A1RH-11A-34R-A169-07

Similarly, the pipeline retrieves the gene expression matrix for the gene IDs obtained from the previous step for the sample "TCGA-E9-A1RH-11A-34R-A169-07" from the TCGA-BRCA dataset. The gene expression matrix is obtained from the provided GDC file (269c35f0-a4f7-4e30-a69f-f1f3b7b5dace.rna_seq.augmented_star_gene_counts.tsv) and saved as matched_lines_TCGA-E9-A1RH-11A-34R-A169-07.txt.
