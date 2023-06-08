params.reference = 'GCA_000001405.15_GRCh38_full_analysis_set.fna'
params.library = 'library.fa'
params.gdc_file1 = '/gdc_download_20230608_115738.970040/bf6e2b60-d65b-4963-8e8e-49ab3f93c3e1'
params.gdc_file2 = 'tca.tsv'

// Step 1: Download reference genome
process downloadReference {
    output:
    file 'reference.fasta'

    script:
    """
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
    gunzip GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
    mv GCA_000001405.15_GRCh38_full_analysis_set.fna reference.fasta
    """
}

// Step 2: Build index with bowtie2
process buildIndex {
    input:
    file reference

    output:
    file 'human_ref.1.bt2'

    script:
    """
    bowtie2-build -f --seed 0 ${reference} human_ref
    """
}

// Step 3: Map library.fa to reference genome
process mapLibrary {
    input:
    file reference, file library

    output:
    file 'alignment.sam'

    script:
    """
    bowtie2 -x ${reference} -f ${library} -S alignment.sam
    """
}

// Step 4: Sort data with samtools
process sortAlignment {
    input:
    file 'alignment.sam'

    output:
    file 'alignment.bam'

    script:
    """
    samtools sort alignment.sam > alignment.bam
    """
}

// Step 5: Get chromosome name, start and end positions, and strand
process extractInfo {
    input:
    file 'alignment.bam'

    output:
    file 'gene_info.txt'

    script:
    """
    samtools view -bS alignment.sam > alignment.bam
    bedtools bamtobed -i alignment.bam | awk -F "\t" '{OFS="\t"}{print $4,$1,$2,$3,$6}' > gene_info.txt
    """
}

// Step 6: Get gene names from fasta file descriptions
process extractGeneNames {
    input:
    file library

    output:
    file 'gene_names.txt'

    script:
    """
    awk -F'|' '/^>/{print \$3}' ${library} > gene_names.txt
    """
}

// Step 7: Get gene names based on mapping results
process mapGeneNames {
    input:
    file 'alignment.sam', file 'gene_names.txt'

    output:
    file 'mapped_genes.txt'

    script:
    """
    samtools view alignment.sam | awk -F'\t' -v OFS='\t' '
        BEGIN {
            while (getline gene < "${gene_names}") {
                gene_names[gene] = 1
            }
            close("${gene_names}")
        }
        {
            if (\$3 != "*") {
                gene_name = (\$0 ~ /XO:i:([0-9]+)/) ? (\$0 ~ /XO:i:([0-9]+)/ && substr(\$0, RSTART+5, RLENGTH-5) in gene_names ? substr(\$0, RSTART+5, RLENGTH-5) : "N/A") : "N/A"
                print \$0, gene_name
            } else {
                print \$0, "N/A"
            }
        }
    ' > mapped_genes.txt
    cut -f 1 mapped_genes.txt > mapped_names.txt
    cut -d "|" -f 3 mapped_names.txt > mapped_genes.txt
    """
}

// Step 8: Check differences in annotation between fasta genes and mapped genes
process checkAnnotation {
    input:
    file 'gene_names.txt', file 'mapped_genes.txt'

    script:
    """
    awk 'NR==FNR { lines[NR]=\$0; next } { if (\$0 == lines[FNR]) print; else print "" }' gene_names.txt mapped_genes.txt
    """
}

// Step 9: Get gene matrix for genes based on GDC file 1
process getGeneMatrix1 {
    input:
    file 'gene_names.txt'

    output:
    file 'matched_lines1.txt'

    script:
    """
    uniq gene_names.txt > uniq_genes.txt
    awk 'FNR==NR { values[\$1]; next } \$2 in values' uniq_genes.txt ${params.gdc_file1} > matched_lines1.txt
    """
}

// Step 10: Get gene matrix for genes based on GDC file 2
process getGeneMatrix2 {
    input:
    file 'gene_names.txt'

    output:
    file 'matched_lines2.txt'

    script:
    """
    awk 'FNR==NR { values[\$1]; next } \$2 in values' uniq_genes.txt ${params.gdc_file2} > matched_lines2.txt
    """
}

// Define workflow
workflow {
    reference = downloadReference()
    index = buildIndex(reference)
    alignment = mapLibrary(reference, params.library)
    sortedAlignment = sortAlignment(alignment)
    geneInfo = extractInfo(sortedAlignment)
    geneNames = extractGeneNames(params.library)
    mappedGenes = mapGeneNames(alignment, geneNames)
    annotationDiff = checkAnnotation(geneNames, mappedGenes)
    geneMatrix1 = getGeneMatrix1(geneNames)
    geneMatrix2 = getGeneMatrix2(geneNames)

    // Specify output dependencies if needed
    sortedAlignment.outpath = file('alignment.bam')
    geneInfo.outpath = file('gene_info.txt')
    geneNames.outpath = file('gene_names.txt')
    mappedGenes.outpath = file('mapped_genes.txt')
    annotationDiff.outpath = file('annotation_diff.txt')
    geneMatrix1.outpath = file('matched_lines1.txt')
    geneMatrix2.outpath = file('matched_lines2.txt')
}
