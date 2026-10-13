#!/usr/bin/env nextflow

process CAT {
    tag "${meta.id}"
    label 'process_high'
    conda "bioconda::cat=5.3"
    container 'quay.io/biocontainers/cat:5.3--hdfd78af_0'
    publishDir path: { "${params.outdir}/coassembly_based/CAT/${meta.species}/${meta.id}/" }, mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path ("${meta.id}.contig2classification.txt"), emit: contig2classification
    tuple val(meta), path ("${meta.id}.ORF2LCA.txt"), emit: orf2lca
    tuple val(meta), path("*.log"), emit: log

    script:
    """
    CAT contigs -c ${reads} -d ${params.cat_db}/db -t ${params.cat_db}/tax/ -o ${meta.id} -n ${task.cpus} --compress
    """
    stub:
    """
    touch ${meta.id}.contig2classification.txt
    touch ${meta.id}.ORF2LCA.txt
    touch ${meta.id}.log
    """
}

process CAT_ADD_NAMES {
    tag "${meta.id}"
    errorStrategy  { task.attempt <= 3 ? 'retry' : 'ignore' }
    label 'min_mem'
    conda "bioconda::cat=5.3"
    container 'quay.io/biocontainers/cat:5.3--hdfd78af_0'
    publishDir path: { "${params.outdir}/coassembly_based/CAT/${meta.species}/${meta.id}/" }, mode: 'symlink'

    input:
    tuple val(meta), path(report)

    output:
    // each task receives one report, so only one of the two outputs exists
    tuple val(meta), path ("*.contig2classification*.txt"), optional: true, emit: contig2classification_taxname
    tuple val(meta), path ("*.ORF2LCA*.txt"), optional: true, emit: orf2lca_taxname


    script:
    """
    CAT add_names -i  ${report} -o ${report}.tax_name.txt  -t ${params.cat_db}/tax --only_official --exclude_scores
    """
    stub:
    """
    touch ${report}.tax_name.txt
    """
}

process CAT_FILTER_BACTERIA {
    tag "${meta.id}"
    label 'min_mem'
    conda "bioconda::seqkit=2.8.2"
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_1'
    publishDir path: { "${params.outdir}/coassembly_based/CAT/${meta.species}/${meta.id}/bacteria/" }, mode: 'symlink'

    input:
    tuple val(meta), path(contigs), path(c2c_names)

    output:
    tuple val(meta), path("${meta.id}.bacteria.fa"), emit: bacteria_fa
    tuple val(meta), path("${meta.id}.contig2classification.bacteria.txt"), emit: bacteria_class

    script:
    """
    awk -F'\\t' '
        NR == 1 { for (i = 1; i <= NF; i++) if (\$i == "superkingdom") col = i
                  if (!col) { print "no superkingdom column in ${c2c_names}" > "/dev/stderr"; exit 1 }
                  print; next }
        \$col == "Bacteria"' ${c2c_names} > ${meta.id}.contig2classification.bacteria.txt
    awk -F'\\t' 'NR > 1 { print \$1 }' ${meta.id}.contig2classification.bacteria.txt > ${meta.id}.contig.list
    seqkit grep -f ${meta.id}.contig.list ${contigs} > ${meta.id}.bacteria.fa
    """
    stub:
    """
    touch ${meta.id}.bacteria.fa
    touch ${meta.id}.contig2classification.bacteria.txt
    """
}

process CAT_PREP_FAMILY {
    tag "${meta.id}"
    label 'min_mem'
    conda "conda-forge::r-base=4.3 conda-forge::r-data.table"
    container 'rocker/tidyverse:4.3.3'

    input:
    tuple val(meta), path(bacteria_class)

    output:
    path("${meta.id}.family.tsv"), emit: family

    script:
    """
    prep_cat.R ${bacteria_class} ${meta.id} ${meta.id}.family.tsv
    """
    stub:
    """
    touch ${meta.id}.family.tsv
    """
}

process CAT_MERGE_TABLES {
    tag "${mode}"
    label 'min_mem'
    conda "conda-forge::r-base=4.3 conda-forge::r-data.table"
    container 'rocker/tidyverse:4.3.3'
    publishDir path: { "${params.outdir}/coassembly_based/quantification/${mode}/" }, mode: 'copy'

    input:
    tuple val(mode), path(coverm_files)
    path cat_files
    path samplesheet

    output:
    path("contig_table.${mode}.tsv.gz"), emit: contig_table
    path("family_table.${mode}.tsv"), emit: family_table
    path("mapping_summary.${mode}.tsv"), emit: mapping_summary

    script:
    """
    merge_tables.R ${mode} ${samplesheet} "${coverm_files}" "${cat_files}"
    """
    stub:
    """
    touch contig_table.${mode}.tsv.gz family_table.${mode}.tsv mapping_summary.${mode}.tsv
    """
}
