#!/usr/bin/env nextflow

process CAT {
    tag "${meta.id}"
    label 'process_high'
    conda "bioconda::cat=5.3"
    container 'quay.io/biocontainers/cat:5.3--hdfd78af_0'
    publishDir path: { "${params.outdir}/CAT/${meta.species}/${meta.id}/" }, mode: 'symlink'

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
    publishDir path: { "${params.outdir}/CAT/${meta.species}/${meta.id}/" }, mode: 'symlink'

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path ("*.contig2classification*.txt"), emit: contig2classification_taxname
    tuple val(meta), path ("*.ORF2LCA*.txt"), emit: orf2lca_taxname


    script:
    """
    CAT add_names -i  ${report} -o ${report}.tax_name.txt  -t ${params.cat_db}/tax
    """
    stub:
    """
    touch ${meta.id}.contig2classification.tax_name.txt
    touch ${meta.id}.ORF2LCA.tax_name.txt
    """
}
