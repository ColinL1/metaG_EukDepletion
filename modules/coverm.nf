#!/usr/bin/env nextflow

process COVERM_CONTIG {
    tag "${meta.id} / ${mode}"
    label 'med_mem'
    conda "bioconda::coverm=0.7.0"
    container 'quay.io/biocontainers/coverm:0.7.0--hcb7b614_4'
    publishDir path: { "${params.outdir}/coassembly_based/quantification/${mode}/coverm/" }, mode: 'copy'

    input:
    tuple val(meta), val(mode), path(bam), path(bai)

    output:
    tuple val(mode), path("${meta.id}.${mode}.coverm.tsv"), emit: coverage

    script:
    """
    coverm contig \\
        --bam-files ${bam} \\
        --methods length count mean trimmed_mean covered_fraction tpm \\
        --min-read-percent-identity ${params.min_read_id} \\
        --min-read-aligned-percent  ${params.min_read_aln} \\
        --threads ${task.cpus} \\
        --output-file ${meta.id}.${mode}.coverm.tsv
    """
    stub:
    """
    touch ${meta.id}.${mode}.coverm.tsv
    """
}
