#!/usr/bin/env nextflow

process MINIMAP2_INDEX {
    tag "${meta.id} / ${mode}"
    label 'process_high'
    conda "bioconda::minimap2"
    container 'quay.io/biocontainers/minimap2:2.31--h118bc1c_0'
    publishDir path: { "${params.outdir}/coassembly_based/quantification/mapping_indexes/" }, mode: 'symlink'

    input:
    tuple val(meta), val(mode), path(fasta)

    output:
    tuple val(meta), val(mode), path("${meta.id}.${mode}.mmi"), emit: index

    script:
    """
    minimap2 -x sr -I 8G -t ${task.cpus} -d ${meta.id}.${mode}.mmi ${fasta}
    """
    stub:
    """
    touch ${meta.id}.${mode}.mmi
    """
}
