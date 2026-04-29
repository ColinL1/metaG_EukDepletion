#!/usr/bin/env nextflow

process MEGAHIT_PE {
    tag "${meta.id}"
    label "big_mem"
    conda "bioconda::megahit"
    publishDir "${params.outdir}/assembly/${meta.id}", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path ("${meta.id}/final.contigs.fa"), emit: contigs_fa
    tuple val(meta), path("${meta.id}/log"), emit: log


    script:
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} --out-dir ${meta.id} --k-min ${params.k_min ?: 27} --k-max ${params.k_max ?: 127} --k-step 10 --num-cpu-threads ${task.cpus}
    """
    stub:
    """
    mkdir ${meta.id}
    touch ${meta.id}/final.contigs.fa
    touch ${meta.id}/log
    """
}
