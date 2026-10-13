#!/usr/bin/env nextflow

process MEGAHIT_PE {
    tag "${meta.id}"
    label "big_mem"
    conda "bioconda::megahit"
    container 'quay.io/biocontainers/megahit:1.2.9--haf24da9_8'
    publishDir path: { "${params.outdir}/${task.process.contains('COASSEMBLY') ? 'coassembly_based' : 'assembly_based'}/assembly/${meta.id}" }, mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path ("${meta.id}/${meta.id}.contigs.fa"), emit: contigs_fa
    tuple val(meta), path("${meta.id}/log"), emit: log


    script:
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} --out-dir ${meta.id} --k-min ${params.k_min ?: 27} --k-max ${params.k_max ?: 127} --k-step 10 --num-cpu-threads ${task.cpus}
    mv ${meta.id}/final.contigs.fa ${meta.id}/${meta.id}.contigs.fa
    """
    stub:
    """
    mkdir ${meta.id}
    touch ${meta.id}/${meta.id}.contigs.fa
    touch ${meta.id}/log
    """
}
