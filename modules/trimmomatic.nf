#!/usr/bin/env nextflow
// Trimmomatic illumina reads trimming

process TRIMMOMATIC {
    tag "${meta.id}"
    label "med_mem"
    conda "bioconda::trimmomatic"
    publishDir "${params.outdir}/trimmed_fastq/${meta.species}/${meta.id}/reads/", mode: 'symlink'

    input: 
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_TRIM_paried*"), emit: trimmed_fq_out
    tuple val(meta), path("${meta.id}_TRIM_unparied*"), emit: singleton_fq_out
    tuple val(meta), path("${meta.id}_trimmomatic.log"), emit: log_trim

    script:
    """
    trimmomatic PE -threads ${task.cpus} \\
            ${reads[0]} ${reads[1]} \\
            ${meta.id}_TRIM_paried_R1.fq.gz ${meta.id}_TRIM_unparied_R1.fq.gz \\
            ${meta.id}_TRIM_paried_R2.fq.gz ${meta.id}_TRIM_unparied_R2.fq.gz \\
            ILLUMINACLIP:${params.adapters}:${params.illuminaclip} \\
            HEADCROP:${params.headcrop} LEADING:${params.leading} TRAILING:${params.trailing} SLIDINGWINDOW:${params.slidingwindow} MINLEN:${params.minlen} \\
            2>> ${meta.id}_trimmomatic.log 1>>file.out
    """
    stub:
    """
    touch ${meta.id}_TRIM_paried_R1.fq.gz
    touch ${meta.id}_TRIM_unparied_R1.fq.gz
    touch ${meta.id}_TRIM_paried_R2.fq.gz
    touch ${meta.id}_TRIM_unparied_R2.fq.gz
    touch ${meta.id}_trimmomatic.log
    """
}
