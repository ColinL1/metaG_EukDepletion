#! /usr/bin/env nextflow
// Process to concatenate the fastq files
process CONCAT_FASTQ {
    tag "${meta.id}"
    label 'min_mem'
    publishDir "${params.outdir}/concatenated_fastq/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.fastq.gz"), emit: concat_reads
    
    script:
    def r1_files = reads[0].collect { it.toString() }.join(' ')
    def r2_files = reads[1].collect { it.toString() }.join(' ')
    """
    cat ${r1_files} > ${meta.id}_R1.fastq.gz
    cat ${r2_files} > ${meta.id}_R2.fastq.gz
    """
    stub:
    
    """
    touch ${meta.id}_R1.fastq.gz
    touch ${meta.id}_R2.fastq.gz
    """
}

