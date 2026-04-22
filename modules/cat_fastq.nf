#! /usr/bin/env nextflow
// Process to concatenate the fastq files
process CONCAT_FASTQ {
    tag "${meta.id}"
    publishDir "${params.outdir}/concatenated_fastq/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), val(reads)

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.fastq.gz"), emit: concat_reads
    
    script:
    """
    cat ${reads[0].join(' ').replaceAll('[\\[\\],]', '')} > ${meta.id}_R1.fastq.gz
    cat ${reads[1].join(' ').replaceAll('[\\[\\],]', '')} > ${meta.id}_R2.fastq.gz
    """
    
    stub:
    
    """
    touch ${meta.id}_R1.fastq.gz
    touch ${meta.id}_R2.fastq.gz
    """
}

