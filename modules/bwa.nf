#!/usr/bin/env nextflow

process BWA_PE {
    tag "${sample_name}"
    // cpus "${params.cpusHigh}"
    publishDir "$params.outdir/bam/${ref_name}/", mode: 'symlink'
    // publishDir "$params.outdir/reads/${ref_name}/", mode: 'symlink'

    input: 
    tuple val(sample_name), path(reads)
    tuple val(ref_name), val (reference)

    output:
    tuple val ("${sample_name}.${ref_name}"), path ("${sample_name}.${ref_name}.bam"), emit: mapp_file
    // tuple val ("${sample_name}.${ref_name}.unmapped"), path ("${sample_name}.${ref_name}.unmapped.fq.{1,2}.gz"), emit: unmapped_reads
    // tuple val ("${sample_name}.${ref_name}.mapped"), path ("${sample_name}.${ref_name}.mapped.fq.{1,2}.gz"), emit: mapped_reads


    script:
    """
    bwa mem -t ${task.cpus} ${reference} ${reads[0]} ${reads[1]} | samtools view -S -b > ${sample_name}.${ref_name}.bam        
    """
    stub:
    """
    touch ${sample_name}.${ref_name}.bam
    """
}

process BWA_ONT {
    tag "${sample_name}"
    // cpus "${params.cpusHigh}"
    publishDir "$params.outdir/bam/${ref_name}/", mode: 'symlink'
    // publishDir "$params.outdir/reads/${ref_name}/", mode: 'symlink'

    input: 
    tuple val(sample_name), path(reads)
    tuple val(ref_name), val (reference)

    output:
    tuple val ("${sample_name}.${ref_name}"), path ("${sample_name}.${ref_name}.bam"), emit: mapp_file
    // tuple val ("${sample_name}.${ref_name}.unmapped"), path ("${sample_name}.${ref_name}.unmapped.fq.{1,2}.gz"), emit: unmapped_reads
    // tuple val ("${sample_name}.${ref_name}.mapped"), path ("${sample_name}.${ref_name}.mapped.fq.{1,2}.gz"), emit: mapped_reads


    script:
    """
    bwa mem -t ${task.cpus} -x ont2d ${reference} ${reads} | samtools view -S -b > ${sample_name}.${ref_name}.bam
    """
    stub:
    """
    touch ${sample_name}.${ref_name}.bam
    """
}
