#!/usr/bin/env nextflow

process map2ref {
    tag "${sample_name}"
    cpus "${params.cpusHigh}"
    publishDir "$params.outdir/bam/${ref_name}", mode: 'symlink'

    input: 
    tuple val(sample_name), path(reads)
    tuple val(ref_name), path (reference)

    output:
    tuple val ("${sample_name}.${ref_name}"), path ("${sample_name}.${ref_name}.bam"), emit: mapp_file

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont ${reference} ${reads} --split-prefix=tmp | samtools view -S -b > ${sample_name}.${ref_name}.bam
    """
}