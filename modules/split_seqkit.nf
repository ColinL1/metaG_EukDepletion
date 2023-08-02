#!/usr/bin/env nextflow

process split_bac {
    tag "${sample_name}"
    // cpus "${params.cpusMin}"
    publishDir "$params.outdir/fastq_split/bacteria/${sample_name}", mode: 'symlink'

    input: 
    tuple val(sample_name), path(kaiju_out), path(reads)

    output:
    tuple val ("${sample_name}.bacteria"), path ("${sample_name}.bacteria.f*.gz"), emit: bacteria_reads
    tuple val ("${sample_name}.non-bacteria"), path ("${sample_name}.non-bacteria.f*.gz"), emit: non_bacteria_reads


    script:
    if ( params.mode == 'contigs' )
    """
    grep "C" ${kaiju_out} | cut -d\$'\t' -f2 > list_bacteria.txt
    grep "U" ${kaiju_out} | cut -d\$'\t' -f2 > list_not-bacteria.txt
    seqkit grep -j ${task.cpus} -f list_bacteria.txt -i ${reads} > ${sample_name}.bacteria.fa
    seqkit grep -j ${task.cpus} -f list_not-bacteria.txt -i ${reads} > ${sample_name}.non-bacteria.fa
    pigz -p ${task.cpus} *.fa
    """
    else
    """
    grep "C" ${kaiju_out} | cut -d\$'\t' -f2 > list_bacteria.txt
    grep "U" ${kaiju_out} | cut -d\$'\t' -f2 > list_not-bacteria.txt
    seqkit grep -j ${task.cpus} -f list_bacteria.txt -i ${reads} > ${sample_name}.bacteria.fq
    seqkit grep -j ${task.cpus} -f list_not-bacteria.txt -i ${reads} > ${sample_name}.non-bacteria.fq
    pigz -p ${task.cpus} *.fq
    """
}