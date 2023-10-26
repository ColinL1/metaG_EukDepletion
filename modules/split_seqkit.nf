#!/usr/bin/env nextflow

process split_bac {
    tag "${sample}"
    // cpus "${params.cpusMin}"
    // errorStrategy 'ignore'
    publishDir "$params.outdir/reads/bacteria/${sample}", mode: 'symlink'

    input: 
    tuple val(sample), path(kaiju_out), path(reads)

    output:
    tuple val ("${sample}.bacteria"), path ("${sample}.bacteria.f*.gz"), emit: bacteria_reads
    tuple val ("${sample}.non-bacteria"), path ("${sample}.non-bacteria.f*.gz"), emit: non_bacteria_reads


    script:
    if ( params.mode == 'contigs' )
    """
    grep  -w "C" ${kaiju_out} | cut -d\$'\t' -f2 > list_bacteria.txt
    grep  -w "U" ${kaiju_out} | cut -d\$'\t' -f2 > list_not-bacteria.txt
    seqkit grep -j ${task.cpus} -f list_bacteria.txt -i ${reads} > ${sample}.bacteria.fa
    seqkit grep -j ${task.cpus} -f list_not-bacteria.txt -i ${reads} > ${sample}.non-bacteria.fa
    pigz -p ${task.cpus} ${sample}.non-bacteria.fa ${sample}.bacteria.fa
    rm -rf list_bacteria.txt list_not-bacteria.txt
    """
    else
    """
    grep -w "C" ${kaiju_out} | cut -d\$'\t' -f2 > list_bacteria.txt
    grep -w "U" ${kaiju_out} | cut -d\$'\t' -f2 > list_not-bacteria.txt
    seqkit grep -j ${task.cpus} -f list_bacteria.txt -i ${reads} > ${sample}.bacteria.fq
    seqkit grep -j ${task.cpus} -f list_not-bacteria.txt -i ${reads} > ${sample}.non-bacteria.fq
    pigz -p ${task.cpus} ${sample}.bacteria.fq ${sample}.non-bacteria.fq
    rm -rf list_bacteria.txt list_not-bacteria.txt
    """
}

process split_bac_pe {
    tag "${sample}"
    // cpus "${params.cpusMin}"
    publishDir "$params.outdir/fastq_split/bacteria/${sample}", mode: 'symlink'

    input: 
    tuple val(sample), path(kaiju_out), path(reads)

    output:
    tuple val ("${sample}.bacteria"), path ("${sample}_PE_*.bacteria.fq.gz"), emit: bacteria_reads
    tuple val ("${sample}.non-bacteria"), path ("${sample}_PE_*.non-bacteria.fq.gz"), emit: non_bacteria_reads


    script:
    """
    grep -w "C" ${kaiju_out} | cut -d\$'\t' -f2 > list_bacteria.txt
    grep -w "U" ${kaiju_out} | cut -d\$'\t' -f2 > list_not-bacteria.txt
    seqkit grep -j ${task.cpus} -f list_bacteria.txt -i ${reads[0]} > ${sample}_PE_1.bacteria.fq
    seqkit grep -j ${task.cpus} -f list_bacteria.txt -i ${reads[1]} > ${sample}_PE_2.bacteria.fq
    seqkit grep -j ${task.cpus} -f list_not-bacteria.txt -i ${reads[0]} > ${sample}_PE_1.non-bacteria.fq
    seqkit grep -j ${task.cpus} -f list_not-bacteria.txt -i ${reads[1]} > ${sample}_PE_2.non-bacteria.fq
    pigz -p ${task.cpus} *.fq
    """
}