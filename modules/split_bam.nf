#!/usr/bin/env nextflow

process split_bam {
    tag "${sample_name}"
    // cpus "${params.cpusMin}"
    publishDir "$params.outdir/fastq_split/${ref_name}/${sample_name}", mode: 'symlink'

    input: 
    tuple val(sample_name), path(bam)
    tuple val(ref_name), path (reference)

    output:
    tuple val ("${sample_name}.unmapped"), path ("${sample_name}.unmapped.fq.gz"), emit: unmapped_reads
    tuple val ("${sample_name}.mapped"), path ("${sample_name}.mapped.fq.gz"), emit: mapped_reads

    script:
    """
    samtools fastq -n -f 4 ${bam} --threads ${task.cpus} | pigz -p ${task.cpus} > ${sample_name}.unmapped.fq.gz
    samtools fastq -n -F 4 ${bam} --threads ${task.cpus} | pigz -p ${task.cpus} > ${sample_name}.mapped.fq.gz
    """
}

