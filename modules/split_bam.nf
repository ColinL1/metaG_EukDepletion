#!/usr/bin/env nextflow

process SPLIT_BAM {
    tag "${sample}"
    // cpus "${params.cpusMin}"
    publishDir "$params.outdir/${seq_type}/reads/${ref_name}/${base_name}", mode: 'symlink'

    input: 
    // tuple val(base_name), path(bam)
    // tuple val(ref_name), path (reference)
    tuple val(sample), val(base_name), path(bam), val (seq_type), val(ref_name), val (reference)


    output:
    tuple val(sample), val ("${base_name}.unmapped"), path ("${base_name}.unmapped.fq.gz"), val (seq_type), emit: unmapped_reads
    tuple val(sample), val ("${base_name}.mapped"), path ("${base_name}.mapped.fq.gz"), val (seq_type), emit: mapped_reads

    script:
    """
    samtools fastq -n -f 4 ${bam} --threads ${task.cpus} | pigz -p ${task.cpus} > ${base_name}.unmapped.fq.gz
    samtools fastq -n -F 4 ${bam} --threads ${task.cpus} | pigz -p ${task.cpus} > ${base_name}.mapped.fq.gz
    """
    stub:
    """
    touch ${base_name}.unmapped.fq.gz
    touch ${base_name}.mapped.fq.gz
    """
}

process SPLIT_BAM_BWA {
    tag "${sample_name}"
    // cpus "${params.cpusMin}"
    publishDir "$params.outdir/${seq_type}/reads/${ref_name}/${sample_name}", mode: 'symlink'

    input: 
    tuple val(sample_name), path(bam)
    tuple val(ref_name), path (reference)

    output:
    tuple val ("${sample_name}.unmapped"), path ("${sample_name}.unmapped_{1,2}.fq.gz"), emit: unmapped_reads
    tuple val ("${sample_name}.mapped"), path ("${sample_name}.mapped_{1,2}.fq.gz"), emit: mapped_reads
    tuple val ("${sample_name}.singleton"), path ("${sample_name}*_singleton.fq.gz"), emit: singleton_reads
    script:
    """
    samtools fastq -f 4 --threads ${task.cpus} -n ${bam}  -1 ${sample_name}.unmapped_1.fq -2 ${sample_name}.unmapped_2.fq -0 /dev/null -s ${sample_name}.unmapped_singleton.fq
    samtools fastq -F 4 --threads ${task.cpus} -n ${bam}  -1 ${sample_name}.mapped_1.fq -2 ${sample_name}.mapped_2.fq -0 /dev/null -s ${sample_name}.mapped_singleton.fq
    pigz -p ${task.cpus} *.fq
    """
    stub:
    """
    touch ${base_name}.unmapped_1.fq.gz
    touch ${base_name}.unmapped_2.fq.gz
    touch ${base_name}.mapped_1.fq.gz
    touch ${base_name}.mapped_2.fq.gz
    touch ${sample_name}1_singleton.fq.gz
    touch ${sample_name}2_singleton.fq.gz
    """
}
