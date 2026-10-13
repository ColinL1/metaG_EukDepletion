#!/usr/bin/env nextflow

process SPLIT_BAM_HOST {
    tag "${meta.id}"
    label "med_mem"
    conda "bioconda::samtools"
    container 'quay.io/biocontainers/samtools:1.23--h96c455f_0'
    publishDir path: { "${params.outdir}/assembly_based/mapping/${meta.species}/${meta.id}/" }, mode: 'symlink'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.unmapped.fa.gz"), emit: unmapped_reads
    tuple val(meta), path("*.mapped.fa.gz"), emit: mapped_reads

    script:
    """
    samtools view -F 0x900 -F 0x4 ${bam} --threads ${task.cpus} -o ${meta.id}.host.mapped.fa
    samtools view -f 0x4 -F 0x900 ${bam} --threads ${task.cpus} -o ${meta.id}.host.unmapped.fa
    gzip ${meta.id}.host.mapped.fa ${meta.id}.host.unmapped.fa
    """
    stub:
    """
    touch ${meta.id}.host.unmapped.fa.gz
    touch ${meta.id}.host.mapped.fa.gz
    """
}

process SPLIT_BAM_SYM {
    tag "${meta.id}"
    label "med_mem"
    conda "bioconda::samtools"
    container 'quay.io/biocontainers/samtools:1.23--h96c455f_0'
    publishDir path: { "${params.outdir}/assembly_based/mapping/symbiodiniaceae/${meta.id}/" }, mode: 'symlink'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.unmapped.fa.gz"), emit: unmapped_reads
    tuple val(meta), path("*.mapped.fa.gz"), emit: mapped_reads

    script:
    """
    samtools view -F 0x900 -F 0x4 ${bam} --threads ${task.cpus} -o ${meta.id}.symbiodiniaceae.mapped.fa
    samtools view -f 0x4 -F 0x900 ${bam} --threads ${task.cpus} -o ${meta.id}.symbiodiniaceae.unmapped.fa
    gzip ${meta.id}.symbiodiniaceae.mapped.fa ${meta.id}.symbiodiniaceae.unmapped.fa
    """
    stub:
    """
    touch ${meta.id}.symbiodiniaceae.unmapped.fa.gz
    touch ${meta.id}.symbiodiniaceae.mapped.fa.gz
    """
}