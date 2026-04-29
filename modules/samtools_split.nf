#!/usr/bin/env nextflow

process SPLIT_BAM_HOST {
    tag "${meta.id}"
    label "med_mem"
    conda "bioconda::samtools" 
    publishDir "${params.outdir}/mapping/${meta.species}/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.unmapped.fa.gz"), emit: unmapped_reads
    tuple val(meta), path("*.mapped.fa.gz"), emit: mapped_reads

    script:
    """
    samtools view -F 0x900 -F 0x4 ${bam} --threads ${task.cpus} -o ${meta.id}.scleractina.mapped.fa
    samtools view -f 0x4 -F 0x900 ${bam} --threads ${task.cpus} -o ${meta.id}.scleractina.unmapped.fa
    gzip ${meta.id}.scleractina.mapped.fa ${meta.id}.scleractina.unmapped.fa
    """
    stub:
    """
    touch ${meta.id}.scleractina.unmapped.fa.gz
    touch ${meta.id}.scleractina.mapped.fa.gz
    """
}

process SPLIT_BAM_SYM {
    tag "${meta.id}"
    label "med_mem"
    conda "bioconda::samtools"
    publishDir "${params.outdir}/mapping/Symbiodiniaceae/${meta.id}/", mode: 'symlink'

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