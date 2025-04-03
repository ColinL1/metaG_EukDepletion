#!/usr/bin/env nextflow

process SPLIT_BAM_HOST {
    tag "${meta.id}"
    label "med_mem"
    publishDir "${
    meta.species == 'H2' ? "$baseDir/results/mapping/Aiptasia/${meta.id}/reads/" :
    meta.species == 'F003' ? "$baseDir/results/mapping/Aiptasia/${meta.id}/reads/" :
    meta.species == 'Porites' ? "$baseDir/results/mapping/Corals/${meta.id}/reads/" :
    meta.species == 'Acropora' ? "$baseDir/results/mapping/Corals/${meta.id}/reads/" :
    meta.species == 'Pocillopora' ? "$baseDir/results/mapping/Corals/${meta.id}/reads/" :
    meta.species == 'unknown' ? 'unknown_samples' :
    'other_samples'}", mode: 'symlink'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.unmapped.fa.gz"), emit: unmapped_reads
    tuple val(meta), path("*.mapped.fa.gz"), emit: mapped_reads

    script:
    """
    samtools view -F 0x900 -F 0x4 ${bam} --threads ${task.cpus} -o ${meta.id}.scleractina.mapped.fa
    samtools view -f 0x4 -F 0x900 ${bam} --threads ${task.cpus} -o ${meta.id}.scleractina.unmapped.fa
    pigz -p ${task.cpus} ${meta.id}.scleractina.mapped.fa ${meta.id}.scleractina.unmapped.fa
    """
    stub:
    """
    touch ${meta.id}.unmapped.fa.gz
    touch ${meta.id}.mapped.fa.gz
    """
}


process SPLIT_BAM_SYM {
    tag "${meta.id}"
    label "med_mem"
    publishDir "$baseDir/results/mapping/Symbiodiniaceae/${meta.id}/reads/", mode: 'symlink'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.unmapped.fa.gz"), emit: unmapped_reads
    tuple val(meta), path("*.mapped.fa.gz"), emit: mapped_reads

    script:
    """
    samtools view -F 0x900 -F 0x4 ${bam} --threads ${task.cpus} -o ${meta.id}.symbiodiniaceae.mapped.fa
    samtools view -f 0x4 -F 0x900 ${bam} --threads ${task.cpus} -o ${meta.id}.symbiodiniaceae.unmapped.fa
    pigz -p ${task.cpus} ${meta.id}.symbiodiniaceae.mapped.fa ${meta.id}.symbiodiniaceae.unmapped.fa
    """
    stub:
    """
    touch ${meta.id}.unmapped.fa.gz
    touch ${meta.id}.mapped.fa.gz
    """
}