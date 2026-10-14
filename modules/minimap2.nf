#!/usr/bin/env nextflow

process MINIMAP2_MAP_HOST {
    tag "${meta.id}"
    label 'process_high'
    conda "bioconda::minimap2 bioconda::samtools"
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:b411340b52d82a9c276d87c7a3dcffc880be762f-0'
    publishDir path: { "${params.outdir}/assembly_based/mapping/${meta.species}/${meta.id}/" }, mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: mapp_file

    script:
    if( "${meta.species}" == 'H2' )
        """
        minimap2 -t ${task.cpus} -ax asm5 ${params.ref_mmi_aip} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.host.bam
        """
    else if( "${meta.species}" == 'F003' )
        """
        minimap2 -t ${task.cpus} -ax asm5 ${params.ref_mmi_aip} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.host.bam
        """
    else if( "${meta.species}" == 'Plob' )
        """
	    minimap2 -t ${task.cpus} -ax asm5 ${params.ref_mmi_scl} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.host.bam

        """
    else if ( "${meta.species}" == 'Ahemp' )
        """
	    minimap2 -t ${task.cpus} -ax asm5 ${params.ref_mmi_scl} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.host.bam
        """
    else if ( "${meta.species}" == 'Pver' )
        """
	    minimap2 -t ${task.cpus} -ax asm5 ${params.ref_mmi_scl} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.host.bam
        """
    else
        error "Invalid alignment mode: ${meta.species}"

    stub:
    if( "${meta.species}" == 'H2' )
        """
        touch ${meta.id}.host.bam
        touch ${meta.id}
        """
    else if( "${meta.species}" == 'F003' )
        """
        touch ${meta.id}.host.bam
        touch ${meta.id}
        """
    else if( "${meta.species}" == 'Plob' )
        """
        touch ${meta.id}.host.bam
        touch ${meta.id}
        """
    else if ( "${meta.species}" == 'Ahemp' )
        """
        touch ${meta.id}.host.bam
        touch ${meta.id}
        """
    else if ( "${meta.species}" == 'Pver' )
        """
        touch ${meta.id}.host.bam
        touch ${meta.id}
        """
    else
        error "Invalid alignment mode: ${meta.species}"
}

process MINIMAP2_MAP_SYM {
    tag "${meta.id}"
    label "process_high"
    conda "bioconda::minimap2 bioconda::samtools"
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:b411340b52d82a9c276d87c7a3dcffc880be762f-0'
    publishDir path: { "${params.outdir}/assembly_based/mapping/symbiodiniaceae/${meta.id}/" }, mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: mapp_file

    script:
	"""
	minimap2 -t ${task.cpus} -ax asm5 ${params.ref_mmi_sym} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.symbiodiniaceae.bam
	"""
    stub:
    """
    touch ${meta.id}.symbiodiniaceae.bam
    touch ${meta.id}
    """
}

process MINIMAP2_MAP_CATALOGUE {
    tag "${meta.id} / ${mode}"
    label 'process_high'
    conda "bioconda::minimap2 bioconda::samtools"
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:b411340b52d82a9c276d87c7a3dcffc880be762f-0'
    publishDir path: { "${params.outdir}/coassembly_based/quantification/${mode}/flagstat/" }, mode: 'copy', pattern: '*.flagstat'

    input:
    tuple val(meta), val(mode), path(reads), path(mmi)

    output:
    tuple val(meta), val(mode), path("${meta.id}.${mode}.bam"), path("${meta.id}.${mode}.bam.bai"), emit: bam
    tuple val(meta), val(mode), path("${meta.id}.${mode}.flagstat"), emit: flagstat

    script:
    """
    minimap2 -ax sr --secondary=no -t ${task.cpus} ${mmi} ${reads[0]} ${reads[1]} \\
      | samtools sort -@ 4 -o ${meta.id}.${mode}.bam -
    samtools index ${meta.id}.${mode}.bam
    samtools flagstat ${meta.id}.${mode}.bam > ${meta.id}.${mode}.flagstat
    """
    stub:
    """
    touch ${meta.id}.${mode}.bam ${meta.id}.${mode}.bam.bai ${meta.id}.${mode}.flagstat
    """
}
