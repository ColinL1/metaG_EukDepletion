#!/usr/bin/env nextflow

process MINIMAP2_MAP_HOST {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/mapping/${meta.species}/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: mapp_file

    script:
    if( "${meta.species}" == 'H2' )
        """
        minimap2 -t ${task.cpus} -ax asm5 ${params.ref_file_aip} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.aiptasiidae.bam
        """
    else if( "${meta.species}" == 'F003' )
        """
        minimap2 -t ${task.cpus} -ax asm5 ${params.ref_file_aip} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.aiptasiidae.bam
        """
    else if( "${meta.species}" == 'Porites' )
        """
	    minimap2 -t ${task.cpus} -ax asm5 ${params.ref_file_scl} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.scleractina.bam
        
        """
    else if ( "${meta.species}" == 'Acropora' )
        """
	    minimap2 -t ${task.cpus} -ax asm5 ${params.ref_file_scl} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.scleractina.bam
        """
    else if ( "${meta.species}" == 'Pocillopora' )
        """
	    minimap2 -t ${task.cpus} -ax asm5 ${params.ref_file_scl} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.scleractina.bam
        """
    else
    // error "Invalid alignment mode: ${meta.species}"
        """
        minimap2 -t ${task.cpus} -ax map-ont ${params.ref_file_scl_ont} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.scleractina.bam
        """

    stub:
    if( "${meta.species}" == 'H2' )
        """
        touch ${meta.id}.${meta.species}-map-H2.bam
        touch ${meta.id}
        """
    else if( "${meta.species}" == 'F003' )
        """
        touch ${meta.id}.${meta.species}-map-F003.bam
        touch ${meta.id}
        """
    else if( "${meta.species}" == 'Porites' )
        """
        touch ${meta.id}.${meta.species}-map-Porites.bam
        touch ${meta.id}
        """
    else if ( "${meta.species}" == 'Acropora' )
        """
        touch ${meta.id}.${meta.species}-map-Acropora.bam
        touch ${meta.id}
        """
    else if ( "${meta.species}" == 'Pocillopora' )
        """
        touch ${meta.id}.${meta.species}-map-Pocillopora.bam
        touch ${meta.id}
        """
    else
    // error "Invalid alignment mode: ${meta.species}"
        """
        touch ${meta.id}.${meta.species}-map-ONT.bam
        touch ${meta.id}
        """
}

process MINIMAP2_MAP_SYM {
    tag "${meta.id}"
    label "process_high"
    publishDir "${params.outdir}/mapping/Symbiodiniaceae/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: mapp_file

    script:
	"""
	minimap2 -t ${task.cpus} -ax map-ont ${params.ref_file_sym_ont} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.symbiodiniaceae.bam
	"""
    stub:
    """
    touch ${meta.id}.symbiodiniaceae.bam
    touch ${meta.id}
    """
}



