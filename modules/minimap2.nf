#!/usr/bin/env nextflow

//TODO: check if possible to inherit conditional directly from queue channel  
process MINIMAP2_MAP_HOST {
    tag "${meta.id}"
    label "process_high"
    publishDir "${
        meta.species == 'H2' ? "$baseDir/results/mapping/Aiptasia/${meta.id}/bam/" :
        meta.species == 'F003' ? "$baseDir/results/mapping/Aiptasia/${meta.id}/bam/" :
        meta.species == 'Porites' ? "$baseDir/results/mapping/Corals/${meta.id}/bam/" :
        meta.species == 'Acropora' ? "$baseDir/results/mapping/Corals/${meta.id}/bam/" :
        meta.species == 'Pocillopora' ? "$baseDir/results/mapping/Corals/${meta.id}/bam/" :
        meta.species == 'unknown' ? 'unknown_samples' :
        'other_samples'}", mode: 'symlink'


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
    error "Invalid alignment mode: ${meta.species}"

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
    error "Invalid alignment mode: ${meta.species}"

}

params.ref_file_aip = '/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_minimap2_contigs/all_aiptasiidae.mmi'
params.ref_file_scl = '/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_minimap2_contigs/all_scleractinia.mmi'

process MINIMAP2_MAP_SYM {
    tag "${meta.id}"
    label "process_high"
    publishDir "$baseDir/results/mapping/Symbiodiniaceae/${meta.id}/bam/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: mapp_file

    script:
	"""
	minimap2 -t ${task.cpus} -ax asm5 ${params.ref_file_2} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.symbiodiniaceae.bam
	"""
    stub:
    """
    touch ${meta.id}.bam
    touch ${meta.id}
    """
}

params.ref_file_2 = '/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_minimap2_contigs/all_symbiodiniaceae.mmi'

