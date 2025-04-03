#!/usr/bin/env nextflow

//TODO: check if possible to inherit conditional directly from queue channel  
process MINIMAP2_INDEX {
    tag "${sample}"
    publishDir "$params.outdir/mapping_indexes/${taxon}/", mode: 'symlink'

    input: 
    tuple val(sample), val(taxon), path(ref_genome_fasta), val(seq_type)

    output:
    tuple val(sample), val(taxon), path("${taxon}.mmi"), val(seq_type), emit: indexed_reference

    script:
    // params.mode = "${seq_type}"
    if( "${seq_type}" == 'ONT' )
    """
    minimap2 -t ${task.cpus} -x map-ont -d ${taxon}.mmi ${ref_genome_fasta}
    """
    else if( "${seq_type}" == 'contigs_ont' )
	"""
	minimap2 -t ${task.cpus} -x asm5 -d ${taxon}.mmi ${ref_genome_fasta}
	"""
    else if( "${seq_type}" == 'contigs_illumina' )
    """
	minimap2 -t ${task.cpus} -x sr -d ${taxon}.mmi ${ref_genome_fasta}
    """
    else
    error "Invalid alignment mode: ${seq_type}"

    stub:
    if( "${seq_type}" == 'ONT' )
        """
        touch ${taxon}.mmi
        touch ONT_pathway
        """
    else if( "${seq_type}" == 'contigs_ont' )
	"""
        touch ${taxon}.mmi
        touch contigs_ont_pathway
	"""
    else if( "${seq_type}" == 'contigs_illumina' )
        """
        touch ${taxon}.mmi
        touch contigs_illumina_pathway 
        """
    else
    error "Invalid alignment mode: ${seq_type}"
}