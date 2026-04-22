#!/usr/bin/env nextflow

process add_kingdom {
    tag "${sample}"
    cpus "${params.cpusMin}"
    publishDir "$params.outdir/diamond_blastx/Kingdom", mode: 'symlink'
    
    input: 
    tuple val(sample), path(blastx)

    output:
    path ("${sample}_blastx_kingdom.txt"), emit: kingdoms_blastx
    

    script:
    """
    python ${params.py_scripts}/add_kingdoms_blastx.py  -b  ${blastx}  -o ${sample}_blastx_kingdom.txt
    """
}