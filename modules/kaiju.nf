#!/usr/bin/env nextflow

process KAIJU_SE {
    tag "${meta.id}"
    label "process_high"
    publishDir "${params.outdir}/mapping/bacteria-kaiju/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.out"), emit: report_kaiju_out
    
    script:
    """
    kaiju -t ${params.nodes} -f ${params.kaiju_db} -i ${reads} -a mem -z ${task.cpus} -o ${meta.id}.out
    """
    stub:
    """
    touch ${meta.id}.out
    """
}

process KAIJU_PE {
    tag "${meta.id}"
    label "process_high"
    publishDir "${params.outdir}/mapping/bacteria-kaiju/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.out"), emit: report_kaiju_out

    script:
    """
    kaiju -t ${params.nodes} -f ${params.kaiju_db} -i ${reads[0]} -j ${reads[0]} -a mem -z ${task.cpus} -o ${meta.id}.out
    """
    stub:
    """
    touch ${meta.id}.out
    """
}
