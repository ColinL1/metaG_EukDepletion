#!/usr/bin/env nextflow

process bracken {
    tag "${sample}"
    // cpus "${params.cpusHigh}"
    // memory "${params.memMax}"
    publishDir "$params.outdir/bracken_reports/", mode: 'symlink'

    input: 
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.bracken"), emit: bracken_out
    tuple val(sample), path("${sample}.breport"), emit: report_bracken_out
    tuple val(sample), path("${sample}.log"), emit: log_bracken

    script:
    """
    bracken -d ${params.kraken_db} -i ${sample}.k2report -o ${sample}.bracken -w ${sample}.breport -r ${params.read_length} -l S -t 10 2>&1 > ${sample}.log
    """
}
