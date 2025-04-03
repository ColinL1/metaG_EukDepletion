#!/usr/bin/env nextflow

process BRACKEN { 
    tag "${sample}"
    publishDir "$params.outdir/bracken_reports/", mode: 'symlink'

    input: 
    tuple val(sample), val(base_name), path(reads), val (seq_type)

    output:
    tuple val(sample), val(base_name), path("${base_name}.bracken"), val (seq_type), emit: bracken_out
    tuple val(sample), val(base_name), path("${base_name}.breport"), val (seq_type), emit: report_bracken_out
    tuple val(sample), val(base_name), path("${base_name}.log"), val (seq_type), emit: log_bracken

    script:
    """
    bracken -d ${params.kraken_db} -i ${base_name}.k2report -o ${base_name}.bracken -w ${base_name}.breport -r ${params.read_length} -l S -t 10 2>&1 > ${base_name}.log
    """
    stub:
    """
    touch ${base_name}.bracken
    touch ${base_name}.breport
    touch ${base_name}.log
    """

}

process BRACKEN_ONT { //TODO: check if read length average works
    tag "${sample}"
    publishDir "$params.outdir/bracken_reports/", mode: 'symlink'

    input: 
    tuple val(sample), val(base_name), path(reads), val (seq_type)
    tuple val(sample), val(base_name), path (fastp_json_report), val (seq_type)

    output:
    tuple val(sample), val(base_name), path("${base_name}.bracken"), val (seq_type), emit: bracken_out
    tuple val(sample), val(base_name), path("${base_name}.breport"), val (seq_type), emit: report_bracken_out
    tuple val(sample), val(base_name), path("${base_name}.log"), val (seq_type), emit: log_bracken

    script:
    """
    bracken -d ${params.kraken_db} -i ${base_name}.k2report -o ${base_name}.bracken -w ${base_name}.breport -r \$(cat ${fastp_json_report} | python -c "import sys, json; print(json.load(sys.stdin)['summary']['before_filtering']['read1_mean_length'])") -l S -t 10 2>&1 > ${base_name}.log
    """
    stub:
    """
    touch ${base_name}.bracken
    touch ${base_name}.breport
    touch ${base_name}.log
    """

}
