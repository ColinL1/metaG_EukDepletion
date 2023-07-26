#!/usr/bin/env nextflow

process fastp_report {
    tag "${sample_name}"
    cpus "${params.cpusMin}"
    publishDir "$params.outdir/fastp/", mode: 'symlink'

    input: 
    tuple val(sample_name), path(reads)

    output:
    path ("${sample_name}.json"), emit: report_json
    path ("${sample_name}.html"), emit: report_html

    script:
    """
    fastp -i ${reads} --json ${sample_name}.json --html ${sample_name}.html -w 16 -Q -L -A &> /dev/null 
    """
}

