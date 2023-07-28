#!/usr/bin/env nextflow

process fastp_parse {
    tag "${sample_name}"
    // cpus "${params.cpusMin}"
    // maxRetries 5
    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    publishDir "$params.outdir/fastp/", mode: 'symlink'

    input: 
    tuple val(sample_name), path(reads)

    output:
    path ("report.csv"), emit: report_csv

    script:
    template 'parse_json.py'

}

