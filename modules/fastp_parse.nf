#!/usr/bin/env nextflow

process fastp_parse {
    // tag "${sample_name}"
    // cpus "${params.cpusMin}"
    // maxRetries 5
    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    publishDir "$params.outdir/fastp/", mode: 'symlink'

    input: 
    // tuple val(sample_name), path(reads)
    tuple path(reports)

    output:
    path ("report.csv"), emit: report_csv

    script:
    // template 'parse_json.py'
    '''
    python /home/colinl/metaG/Git/metaG_EukDepletion/parse_json_cl.py -j . -o report.csv
    '''

}

