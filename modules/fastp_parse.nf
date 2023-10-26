#!/usr/bin/env nextflow

process fastp_parse {
    // tag "${sample_name}"
    // cpus "${params.cpusMin}"
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    publishDir "$params.outdir/fastp/", mode: 'symlink'
    label 'min_mem'


    input: 
    // tuple val(sample_name), path(reads)
    tuple val(reports), path(reports)

    output:
    path ("report.csv"), emit: report_csv

    script:
    '''
    python /home/colinl/metaG/Git/metaG_EukDepletion/parse_json_cl.py -j . -o report.csv
    '''

}

process fastp_plot {
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    publishDir "$params.outdir/out/", mode: 'symlink'
    label 'min_mem'

    input: 
    tuple val(reports), path(reports)
    // path (sample_metadata_sheet)
    // path(parse_script)
    
    output:
    path ("report.csv"), emit: report_csv
    path ("*.png"), emit: barplot_png
    path ("*.svg"), emit: barplot_svg

    script:
    """
    python ${params.parse_script} -j . -o report.csv
    Rscript ${params.barplot_script} -f report.csv -m ${params.sample_metadata_sheet} 
    """

}