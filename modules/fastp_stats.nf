#!/usr/bin/env nextflow

process fastp_report {
    tag "${sample_name}"
    // cpus "${params.cpusMin}"
    maxRetries 5
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    publishDir "$params.outdir/fastp/${sample_name}", mode: 'symlink'

    input: 
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path ("${sample_name}.json"), emit: report_json
    tuple val(sample_name), path ("${sample_name}.html"), emit: report_html

    script:
    """
    fastp -i ${reads} --json ${sample_name}.json --html ${sample_name}.html -w 16 -Q -L -A &> /dev/null
    """
}

process fastp_fasta_report {
    tag "${sample_name}"
    // cpus "${params.cpusMin}"
    maxRetries 5
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    publishDir "$params.outdir/fastp/${sample_name}", mode: 'symlink'

    input: 
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path ("${sample_name}.json"), emit: report_json
    tuple val(sample_name), path ("${sample_name}.html"), emit: report_html

    script:
    """
    seqtk seq -F '#' ${reads} > ${sample_name}.fq
    fastp -i ${sample_name}.fq --json ${sample_name}.json --html ${sample_name}.html -w 16 -Q -L -A &> /dev/null
    """
}