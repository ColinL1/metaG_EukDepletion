#!/usr/bin/env nextflow

process FASTP_REPORT {
    tag "${sample}"
    // cpus "${params.cpusMin}"
    maxRetries 5
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    publishDir "$params.outdir/${seq_type}/fastp/${base_name}", mode: 'symlink'

    input: 
    tuple val(sample), val(base_name), path(reads), val (seq_type)
    
    output:
    tuple val(sample), val(base_name), path ("${base_name}_fastp.json"), val (seq_type), emit: report_json
    tuple val(sample), val(base_name), path ("${base_name}_fastp.html"), val (seq_type), emit: report_html

    script:
    """
    fastp -i ${reads} --json ${base_name}_fastp.json --html ${base_name}_fastp.html -w 16 -Q -L -A &> /dev/null
    """
    stub:
    """
    touch ${base_name}_fastp.json
    touch ${base_name}_fastp.html
    """
}

process fastp_fasta_report { //TODO: add ${seq_type}
    tag "${sample_name}"
    // cpus "${params.cpusMin}"
    maxRetries 5
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    publishDir "$params.outdir/${seq_type}/fastp/${base_name}", mode: 'symlink'

    input: 
    tuple val(base_name), path(reads)

    output:
    tuple val(base_name), path ("${base_name}_fastp.json"), emit: report_json
    tuple val(base_name), path ("${base_name}_fastp.html"), emit: report_html

    script:
    """
    seqtk seq -F '#' ${reads} > ${base_name}.fq
    fastp -i ${base_name}.fq --json ${base_name}.json --html ${base_name}.html -w 16 -Q -L -A &> /dev/null
    """
    stub:
    """
    touch ${base_name}_fastp.json
    touch ${base_name}_fastp.html
    """
}