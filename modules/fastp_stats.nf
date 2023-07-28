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

//     python /home/colinl/metaG/Git/metaG_EukDepletion/modules/parse_json.py -j /home/colinl/metaG/Git/metaG_EukDepletion/results_full_1/fastp -o PROVA_PARSER_V2_2.csv