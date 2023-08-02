#!/usr/bin/env nextflow

process kaiju {
    tag "${sample}"
    // cpus "${params.cpusHigh}"
    // memory "${params.memMax}"
    publishDir "$params.outdir/Kaiju_reports/", mode: 'symlink'

    input: 
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.out"), emit: report_kaiju_out

    script:
    """
    kaiju -t ${params.nodes} -f ${params.kaiju_db} -i ${reads} -a mem -z ${task.cpus} -o ${sample}.out
    """
}

process kaiju_multi {
    tag "${out_names}"
    label 'big_mem'
    // cpus "${params.cpusVHigh}"
    // memory "${params.memMax}"
    publishDir "$params.outdir/Kaiju/", mode: 'symlink'

    input: 
    tuple val(names), path(reads), val(out_names)

    output:
    path '*_out', emit: kaiju_out

    script:
    """
    kaiju-multi -t ${params.nodes} -f ${params.kaiju_db} -i ${(reads as List).join(',')} -a mem -z ${task.cpus} -o ${(out_names as List).join(',')}
    """
}
    // kaiju-multi -t ${params.nodes} -f ${params.kaiju_db} -i ${(reads as List).join(',')} -a mem -z ${task.cpus} -o ${(out_names as List).join(',')}

            // kaiju-multi -t ${params.nodes} -f ${params.kaiju_db} -i ${reads} -a mem -z ${task.cpus} -o ${out_names}

// process tax_mmseqs_nr {
//     tag "${sample}"
//     cpus "${params.cpusHigh}"
//     memory "${params.memMax}"
//     publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

//     input: 
//     path(reads)

//     output:
//     path("ONT_reads_NR_lca.tsv"), emit: report_lca_tsv
//     path("ONT_reads_NR_report") , emit: report
//     path("ONT_reads_NR_tophit_aln") , emit: report_tophit_aln
//     path("ONT_reads_NR_tophit_report") , emit: report_tophit_report
//     tuple file("run_setting.txt"), file("run_error.txt"), emit: logs


//     script:
//     """
//     mmseqs easy-taxonomy ${reads} --lca-ranks species,genus,family,order,class,phylum,superkingdom --threads ${task.cpus} ${params.mmseqs2_db_nr} ONT_reads_NR ${params.tmp} > run_setting.txt 2> run_error.txt
//     """
// }