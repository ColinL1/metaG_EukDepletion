#!/usr/bin/env nextflow

process KAIJU_SE {
    tag "${sample}"
    // cpus "${params.cpusHigh}"
    // memory "${params.memMax}"
    publishDir "$params.outdir/${seq_type}/Kaiju_reports/", mode: 'symlink'

    input: 
    tuple val(sample), val(base_name), path(reads), val(seq_type)


    output:
    tuple val(sample), val(base_name), path("${base_name}.out"), val(seq_type), emit: report_kaiju_out

    script:
    """
    kaiju -t ${params.nodes} -f ${params.kaiju_db} -i ${reads} -a mem -z ${task.cpus} -o ${base_name}.out
    """
    stub:
    """
    touch ${base_name}.out
    """
}

process KAIJU_PE {
    tag "${sample}"
    // cpus "${params.cpusHigh}"
    // memory "${params.memMax}"
    publishDir "$params.outdir/Kaiju_reports/", mode: 'symlink'

    input: 
    // tuple val(sample), path(reads)
    tuple val(sample), val(base_name), path(reads), val(seq_type)

    output:
    tuple val(sample), val(base_name), path("${base_name}.out"), val(seq_type), emit: report_kaiju_out

    script:
    """
    kaiju -t ${params.nodes} -f ${params.kaiju_db} -i ${reads[0]} -j ${reads[0]} -a mem -z ${task.cpus} -o ${base_name}.out
    """
    stub:
    """
    touch ${base_name}.out
    """
}


//all below is to fix and NOT USED currently 

process kaiju_multi { //to fix NOT USED AS IS 
    tag "${out_names}"
    label 'big_mem'
    // cpus "${params.cpusVHigh}"
    // memory "${params.memMax}"
    publishDir "$params.outdir/Kaiju/", mode: 'symlink'

    input: 
    tuple val(names), path(reads, stageAs: "?/*"), val(out_names)

    output:
    path '*_out', emit: kaiju_out

    script:
    """
    kaiju-multi -t ${params.nodes} -f ${params.kaiju_db} -i ${(reads as List).join(',')} -a mem -z ${task.cpus} -o ${(out_names as List).join(',')}
    """
    stub:
    """
    touch ${base_name}_1.out
    touch ${base_name}_2.out
    touch ${base_name}_3.out
    touch ${base_name}_4.out
    touch ${base_name}_5.out
    """
}

process kaiju_multi_pe { //to fix NOT USED AS IS 
    tag "${out_names}"
    label 'big_mem'
    // cpus 230
    // memory "${params.memMax}"
    publishDir "$params.outdir/Kaiju/", mode: 'symlink'

    input: 
    tuple val(names), path(reads_1), path(reads_2), val(out_names)

    output:
    path '*_out', emit: kaiju_out

    script:
    """
    kaiju-multi -t ${params.nodes} -f ${params.kaiju_db} -i ${(reads_1 as List).join(',')} -j ${(reads_2 as List).join(',')} -a mem -z ${task.cpus} -o ${(out_names as List).join(',')}
    """
    stub:
    """
    touch ${base_name}_1.out
    touch ${base_name}_2.out
    touch ${base_name}_3.out
    touch ${base_name}_4.out
    touch ${base_name}_5.out
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