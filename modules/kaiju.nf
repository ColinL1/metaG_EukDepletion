#!/usr/bin/env nextflow

//params
params.nodes = '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/bacteria_kaiju/nodes.dmp'
params.names = '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/bacteria_kaiju/names.dmp'
params.kaiju_db = '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/bacteria_kaiju/kaiju_db_refseq.fmi'


process KAIJU_SE {
    tag "${meta.id}"
    label "process_high"
    publishDir "$baseDir/results/mapping/bacteria-kaiju/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.out"), emit: report_kaiju_out
    
    script:
    """
    kaiju -t ${params.nodes} -f ${params.kaiju_db} -i ${reads} -a mem -z ${task.cpus} -o ${meta.id}.out
    """
    stub:
    """
    touch ${meta.id}.out
    """
}

process KAIJU_PE {
    tag "${meta.id}"
    label "process_high"
    publishDir "$baseDir/results/mapping/bacteria-kaiju/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.out"), emit: report_kaiju_out

    script:
    """
    kaiju -t ${params.nodes} -f ${params.kaiju_db} -i ${reads[0]} -j ${reads[0]} -a mem -z ${task.cpus} -o ${meta.id}.out
    """
    stub:
    """
    touch ${meta.id}.out
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
    touch ${meta.id}_1.out
    touch ${meta.id}_2.out
    touch ${meta.id}_3.out
    touch ${meta.id}_4.out
    touch ${meta.id}_5.out
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
    touch ${meta.id}_1.out
    touch ${meta.id}_2.out
    touch ${meta.id}_3.out
    touch ${meta.id}_4.out
    touch ${meta.id}_5.out
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