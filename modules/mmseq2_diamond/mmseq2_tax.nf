#!/usr/bin/env nextflow

process tax_mmseqs_nt {
    tag "${sample}"
    cpus "${params.cpusHigh}"
    memory "${params.memMax}"
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input: 
    path(reads)

    output:
    path("ONT_reads_NT_lca.tsv"), emit: report_lca_tsv
    path("ONT_reads_NT_report"), emit: report
    path("ONT_reads_NT_tophit_aln"), emit: report_tophit_aln
    path("ONT_reads_NT_tophit_report"), emit: report_tophit_report
    tuple file("run_setting.txt"), file("run_error.txt"), emit: logs

    script:
    """
    mmseqs easy-taxonomy ${reads} --lca-ranks species,genus,family,order,class,phylum,superkingdom --search-type 3 --threads ${task.cpus} ${params.mmseqs2_db_nt} ONT_reads_NT ${params.tmp} > run_setting.txt 2> run_error.txt
    """
}

process tax_mmseqs_nr {
    tag "${sample}"
    cpus "${params.cpusHigh}"
    memory "${params.memMax}"
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input: 
    path(reads)

    output:
    path("ONT_reads_NR_lca.tsv"), emit: report_lca_tsv
    path("ONT_reads_NR_report") , emit: report
    path("ONT_reads_NR_tophit_aln") , emit: report_tophit_aln
    path("ONT_reads_NR_tophit_report") , emit: report_tophit_report
    tuple file("run_setting.txt"), file("run_error.txt"), emit: logs


    script:
    """
    mmseqs easy-taxonomy ${reads} --lca-ranks species,genus,family,order,class,phylum,superkingdom --threads ${task.cpus} ${params.mmseqs2_db_nr} ONT_reads_NR ${params.tmp} > run_setting.txt 2> run_error.txt
    """
}