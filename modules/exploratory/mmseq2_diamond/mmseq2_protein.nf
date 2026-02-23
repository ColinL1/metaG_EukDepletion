#!/usr/bin/env nextflow

process aa_eggNOG_mmseqs {
    tag "${sample}"
    cpus "${params.cpusHigh}"
    memory "${params.memMax}"
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input:
    path (sample)

    output:
    path("ONT_reads_eggNOG"), emit: eggNOG_search_out
    tuple file("run_setting.txt"), file("run_error.txt"), emit: logs


    script:
    """
    mmseqs easy-search ${sample} --threads ${task.cpus}  ${params.mmseqs2_db_eggNOG} ONT_reads_eggNOG ${params.tmp} > run_setting.txt 2> run_error.txt
    """
}

process pfam_A_mmseqs {
    tag "${sample}"
    cpus "${params.cpusHigh}"
    memory "${params.memMax}"
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input:
    path (sample)

    output:
    path("ONT_reads_pfam_A"), emit: pfam_A_search_out
    tuple file("run_setting.txt"), file("run_error.txt"), emit: logs


    script:
    """
    mmseqs easy-search ${sample} --threads ${task.cpus} ${params.mmseqs2_db_pfamA} ONT_reads_pfam_A ${params.tmp} > run_setting.txt 2> run_error.txt
    """
}

process pfam_B_mmseqs {
    tag "${sample}"
    cpus "${params.cpusHigh}"
    memory "${params.memMax}"
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input:
    path (sample)


    output:
    path("ONT_reads_pfam_B"), emit: pfam_B_search_out
    tuple file("run_setting.txt"), file("run_error.txt"), emit: logs


    script:
    """
    mmseqs easy-search ${sample} --threads ${task.cpus} ${params.mmseqs2_db_pfamB} ONT_reads_pfam_B ${params.tmp} > run_setting.txt 2> run_error.txt
    """
}

process swiss_prot_mmseqs {
    tag "${sample}"
    cpus "${params.cpusHigh}"
    memory "${params.memMax}"
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input:
    path (sample)

    output:
    path("ONT_reads_Swiss_prot"), emit: swiss_prot_search_out
    tuple file("run_setting.txt"), file("run_error.txt"), emit: logs


    script:
    """
    mmseqs easy-search ${sample} --threads ${task.cpus} ${params.mmseqs2_db_SwProt} ONT_reads_Swiss_prot ${params.tmp} > run_setting.txt 2> run_error.txt
    """
}

process trEMBL_mmseqs {
    tag "${sample}"
    cpus "${params.cpusHigh}"
    memory "${params.memMax}"
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input:
    path (sample)


    output:
    path("ONT_reads_TrEMBL"), emit: trEMBL_search_out
    tuple file("run_setting.txt"), file("run_error.txt"), emit: logs


    script:
    """
    mmseqs easy-search ${sample} --threads ${task.cpus} ${params.mmseqs2_db_TrEMBL} ONT_reads_TrEMBL ${params.tmp} > run_setting.txt 2> run_error.txt
    """
}