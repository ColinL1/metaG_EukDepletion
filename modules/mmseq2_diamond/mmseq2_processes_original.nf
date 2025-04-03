#!/usr/bin/env nextflow

process tax_mmseqs_nt {
    tag "${sample}"
    cpus $params.cpus.max
    memory $params.mem.max
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input:
    path (sample)

    output:
    path("ONT_reads_NT_lca.tsv"), emit: report_lca_tsv
    path("ONT_reads_NT_report"), emit: report
    path("ONT_reads_NT_tophit_aln"), emit: report_tophit_aln
    path("ONT_reads_NT_tophit_report"), emit: report_tophit_report
    tuple file("run_setting.txt"), file("run_error.txt"), emit: logs

    script:
    """
    mmseqs easy-taxonomy ${sample} --lca-ranks species,genus,family,order,class,phylum,superkingdom --search-type 3 --threads ${task.cpus} ${params.mmseqs2_db_nt} ONT_reads_NT ${params.tmp} > run_setting.txt 2> run_error.txt
    """
}

process tax_mmseqs_nr {
    tag "${sample}"
    cpus $params.cpus.max
    memory $params.mem.max
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input:
    path (sample) from nr.toList()

    output:
    path("ONT_reads_NR_lca.tsv"), emit: report_lca_tsv
    path("ONT_reads_NR_report") , emit: report_lca_tsv
    path("ONT_reads_NR_tophit_aln") , emit: report_lca_tsv
    path("ONT_reads_NR_tophit_report") , emit: report_lca_tsv
    tuple file("run_setting.txt"), file("run_error.txt"), emit: logs


    script:
    """
    mmseqs easy-taxonomy ${sample} --lca-ranks species,genus,family,order,class,phylum,superkingdom --threads ${task.cpus} ${params.mmseqs2_db_nr} ONT_reads_NR ${params.tmp} > run_setting.txt 2> run_error.txt
    """
}

process aa_eggNOG_mmseqs {
    tag "${sample}"
    cpus $params.cpus.max
    memory $params.mem.max
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input:
    path (sample) from eggNOG.toList()

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
    cpus $params.cpus.max
    memory $params.mem.max
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input:
    path (sample) from pfam_A.toList()

    output:
    path("ONT_reads_pfam_A"), emit: pfam_B_search_out
    tuple file("run_setting.txt"), file("run_error.txt"), emit: logs


    script:
    """
    mmseqs easy-search ${sample} --threads ${task.cpus} ${params.mmseqs2_db_pfamA} ONT_reads_pfam_A ${params.tmp} > run_setting.txt 2> run_error.txt
    """
}

process pfam_B_mmseqs {
    tag "${sample}"
    cpus $params.cpus.max
    memory $params.mem.max
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input:
    path (sample) from pfam_B.toList()


    output:
    path("ONT_reads_pfam_B"), emit: pfam_B_search_out
    tuple file("run_setting.txt"), file("run_error.txt"), emit: logs


    script:
    """
    mmseqs easy-search ${sample} --threads ${task.cpus} ${params.mmseqs2_db_pfamB} ONT_reads_pfam_B ${params.tmp} > run_setting.txt 2> run_error.txt
    """
}

process Swiss_Prot_mmseqs {
    tag "${sample}"
    cpus $params.cpus.max
    memory $params.mem.max
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input:
    path (sample).toList()

    output:
    path("ONT_reads_Swiss_prot"), emit: swiss_prot_search_out
    tuple file("run_setting.txt"), file("run_error.txt"), emit: logs


    script:
    """
    mmseqs easy-search ${sample} --threads ${task.cpus} ${params.mmseqs2_db_SwProt} ONT_reads_Swiss_prot ${params.tmp} > run_setting.txt 2> run_error.txt
    """
}

process TrEMBL_mmseqs {
    tag "${sample}"
    cpus $params.cpus.max
    memory $params.mem.max
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input:
    path (sample).toList()


    output:
    path("ONT_reads_TrEMBL"), emit: trEMBL_search_out
    tuple file("run_setting.txt"), file("run_error.txt"), emit: logs


    script:
    """
    mmseqs easy-search ${sample} --threads ${task.cpus} ${params.mmseqs2_db_TrEMBL} ONT_reads_TrEMBL ${params.tmp} > run_setting.txt 2> run_error.txt
    """
}