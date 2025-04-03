#!/usr/bin/env nextflow

process plot_tax_mmseqs_nt {
    tag "${lca_report}"
    cpus "${params.cpusMin}"
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input:
    path (lca_report)

    output:
    path("ONT_NT_plots"), emit: nt_plots_folder

    script:
    """
    Rscript ${params.tax_plot_LCA_R} -f ${lca_report} -o ONT_NT_plots -d NT
    rm Rplots.pdf
    """
}

process plot_tax_mmseqs_nr {
    tag "${lca_report}"
    cpus "${params.cpusMin}"
    publishDir "$params.outdir/mmseqs2_reports/", mode: 'symlink'

    input:
    path (lca_report)

    output:
    path("ONT_NR_plots"), emit: nr_plots_folder

    script:
    """
    Rscript ${params.tax_plot_LCA_R} -f ${lca_report} -o ONT_NR_plots -d NR
    rm Rplots.pdf
    """
}
