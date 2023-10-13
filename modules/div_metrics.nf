#!/usr/bin/env nextflow
// Trimmomatic illumina reads trimming

process alpha_diversity {
    tag "${sample}"
    // cpus "${params.cpusHigh}"
    // memory "${params.memMax}"
    label "min_mem"
    publishDir "$params.outdir/diversity_reports/", mode: 'symlink'

    input: 
    tuple val(sample), path(bracken_report)

    output:
    tuple val(sample), path("${sample}_alpha_diversity.txt"), emit: alpha_diversity_report

    script:
    """
    alpha_diversity.py -f ${bracken_report} -a Sh >> ${sample}_alpha_diversity.txt
    alpha_diversity.py -f ${bracken_report} -a BP >> ${sample}_alpha_diversity.txt
    alpha_diversity.py -f ${bracken_report} -a Si >> ${sample}_alpha_diversity.txt
    alpha_diversity.py -f ${bracken_report} -a ISi >> ${sample}_alpha_diversity.txt
    alpha_diversity.py -f ${bracken_report} -a F >> ${sample}_alpha_diversity.txt
    """
}

process beta_diversity {
    // tag "${sample}"
    // cpus "${params.cpusHigh}"
    // memory "${params.memMax}"
    label "min_mem"
    publishDir "$params.outdir/diversity_reports/", mode: 'symlink'

    input: 
    path(bracken_report)

    output:
    path("beta_diversity.txt"), emit: beta_diversity_report

    script:
    """
    beta_diversity.py -i ${(bracken_report as List).join(' ')} --type bracken > beta_diversity.txt
    """
}

// // individual tester #TODO:remove once complete
// workflow {
//     Channel.fromPath("/home/colinl/metaG/Git/kraken-protocol/*.bracken")
//         .map {tuple( it.name.split('.bracken')[0], it )}
//         .set { input }
//     // ("/home/colinl/metaG/Git/metaG_EukDepletion/input/test_illumina/*_{1,2}_subsample.fq.gz"). set {input_fq}
//     // input.view()
//     alpha_diversity(input)
//     // input.collect{it[1]}.view()
//     beta_diversity(input.collect{it[1]})
//     alpha_diversity.out.alpha_diversity_report.view()
//     beta_diversity.out.beta_diversity_report.view()
//     // .out.trimmed_fq_out.view()
// }