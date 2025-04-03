#!/usr/bin/env nextflow
// Trimmomatic illumina reads trimming

process ALPHA_DIVERSITY {
    tag "${sample}"
    label "min_mem"
    publishDir "$params.outdir/diversity_reports/", mode: 'symlink'

    input: 
    tuple val(sample), val(base_name), path(bracken_report), val (seq_type)

    output:
    tuple val(sample), val(base_name), path("${base_name}_alpha_diversity.txt"), val (seq_type), emit: alpha_diversity_report

    script:
    """
    alpha_diversity.py -f ${bracken_report} -a Sh >> ${base_name}_alpha_diversity.txt
    alpha_diversity.py -f ${bracken_report} -a BP >> ${base_name}_alpha_diversity.txt
    alpha_diversity.py -f ${bracken_report} -a Si >> ${base_name}_alpha_diversity.txt
    alpha_diversity.py -f ${bracken_report} -a ISi >> ${base_name}_alpha_diversity.txt
    alpha_diversity.py -f ${bracken_report} -a F >> ${base_name}_alpha_diversity.txt
    """ //TODO: script so ${base_name}_alpha_diversity.txt becomes a single tsv table
    stub:
    """
    touch ${base_name}_alpha_diversity.txt
    """
}

process BETA_DIVERSITY {
    // tag "${sample}"
    // cpus "${params.cpusHigh}"
    // memory "${params.memMax}"
    label "min_mem"
    publishDir "$params.outdir/diversity_reports/", mode: 'symlink'

    input: 
    path (bracken_report)

    output:
    path("beta_diversity.txt"), emit: beta_diversity_report

    script:
    """
    beta_diversity.py -i ${(bracken_report as List).join(' ')} --type bracken > beta_diversity.txt
    """
    stub:
    """
    touch beta_diversity.txt
    """
}

// // individual tester #TODO:remove once complete
// workflow {
//     Channel.fromPath("/home/colinl/metaG/Git/kraken-protocol/*.bracken")
//         .map {tuple( it.name.split('.bracken')[0], it )}
//         .set { input }
//     // ("/home/colinl/metaG/Git/metaG_EukDepletion/input/test_illumina/*_{1,2}_subsample.fq.gz"). set {input_fq}
//     // input.view()
//     ALPHA_DIVERSITY(input)
//     // input.collect{it[1]}.view()
//     BETA_DIVERSITY(input.collect{it[1]})
//     ALPHA_DIVERSITY.out.alpha_diversity_report.view()
//     BETA_DIVERSITY.out.beta_diversity_report.view()
//     // .out.trimmed_fq_out.view()
// }