#!/usr/bin/env nextflow
// fastqc on samples

process MULTIQC { //TODO: add ${seq_type} 
    tag "${qc}"
    label 'min_mem'

    // cpus "${params.cpusHigh}"
    // memory "${params.memMax}"
    publishDir "$params.outdir/${seq_type}/QC/multiqc", mode: 'symlink'
    errorStrategy 'ignore'

    input: 
    val(qc)

    output:
    path("multiqc_data"), emit: data //path to folder
    path("multiqc_report.html"), emit:html // path to html

    script:
    """
    multiqc --quiet --interactive ${qc}
    """
    stub:
    """
    mkdir multiqc_data
    touch multiqc_data/stub.txt
    touch multiqc_report.html
    """
}

// individual tester #TODO:remove once complete
// workflow {
//     // Channel.fromFilePairs("/home/colinl/metaG/Git/metaG_EukDepletion/input/test_illumina/*_{1,2}_subsample.fq.gz").collect(it[0]).set {input_fq}
//     Channel.fromPath("/home/colinl/metaG/Git/metaG_EukDepletion/modules/QC").set {input_fq}
//     // input_fq.view()
//     multiqc(input_fq)
//     multiqc.out.html.view()
//     multiqc.out.data.view()
// }