#!/usr/bin/env nextflow
// fastqc on samples

process fastqc {
    tag "${sample}"
    label 'min_mem'
    // cpus "${params.cpusHigh}"
    // memory "${params.memMax}"
    publishDir "$params.outdir/QC/fastqc", mode: 'symlink'

    input: 
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*.zip"), emit: zip //["sample name", [path to zip 1.zip, path to zip 2.zip]]
    tuple val(sample), path("*.html"), emit:html // as above 

    script:
    """
    fastqc -q -t ${task.cpus} ${reads[0]} ${reads[1]} 
    """
}

// // individual tester #TODO:remove once complete
// workflow {
//     Channel.fromFilePairs("/home/colinl/metaG/Git/metaG_EukDepletion/input/test_illumina/*_{1,2}_subsample.fq.gz"). set {input_fq}
//     // input_fq.view()
//     fastqc(input_fq)
//     fastqc.out.zip.view()
// }