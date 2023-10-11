#!/usr/bin/env nextflow
// Trimmomatic illumina reads trimming

process trimmomatic_pe {
    tag "${sample}"
    // cpus "${params.cpusHigh}"
    // memory "${params.memMax}"
    publishDir "$params.outdir/Trimmed_fastq/", mode: 'symlink'

    input: 
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_TRIM_paried*"), emit: trimmed_fq_out
    tuple val(sample), path("${sample}_TRIM_unparied*"), emit: singleton_fq_out
    tuple val(sample), path("${sample}_trimmomatic.log"), emit: log_trim

    script:
    """
    trimmomatic PE -threads ${task.cpus} ${reads[0]} ${reads[1]} ${sample}_TRIM_paried_R1.fq.gz ${sample}_TRIM_unparied_R1.fq.gz ${sample}_TRIM_paried_R2.fq.gz ${sample}_TRIM_unparied_R2.fq.gz ILLUMINACLIP:/home/colinl/metaG/Git/metaG_EukDepletion/modules/adapters/Truseq_V.3_edited.fa:2:30:10:8:TRUE HEADCROP:5 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 2>> ${sample}_trimmomatic.log 1>>file.out
    """
}

// individual tester #TODO:remove once complete
// workflow {
//     Channel.fromFilePairs("/home/colinl/metaG/Git/metaG_EukDepletion/input/test_illumina/*_{1,2}_subsample.fq.gz"). set {input_fq}
//     // input_fq.view()
//     trimmomatic_pe(input_fq)
//     trimmomatic_pe.out.trimmed_fq_out.view()
// }