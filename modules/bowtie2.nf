#!/usr/bin/env nextflow

process map2ref_pe {
    tag "${sample_name}"
    // cpus "${params.cpusHigh}"
    publishDir "$params.outdir/bam/${ref_name}/", mode: 'symlink'
    publishDir "$params.outdir/reads/${ref_name}/", mode: 'symlink'

    input: 
    tuple val(sample_name), path(reads)
    tuple val(ref_name), val (reference)

    output:
    tuple val ("${sample_name}.${ref_name}"), path ("${sample_name}.${ref_name}.bam"), emit: mapp_file
    tuple val ("${sample_name}.${ref_name}.unmapped"), path ("${sample_name}.${ref_name}.unmapped.fq.{1,2}.gz"), emit: unmapped_reads
    tuple val ("${sample_name}.${ref_name}.mapped"), path ("${sample_name}.${ref_name}.mapped.fq.{1,2}.gz"), emit: mapped_reads


    script:
        """
        bowtie2 -x ${reference} -1 ${reads[0]} -2 ${reads[1]} --un-conc-gz ${sample_name}.${ref_name}.unmapped.fq.gz --al-conc-gz ${sample_name}.${ref_name}.mapped.fq.gz -p ${task.cpus} | samtools view -S -b > ${sample_name}.${ref_name}.bam
        """
}
