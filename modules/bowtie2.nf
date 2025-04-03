#!/usr/bin/env nextflow

process BOWTIE2_MAP {
    tag "${sample}"
    // cpus "${params.cpusHigh}"
    publishDir "$params.outdir/${seq_type}/bam/${ref_name}/", mode: 'symlink'
    publishDir "$params.outdir/${seq_type}/reads/${ref_name}/", mode: 'symlink'

    input: 
    // tuple val(sample_name), path(reads)
    tuple val(sample), val(base_name), path(reads), val (seq_type), val(ref_name), val (reference)
    // tuple val(ref_name), val (reference)

    output:
    tuple val(sample), val ("${base_name}.${ref_name}"), path ("${base_name}.${ref_name}.bam"), val (seq_type), emit: mapp_file
    tuple val(sample), val ("${base_name}.${ref_name}.unmapped"), path ("${base_name}.${ref_name}.unmapped_{1,2}.fq.gz"), val (seq_type), emit: unmapped_reads
    tuple val(sample), val ("${base_name}.${ref_name}.mapped"), path ("${base_name}.${ref_name}.mapped_{1,2}.fq.gz"), val (seq_type), emit: mapped_reads


    script:
    """
    bowtie2 -x ${reference} -1 ${reads[0]} -2 ${reads[1]} --un-conc-gz ${base_name}.${ref_name}.unmapped.fq.gz --al-conc-gz ${base_name}.${ref_name}.mapped.fq.gz -p ${task.cpus} | samtools view -S -b > ${base_name}.${ref_name}.bam
    ls *.1.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.1.gz/_1.fq.gz/'g) ; done
    ls *.2.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.2.gz/_2.fq.gz/'g) ; done
    """
    stub:
    """
    touch ${base_name}.${ref_name}.bam
    touch ${base_name}.${ref_name}.unmapped_1.fq.gz
    touch ${base_name}.${ref_name}.unmapped_2.fq.gz
    touch ${base_name}.${ref_name}.mapped_1.fq.gz
    touch ${base_name}.${ref_name}.mapped_2.fq.gz
    """
}

process BOWTIE2_MAP_VS {
    tag "${sample}"
    // cpus "${params.cpusHigh}"
    publishDir "$params.outdir/${seq_type}/bam/${ref_name}/", mode: 'symlink'
    publishDir "$params.outdir/${seq_type}/reads/${ref_name}/", mode: 'symlink'

    input: 
    // tuple val(base_name), path(reads)
    tuple val(sample), val(base_name), path(reads), val (seq_type), val(ref_name), val (reference)

    // tuple val(ref_name), val (reference)

    output:
    tuple val(sample), val ("${base_name}.${ref_name}"), path ("${base_name}.${ref_name}.bam"), val (seq_type), emit: mapp_file
    tuple val(sample), val ("${base_name}.${ref_name}.unmapped"), path ("${base_name}.${ref_name}.unmapped_{1,2}.fq.gz"), val (seq_type), emit: unmapped_reads
    tuple val(sample), val ("${base_name}.${ref_name}.mapped"), path ("${base_name}.${ref_name}.mapped_{1,2}.fq.gz"), val (seq_type), emit: mapped_reads

    script:
        """
        bowtie2 --very-sensitive -x ${reference} -1 ${reads[0]} -2 ${reads[1]} --un-conc-gz ${base_name}.${ref_name}.unmapped.fq.gz --al-conc-gz ${base_name}.${ref_name}.mapped.fq.gz -p ${task.cpus} | samtools view -S -b > ${base_name}.${ref_name}.bam
        """
}
