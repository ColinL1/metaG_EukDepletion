#!/usr/bin/env nextflow

process ont_trim {
    tag "${sample}"
    // cpus "${params.cpusMin}"
    publishDir "$params.outdir/porechop/${sample}", mode: 'symlink'
    
    input: 
    tuple val(sample), val(reads)

    output:
    tuple val ("${sample}"), path ("${sample}_trim.fastq.gz") , emit: trimmed_fastq
    path ("${sample}.log") , emit: log

    script:
    """
    porechop -i ${reads} -o ${sample}_trim.fastq.gz --threads ${task.cpus} --check_reads 200000 > ${sample}.log 
    """
}

process clean_names {
    tag "${sample}"
    cpus "${params.cpusMin}"
    publishDir "$params.outdir/fixed_names/", mode: 'symlink'

    input: 
    tuple val(sample), val(reads)

    output:
    tuple val ("${sample}"), path ("${sample}_clean_name.fastq.gz"), emit: clean_names_reads
    

    script:
    """
    bioawk -v FILE="${sample}" -c fastx '{ print \$name , FILE}' < ${reads} > ${sample}_new_keyfile.txt
    seqkit replace -k ${sample}_new_keyfile.txt -p '^(\\S+)' -r '{kv}_\$1' ${reads} -j ${task.cpus} -o ${sample}_clean_name-temp.fastq >/dev/null 2>/dev/null
    seqtk seq -l 0 ${sample}_clean_name-temp.fastq > ${sample}_clean_name.fastq
    rm ${sample}_clean_name-temp.fastq ${sample}_new_keyfile.txt
    pigz -p ${task.cpus} ${sample}_clean_name.fastq
    """
}