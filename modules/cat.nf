#!/usr/bin/env nextflow

params.cat_db = "/home/colinl/databases/CAT/20231120_CAT_gtdb"
params.cat_ex = "/home/colinl/databases/CAT/CAT_pack-5.3/CAT_pack/CAT"

process CAT {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/CAT/${meta.species}/${meta.id}/", mode: 'symlink'
    publishDir "${params.outdir}/concatenated_fastq/${meta.id}/", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path ("${meta.id}.contig2classification.txt"), emit: contig2classification
    tuple val(meta), path ("${meta.id}.ORF2LCA.txt"), emit: orf2lca
    tuple val(meta), path("*.log"), emit: log

    script:
    """
    ${params.cat_ex} contigs -c ${reads} -d ${params.cat_db}/db -t ${params.cat_db}/tax/ -o ${meta.id} -n ${task.cpus} --compress
    """
    stub:
    """
    touch ${meta.id}.contig2classification.txt
    touch ${meta.id}.ORF2LCA.txt
    touch ${meta.id}.log
    """
}

process CAT_ADD_NAMES {
    tag "${meta.id}"
    errorStrategy  { task.attempt <= 3 ? 'retry' : 'ignore' }
    label 'min_mem'
    publishDir "${params.outdir}/CAT/${meta.species}/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path ("*.contig2classification*.txt"), emit: contig2classification_taxname
    tuple val(meta), path ("*.ORF2LCA*.txt"), emit: orf2lca_taxname


    script:
    """
    ${params.cat_ex} add_names -i  ${report} -o ${report}.tax_name.txt  -t ${params.cat_db}/tax
    """
    stub:
    """
    touch ${meta.id}.contig2classification.tax_name.txt
    touch ${meta.id}.ORF2LCA.tax_name.txt
    """
}
