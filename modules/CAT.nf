#!/usr/bin/env nextflow

process CAT_NR {
    tag "${meta.id}"
    publishDir "$baseDir/results/CAT/NR/${meta.species}/${meta.id}/", mode: 'symlink'
    publishDir "$baseDir/results/concatenated_fastq/${meta.id}/", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path ("${meta.id}.contig2classification.txt"), emit: contig2classification
    tuple val(meta), path ("${meta.id}.ORF2LCA.txt"), emit: orf2lca
    tuple val(meta), path("*.log"), emit: log


    script:
    """
    /home/colinl/databases/CAT/CAT_pack-5.3/CAT_pack/CAT contigs -c ${reads} -d /home/colinl/databases/CAT/20231120_CAT_nr/db -t /home/colinl/databases/CAT/20231120_CAT_nr/tax/ -o ${meta.id} -n ${task.cpus} --compress
    """
    stub:
    """
    touch ${meta.id}.contig2classification.txt
    touch ${meta.id}.ORF2LCA.txt
    touch ${meta.id}.log
    """
}

process CAT_GTDB {
    tag "${meta.id}"
    publishDir "$baseDir/results/CAT/results_without_euk/GTDB/${meta.species}/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path ("${meta.id}.contig2classification.txt"), emit: contig2classification
    tuple val(meta), path ("${meta.id}.ORF2LCA.txt"), emit: orf2lca
    tuple val(meta), path("*.log"), emit: log
    
    script:
    """
    /home/colinl/databases/CAT/CAT_pack-5.3/CAT_pack/CAT contigs -c ${reads} -d /home/colinl/databases/CAT/20231120_CAT_gtdb/db -t /home/colinl/databases/CAT/20231120_CAT_gtdb/tax/ -o ${meta.id} -n ${task.cpus} --compress
    """
    stub:
    """
    touch ${meta.id}.contig2classification.txt
    touch ${meta.id}.ORF2LCA.txt
    touch ${meta.id}.log
    """
}

process CAT_ADD_NAMES_NR {
    tag "${meta.id}"
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    label 'min_mem'
    publishDir "$baseDir/results/CAT/NR/${meta.species}/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path ("*.contig2classification*.txt"), emit: contig2classification_taxname
    tuple val(meta), path ("*.ORF2LCA*.txt"), emit: orf2lca_taxname


    script:
    """
    /home/colinl/databases/CAT/CAT_pack-5.3/CAT_pack/CAT add_names -i  ${report} -o ${report}.tax_name.txt  -t /home/colinl/databases/CAT/20231120_CAT_nr/tax
    """
    stub:
    """
    touch ${meta.id}.contig2classification.tax_name.txt
    touch ${meta.id}.ORF2LCA.tax_name.txt
    """
}

process CAT_ADD_NAMES_GTDB {
    tag "${meta.id}"
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' } 
    label 'min_mem'
    publishDir "$baseDir/results/CAT/results_without_euk/GTDB/${meta.species}/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path ("*.contig2classification*.txt"), emit: contig2classification_taxname
    tuple val(meta), path ("*.ORF2LCA*.txt"), emit: orf2lca_taxname

    script:
    """
    /home/colinl/databases/CAT/CAT_pack-5.3/CAT_pack/CAT add_names -i ${report} -o ${report}.tax_name.txt  -t /home/colinl/databases/CAT/20231120_CAT_gtdb/tax
    """
    stub:
    """
    touch ${meta.id}.contig2classification.tax_name.txt
    touch ${meta.id}.ORF2LCA.tax_name.txt
    """
}