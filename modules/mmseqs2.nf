#!/usr/bin/env nextflow

process POOL_CONTIGS {
    tag "${species}"
    label 'min_mem'

    input:
    tuple val(species), path(fastas)

    output:
    tuple val(species), path("${species}.pooled.fa"), emit: pooled

    script:
    """
    cat ${fastas} > ${species}.pooled.fa
    """
    stub:
    """
    touch ${species}.pooled.fa
    """
}

process MMSEQS_LINCLUST {
    tag "${species}"
    label 'process_high'
    // easy-linclust can crash intermittently (Sls::error); a retry clears it
    errorStrategy 'retry'
    maxRetries 2
    conda "bioconda::mmseqs2=18.8cc5c"
    container 'quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0'
    publishDir path: { "${params.outdir}/coassembly_based/quantification/pooled/catalogue/" }, mode: 'copy'

    input:
    tuple val(species), path(pooled)

    output:
    tuple val(species), path("${species}.rep.fa"), path("${species}.clusters.tsv"), emit: catalogue

    script:
    if( params.derep )
        """
        mmseqs easy-linclust ${pooled} ${species} tmp \\
            --min-seq-id ${params.derep_id} -c ${params.derep_cov} \\
            --cov-mode 1 --threads ${task.cpus} >/dev/null
        mv ${species}_rep_seq.fasta ${species}.rep.fa
        mv ${species}_cluster.tsv   ${species}.clusters.tsv
        """
    else
        """
        cp ${pooled} ${species}.rep.fa
        awk '/^>/{n=substr(\$1,2); print n "\\t" n}' ${pooled} > ${species}.clusters.tsv
        """
    stub:
    """
    touch ${species}.rep.fa ${species}.clusters.tsv
    """
}
