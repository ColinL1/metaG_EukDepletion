#!/usr/bin/env nextflow

process K2_SE {
    tag "${meta.id}"
    label "med_mem"
    conda "bioconda::kraken2"
    publishDir "${params.outdir}/mapping/unmapped-K2/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.out.gz"), emit: k2_out
    tuple val(meta), path("*.report.txt"), emit: report_k2_out

    script:
        """
    kraken2 \\
		--db /dev/shm/${params.kraken_dbName} \\
        --threads ${task.cpus} \\
        --memory-mapping \\
        --report ${meta.id}.report.txt \\
        --output ${meta.id}.out \\
        ${reads}
    pigz ${meta.id}.out
    """
    stub:
    """
    touch ${meta.id}.out.gz
    touch ${meta.id}.report.txt
    """
}

process K2_PE {
    tag "${meta.id}"
    label "med_mem"
    conda "bioconda::kraken2"
    publishDir "${params.outdir}/mapping/unmapped-K2/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.out.gz"), emit: k2_out
    tuple val(meta), path("*.report.txt"), emit: report_k2_out

    script:
    """
    kraken2 \\
        --db /dev/shm/${params.kraken_dbName} \\
        --paired \\
        --threads ${task.cpus} \\
        --memory-mapping \\
        --report ${meta.id}.report.txt \\
        --output ${meta.id}.out \\
        ${reads[0]} ${reads[1]}
    pigz ${meta.id}.out
    """
    stub:
    """
    touch ${meta.id}.out.gz
    touch ${meta.id}.report.txt
    """
}
