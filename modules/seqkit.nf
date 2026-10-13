#!/usr/bin/env nextflow

// Prefix contig names with the sample id so they stay unique once pooled
process SEQKIT_PREFIX_CONTIGS {
    tag "${meta.id}"
    label 'min_mem'
    conda "bioconda::seqkit=2.8.2"
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_1'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}.renamed.fa"), emit: renamed

    script:
    """
    seqkit seq -i ${contigs} \\
      | awk -v s="${meta.id}" '/^>/{print ">" s "__" substr(\$1,2); next}1' \\
      > ${meta.id}.renamed.fa
    """
    stub:
    """
    touch ${meta.id}.renamed.fa
    """
}

process SEQKIT_SPLIT_KAIJU_PE {
    tag "${meta.id}"
    label 'min_mem'
    conda "bioconda::seqkit=2.8.2"
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_1'
    publishDir path: { "${params.outdir}/reads_based/mapping/bacteria-kaiju/${meta.id}/" }, mode: 'symlink'

    input:
    tuple val(meta), path(reads), path(ids)

    output:
    tuple val(meta), path("${meta.id}.bacteria.mapped_{1,2}.fq.gz"), emit: bacteria
    tuple val(meta), path("${meta.id}.archaea.mapped_{1,2}.fq.gz"), emit: archaea
    tuple val(meta), path("${meta.id}.viruses.mapped_{1,2}.fq.gz"), emit: viruses
    tuple val(meta), path("${meta.id}.non-mapped.unmapped_{1,2}.fq.gz"), emit: non_mapped

    script:
    """
    for group in bacteria archaea viruses; do
        seqkit grep -j ${task.cpus} -f ${meta.id}.\$group.ids ${reads[0]} -o ${meta.id}.\$group.mapped_1.fq.gz
        seqkit grep -j ${task.cpus} -f ${meta.id}.\$group.ids ${reads[1]} -o ${meta.id}.\$group.mapped_2.fq.gz
    done
    seqkit grep -j ${task.cpus} -f ${meta.id}.non-mapped.ids ${reads[0]} -o ${meta.id}.non-mapped.unmapped_1.fq.gz
    seqkit grep -j ${task.cpus} -f ${meta.id}.non-mapped.ids ${reads[1]} -o ${meta.id}.non-mapped.unmapped_2.fq.gz
    """
    stub:
    """
    for group in bacteria archaea viruses; do
        touch ${meta.id}.\$group.mapped_1.fq.gz ${meta.id}.\$group.mapped_2.fq.gz
    done
    touch ${meta.id}.non-mapped.unmapped_1.fq.gz ${meta.id}.non-mapped.unmapped_2.fq.gz
    """
}

process SEQKIT_SPLIT_KAIJU_SE {
    tag "${meta.id}"
    label 'min_mem'
    conda "bioconda::seqkit=2.8.2"
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_1'
    publishDir path: { "${params.outdir}/assembly_based/mapping/bacteria-kaiju/${meta.id}/" }, mode: 'symlink'

    input:
    tuple val(meta), path(contigs), path(ids)

    output:
    tuple val(meta), path("${meta.id}.bacteria.fa.gz"), emit: bacteria
    tuple val(meta), path("${meta.id}.archaea.fa.gz"), emit: archaea
    tuple val(meta), path("${meta.id}.viruses.fa.gz"), emit: viruses
    tuple val(meta), path("${meta.id}.non-mapped.fa.gz"), emit: non_mapped

    script:
    """
    for group in bacteria archaea viruses non-mapped; do
        seqkit grep -j ${task.cpus} -f ${meta.id}.\$group.ids ${contigs} -o ${meta.id}.\$group.fa.gz
    done
    """
    stub:
    """
    for group in bacteria archaea viruses non-mapped; do
        touch ${meta.id}.\$group.fa.gz
    done
    """
}

process SEQKIT_STATS {
    tag "${meta.id}"
    label 'min_mem'
    conda "bioconda::seqkit=2.8.2"
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_1'

    input:
    tuple val(meta), path(files)

    output:
    path("${meta.id}.stats.tsv"), emit: stats

    script:
    """
    seqkit stats -T -j ${task.cpus} ${files} \\
      | awk -v s="${meta.id}" 'BEGIN{OFS="\\t"} NR==1{print "sample", \$0; next} {print s, \$0}' \\
      > ${meta.id}.stats.tsv
    """
    stub:
    """
    printf 'sample\\tfile\\tformat\\ttype\\tnum_seqs\\tsum_len\\tmin_len\\tavg_len\\tmax_len\\n' > ${meta.id}.stats.tsv
    for f in ${files}; do printf '${meta.id}\\t%s\\tNA\\tNA\\t0\\t0\\t0\\t0\\t0\\n' "\$f" >> ${meta.id}.stats.tsv; done
    """
}
