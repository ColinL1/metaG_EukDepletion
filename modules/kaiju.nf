#!/usr/bin/env nextflow

process KAIJU_SE {
    tag "${meta.id}"
    label "big_mem"
    conda "bioconda::kaiju"
    container 'quay.io/biocontainers/kaiju:1.10.1--h5ca1c30_3'
    publishDir path: { "${params.outdir}/assembly_based/mapping/bacteria-kaiju/${meta.id}/" }, mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.out"), emit: report_kaiju_out
    
    script:
    """
    kaiju -t ${params.nodes} -f ${params.kaiju_db} -i ${reads} -a mem -z ${task.cpus} -o ${meta.id}.out
    """
    stub:
    """
    touch ${meta.id}.out
    """
}

process KAIJU_PE {
    tag "${meta.id}"
    label "big_mem"
    conda "bioconda::kaiju"
    container 'quay.io/biocontainers/kaiju:1.10.1--h5ca1c30_3'
    publishDir path: { "${params.outdir}/reads_based/mapping/bacteria-kaiju/${meta.id}/" }, mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.out"), emit: report_kaiju_out

    script:
    """
    kaiju -t ${params.nodes} -f ${params.kaiju_db} -i ${reads[0]} -j ${reads[1]} -a mem -z ${task.cpus} -o ${meta.id}.out
    """
    stub:
    """
    touch ${meta.id}.out
    """
}

process KAIJU_ADD_TAXON_NAMES {
    tag "${meta.id}"
    label "med_mem"
    conda "bioconda::kaiju"
    container 'quay.io/biocontainers/kaiju:1.10.1--h5ca1c30_3'
    publishDir path: { "${params.outdir}/${task.process.contains('MAP_CONSE_PE') ? 'reads_based' : 'assembly_based'}/mapping/bacteria-kaiju/${meta.id}/" }, mode: 'symlink'

    input:
    tuple val(meta), path(kaiju_out)

    output:
    tuple val(meta), path("${meta.id}.taxon.out"), emit: taxon

    script:
    """
    kaiju-addTaxonNames -t ${params.nodes} -n ${params.names} -i ${kaiju_out} -o ${meta.id}.taxon.out -r ${params.kaiju_ranks}
    """
    stub:
    """
    touch ${meta.id}.taxon.out
    """
}

process KAIJU_SUMMARISE {
    tag "${meta.id}"
    label "min_mem"
    publishDir path: { "${params.outdir}/${task.process.contains('MAP_CONSE_PE') ? 'reads_based' : 'assembly_based'}/mapping/bacteria-kaiju/${meta.id}/" }, mode: 'copy', pattern: '*.summary.tsv'

    input:
    tuple val(meta), path(taxon)

    output:
    tuple val(meta), path("${meta.id}.{bacteria,archaea,viruses,non-mapped}.ids"), emit: ids
    path("${meta.id}.kaiju.summary.tsv"), emit: summary

    script:
    """
    touch ${meta.id}.bacteria.ids ${meta.id}.archaea.ids ${meta.id}.viruses.ids ${meta.id}.non-mapped.ids
    awk -F'\\t' -v s="${meta.id}" '
        { total++ }
        \$1 == "C" && \$4 ~ /^Bacteria;/ { bacteria++; print \$2 > (s ".bacteria.ids"); next }
        \$1 == "C" && \$4 ~ /^Archaea;/  { archaea++;  print \$2 > (s ".archaea.ids");  next }
        \$1 == "C" && \$4 ~ /^Viruses;/  { viruses++;  print \$2 > (s ".viruses.ids");  next }
        { if (\$1 == "U") unclassified++; else root_level++
          print \$2 > (s ".non-mapped.ids") }
        END {
            print "sample\\ttotal\\tbacteria\\tarchaea\\tviruses\\tnon_mapped\\tnon_mapped_unclassified\\tnon_mapped_root_level"
            printf "%s\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\n", s, total, bacteria, archaea, viruses, unclassified + root_level, unclassified, root_level
        }' ${taxon} > ${meta.id}.kaiju.summary.tsv
    """
    stub:
    """
    touch ${meta.id}.bacteria.ids ${meta.id}.archaea.ids ${meta.id}.viruses.ids ${meta.id}.non-mapped.ids ${meta.id}.kaiju.summary.tsv
    """
}
