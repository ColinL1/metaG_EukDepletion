#!/usr/bin/env nextflow

// TODO: check if worth adding SE option
process KRAKEN2_PE {
    tag "${meta.id}"
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' } 
    label 'med_mem'
    publishDir "$baseDir/results/kraken2_reports/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path ("*.tsv"), emit: report_kraken_out
    tuple val(meta), path ("*.k2.report"), emit: kraken_out
    tuple val(meta), path ("*.log"), emit: log_kraken

    script:
    """
    kraken2 --db ${params.kraken_db} ${params.mem_mapping} --threads ${task.cpus} --report ${meta.id}.k2.report --report-minimizer-data --minimum-hit-groups 3 --output ${meta.id}.tsv --paired ${reads[0]} ${reads[1]} 2> ${meta.id}.log
    """
    stub:
    """
    touch ${meta.id}.tsv
    touch ${meta.id}.k2.report
    touch ${meta.id}.log
    """
}
params.kraken_db = '/dev/shm/k2_nt'
params.mem_mapping = "--memory-mapping"

// process KRAKEN2_SE {
//     tag "${sample}"
//     // cpus "${params.cpusHigh}"
//     // memory "${params.memMax}"
//     publishDir "$params.outdir/kraken2_reports/", mode: 'symlink'

//     input: 
//     tuple val(sample), val(base_name), path(reads), val (seq_type)

//     output:
//     tuple val(sample), val(base_name), path("${base_name}.tsv"), val (seq_type), emit: report_kraken_out
//     tuple val(sample), val(base_name), path("${base_name}.k2report"), val (seq_type), emit: kraken_out
//     tuple val(sample), val(base_name), path("${base_name}.log"), val (seq_type), emit: log_kraken

//     script:
//     """
//     kraken2 --db ${params.kraken_db} --threads ${task.cpus} --report ${base_name}.k2report --report-minimizer-data --minimum-hit-groups 3 --output ${base_name}.tsv ${base_name} 2> ${base_name}.log
//     """
//     stub:
//     """
//     touch ${base_name}.tsv
//     touch ${base_name}.k2report
//     touch ${base_name}.log
//     """
// }

process TAXPASTA {
    tag "${meta.id}"
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' } 
    label 'min_mem'
    publishDir "$baseDir/reports/", mode: 'symlink'

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path ("*.tax_named.tsv"), emit: taxpasta_out

    script:
    """
    taxpasta standardise -p kraken2 -o ${meta.id}.tax.tsv ${report}
    python /home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/mapping/bacteria-kaiju/reads/kraken_exploration/add_taxa.py ${meta.id}.tax.tsv
    """
    stub:
    """
    touch ${meta.id}.tax_named.tsv
    """
}

params.input = "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/mapping/bacteria-kaiju/reads/kraken_exploration/input/*_{1,2}.non-bacteria.fq.gz"

Channel
    .fromFilePairs(params.input)
    .map { id, reads ->
        (species, replicate, method, buffer, other, unspecified) = id.tokenize("_")
        meta = [
            id:id,
            species:species,
            replicate:replicate,
            method:method,
            buffer:buffer,
            // other:other,
            // unspecified:unspecified,
        ]
        [meta, reads]
    }
    .map { meta, reads ->
        def newmap = [
            species: meta.species == "Po" ? "Porites" :
                    meta.species == "Por" ? "Porites" :
                    meta.species == "Ac" ? "Acropora" :
                    meta.species == "Acro" ? "Acropora" :
                    meta.species == "F003" ? "F003" :
                    meta.species == "F3" ? "F003" :
                    meta.species == "H2" ? "H2" :
                    meta.species == "Poci" ? "Pocillopora" :
                    meta.species == "Pr" ? "Pocillopora" :
                    "Unknown"
        ]
        [meta + newmap, reads]
    }
    .set { samples_coral }
    
workflow {
    // samples_coral.view()
    KRAKEN2_PE(samples_coral)
        TAXPASTA(KRAKEN2_PE.out.kraken_out)
}