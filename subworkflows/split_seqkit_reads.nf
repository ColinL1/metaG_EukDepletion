#!/usr/bin/env nextflow

process SEQKIT_SPLIT_BAC {
    tag "${sample}"
    // cpus "${params.cpusMin}"
    // errorStrategy 'ignore'
    publishDir "$params.outdir/results/${seq_type}/reads/bacteria/${sample}", mode: 'symlink'

    input: 
    tuple val(sample), val(base_name), path(kaiju_out), path(reads), val (seq_type)
    // tuple val(sample), , path(reads)


    output:
    tuple val(sample), val ("${base_name}.bacteria"), path ("${base_name}.bacteria.f*.gz"), val (seq_type), emit: bacteria_reads
    tuple val(sample), val ("${base_name}.non-bacteria"), path ("${base_name}.non-bacteria.f*.gz"), val (seq_type), emit: non_bacteria_reads


    script:
    if ( params.mode == 'contigs' )
    """
    grep  -w "C" ${kaiju_out} | cut -d\$'\t' -f2 > list_bacteria.txt
    grep  -w "U" ${kaiju_out} | cut -d\$'\t' -f2 > list_not-bacteria.txt
    seqkit grep -j ${task.cpus} -f list_bacteria.txt -i ${reads} > ${sample}.bacteria.fa
    seqkit grep -j ${task.cpus} -f list_not-bacteria.txt -i ${reads} > ${sample}.non-bacteria.fa
    pigz -p ${task.cpus} ${sample}.non-bacteria.fa ${sample}.bacteria.fa
    rm -rf list_bacteria.txt list_not-bacteria.txt
    """
    else
    """
    grep -w "C" ${kaiju_out} | cut -d\$'\t' -f2 > list_bacteria.txt
    grep -w "U" ${kaiju_out} | cut -d\$'\t' -f2 > list_not-bacteria.txt
    seqkit grep -j ${task.cpus} -f list_bacteria.txt -i ${reads} > ${base_name}.bacteria.fq
    seqkit grep -j ${task.cpus} -f list_not-bacteria.txt -i ${reads} > ${base_name}.non-bacteria.fq
    pigz -p ${task.cpus} ${base_name}.bacteria.fq ${base_name}.non-bacteria.fq
    rm -rf list_bacteria.txt list_not-bacteria.txt
    """
    stub: 
    if ( params.mode == 'contigs' )
    """
    touch ${base_name}.bacteria.fa.gz
    touch ${base_name}.non-bacteria.fa.gz
    """
    else
    """
    touch ${base_name}.bacteria.fq.gz
    touch ${base_name}.non-bacteria.fq.gz
    """
}

process SEQKIT_SPLIT_BAC_PE {
    tag "${sample}"
    // cpus "${params.cpusMin}"
    publishDir "$params.outdir/results/${seq_type}/reads/bacteria/${sample}", mode: 'symlink'

    input: 
    // tuple val(sample), path(kaiju_out), path(reads)
    tuple val(sample), val(base_name), path(kaiju_out), path(reads), val (seq_type)

    output:
    tuple val(sample), val ("${base_name}.bacteria"), path ("${base_name}_PE_*.bacteria.fq.gz"), val (seq_type), emit: bacteria_reads
    tuple val(sample), val ("${base_name}.non-bacteria"), path ("${base_name}_PE_*.non-bacteria.fq.gz"), val (seq_type), emit: non_bacteria_reads


    script:
    """
    grep -w "C" ${kaiju_out} | cut -d\$'\t' -f2 > list_bacteria.txt
    grep -w "U" ${kaiju_out} | cut -d\$'\t' -f2 > list_not-bacteria.txt
    seqkit grep -j ${task.cpus} -f list_bacteria.txt -i ${reads[0]} > ${base_name}_PE_1.bacteria.fq
    seqkit grep -j ${task.cpus} -f list_bacteria.txt -i ${reads[1]} > ${base_name}_PE_2.bacteria.fq
    seqkit grep -j ${task.cpus} -f list_not-bacteria.txt -i ${reads[0]} > ${base_name}_PE_1.non-bacteria.fq
    seqkit grep -j ${task.cpus} -f list_not-bacteria.txt -i ${reads[1]} > ${base_name}_PE_2.non-bacteria.fq
    pigz -p ${task.cpus} *.fq
    rm -rf list_bacteria.txt list_not-bacteria.txt
    """
    stub: 
    """
    touch ${base_name}_PE_1.bacteria.fq.gz
    touch ${base_name}_PE_2.bacteria.fq.gz
    touch ${base_name}_PE_1.non-bacteria.fq.gz
    touch ${base_name}_PE_2.non-bacteria.fq.gz
    """
}

// fixed proces start here.

process SEQKIT_SPLIT_KAIJU_PE {
    tag "${meta.id}"
    label "min_mem"
    publishDir "$baseDir/results/mapping/bacteria-kaiju/reads/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(kaiju_out)

    output:
    tuple val(meta), path ("${meta.id}_{1,2}.bacteria.fq.gz"), emit: bacteria_reads
    tuple val(meta), path ("${meta.id}_{1,2}.non-bacteria.fq.gz"), emit: non_bacteria_reads

    script:
    """
    grep  -w "C" ${kaiju_out} | cut -d\$'\t' -f2 > list_bacteria.txt
    grep  -w "U" ${kaiju_out} | cut -d\$'\t' -f2 > list_not-bacteria.txt
    seqkit grep -j ${task.cpus} -f list_bacteria.txt -i ${reads[0]} > ${meta.id}_1.bacteria.fq
    seqkit grep -j ${task.cpus} -f list_bacteria.txt -i ${reads[1]} > ${meta.id}_2.bacteria.fq
    seqkit grep -j ${task.cpus} -f list_not-bacteria.txt -i ${reads[0]} > ${meta.id}_1.non-bacteria.fq
    seqkit grep -j ${task.cpus} -f list_not-bacteria.txt -i ${reads[1]} > ${meta.id}_2.non-bacteria.fq
    pigz -p ${task.cpus} ${meta.id}*.fq
    rm -rf list_bacteria.txt list_not-bacteria.txt
    """
    stub: 
    """
    touch ${meta.id}_1.bacteria.fq.gz
    touch ${meta.id}_2.bacteria.fq.gz
    touch ${meta.id}_1.non-bacteria.fq.gz
    touch ${meta.id}_2.non-bacteria.fq.gz
    """
}

params.input = "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/mapping/symbiodiniaceae/*/*.symbiodiniaceae.unmapped_{1,2}.fq.gz"

Channel
    .fromFilePairs(params.input)
    .map { id, reads ->
        (species, replicate, method, buffer, other, unspecified) = id.minus(~/_TRIM_paried/).tokenize("_")
        meta = [
            id:id.minus(~/_TRIM_paried/),
            species:species,
            replicate:replicate,
            method:method,
            buffer:buffer.minus(~/_TRIM_paried/),
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

params.input_kaiju = "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/mapping/bacteria-kaiju/*/*.out"

Channel
    .fromPath(params.input_kaiju)
    .map {reads ->
        (species, replicate, method, buffer) = reads.getParent().getName().tokenize("_")
        meta = [
            id:reads.getParent().getName().minus(~/_TRIM_paried/),
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
    .set { samples_kaiju_reports }

workflow {
    // samples_coral.view()
    // samples_kaiju_reports.view()
        SEQKIT_SPLIT_KAIJU_PE(samples_coral, samples_kaiju_reports)
}
