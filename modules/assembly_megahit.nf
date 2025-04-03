#!/usr/bin/env nextflow

process MEGAHIT_PE {
    tag "${meta.id}"
    //label "process_high"

    publishDir "$baseDir/results/assembly/${meta.id}", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path ("${meta.id}/final.contigs.fa"), emit: contigs_fa
    tuple val(meta), path("${meta.id}/log"), emit: log


    script:
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} --out-dir ${meta.id} --k-min 27 --k-max 127 --k-step 10 --num-cpu-threads ${task.cpus}
    """
    stub:
    """
    mkdir ${meta.id}
    touch ${meta.id}/final.contigs.fa
    touch ${meta.id}/log
    """
}

// process MEGAHIT_SE {
//     tag "${meta.id}"
//     publishDir "$baseDir/${seq_type}_contigs", mode: 'symlink'

//     input:
//     tuple val(sample), val(base_name), path(reads), val (seq_type)

//     output:
//     tuple val(sample), val("${base_name}.contigs"), path ("${base_name}/final.contigs.fa"), val("contigs_ont"), emit: contigs_fa

//     script:
//     """
//     megahit -r ${reads}  --out-dir ${base_name} --k-min 27 --k-max 127 --k-step 10 --num-cpu-threads ${task.cpus}
//     """
//     stub:
//     """
//     mkdir ${base_name}
//     touch ${base_name}/final.contigs.fa
//     """
// }

workflow {
    Channel.fromFilePairs("/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/assembly/input/*_{1,2}.fq.gz")
        | map { id, reads ->
            (species, replicate, method, buffer, other, unspecified) = id.tokenize("_")
            meta = [
                id:id,
                // single_end:'PE',
                species:species,
                replicate:replicate,
                method:method,
                buffer:buffer,
                other:other,
                unspecified:unspecified,
            ]
            [meta, reads]
    }  // | view()
    | MEGAHIT_PE
    // MEGAHIT_PE(input_ch_pe)
}