// process concatenateFastq {
//     tag "${meta.species}_${meta.method}_${meta.buffer}"

//     input:
//     tuple val(meta), path(reads)

//     output:
//     tuple val(meta), path("*.fq.gz")

//     script:
//     def outputR1 = "${meta.species}_${meta.method}_${meta.buffer}_R1.fq.gz"
//     def outputR2 = "${meta.species}_${meta.method}_${meta.buffer}_R2.fq.gz"

//     """
//     # Extract R1 and R2 reads
//     cat ${reads[0]} >> ${outputR1}
//     cat ${reads[1]} >> ${outputR2}

//     # Ensure output files are created
//     ls ${outputR1} ${outputR2}
//     """
//     stub:
//         script:
//     """
//     # Stub: Simulate concatenation
//     touch ${meta.species}_${meta.method}_${meta.buffer}_R1.fq.gz
//     touch ${meta.species}_${meta.method}_${meta.buffer}_R2.fq.gz
//     """
// }


// params.input_fq = "/home/colinl/metaG/Git/metaG_EukDepletion/input/test/*_R{1,2}.fq.gz"

// Channel
//     .fromFilePairs(params.input_fq)
//     .map { id, reads ->
//         (species, replicate, method, buffer) = id.tokenize("_") //, other, unspecified
//         meta = [
//             coassembly_grouping: "${species}_${method}_${buffer}",
//             id:id,
//             species:species,
//             replicate:replicate,
//             method:method,
//             buffer:buffer,
//             // other:other,
//             // unspecified:unspecified,
//         ]
//         [meta, reads]
//     }
//     .flatMap{ meta -> [ ] } // Collect all tuples into a list
//     .groupTuple() // Group by coassembly_grouping from the meta part of the tuple
//     // .map { it, reads ->
//     //     def group_key = coassembly_grouping
//     //     def group_reads = reads.collect { it[1] } // Collect the reads from the grouped tuples
//     //     def meta = [
//     //         id: coassembly_grouping,
//     //         species: group_key.tokenize("_")[0],
//     //         method: group_key.tokenize("_")[1],
//     //         buffer: group_key.tokenize("_")[2]
//     //     ]
//     //     [meta, group_reads.flatten()]
//     // }
//     .set { grouped_reads } // Create a new channel for grouped reads


// Channel
//     .fromFilePairs(params.input_fq)
//     .map { id, reads ->
//         (species, replicate, method, buffer) = id.tokenize("_") //, other, unspecified
//         meta = [
//             id:id,
//             species:species,
//             replicate:replicate,
//             method:method,
//             buffer:buffer,
//             // other:other,
//             // unspecified:unspecified,
//         ]
//         [meta, reads]
//     }
//     .map { meta, reads ->
//         def newmap = [
//             species: meta.species == "Po" ? "Porites" :
//                     meta.species == "Por" ? "Porites" :
//                     meta.species == "Ac" ? "Acropora" :
//                     meta.species == "Acro" ? "Acropora" :
//                     meta.species == "F003" ? "F003" :
//                     meta.species == "F3" ? "F003" :
//                     meta.species == "H2" ? "H2" :
//                     meta.species == "Poci" ? "Pocillopora" :
//                     meta.species == "Pr" ? "Pocillopora" :
//                     "Unknown"
//         ]
//         [meta + newmap, reads]
//     }
//     .set { samples_coral }

// Group the files by species, method, and buffer
// samples_coral
//     .map { meta, reads ->
//         def key = [meta.species, meta.method, meta.buffer ]
//         tuple( key, reads[0], reads[1] )
//     }
//     .groupTuple()
//     .map { key, reads_0, reads_1 ->
//         def meta = [
//             id: key.join('_'),
//             species: key[0],
//             method: key[1],
//             buffer: key[2]
//         ]
//         def reads = [
//             reads_0.collect { it },
//             reads_1.collect { it }
//         ]
//         [meta, reads]
//     }
//     .set { grouped_samples }

// Process to concatenate the fastq files
process CONCAT_FASTQ {
    tag "$meta.id"
    publishDir "$baseDir/results/concatenated_fastq/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), val(reads)

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.fastq.gz"), emit: concat_reads
    
    script:
    """
    cat ${reads[0].join(' ').replaceAll('[\\[\\],]', '')} > ${meta.id}_R1.fastq.gz
    cat ${reads[1].join(' ').replaceAll('[\\[\\],]', '')} > ${meta.id}_R2.fastq.gz
    """
    
    stub:
    
    """
    touch ${meta.id}_R1.fastq.gz
    touch ${meta.id}_R2.fastq.gz
    """
}

// workflow {
//     // grouped_samples.view()
//     // CONCAT_FASTQ(grouped_samples) // Use the new grouped_reads channel
// }

