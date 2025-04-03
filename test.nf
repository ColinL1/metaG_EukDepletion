#!/usr/bin/env nextflow

include { MAP_CONSE_PE } from './subworkflows/map_conse_pe.nf'
include { MAP_CONSE_ASSEMBLY } from './subworkflows/map_conse_assembly.nf'
include { CAT } from './subworkflows/CAT.nf'

params.input_fq = "/home/colinl/metaG/Git/metaG_EukDepletion/input/*_{1,2}.fq.gz"

Channel
    .fromFilePairs(params.input_fq)
    .map { id, reads ->
        (species, replicate, method, buffer) = id.tokenize("_") //, other, unspecified
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



// params.input = "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/assembly/assembly/*/*/final.contigs.fa"

// Channel
//     .fromPath(params.input)
//     .map { reads ->
//         (species, replicate, method, buffer, other, unspecified) = reads.getParent().getName().tokenize("_")
//         meta = [
//             id:reads.getParent().getName().minus(~/_TRIM_paried/),
//             // single_end:'PE',
//             species:species,
//             replicate:replicate,
//             method:method,
//             buffer:buffer,
//             other:other,
//             unspecified:unspecified,
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
//     .set { samples_coral_ass }
samples_coral
.map { meta, reads ->
    def key = [meta.species, meta.method, meta.buffer ]
    tuple( key, reads[0], reads[1] )
}
.groupTuple()
.map { key, reads_0, reads_1 ->
    def meta = [
        id: key.join('_'),
        species: key[0],
        method: key[1],
        buffer: key[2]
    ]
    def reads = [
        reads_0.collect { it },
        reads_1.collect { it }
    ]
    [meta, reads]
}
.set { grouped_samples }

workflow {
    MAP_CONSE_PE(samples_coral)
    MAP_CONSE_ASSEMBLY(samples_coral)
    CAT(grouped_samples)
    // samples_coral.view()
}