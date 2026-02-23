#!/usr/bin/env nextflow

include { TRIMMOMATIC } from './modules/trimmomatic.nf'
include { MAP_CONSE_PE } from './subworkflows/map_conse_pe.nf'
include { MAP_CONSE_ASSEMBLY } from './subworkflows/map_conse_assembly.nf'
include { CAT_WORKFLOW } from './subworkflows/CAT.nf'

workflow {
Channel
    .fromFilePairs(params.input_fq)
    .map { id, reads ->
        def (species, replicate, method, buffer) = id.tokenize("_") //, other, unspecified
        def meta = [
            id:id,
            species:species,
            replicate:replicate,
            method:method,
            buffer:buffer,
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

    

    TRIMMOMATIC(samples_coral)
        MAP_CONSE_PE(TRIMMOMATIC.out.trimmed_fq_out)
        MAP_CONSE_ASSEMBLY(TRIMMOMATIC.out.trimmed_fq_out)
        
        TRIMMOMATIC.out.trimmed_fq_out
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

        CAT_WORKFLOW(grouped_samples)

}