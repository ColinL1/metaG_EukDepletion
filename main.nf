#!/usr/bin/env nextflow

include { TRIMMOMATIC } from './modules/trimmomatic.nf'
include { MAP_CONSE_PE } from './subworkflows/map_conse_pe.nf'
include { MAP_CONSE_ASSEMBLY } from './subworkflows/map_conse_assembly.nf'
include { CAT_WORKFLOW } from './subworkflows/CAT.nf'

workflow {
channel //TODO: swtich to sample sheet input for cleaner more flexible handling. Low priority. 
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
        [meta, reads_0, reads_1]
    }
    .set { grouped_samples }

        CAT_WORKFLOW(grouped_samples)

}