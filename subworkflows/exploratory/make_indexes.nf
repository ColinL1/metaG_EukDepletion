#!/usr/bin/env nextflow

/*
========================================================================================
    Run kaiju and extract bacteria matching reads
========================================================================================
*/

params.nodes = "/share/databases/kaiju/nodes.dmp"
params.kaiju_db = "/share/databases/kaiju/refseq/kaiju_db_refseq.fmi"

/*
========================================================================================
    Include Modules
========================================================================================
*/

include { NCBI_DOWNLOAD } from '../modules/NCBI_Datasets.nf' 
include { MINIMAP2_INDEX } from '../modules/minimap2_index.nf' 

/*
========================================================================================
    Workflow EXTRACT_BACTERIA
========================================================================================
*/


//TODO: find a way to clean up after use
workflow BUILD_REF_INDEX {
    take: 
        file
    main:
        NCBI_DOWNLOAD(file)
            MINIMAP2_INDEX(NCBI_DOWNLOAD.out.reference_genomes)

    emit:
        
}
