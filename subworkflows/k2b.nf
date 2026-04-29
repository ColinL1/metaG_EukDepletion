#!/usr/bin/env nextflow

/*
========================================================================================
    Set params for reads input and output
    This subworkflow is designed to run Kraken2 and Bracken with a shared DB loaded in RAM.
    This is very fast and efficient in a high memory workstation local environment, but not advised for a HPC with a queue manager. 
========================================================================================
*/

//K2 db and db name # Included in nextflow.config normally, but can be overridden here if needed.
params.kraken_db = '/share/databases/k2_nt'
params.kraken_dbName = 'k2_nt'

// Bracken (uses same DB directory as K2; bracken adds .kmer_distrib files there)
// read_len should match the value used when building the bracken DB (e.g. 150 for Illumina)
params.bracken_read_len_pe  = 150   // Illumina PE read length
params.bracken_read_len_ont = 300   // ONT: bracken is short-read oriented; adjust if DB was built differently
params.bracken_tax_level    = 'F'   // S=Species, G=Genus, F=Family, etc.
params.bracken_threshold    = 10  // Threshold for Bracken classification


include { load_in_shm; clean_shm } from '../modules/load_k2_db_shm.nf'
include { K2_SE ; K2_PE } from '../modules/K2.nf'
include { BRACKEN_PE ; BRACKEN_ONT } from '../modules/bracken.nf'


workflow K2_BRACKEN {
    take:
        reads_ch
        reads_ont_ch

    main:
    // ── Phase 1: Load DB ────────────────────────────────────────────────────
    // This can be change to a collect step to make the trigger be the end of all other processes
    // This way there is no conflict between loading the DB and other processes that may need RAM. 
    load_in_shm(Channel.value(true))
    db_signal = load_in_shm.out.db_loaded

    // ── Phase 2: Classification (K2 + Bracken) ─────────────────────────────
    // .combine(db_signal) ensures each K2 job only starts once the DB is loaded.
    K2_PE(
        reads_ch 
        .combine(db_signal)
        .map { meta, reads, _loaded -> [meta, reads] }
    )
    BRACKEN_PE(K2_PE.out.report_k2_out)

    K2_SE(
        reads_ont_ch
        .combine(db_signal)
        .map { meta, reads, _loaded -> [meta, reads] }
    )
    BRACKEN_ONT(K2_SE.out.report_k2_out)

    // ── Phase 3: Cleanup ────────────────────────────────────────────────────
    // Wait for ALL Bracken jobs (PE + ONT) before removing the DB from /dev/shm.
    // Gating on Bracken (not K2) is safer if Bracken is ever pointed at the shm path.
    all_k2_done = BRACKEN_PE.out.bracken
        .mix(BRACKEN_ONT.out.bracken)
        .collect()
    clean_shm(all_k2_done)

}