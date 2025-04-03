#!/usr/bin/env nextflow

/*
========================================================================================
    Run kraken and extract BRACKEN
========================================================================================
*/
params.kraken_db = "/share/databases/k2_nt"
params.read_length = "150"
/*
========================================================================================
    Include Modules
========================================================================================
*/

include { KRAKEN2_PE; KRAKEN2_SE } from '../modules/kraken2.nf'
include { BRACKEN; BRACKEN_ONT } from '../modules/bracken.nf'
include { ALPHA_DIVERSITY; BETA_DIVERSITY } from '../modules/div_metrics.nf'
// include {FASTP_REPORT} from './modules/fastp_stats.nf'

/*
========================================================================================
    Workflow KRAKEN_BRACKEN
========================================================================================
*/

workflow KRAKEN_BRACKEN {
    take: 
        file
        fastp_report_out
    main:
        KRAKEN2_SE(file)
            BRACKEN_ONT(KRAKEN2_SE.out.kraken_out, fastp_report_out)
                ALPHA_DIVERSITY(BRACKEN_ONT.out.bracken_out)  
                BETA_DIVERSITY(BRACKEN_ONT.out.bracken_out.collect{it[2]})
    emit:
        kraken2_report = KRAKEN2_SE.out.report_kraken_out
        bracken_report = BRACKEN_ONT.out.report_bracken_out 
        kraken_log = KRAKEN2_SE.out.log_kraken
        braken_log = BRACKEN_ONT.out.log_bracken
        BETA_DIVERSITY = BETA_DIVERSITY.out.beta_diversity_report
        ALPHA_DIVERSITY = ALPHA_DIVERSITY.out.alpha_diversity_report

}

workflow KRAKEN_BRACKEN_PE {
    take: 
        file
    main:
    // if (exists("/dev/shm/k2_nt") == true)
    //         KRAKEN2_PE(file)
    //             BRACKEN(KRAKEN2_PE.out.kraken_out)
    //                 ALPHA_DIVERSITY(BRACKEN.out.bracken_out)  
    //                 BETA_DIVERSITY(BRACKEN.out.bracken_out.collect{it[1]})
    //     shm_k2db = files("/dev/shm/*/*")
    //     result = shm_k2db.delete()
    //     println result ? "OK" : "Cannot delete: $myFile"
    // else
        // file("${params.kraken_db}/").copyTo("/dev/shm/") //, overwrit checkIfExists = false)
            KRAKEN2_PE(file)
                BRACKEN(KRAKEN2_PE.out.kraken_out)
                    ALPHA_DIVERSITY(BRACKEN.out.bracken_out)  
                    BETA_DIVERSITY(BRACKEN.out.bracken_out.collect{it[2]})
        // shm_k2db = files("/dev/shm/*/*")
        // result = shm_k2db.deleteDir()
        // println result ? "OK" : "Cannot delete: $myFile"
    emit:
        kraken2_report = KRAKEN2_PE.out.report_kraken_out
        bracken_report = BRACKEN.out.report_bracken_out
        BETA_DIVERSITY = BETA_DIVERSITY.out.beta_diversity_report
        ALPHA_DIVERSITY = ALPHA_DIVERSITY.out.alpha_diversity_report
}