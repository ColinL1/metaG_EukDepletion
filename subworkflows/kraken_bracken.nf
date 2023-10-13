#!/usr/bin/env nextflow

/*
========================================================================================
    Run kraken and extract bracken
========================================================================================
*/
params.kraken_db = "/share/databases/k2_nt"
params.read_length = "150"
/*
========================================================================================
    Include Modules
========================================================================================
*/

include { kraken2_pe; kraken2 } from '../modules/kraken2.nf'
include { bracken } from '../modules/bracken.nf'
include { alpha_diversity; beta_diversity } from '../modules/div_metrics.nf'

/*
========================================================================================
    Workflow KRAKEN_BRACKEN
========================================================================================
*/

workflow KRAKEN_BRACKEN {
    take: 
        file
    main:
        kraken2(file)
            bracken(kraken2_pe.out.kraken_out)
                alpha_diversity( bracken.out.bracken_out)  
                beta_diversity(bracken.out.bracken_out.collect{it[1]})
    emit:
        kraken2_report = kraken2.out.report_kraken_out
        bracken_report = bracken.out.report_bracken_out 
        kraken_log = kraken.out.log_kraken
        braken_log = bracken.out.log_bracken
        beta_diversity = beta_diversity.out.beta_diversity_report
        alpha_diversity = alpha_diversity.out.alpha_diversity_report

}

workflow KRAKEN_BRACKEN_PE {
    take: 
        file
    main:
    // if (exists("/dev/shm/k2_nt") == true)
    //         kraken2_pe(file)
    //             bracken(kraken2_pe.out.kraken_out)
    //                 alpha_diversity(bracken.out.bracken_out)  
    //                 beta_diversity(bracken.out.bracken_out.collect{it[1]})
    //     shm_k2db = files("/dev/shm/*/*")
    //     result = shm_k2db.delete()
    //     println result ? "OK" : "Cannot delete: $myFile"
    // else
        // file("${params.kraken_db}/").copyTo("/dev/shm/") //, overwrit checkIfExists = false)
            kraken2_pe(file)
                bracken(kraken2_pe.out.kraken_out)
                    alpha_diversity(bracken.out.bracken_out)  
                    beta_diversity(bracken.out.bracken_out.collect{it[1]})
        // shm_k2db = files("/dev/shm/*/*")
        // result = shm_k2db.deleteDir()
        // println result ? "OK" : "Cannot delete: $myFile"
    emit:
        kraken2_report = kraken2_pe.out.report_kraken_out
        bracken_report = bracken.out.report_bracken_out
        beta_diversity = beta_diversity.out.beta_diversity_report
        alpha_diversity = alpha_diversity.out.alpha_diversity_report
}

// // individual tester #TODO:remove once complete
// workflow {
//     Channel.fromFilePairs("/home/colinl/metaG/Git/metaG_EukDepletion/input/test_illumina/*_{1,2}_subsample.fq.gz"). set {input_fq}
//     // input_fq.view()
//     KRAKEN_BRACKEN_PE(input_fq)
// }