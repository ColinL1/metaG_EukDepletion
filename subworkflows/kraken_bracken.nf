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
include { bracken } from '../modules/braken.nf'

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
    emit:
        kraken2_report = kraken2.out.report_kraken_out
        bracken_report = bracken.out.report_bracken_out 
        kraken_log = kraken.out.log_kraken
        braken_log = bracken.out.log_bracken
}

workflow KRAKEN_BRACKEN_PE {
    take: 
        file
    main:
        kraken2_pe(file)
            bracken(kraken2_pe.out.kraken_out)

    emit:
        kraken2_report = kraken2_pe.out.report_kraken_out
        bracken_report = bracken.out.report_bracken_out
}