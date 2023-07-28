#!/usr/bin/env nextflow

/*
========================================================================================
    MULTIQC module
========================================================================================
    Website: https://multiqc.info/
========================================================================================
*/

/*
========================================================================================
    Map to reference and split fastq
========================================================================================
*/

/*
========================================================================================
    Include Modules
========================================================================================
*/

include { kaiju } from '../modules/kaiju_multi.nf'
include { split_bac } from '../modules/split_seqkit.nf'
include { fastp_report } from '../modules/fastp_stats.nf'


/*
========================================================================================
    Workflow READ_QC
========================================================================================
*/

workflow EXTRACT_BACTERIA {
    take: 
        reads
    main:
        kaiju(reads)
            split_bac(kaiju.out.report_kaiju_out.join(reads))
                fastp_report(split_bac.out.concat())
    
    emit:
        bacteria_reads = split_bac.out.bacteria_reads
        non_bacteria_reads = split_bac.out.non_bacteria_reads
        report_json = fastp_report.out.report_json
}