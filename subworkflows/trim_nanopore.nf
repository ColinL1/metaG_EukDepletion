#!/usr/bin/env nextflow

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
include { ont_trim } from '../modules/ONT_prep.nf'
include { fastp_report } from '../modules/fastp_stats.nf'
/*
========================================================================================
    Workflow READ_QC
========================================================================================
*/

workflow TRIM_NANOPORE {
    take:
        reads

    main:
        ont_trim(reads)
            fastp_report(ont_trim.out.trimmed_fastq)
    emit:
        trimmed_reads = ont_trim.out.trimmed_fastq
        report_json = fastp_report.out.report_json
}