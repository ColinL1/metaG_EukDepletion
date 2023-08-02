#!/usr/bin/env nextflow

/*
========================================================================================
    Map to reference file and split fastq and generate fastp report
========================================================================================
*/

/*
========================================================================================
    Include Modules
========================================================================================
*/

include { map2ref } from '../modules/minimap2.nf' 
include { split_bam } from '../modules/split_bam.nf'
include { fastp_report } from '../modules/fastp_stats.nf'

/*
========================================================================================
    Workflow MAP2REF_SWF
========================================================================================
*/

workflow MAP2REF_SWF {
    take:
        reads
        reference_ch

    main:    
        map2ref(reads, reference_ch)
            split_bam(map2ref.out.mapp_file, reference_ch)
                fastp_report(split_bam.out.mapped_reads.concat(split_bam.out.unmapped_reads))

    emit:
        mapped_reads = split_bam.out.mapped_reads
        non_mapped_reads = split_bam.out.unmapped_reads
        report_json = fastp_report.out.report_json
}