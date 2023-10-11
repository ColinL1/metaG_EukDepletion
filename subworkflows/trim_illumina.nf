#!/usr/bin/env nextflow

/*
========================================================================================
    Trim illumina long read data and generate fastp and multiqc report 
========================================================================================
*/

/*
========================================================================================
    Include Modules
========================================================================================
*/ 
include { trimmomatic_pe } from '../modules/trimmomatic.nf' 
include { fastqc } from '../modules/fastqc.nf' 
include { multiqc } from '../modules/multiqc.nf' 

// /home/colinl/metaG/Git/metaG_EukDepletion/modules/trimmomatic.nf 
// /home/colinl/metaG/Git/metaG_EukDepletion/modules/multiqc.nf
// /home/colinl/metaG/Git/metaG_EukDepletion/modules/fastqc.nf

include { fastp_report } from '../modules/fastp_stats.nf'
/*
========================================================================================
    Workflow TRIM_PE
========================================================================================
*/
// #TODO: integrate multiqc in this or main workflow. 
workflow TRIM_PE {
    take:
        reads
    main:
        trimmomatic_pe(reads)
            fastqc(trimmomatic_pe.out.trimmed_fq_out)
            fastp_report(trimmomatic_pe.out.trimmed_fq_out)
        
    emit:
        trimmed_reads = trimmomatic_pe.out.trimmed_fq_out
        report_json = fastp_report.out.report_json
        fastqc = fastqc.out.zip
        trim_log = trimmomatic_pe.out.log_trim
}

// // individual tester #TODO:remove once complete
// workflow {
//     // Channel.fromFilePairs("/home/colinl/metaG/Git/metaG_EukDepletion/input/test_illumina/*_{1,2}_subsample.fq.gz").collect(it[0]).set {input_fq}
//     Channel.fromFilePairs("/home/colinl/metaG/Git/metaG_EukDepletion/input/test_illumina/*_{1,2}_subsample.fq.gz"). set {input_fq}
//     TRIM_PE(input_fq)
//     // TRIM_PE.out.trimmed_reads.view()
//     // TRIM_PE.out.qc.view()
// }