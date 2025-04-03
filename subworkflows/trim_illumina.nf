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
include { TRIMMOMATIC } from '../modules/trimmomatic.nf' 
include { FASTQC } from '../modules/fastqc.nf' 
// include { MULTIQC } from '../modules/multiqc.nf' 

// /home/colinl/metaG/Git/metaG_EukDepletion/modules/trimmomatic.nf 
// /home/colinl/metaG/Git/metaG_EukDepletion/modules/multiqc.nf
// /home/colinl/metaG/Git/metaG_EukDepletion/modules/FASTQC.nf

include { FASTP_REPORT } from '../modules/fastp_stats.nf'
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
        TRIMMOMATIC(reads)
            FASTQC(TRIMMOMATIC.out.trimmed_fq_out)
            FASTP_REPORT(TRIMMOMATIC.out.trimmed_fq_out)
        
    emit:
        trimmed_reads = TRIMMOMATIC.out.trimmed_fq_out
        report_json = FASTP_REPORT.out.report_json
        FASTQC = FASTQC.out.zip
        trim_log = TRIMMOMATIC.out.log_trim
}

// // individual tester #TODO:remove once complete
// workflow {
//     // Channel.fromFilePairs("/home/colinl/metaG/Git/metaG_EukDepletion/input/test_illumina/*_{1,2}_subsample.fq.gz").collect(it[0]).set {input_fq}
//     Channel.fromFilePairs("/home/colinl/metaG/Git/metaG_EukDepletion/input/test_illumina/*_{1,2}_subsample.fq.gz"). set {input_fq}
//     TRIM_PE(input_fq)
//     // TRIM_PE.out.trimmed_reads.view()
//     // TRIM_PE.out.qc.view()
// }