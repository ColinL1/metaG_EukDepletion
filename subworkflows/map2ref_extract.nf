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

include { MINIMAP2_MAP } from '../modules/minimap2.nf' 
include { SPLIT_BAM; SPLIT_BAM_BWA } from '../modules/split_bam.nf'
include { FASTP_REPORT } from '../modules/fastp_stats.nf'
// include { BOWTIE2_MAP; map2ref_pe_vsensitive } from '../modules/bowtie2.nf' 
include { BOWTIE2_MAP } from '../modules/bowtie2.nf' 
// include { BWA_PE; BWA_ONT } from '../modules/bwa.nf' 

/*
========================================================================================
    Workflow MAP2REF_SWF
========================================================================================
*/

workflow MAP2REF_SWF {
    take:
        reads

    main:    
        MINIMAP2_MAP(reads)
            SPLIT_BAM(MINIMAP2_MAP.out.mapp_file)
                FASTP_REPORT(SPLIT_BAM.out.mapped_reads.concat(SPLIT_BAM.out.unmapped_reads))

    emit:
        mapped_reads = SPLIT_BAM.out.mapped_reads
        non_mapped_reads = SPLIT_BAM.out.unmapped_reads
        report_json = FASTP_REPORT.out.report_json
}

workflow MAP2REF_SWF_PE {
    take:
        data

    main:    
        BOWTIE2_MAP(data)
            FASTP_REPORT(BOWTIE2_MAP.out.mapped_reads.concat(BOWTIE2_MAP.out.unmapped_reads))

    emit:
        mapped_reads = BOWTIE2_MAP.out.mapped_reads
        non_mapped_reads = BOWTIE2_MAP.out.unmapped_reads
        report_json = FASTP_REPORT.out.report_json
}
workflow MAP2REF_S_SWF_PE {
    take:
        reads
        reference_ch

    main:    
        map2ref_pe_vsensitive(reads, reference_ch)
            FASTP_REPORT(map2ref_pe_vsensitive.out.mapped_reads.concat(map2ref_pe_vsensitive.out.unmapped_reads))

    emit:
        mapped_reads = map2ref_pe_vsensitive.out.mapped_reads
        non_mapped_reads = map2ref_pe_vsensitive.out.unmapped_reads
        report_json = FASTP_REPORT.out.report_json
}

workflow MAP2REF_BWA_PE {
    take:
        reads
        reference_ch

    main:    
        BWA_PE(reads, reference_ch)
            SPLIT_BAM_BWA(BWA_PE.out.mapp_file, reference_ch)
                FASTP_REPORT(SPLIT_BAM_BWA.out.mapped_reads.concat(SPLIT_BAM_BWA.out.unmapped_reads))

    emit:
        mapped_reads = SPLIT_BAM_BWA.out.mapped_reads
        non_mapped_reads = SPLIT_BAM_BWA.out.unmapped_reads
        report_json = FASTP_REPORT.out.report_json
}

workflow MAP2REF_BWA_ONT {
    take:
        reads
        reference_ch

    main:    
        BWA_ONT(reads, reference_ch)
            SPLIT_BAM(BWA_PE.out.mapp_file, reference_ch)
                FASTP_REPORT(SPLIT_BAM.out.mapped_reads.concat(SPLIT_BAM.out.unmapped_reads))

    emit:
        mapped_reads = SPLIT_BAM.out.mapped_reads
        non_mapped_reads = SPLIT_BAM.out.unmapped_reads
        report_json = FASTP_REPORT.out.report_json
}