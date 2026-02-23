#!/usr/bin/env nextflow

include { BOWTIE2_MAP_HOST ; BOWTIE2_MAP_SYM } from '../modules/bowtie2.nf'
include {KAIJU_PE} from '../modules/kaiju.nf'

workflow MAP_CONSE_PE {
    take: 
        reads
    main:
        BOWTIE2_MAP_HOST(reads)
            BOWTIE2_MAP_SYM(BOWTIE2_MAP_HOST.out.unmapped_reads)
                KAIJU_PE(BOWTIE2_MAP_SYM.out.unmapped_reads)
    emit:
        kaiju_report = KAIJU_PE.out.report_kaiju_out
        
        bowtie2_host_bam = BOWTIE2_MAP_HOST.out.mapp_file
        bowtie2_host_map_fq = BOWTIE2_MAP_HOST.out.mapped_reads
        bowtie2_host_unmap_fq = BOWTIE2_MAP_HOST.out.unmapped_reads

        bowtie2_sym_bam = BOWTIE2_MAP_SYM.out.mapp_file
        bowtie2_sym_map_fq = BOWTIE2_MAP_SYM.out.mapped_reads
        bowtie2_sym_unmap_fq = BOWTIE2_MAP_SYM.out.unmapped_reads
}
