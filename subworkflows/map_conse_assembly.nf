#!/usr/bin/env nextflow

include { MINIMAP2_MAP_HOST ; MINIMAP2_MAP_SYM } from '../modules/minimap2.nf'
include { SPLIT_BAM_HOST ; SPLIT_BAM_SYM } from '../modules/samtools_split.nf'
include { KAIJU_SE } from '../modules/kaiju.nf'
include { MEGAHIT_PE } from '../modules/assembly_megahit.nf'

workflow MAP_CONSE_ASSEMBLY {
    take: 
        reads
    main:
    MEGAHIT_PE(reads)
        MINIMAP2_MAP_HOST(MEGAHIT_PE.out.contigs_fa) 
            SPLIT_BAM_HOST(MINIMAP2_MAP_HOST.out.mapp_file)
                MINIMAP2_MAP_SYM(SPLIT_BAM_HOST.out.unmapped_reads) // was run against mapped reads. not good. 
                    SPLIT_BAM_SYM(MINIMAP2_MAP_SYM.out.mapp_file)
                        KAIJU_SE(SPLIT_BAM_SYM.out.unmapped_reads)
    emit:
    kaiju_report = KAIJU_SE.out.report_kaiju_out

    bowtie2_host_bam = MINIMAP2_MAP_HOST.out.mapp_file
    bowtie2_host_map_fq = SPLIT_BAM_HOST.out.mapped_reads
    bowtie2_host_unmap_fq = SPLIT_BAM_HOST.out.unmapped_reads
    
    bowtie2_sym_bam = MINIMAP2_MAP_SYM.out.mapp_file
    bowtie2_sym_map_fq = SPLIT_BAM_SYM.out.mapped_reads
    bowtie2_sym_unmap_fq = SPLIT_BAM_SYM.out.unmapped_reads
}
