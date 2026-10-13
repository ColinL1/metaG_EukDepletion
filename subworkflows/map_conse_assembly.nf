#!/usr/bin/env nextflow

include { MINIMAP2_MAP_HOST ; MINIMAP2_MAP_SYM } from '../modules/minimap2.nf'
include { SPLIT_BAM_HOST ; SPLIT_BAM_SYM } from '../modules/samtools_split.nf'
include { KAIJU_SE ; KAIJU_ADD_TAXON_NAMES ; KAIJU_SUMMARISE } from '../modules/kaiju.nf'
include { SEQKIT_SPLIT_KAIJU_SE ; SEQKIT_STATS } from '../modules/seqkit.nf'
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
                            KAIJU_ADD_TAXON_NAMES(KAIJU_SE.out.report_kaiju_out)
                                KAIJU_SUMMARISE(KAIJU_ADD_TAXON_NAMES.out.taxon)
                                    SEQKIT_SPLIT_KAIJU_SE(SPLIT_BAM_SYM.out.unmapped_reads.join(KAIJU_SUMMARISE.out.ids))

    SEQKIT_STATS(
        MEGAHIT_PE.out.contigs_fa
            .join(SPLIT_BAM_HOST.out.mapped_reads)
            .join(SPLIT_BAM_SYM.out.mapped_reads)
            .join(SEQKIT_SPLIT_KAIJU_SE.out.bacteria)
            .join(SEQKIT_SPLIT_KAIJU_SE.out.archaea)
            .join(SEQKIT_SPLIT_KAIJU_SE.out.viruses)
            .join(SEQKIT_SPLIT_KAIJU_SE.out.non_mapped)
            .map { row -> tuple(row[0], row[1..-1].flatten()) }
    )
    stats = SEQKIT_STATS.out.stats
        .collectFile(name: 'stats_assembly.tsv', keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/assembly_based/stats")
    kaiju_summary = KAIJU_SUMMARISE.out.summary
        .collectFile(name: 'kaiju_summary_assembly.tsv', keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/assembly_based/stats")
    emit:
    kaiju_report = KAIJU_SE.out.report_kaiju_out
    kaiju_summary = kaiju_summary
    bacteria_fa = SEQKIT_SPLIT_KAIJU_SE.out.bacteria
    archaea_fa = SEQKIT_SPLIT_KAIJU_SE.out.archaea
    viruses_fa = SEQKIT_SPLIT_KAIJU_SE.out.viruses
    non_mapped_fa = SEQKIT_SPLIT_KAIJU_SE.out.non_mapped
    stats = stats

    bowtie2_host_bam = MINIMAP2_MAP_HOST.out.mapp_file
    bowtie2_host_map_fq = SPLIT_BAM_HOST.out.mapped_reads
    bowtie2_host_unmap_fq = SPLIT_BAM_HOST.out.unmapped_reads
    
    bowtie2_sym_bam = MINIMAP2_MAP_SYM.out.mapp_file
    bowtie2_sym_map_fq = SPLIT_BAM_SYM.out.mapped_reads
    bowtie2_sym_unmap_fq = SPLIT_BAM_SYM.out.unmapped_reads
}
