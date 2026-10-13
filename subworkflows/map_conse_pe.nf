#!/usr/bin/env nextflow

include { BOWTIE2_MAP_HOST ; BOWTIE2_MAP_SYM } from '../modules/bowtie2.nf'
include { KAIJU_PE ; KAIJU_ADD_TAXON_NAMES ; KAIJU_SUMMARISE } from '../modules/kaiju.nf'
include { SEQKIT_SPLIT_KAIJU_PE ; SEQKIT_STATS } from '../modules/seqkit.nf'

workflow MAP_CONSE_PE {
    take: 
        reads
    main:
        BOWTIE2_MAP_HOST(reads)
            BOWTIE2_MAP_SYM(BOWTIE2_MAP_HOST.out.unmapped_reads)
                KAIJU_PE(BOWTIE2_MAP_SYM.out.unmapped_reads)
                    KAIJU_ADD_TAXON_NAMES(KAIJU_PE.out.report_kaiju_out)
                        KAIJU_SUMMARISE(KAIJU_ADD_TAXON_NAMES.out.taxon)
                            SEQKIT_SPLIT_KAIJU_PE(BOWTIE2_MAP_SYM.out.unmapped_reads.join(KAIJU_SUMMARISE.out.ids))

        SEQKIT_STATS(
            reads
                .join(BOWTIE2_MAP_HOST.out.mapped_reads)
                .join(BOWTIE2_MAP_SYM.out.mapped_reads)
                .join(SEQKIT_SPLIT_KAIJU_PE.out.bacteria)
                .join(SEQKIT_SPLIT_KAIJU_PE.out.archaea)
                .join(SEQKIT_SPLIT_KAIJU_PE.out.viruses)
                .join(SEQKIT_SPLIT_KAIJU_PE.out.non_mapped)
                .map { row -> tuple(row[0], row[1..-1].flatten()) }
        )
        stats = SEQKIT_STATS.out.stats
            .collectFile(name: 'stats_reads.tsv', keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/reads_based/stats")
        kaiju_summary = KAIJU_SUMMARISE.out.summary
            .collectFile(name: 'kaiju_summary_reads.tsv', keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/reads_based/stats")
    emit:
        kaiju_report = KAIJU_PE.out.report_kaiju_out
        kaiju_summary = kaiju_summary
        bacteria_fq = SEQKIT_SPLIT_KAIJU_PE.out.bacteria
        archaea_fq = SEQKIT_SPLIT_KAIJU_PE.out.archaea
        viruses_fq = SEQKIT_SPLIT_KAIJU_PE.out.viruses
        non_mapped_fq = SEQKIT_SPLIT_KAIJU_PE.out.non_mapped
        stats = stats
        
        bowtie2_host_bam = BOWTIE2_MAP_HOST.out.mapp_file
        bowtie2_host_map_fq = BOWTIE2_MAP_HOST.out.mapped_reads
        bowtie2_host_unmap_fq = BOWTIE2_MAP_HOST.out.unmapped_reads

        bowtie2_sym_bam = BOWTIE2_MAP_SYM.out.mapp_file
        bowtie2_sym_map_fq = BOWTIE2_MAP_SYM.out.mapped_reads
        bowtie2_sym_unmap_fq = BOWTIE2_MAP_SYM.out.unmapped_reads
}
