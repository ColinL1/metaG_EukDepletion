#!/usr/bin/env nextflow

include { CAT; CAT_ADD_NAMES; CAT_FILTER_BACTERIA; CAT_PREP_FAMILY; CAT_MERGE_TABLES } from '../modules/cat.nf'
include { CONCAT_FASTQ } from '../modules/concat_fastq.nf'
include { MEGAHIT_PE as MEGAHIT_COASSEMBLY } from '../modules/assembly_megahit.nf'
include { SEQKIT_PREFIX_CONTIGS } from '../modules/seqkit.nf'
include { POOL_CONTIGS; MMSEQS_LINCLUST } from '../modules/mmseqs2.nf'
include { MINIMAP2_INDEX } from '../modules/minimap2_index.nf'
include { MINIMAP2_MAP_CATALOGUE } from '../modules/minimap2.nf'
include { COVERM_CONTIG } from '../modules/coverm.nf'

workflow CAT_WORKFLOW {
    take:
        reads
    main:
        def quant_mode = params.cat_quant_mode
        if( !(quant_mode in ['persample', 'pooled', 'both']) ) {
            error "Invalid cat_quant_mode: ${quant_mode} (persample | pooled | both)"
        }

        CONCAT_FASTQ(reads)
            MEGAHIT_COASSEMBLY(CONCAT_FASTQ.out.concat_reads)
                CAT(MEGAHIT_COASSEMBLY.out.contigs_fa)
                    CAT_ADD_NAMES(CAT.out.contig2classification.concat(CAT.out.orf2lca))

        CAT_FILTER_BACTERIA(MEGAHIT_COASSEMBLY.out.contigs_fa.join(CAT_ADD_NAMES.out.contig2classification_taxname))
            CAT_PREP_FAMILY(CAT_FILTER_BACTERIA.out.bacteria_class)
            SEQKIT_PREFIX_CONTIGS(CAT_FILTER_BACTERIA.out.bacteria_fa)

        // MINIMAP2_INDEX can be called once only, so both reference types
        // (persample / pooled) go through one mixed channel
        index_in = channel.empty()

        if( quant_mode in ['persample', 'both'] ) {
            index_in = index_in.mix(
                SEQKIT_PREFIX_CONTIGS.out.renamed
                    .map { meta, fa -> tuple(meta, 'persample', fa) }
            )
        }

        if( quant_mode in ['pooled', 'both'] ) {
            POOL_CONTIGS(
                SEQKIT_PREFIX_CONTIGS.out.renamed
                    .map { meta, fa -> tuple(meta.species, fa) }
                    .groupTuple()
            )
                MMSEQS_LINCLUST(POOL_CONTIGS.out.pooled)
            index_in = index_in.mix(
                MMSEQS_LINCLUST.out.catalogue
                    .map { species, fa, _clusters -> tuple([id: species, species: species], 'pooled', fa) }
            )
        }

        MINIMAP2_INDEX(index_in)

        map_in = channel.empty()

        if( quant_mode in ['persample', 'both'] ) {
            persample_idx = MINIMAP2_INDEX.out.index
                .filter { _meta, mode, _mmi -> mode == 'persample' }
                .map { meta, _mode, mmi -> tuple(meta.id, mmi) }
            map_in = map_in.mix(
                CONCAT_FASTQ.out.concat_reads
                    .map { meta, fq -> tuple(meta.id, meta, fq) }
                    .join(persample_idx)
                    .map { _key, meta, fq, mmi -> tuple(meta, 'persample', fq, mmi) }
            )
        }

        if( quant_mode in ['pooled', 'both'] ) {
            pooled_idx = MINIMAP2_INDEX.out.index
                .filter { _meta, mode, _mmi -> mode == 'pooled' }
                .map { meta, _mode, mmi -> tuple(meta.species, mmi) }
            map_in = map_in.mix(
                CONCAT_FASTQ.out.concat_reads
                    .map { meta, fq -> tuple(meta.species, meta, fq) }
                    .combine(pooled_idx, by: 0)
                    .map { _key, meta, fq, mmi -> tuple(meta, 'pooled', fq, mmi) }
            )
        }

        MINIMAP2_MAP_CATALOGUE(map_in)
            COVERM_CONTIG(MINIMAP2_MAP_CATALOGUE.out.bam)

        // sample metadata table expected by merge_tables.R
        sample_table = CONCAT_FASTQ.out.concat_reads
            .map { meta, _fq -> "${meta.id},${meta.species},${meta.method},${meta.buffer}" }
            .collect()
            .map { rows -> (['sample,species,treatment,buffer'] + rows.sort()).join('\n') + '\n' }
            .collectFile(name: 'cat_quant_samples.csv')
            .first()

        CAT_MERGE_TABLES(
            COVERM_CONTIG.out.coverage.groupTuple(),
            CAT_PREP_FAMILY.out.family.collect(),
            sample_table
        )

    emit:
        contig2classification_taxname = CAT_ADD_NAMES.out.contig2classification_taxname
        orf2lca_taxname = CAT_ADD_NAMES.out.orf2lca_taxname
        cat_log = CAT.out.log
        coverm = COVERM_CONTIG.out.coverage
        contig_table = CAT_MERGE_TABLES.out.contig_table
        family_table = CAT_MERGE_TABLES.out.family_table
        mapping_summary = CAT_MERGE_TABLES.out.mapping_summary

}
