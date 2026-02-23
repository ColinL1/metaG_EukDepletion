#!/usr/bin/env nextflow

include {CAT; CAT_ADD_NAMES} from '../modules/cat.nf'
include {CONCAT_FASTQ} from '../modules/cat_fastq.nf'

workflow CAT_WORKFLOW {
    take:
        reads
    main:
        CONCAT_FASTQ(reads)
            CAT(CONCAT_FASTQ.out.concat_reads)
                CAT_ADD_NAMES(CAT.out.contig2classification.concat(CAT.out.orf2lca))

        emit:
        contig2classification_taxname = CAT_ADD_NAMES.out.contig2classification_taxname
        orf2lca_taxname = CAT_ADD_NAMES.out.orf2lca_taxname
        cat_log = CAT.out.log

}