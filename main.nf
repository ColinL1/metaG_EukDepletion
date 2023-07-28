#!/usr/bin/env nextflow

// Set params for reads and power 
// params.reads = "$baseDir/input/*.fastq.gz"
params.reads = "$baseDir/input/ONT/subset/*.fastq.gz"
params.contigs = "$baseDir/input/*.fasta"
params.outdir = "$baseDir/results_ont/"

//temporary fix to multiply the reference channel to match the number of samples
params.n_samples = 120

// params for mmeseqs2

// params for kaiju
params.nodes = "/home/colinl/database/kaiju/nodes.dmp"
params.kaiju_db = "/home/colinl/database/kaiju/refseq/kaiju_db_refseq.fmi"

// parms for minimap2
params.ref_coral = "/home/colinl/metaG/Git/metaG_EukDepletion/kaiju/mem/split/ref_genomes/corals.mmi"
params.ref_sym = "/home/colinl/metaG/Git/metaG_EukDepletion/kaiju/mem/split/ref_genomes/symbiodiniaceae.mmi"

log.info """\

    metaG kingdom taxonomy - NF   PIPELINE
    =======================================
    Diamond_db : ${params.kaiju_db}
    reads      : ${params.reads}
    outdir     : ${params.outdir}

    """

Channel.fromPath(params.reads).map {tuple( it.name.split('.fastq.gz')[0], it )}.set { seq_reads_ch }
// Channel.fromPath(params.contigs).map {tuple( (it.name.split('.fasta')[0]), it )}.set { seq_contigs_ch }
// Channel.fromPath(params.contigs).map {tuple( (it.name.split('.fasta')[0]).split('_')[0..3].join("_"), it )}.set { seq_contigs_ch } // version that cleans the names (specific to my use case)

include {MAP2REF_SWF as MAP2REF_SWF_CORAL; MAP2REF_SWF as MAP2REF_SWF_SYM} from './subworkflows/map2ref_extract.nf'
include {EXTRACT_BACTERIA } from './subworkflows/extract_bacteria.nf'
include {TRIM_NANOPORE } from './subworkflows/trim_nanopore.nf'

include { fastp_parse } from './modules/fastp_parse.nf'


Channel.fromPath(params.ref_coral) 
            .map {tuple( it.name.split('.mmi')[0], it )}
            .flatMap {it * 120 }
            .collate(2)
            .set { ref_coral_ch }

Channel.fromPath(params.ref_sym) 
            .map {tuple( it.name.split('.mmi')[0], it )}
            .flatMap {it * 120 }
            .collate(2)
            .set { ref_sym_ch }

workflow ONT_VF {
        params.mode = 'ONT'

    TRIM_NANOPORE(seq_reads_ch)
        EXTRACT_BACTERIA(TRIM_NANOPORE.out.trimmed_reads)
            MAP2REF_SWF_CORAL(EXTRACT_BACTERIA.out.non_bacteria_reads, ref_coral_ch)
                MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        

        //     fastp_parse_ch = TRIM_NANOPORE.out.report_json.concat(
        //     EXTRACT_BACTERIA.out.report_json,
        //     MAP2REF_SWF_CORAL.out.report_json,
        //     MAP2REF_SWF_SYM.out.report_json)
        // fastp_parse(fastp_parse_ch)
}

// workflow ONT_pipeline {
//     main:
//         mode = 'ONT'
//         extract_sym()
        
//         fastp_parse(extract_sym.out.report_json.concat(
//             extract_sym.out.report_json_co,
//             extract_sym.out.report_json_bc,
//             extract_sym.out.report_json_og))
// }


workflow.onComplete {    
    log.info """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

	log.info ( workflow.success ? "\nDone! --> $params.outdir\n" : "Oops .. something went wrong" )
}