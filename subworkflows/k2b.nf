#!/usr/bin/env nextflow

/*
========================================================================================
    Set params for reads input and output
========================================================================================
*/

params.input = "$baseDir/input/illumina/*_{1,2}.fastq.gz"
// params.contigs = "$baseDir/input/coral"
params.outdir = "$baseDir/results/illumina"
input_pe = "/home/colinl/metaG/Git/metaG_EukDepletion/input/illumina/*_{1,2}.fastq.gz"

/*
========================================================================================
    Set params for kaiju database and minimap2 index references files
========================================================================================
*/
// params.nodes = "/home/colinl/database/kaiju/nodes.dmp"
// params.kaiju_db = "/home/colinl/database/kaiju/refseq/kaiju_db_refseq.fmi"

params.kraken_db = "/share/databases/k2_nt"
params.read_length = "150"

// params.ref_host = "/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_bowtie2/all_scleractina" //"/home/colinl/metaG/Git/metaG_EukDepletion/kaiju/mem/split/ref_genomes/aiptasiidae.mmi" //"/home/colinl/metaG/Git/metaG_EukDepletion/kaiju/mem/split/ref_genomes/corals.mmi"
// params.ref_sym = "/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_bowtie2/all_symbiodiniaceae"

//temporary fix to multiply the reference channel to match the number of samples
// params.n_samples = 120

log.info """\

    metaG kingdom taxonomy - NF  PIPELINE
    ======================================
    kraken_db  : ${params.kraken_db}
    reads      : ${params.input}
    outdir     : ${params.outdir}

    """

// Channel.fromPath(params.input+'/*.fastq.gz').map {tuple( it.name.split('.fastq.gz')[0], it )}.set { seq_reads_ch }
// Channel.fromPath(params.contigs).map {tuple( (it.name.split('.fasta')[0]), it )}.set { seq_contigs_ch }
// Channel.fromPath(params.input+'/*.fasta').map {tuple( (it.name.split('.fasta')[0]).split('_')[0..3].join("_"), it )}.set { seq_contigs_ch } // version that cleans the names (specific to my use case)

Channel.fromFilePairs(input_pe).set { seq_reads_pe_ch }

// include {MAP2REF_SWF as MAP2REF_SWF_CORAL; MAP2REF_SWF as MAP2REF_SWF_SYM} from './subworkflows/map2ref_extract.nf'
// include {MAP2REF_SWF_PE as MAP2REF_SWF_CORAL_PE; MAP2REF_SWF_PE as MAP2REF_SWF_SYM_PE} from './subworkflows/map2ref_extract.nf'
// include {EXTRACT_BACTERIA; EXTRACT_BACTERIA_FA; EXTRACT_BACTERIA_PE } from './subworkflows/extract_bacteria.nf'
// include {TRIM_NANOPORE } from './subworkflows/trim_nanopore.nf'

include {KRAKEN_BRACKEN; KRAKEN_BRACKEN_PE} from './subworkflows/kraken_bracken.nf'

// include { megahit_se } from './modules/assembly_megahit.nf'
// include { fastp_report } from './modules/fastp_stats.nf'
// include { fastp_parse } from './modules/fastp_parse.nf'
// include { fasta2fastq } from './modules/ONT_prep.nf'

// Channel.fromPath(params.ref_host) 
//             .map {tuple( it.name.split('.mmi')[0], it )}
//             .flatMap {it * 120 }
//             .collate(2)
//             .set { ref_host_ch }

// Channel.fromPath(params.ref_sym) 
//             .map {tuple( it.name.split('.mmi')[0], it )}
//             .flatMap {it * 120 }
//             .collate(2)
//             .set { ref_sym_ch }

// workflow ONT_VF {
//         params.mode = 'ONT'

//     TRIM_NANOPORE(seq_reads_ch)
//         EXTRACT_BACTERIA(TRIM_NANOPORE.out.trimmed_reads)
//             MAP2REF_SWF_CORAL(EXTRACT_BACTERIA.out.non_bacteria_reads, ref_host_ch)
//                 MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        

//             fastp_parse_ch = TRIM_NANOPORE.out.report_json.concat(
//             EXTRACT_BACTERIA.out.report_json,
//             MAP2REF_SWF_CORAL.out.report_json,
//             MAP2REF_SWF_SYM.out.report_json)
//                 // fastp_parse_ch.collect{it[1]}.toList().view()
//                 fastp_parse(fastp_parse_ch.collect{it[1]}.toList())
// }

// workflow CONTIGS {
//         params.mode = 'contigs'

//     fastp_fasta_report(seq_contigs_ch)
//     EXTRACT_BACTERIA_FA(seq_contigs_ch)
//         MAP2REF_SWF_CORAL(EXTRACT_BACTERIA_FA.out.non_bacteria_reads, ref_host_ch)
//             MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        

//             fastp_parse_ch = fastp_fasta_report.out.report_json.concat(
//             EXTRACT_BACTERIA_FA.out.report_json,
//             MAP2REF_SWF_CORAL.out.report_json,
//             MAP2REF_SWF_SYM.out.report_json)
//                 // fastp_parse_ch.collect{it[1]}.toList().view()
//                 fastp_parse(fastp_parse_ch.collect{it[1]}.toList())
// }

workflow  {
    KRAKEN_BRACKEN_PE(seq_reads_pe_ch)
}

workflow TEST {
    seq_reads_pe_ch.view()
    // KRAKEN_BRACKEN(seq_reads_pe_ch)
}

// workflow NON_TRIMMED_assembly  {
    
//     megahit_se(seq_reads_ch)
//     fasta2fastq(seq_contigs_ch.concat(megahit_se.out.contigs_fa))
//         fastp_report(fasta2fastq.out.converted_fasta)

//         EXTRACT_BACTERIA(fasta2fastq.out.converted_fasta)
//             MAP2REF_SWF_CORAL(EXTRACT_BACTERIA.out.non_bacteria_reads, ref_host_ch)
//                 MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        

//             fastp_parse_ch = fastp_report.out.report_json.concat(
//             EXTRACT_BACTERIA.out.report_json,
//             MAP2REF_SWF_CORAL.out.report_json,
//             MAP2REF_SWF_SYM.out.report_json)
//                 // fastp_parse_ch.collect{it[1]}.toList().view()
//                 fastp_parse(fastp_parse_ch.collect{it[1]}.toList())
// }

// workflow NON_TRIMMED  {
    
//     fasta2fastq(seq_contigs_ch)
//         fastp_report(fasta2fastq.out.converted_fasta.concat(seq_reads_ch))

//         EXTRACT_BACTERIA(fasta2fastq.out.converted_fasta.concat(seq_reads_ch))
//             MAP2REF_SWF_CORAL(EXTRACT_BACTERIA.out.non_bacteria_reads, ref_host_ch)
//                 MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        

//             fastp_parse_ch = fastp_report.out.report_json.concat(
//             EXTRACT_BACTERIA.out.report_json,
//             MAP2REF_SWF_CORAL.out.report_json,
//             MAP2REF_SWF_SYM.out.report_json)
//                 // fastp_parse_ch.collect{it[1]}.toList().view()
//                 fastp_parse(fastp_parse_ch.collect{it[1]}.toList())
// }


// workflow PE_ILLUMINA  {
//     // seq_reads_pe_ch.view()
//     fastp_report(seq_reads_pe_ch)
//     EXTRACT_BACTERIA_PE(seq_reads_pe_ch)
//         MAP2REF_SWF_CORAL_PE(EXTRACT_BACTERIA_PE.out.non_bacteria_reads, ref_host_ch)
//             MAP2REF_SWF_SYM_PE(MAP2REF_SWF_CORAL_PE.out.non_mapped_reads, ref_sym_ch)        

//             fastp_parse_ch = fastp_report.out.report_json.concat(
//             EXTRACT_BACTERIA_PE.out.report_json,
//             MAP2REF_SWF_CORAL_PE.out.report_json,
//             MAP2REF_SWF_SYM_PE.out.report_json)
//                 fastp_parse_ch.collect{it[1]}.toList().view()
//                 // fastp_parse(fastp_parse_ch.collect{it[1]}.toList())
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