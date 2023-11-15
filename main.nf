#!/usr/bin/env nextflow
//TODO: general cleanup and reestablishment workflows and subworkflows!!
//TODO: add stub
/*
========================================================================================
    Set params for reads input and output
========================================================================================
*/
params.input = "$baseDir/input/coral_pe_data/*_{1,2}.fq.gz"
params.outdir = "$baseDir/results/illumina"
// input_pe = '/home/colinl/metaG/Git/metaG_EukDepletion/input/illumina/corals/*_{1,2}.fastq.gz'

Channel.fromFilePairs(params.input).set { seq_reads_pe_ch }
Channel.fromPath(params.input).map {tuple( it.name.split('.fastq.gz')[0], it )}.set { seq_reads_ch_ont }

/*
========================================================================================
    Set params for kaiju database and minimap2 index references files
========================================================================================
*/
// import kaiju subworkflow 
include { MAP2REF_SWF as MAP2REF_SWF_CORAL; MAP2REF_SWF as MAP2REF_SWF_SYM } from './subworkflows/map2ref_extract.nf'
include { MAP2REF_SWF_PE as MAP2REF_SWF_CORAL_PE; MAP2REF_SWF_PE as MAP2REF_SWF_SYM_PE } from './subworkflows/map2ref_extract.nf'
include { EXTRACT_BACTERIA; EXTRACT_BACTERIA_FA; EXTRACT_BACTERIA_PE } from './subworkflows/extract_bacteria.nf'
include { TRIM_PE } from './subworkflows/trim_illumina.nf'

params.nodes = "/share/databases/kaiju/nodes.dmp"
params.kaiju_db = "/share/databases/kaiju/refseq/kaiju_db_refseq.fmi"
params.barplot_script = "/home/colinl/metaG/Git/metaG_EukDepletion/R-scripts/Bar_plot_kaiju.r"
params.parse_script = "/home/colinl/metaG/Git/metaG_EukDepletion/parse_json_cl.py"
params.sample_metadata_sheet = "/home/colinl/metaG/Git/metaG_EukDepletion/sample_metadata_sheet.csv"

params.ref_host = "/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_bowtie2/all_scleractina" //"/home/colinl/metaG/Git/metaG_EukDepletion/kaiju/mem/split/ref_genomes/aiptasiidae.mmi" //"/home/colinl/metaG/Git/metaG_EukDepletion/kaiju/mem/split/ref_genomes/corals.mmi"
params.ref_sym = "/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_bowtie2/all_symbiodiniaceae"

// TODO: implemnt the "untill" operator to match number of repeats time the number of fastq input
Channel.fromPath(params.ref_host) 
            .map {tuple( it.name.split('.mmi')[0], it )}
            .flatMap { it * 120 }
            .collate(2)
            .set { ref_host_ch }

Channel.fromPath(params.ref_sym) 
            .map {tuple( it.name.split('.mmi')[0], it )}
            .flatMap { it * 120 }
            .collate(2)
            .set { ref_sym_ch }
//temporary fix to multiply the reference channel to match the number of samples
params.n_samples = 120

/*
========================================================================================
    Set params for kaiju database and minimap2 index references files
========================================================================================
*/
//import kraken bracken subworkflows
include {KRAKEN_BRACKEN; KRAKEN_BRACKEN_PE} from './subworkflows/kraken_bracken.nf'

params.kraken_db = "/share/databases/kraken2_custom/metaG"
params.read_length = "150"

// import qc modules and workflow
include { fastp_report } from './modules/fastp_stats.nf'
include { fastp_parse } from './modules/fastp_parse.nf'
include { fastp_plot } from './modules/fastp_parse.nf'
include { fasta2fastq } from './modules/ONT_prep.nf'

// TODO import module for assembly variant (currently unused.) 
include { megahit_se; megahit_pe } from './modules/assembly_megahit.nf'

// TODO check for importance of this for ONT workflow
// Channel.fromPath(params.input+'/*.fastq.gz').map {tuple( it.name.split('.fastq.gz')[0], it )}.set { seq_reads_ch }
// // Channel.fromPath(params.contigs).map {tuple( (it.name.split('.fasta')[0]), it )}.set { seq_contigs_ch }
// Channel.fromPath(params.input+'/*.fasta').map {tuple( (it.name.split('.fasta')[0]).split('_')[0..3].join("_"), it )}.set { seq_contigs_ch } // version that cleans the names (specific to my use case)

// include { MAP2REF_BWA_PE  } from './subworkflows/extract_bacteria.nf'
// include { MAP2REF_BWA_PE as MAP2REF_BWA_PE_CORAL; MAP2REF_BWA_PE as MAP2REF_BWA_PE_SYM } from './subworkflows/map2ref_extract.nf'
// include { MAP2REF_S_SWF_PE as MAP2REF_S_SWF_CORAL_PE; MAP2REF_S_SWF_PE as MAP2REF_S_SWF_SYM_PE } from './subworkflows/map2ref_extract.nf'

workflow test_ch {
    // ref_sym_ch.view()
    // ref_host_ch.view()
    seq_reads_pe_ch.view()
}

// TODO: add switch for kaiju vs kraken bracken
log.info """\

    metaG kingdom taxonomy - NF   PIPELINE
    =======================================
    kaiju_db   : ${params.kaiju_db}
    kraken2    : ${params.kraken_db}
    reads      : ${params.input}
    outdir     : ${params.outdir}

    """

//TODO : WORKING AS IS. TO BE optimised. (adding database mv to shared memory)
workflow KRBR_PE_ILLUMINA{
    TRIM_PE(seq_reads_pe_ch)
        fastp_report(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))
            KRAKEN_BRACKEN_PE(TRIM_PE.out.trimmed_reads)
            //TODO add qc
}

//TODO : WORKING AS IS. TO BE optimised. 
workflow KAIJU_PE_ILLUMINA {
    // main:
        kaiju_nodes = file(params.nodes)
        kaiju_db = file(params.kaiju_db)
// TODO check integration with fastqc + multiqc and optional make contigs assembly version. and fastp
    // seq_reads_pe_ch.view()
        TRIM_PE(seq_reads_pe_ch)
            // fastp_report(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))
            fastp_report(TRIM_PE.out.trimmed_reads)
        EXTRACT_BACTERIA_PE(TRIM_PE.out.trimmed_reads)
            MAP2REF_SWF_CORAL_PE(EXTRACT_BACTERIA_PE.out.non_bacteria_reads, ref_host_ch)
                MAP2REF_SWF_SYM_PE(MAP2REF_SWF_CORAL_PE.out.non_mapped_reads, ref_sym_ch)        

                fastp_parse_ch = fastp_report.out.report_json.concat(
                EXTRACT_BACTERIA_PE.out.report_json,
                MAP2REF_SWF_CORAL_PE.out.report_json,
                MAP2REF_SWF_SYM_PE.out.report_json)
                    // fastp_parse_ch.collect{it[1]}.toList().view()
                names = fastp_parse_ch.collect{it[0]}.toList()   
                paths = fastp_parse_ch.collect{it[1]}.toList()
                names.merge(paths).set { fastp_parse_ch_2 }
                    fastp_plot(fastp_parse_ch_2)
    // emit:
    //     EXTRACT_BACTERIA_PE.out.bacteria_reads
    //     MAP2REF_SWF_CORAL_PE.out.mapped_reads
    //     MAP2REF_SWF_SYM_PE.out.mapped_reads
    //     EXTRACT_BACTERIA_PE.out.report_json
    //     MAP2REF_SWF_CORAL_PE.out.report_json
    //     MAP2REF_SWF_SYM_PE.out.report_json
    //     fastp_parse.out.report_csv
}

//TODO : WORKING AS IS. TO BE optimised. 
workflow EXTRACT_CORAL_SYMB {
    // main:
        kaiju_nodes = file(params.nodes)
        kaiju_db = file(params.kaiju_db)
// TODO check integration with fastqc + multiqc and optional make contigs assembly version. and fastp
    // seq_reads_pe_ch.view()
        TRIM_PE(seq_reads_pe_ch)
            // fastp_report(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))
            fastp_report(TRIM_PE.out.trimmed_reads)
        // EXTRACT_BACTERIA_PE(TRIM_PE.out.trimmed_reads)
            MAP2REF_S_SWF_CORAL_PE(TRIM_PE.out.trimmed_reads, ref_host_ch)
                MAP2REF_S_SWF_SYM_PE(MAP2REF_S_SWF_CORAL_PE.out.non_mapped_reads, ref_sym_ch)        

                fastp_parse_ch = fastp_report.out.report_json.concat(
                // EXTRACT_BACTERIA_PE.out.report_json,
                MAP2REF_S_SWF_CORAL_PE.out.report_json,
                MAP2REF_S_SWF_SYM_PE.out.report_json)
                    // fastp_parse_ch.collect{it[1]}.toList().view()
                names = fastp_parse_ch.collect{it[0]}.toList()   
                paths = fastp_parse_ch.collect{it[1]}.toList()
                names.merge(paths).set { fastp_parse_ch_2 }
                    fastp_plot(fastp_parse_ch_2)
}

//TODO : WORKING AS IS. TO BE optimised. 
workflow EXTRACT_CORAL_SYMB_BWA {
    // main:
        kaiju_nodes = file(params.nodes)
        kaiju_db = file(params.kaiju_db)
// TODO check integration with fastqc + multiqc and optional make contigs assembly version. and fastp
    // seq_reads_pe_ch.view()
        TRIM_PE(seq_reads_pe_ch)
            // fastp_report(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))
            fastp_report(TRIM_PE.out.trimmed_reads)
        // EXTRACT_BACTERIA_PE(TRIM_PE.out.trimmed_reads)
            MAP2REF_BWA_PE_CORAL(TRIM_PE.out.trimmed_reads, ref_host_ch)
                MAP2REF_BWA_PE_SYM(MAP2REF_BWA_PE_CORAL.out.non_mapped_reads, ref_sym_ch)        

                fastp_parse_ch = fastp_report.out.report_json.concat(
                // EXTRACT_BACTERIA_PE.out.report_json,
                MAP2REF_BWA_PE_CORAL.out.report_json,
                MAP2REF_BWA_PE_SYM.out.report_json)
                    // fastp_parse_ch.collect{it[1]}.toList().view()
                names = fastp_parse_ch.collect{it[0]}.toList()   
                paths = fastp_parse_ch.collect{it[1]}.toList()
                names.merge(paths).set { fastp_parse_ch_2 }
                    fastp_plot(fastp_parse_ch_2)
}

//TODO : WORKING AS IS. TO BE optimised. 
workflow KAIJU_CONTIGS_PE {
// TODO check integration with fastqc + multiqc and optional make contigs assembly version. and fastp
    // seq_reads_pe_ch.view()
        params.mode == 'contigs'

        TRIM_PE(seq_reads_pe_ch)
            // fastp_report(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))
            fastp_report(TRIM_PE.out.trimmed_reads)
        megahit_pe(TRIM_PE.out.trimmed_reads)
        EXTRACT_BACTERIA_FA(megahit_pe.out.contigs_fa)
            MAP2REF_SWF_CORAL(EXTRACT_BACTERIA_FA.out.non_bacteria_reads, ref_host_ch)
                MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        
                
                fastp_parse_ch = fastp_report.out.report_json.concat(
                EXTRACT_BACTERIA_FA.out.report_json,
                MAP2REF_SWF_CORAL.out.report_json,
                MAP2REF_SWF_SYM.out.report_json)
                    // fastp_parse_ch.collect{it[1]}.toList().view()
                    // fastp_parse(fastp_parse_ch.collect{it[1]}.toList())
                names = fastp_parse_ch.collect{it[0]}.toList() 
                paths = fastp_parse_ch.collect{it[1]}.toList()
                names.merge(paths).set { fastp_parse_ch_2 }
                    fastp_plot(fastp_parse_ch_2)
}

workflow ONT_ANNA_SYM_CORAL {

// TODO check integration with fastqc + multiqc and optional make contigs assembly version. and fastp
    // seq_reads_pe_ch.view()
        params.mode == 'ONT'
        // TRIM_PE(seq_reads_pe_ch)
            // fastp_report(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))
            fastp_report(seq_reads_ch_ont)
        // megahit_pe(TRIM_PE.out.trimmed_reads)
        // EXTRACT_BACTERIA_FA(seq_reads_ch)
            MAP2REF_SWF_CORAL(seq_reads_ch_ont, ref_host_ch)
                MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        
                fastp_parse_ch = fastp_report.out.report_json.concat(
                // EXTRACT_BACTERIA_FA.out.report_json,
                MAP2REF_SWF_CORAL.out.report_json,
                MAP2REF_SWF_SYM.out.report_json)
                    // fastp_parse_ch.collect{it[1]}.toList().view()
                    // fastp_parse(fastp_parse_ch.collect{it[1]}.toList())
                names = fastp_parse_ch.collect{it[0]}.toList()   
                paths = fastp_parse_ch.collect{it[1]}.toList()
                names.merge(paths).set { fastp_parse_ch_2 }
                    fastp_plot(fastp_parse_ch_2)
}

workflow {
    // TODO: add R or python plots (long term add)
    KAIJU_PE_ILLUMINA()
    KRBR_PE_ILLUMINA()

}
//TODO: cleanup and add same as above for ONT reads.
// workflow {
// //     fasta2fastq(seq_contigs_ch)
// //         fastp_report(seq_reads_pe_ch.out.converted_fasta)
//     // seq_reads_pe_ch.view()
//     TRIM_PE(seq_reads_pe_ch)
//         EXTRACT_BACTERIA_PE(TRIM_PE.out.trimmed_reads)
//             MAP2REF_SWF_CORAL(EXTRACT_BACTERIA_PE.out.non_bacteria_reads, ref_host_ch)
//                 MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        

//             fastp_parse_ch = TRIM_PE.out.report_json.concat(
//             EXTRACT_BACTERIA_PE.out.report_json,
//             MAP2REF_SWF_CORAL.out.report_json,
//             MAP2REF_SWF_SYM.out.report_json)
//                 // fastp_parse_ch.collect{it[1]}.toList().view()
//                 fastp_parse(fastp_parse_ch.collect{it[1]}.toList())
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

// workflow  {
    
//     fasta2fastq(seq_contigs_ch)
//         fastp_report(seq_reads_pe_ch.out.converted_fasta)

//     TRIM_NANOPORE(seq_reads_ch)
//         EXTRACT_BACTERIA(TRIM_NANOPORE.out.trimmed_reads.concat(fasta2fastq.out.converted_fasta))
//             MAP2REF_SWF_CORAL(EXTRACT_BACTERIA.out.non_bacteria_reads, ref_host_ch)
//                 MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        

//             fastp_parse_ch = TRIM_NANOPORE.out.report_json.concat(
//             fastp_report.out.report_json,
//             EXTRACT_BACTERIA.out.report_json,
//             MAP2REF_SWF_CORAL.out.report_json,
//             MAP2REF_SWF_SYM.out.report_json)
//                 // fastp_parse_ch.collect{it[1]}.toList().view()
//                 fastp_parse(fastp_parse_ch.collect{it[1]}.toList())
// }

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
// // }
// workflow.onComplete {    
//     log.info """\
//         Pipeline execution summary
//         ---------------------------
//         Completed at: ${workflow.complete}
//         Duration    : ${workflow.duration}
//         Success     : ${workflow.success}
//         workDir     : ${workflow.workDir}
//         exit status : ${workflow.exitStatus}
//         """
//         .stripIndent()

// 	log.info ( workflow.success ? "\nDone! --> $params.outdir\n" : "Oops .. something went wrong" )
// }