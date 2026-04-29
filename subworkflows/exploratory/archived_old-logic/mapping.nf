#!/usr/bin/env nextflow
//TODO: general cleanup and reestablishment workflows and subworkflows!!
//TODO: add stub
/*
/*
========================================================================================
    Set params for kaiju database and mapping (bowtie2 minimap2) index references files
========================================================================================
*/
// import kaiju subworkflow 
include { MAP2REF_SWF as MAP2REF_SWF_CORAL; MAP2REF_SWF as MAP2REF_SWF_SYM } from '../subworkflows/map2ref_extract.nf'
include { MAP2REF_SWF_PE as MAP2REF_SWF_CORAL_PE; MAP2REF_SWF_PE as MAP2REF_SWF_SYM_PE } from '../subworkflows/map2ref_extract.nf'
include { EXTRACT_BACTERIA_ONT; EXTRACT_BACTERIA_PE } from '../subworkflows/extract_bacteria.nf' //EXTRACT_BACTERIA_FA EXTRACT_BACTERIA
include { TRIM_PE } from '../subworkflows/trim_illumina.nf'

params.nodes = "/share/databases/kaiju/nodes.dmp"
params.kaiju_db = "/share/databases/kaiju/refseq/kaiju_db_refseq.fmi"
params.barplot_script = "/home/colinl/metaG/Git/metaG_EukDepletion/R-scripts/Bar_plot_kaiju.r"
params.parse_script = "/home/colinl/metaG/Git/metaG_EukDepletion/parse_json_cl_3f.py"
params.sample_metadata_sheet = "/home/colinl/metaG/Git/metaG_EukDepletion/sample_metadata_sheet.csv"

params.ref_host = "/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_bowtie2/all_scleractina" //"/home/colinl/metaG/Git/metaG_EukDepletion/kaiju/mem/split/ref_genomes/aiptasiidae.mmi" //"/home/colinl/metaG/Git/metaG_EukDepletion/kaiju/mem/split/ref_genomes/corals.mmi"
params.ref_sym = "/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_bowtie2/all_symbiodiniaceae"

/*
========================================================================================
    Set params for kaiju database and minimap2 index references files
========================================================================================
*/

// import qc modules and workflow
include { FASTP_REPORT as FASTP_REPORT_PE; FASTP_REPORT as FASTP_REPORT_ONT } from '../modules/fastp_stats.nf'
// include { FASTP_PARSE } from '../modules/FASTP_PARSE.nf'
include { FASTP_PLOT as FASTP_PLOT_PE; FASTP_PLOT_CONTIGS as FASTP_PLOT_ONT } from '../modules/fastp_parse.nf'
// include { fasta2fastq } from '../modules/ONT_prep.nf'

        Channel.fromPath(params.input)
            .splitCsv(header:true, quote: '\"')
            .map {row -> tuple(row.Sample_name, row.Sample_name, tuple((row.R1), (row.R2)), (row.Seqencing_type), (row.Entry))}
            .filter { it[3] == 'Illumina' }
            .filter { it[4] == 'TEST_MAPPING_1' }
            .map { tuple(it[0], it[1], it[2], it[3]) }
            .set{ seq_reads_pe_ch }

        Channel.fromPath(params.input)
            .splitCsv(header:true, quote: '\"')
            .map {row -> tuple(row.Sample_name, row.Sample_name, (row.R), (row.Seqencing_type), (row.Entry))}
            .filter { it[3] == 'ONT' }
            .filter { it[4] == 'TEST_MAPPING_1' }
            .map { tuple(it[0], it[1], it[2], it[3]) }
            .set{ seq_reads_ont_ch }

        Channel.fromPath(params.input)
            .splitCsv(header:true, quote: '\"')
            .map {row -> tuple(row.Sample_name, (row.reference_host_path))}
            .map {tuple (it[0], it[1].split('/')[-1], it[1].split('.mmi')[0])}
            .set{ch_host}

        Channel.fromPath(params.input)
            .splitCsv(header:true, quote: '\"')
            .map {row -> tuple(row.Sample_name, (row.reference_Sym_path))}
            .map {tuple (it[0], it[1].split('/')[-1], it[1].split('.mmi')[0])}
            .set{ch_sym}


workflow TEST_PARSE_PE {
    // take: 
    //     seq_reads_pe_ch
    //     ch_host
    //     ch_sym
    // main:
            kaiju_nodes = file(params.nodes)
            kaiju_db = file(params.kaiju_db)

            // TRIM_PE(seq_reads_pe_ch)// provide trimmed reads directly from main workflow. 
            MAP2REF_SWF_CORAL_PE(seq_reads_pe_ch.join(ch_host))
                MAP2REF_SWF_SYM_PE(MAP2REF_SWF_CORAL_PE.out.non_mapped_reads.join(ch_sym))
                    EXTRACT_BACTERIA_PE(MAP2REF_SWF_SYM_PE.out.non_mapped_reads)

                FASTP_REPORT_PE(seq_reads_pe_ch)
                fastp_parse_pe_ch = FASTP_REPORT_PE.out.report_json.concat(
                    EXTRACT_BACTERIA_PE.out.report_json,
                    MAP2REF_SWF_CORAL_PE.out.report_json,
                    MAP2REF_SWF_SYM_PE.out.report_json)
                
                // FASTP_REPORT_PE.out.report_json.view()
                    names_pe = fastp_parse_pe_ch.collect{it[1]}.toList()   
                    paths_pe = fastp_parse_pe_ch.collect{it[2]}.toList()
                    seq_typ_pe = fastp_parse_pe_ch.collect{it[3]}.toList()
                    names_pe.merge(paths_pe, seq_typ_pe).set { fastp_parse_ch_pe }
                fastp_parse_ch_pe.view()
                    // FASTP_PLOT_PE(fastp_parse_ch_pe)
}

workflow MAPPING_CONSECUTIVE_PE {
    take: 
        seq_reads_pe_ch
        ch_host
        ch_sym
        untrimmed_reads_pe_ch //add channel for untrimmed reads fastp json files from TRIM_PE(?)
    main:
            kaiju_nodes = file(params.nodes)
            kaiju_db = file(params.kaiju_db)

            // TRIM_PE(seq_reads_pe_ch)// provide trimmed reads directly from main workflow. 
            MAP2REF_SWF_CORAL_PE(seq_reads_pe_ch.join(ch_host))
                MAP2REF_SWF_SYM_PE(MAP2REF_SWF_CORAL_PE.out.non_mapped_reads.join(ch_sym))
                    EXTRACT_BACTERIA_PE(MAP2REF_SWF_SYM_PE.out.non_mapped_reads)

                FASTP_REPORT_PE(seq_reads_pe_ch.concat(untrimmed_reads_pe_ch))
                fastp_parse_pe_ch = FASTP_REPORT_PE.out.report_json.concat(
                    EXTRACT_BACTERIA_PE.out.report_json,
                    MAP2REF_SWF_CORAL_PE.out.report_json,
                    MAP2REF_SWF_SYM_PE.out.report_json)

                    names_pe = fastp_parse_pe_ch.collect{it[1]}.toList()   
                    paths_pe = fastp_parse_pe_ch.collect{it[2]}.toList()
                    seq_typ_pe = fastp_parse_pe_ch.collect{it[3]}.toList()
                    names_pe.merge(paths_pe, seq_typ_pe).set { fastp_parse_ch_pe }
                    FASTP_PLOT_PE(fastp_parse_ch_pe)
}

workflow MAPPING_CONSECUTIVE {
    take: 
        seq_reads_ont_ch
        ch_host
        ch_sym
    main:
            kaiju_nodes = file(params.nodes)
            kaiju_db = file(params.kaiju_db)

            MAP2REF_SWF_CORAL(seq_reads_ont_ch.join(ch_host))
                MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads.join(ch_sym))        
                    EXTRACT_BACTERIA_ONT(MAP2REF_SWF_SYM.out.non_mapped_reads)
                FASTP_REPORT_ONT(seq_reads_ont_ch)
                fastp_parse_ont_ch = FASTP_REPORT_ONT.out.report_json.concat(
                    MAP2REF_SWF_CORAL.out.report_json,
                    MAP2REF_SWF_SYM.out.report_json,
                    EXTRACT_BACTERIA_ONT.out.report_json)

                    names_ont = fastp_parse_ont_ch.collect{it[1]}.toList()   
                    paths_ont = fastp_parse_ont_ch.collect{it[2]}.toList()
                    seq_typ_ont = fastp_parse_ont_ch.collect{it[3]}.toList()
                    names_ont.merge(paths_ont, seq_typ_ont).set { fastp_plot_ont_ch }
                    FASTP_PLOT_ONT(fastp_plot_ont_ch)
}

workflow MAPPING_PARALLEL_PE {
    take: 
        seq_reads_pe_ch
        ch_host
        ch_sym

    main:
            kaiju_nodes = file(params.nodes)
            kaiju_db = file(params.kaiju_db)

            // TRIM_PE(seq_reads_pe_ch) // provide trimmed reads directly from main workflow. 
                FASTP_REPORT_PE(seq_reads_pe_ch)

            MAP2REF_SWF_CORAL_PE(seq_reads_pe_ch.join(ch_host))
            
            MAP2REF_SWF_SYM_PE(seq_reads_pe_ch.join(ch_sym))        
            
            EXTRACT_BACTERIA_PE(seq_reads_pe_ch)

            fastp_parse_pe_ch = FASTP_REPORT_PE.out.report_json.concat(
                    MAP2REF_SWF_CORAL_PE.out.report_json,
                    MAP2REF_SWF_SYM_PE.out.report_json,
                    EXTRACT_BACTERIA_PE.out.report_json)

                    names_pe = fastp_parse_pe_ch.collect{it[1]}.toList()   
                    paths_pe = fastp_parse_pe_ch.collect{it[2]}.toList()
                    seq_typ_pe = fastp_parse_pe_ch.collect{it[3]}.toList()
                    names_pe.merge(paths_pe, seq_typ_pe).set { fastp_parse_ch_pe }
                    FASTP_PLOT_PE(fastp_parse_ch_pe)
}

workflow MAPPING_PARALLEL {
    take: 
        seq_reads_ont_ch
        ch_host
        ch_sym

    main:
            kaiju_nodes = file(params.nodes)
            kaiju_db = file(params.kaiju_db)

            MAP2REF_SWF_CORAL(seq_reads_ont_ch.join(ch_host))
            
            MAP2REF_SWF_SYM(seq_reads_ont_ch.join(ch_sym))        
            
            EXTRACT_BACTERIA_ONT(seq_reads_ont_ch)

                FASTP_REPORT_ONT(seq_reads_ont_ch)
                fastp_parse_ont_ch = FASTP_REPORT_ONT.out.report_json.concat(
                    MAP2REF_SWF_CORAL.out.report_json,
                    MAP2REF_SWF_SYM.out.report_json,
                    EXTRACT_BACTERIA_ONT.out.report_json)

                    names_ont = fastp_parse_ont_ch.collect{it[1]}.toList()   
                    paths_ont = fastp_parse_ont_ch.collect{it[2]}.toList()
                    seq_typ_ont = fastp_parse_ont_ch.collect{it[3]}.toList()
                    names_ont.merge(paths_ont, seq_typ_ont).set { fastp_plot_ont_ch }
                    FASTP_PLOT_ONT(fastp_plot_ont_ch)
}

// workflow MAPPING_PARALLEL_ORIGINAL {
//     take: 
//         seq_reads_pe_ch
//         seq_reads_ont_ch
//         ch_host
//         ch_sym

//     main:
//             kaiju_nodes = file(params.nodes)
//             kaiju_db = file(params.kaiju_db)

//             TRIM_PE(seq_reads_pe_ch)
//                 FASTP_REPORT_PE(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))

//             MAP2REF_SWF_CORAL_PE(TRIM_PE.out.trimmed_reads.join(ch_host))
            
//             MAP2REF_SWF_SYM_PE(TRIM_PE.out.trimmed_reads.join(ch_sym))        
            
//             EXTRACT_BACTERIA_PE(TRIM_PE.out.trimmed_reads)

//             fastp_parse_pe_ch = FASTP_REPORT_PE.out.report_json.concat(
//                     MAP2REF_SWF_CORAL_PE.out.report_json,
//                     MAP2REF_SWF_SYM_PE.out.report_json,
//                     EXTRACT_BACTERIA_PE.out.report_json)

//                     names_pe = fastp_parse_pe_ch.collect{it[0]}.toList()   
//                     paths_pe = fastp_parse_pe_ch.collect{it[1]}.toList()
//                     names_pe.merge(paths_pe).set { fastp_parse_ch_pe }

//             MAP2REF_SWF_CORAL(seq_reads_ont_ch.join(ch_host))
            
//             MAP2REF_SWF_SYM(seq_reads_ont_ch.join(ch_sym))        
            
//             EXTRACT_BACTERIA_ONT(seq_reads_ont_ch)

//                 FASTP_REPORT_ONT(seq_reads_ont_ch)
//                 fastp_parse_ont_ch = FASTP_REPORT_ONT.out.report_json.concat(
//                     MAP2REF_SWF_CORAL.out.report_json,
//                     MAP2REF_SWF_SYM.out.report_json,
//                     EXTRACT_BACTERIA_ONT.out.report_json)

//                     names_ont = fastp_parse_ont_ch.collect{it[0]}.toList()   
//                     paths_ont = fastp_parse_ont_ch.collect{it[1]}.toList()
//                     names_ont.merge(paths_ont).set { fastp_parse_ch_ont }

//                     FASTP_PLOT_PE(fastp_parse_ch_pe)
//                     FASTP_PLOT_ONT(fastp_parse_ch_ont)
// }

// workflow MAPPING_CONSECUTIVe ORIGINAL {
//     take: 
//         seq_reads_pe_ch
//         seq_reads_ont_ch
//         ch_host
//         ch_sym

//     main:
//             kaiju_nodes = file(params.nodes)
//             kaiju_db = file(params.kaiju_db)

//             TRIM_PE(seq_reads_pe_ch)
//             MAP2REF_SWF_CORAL_PE(TRIM_PE.out.trimmed_reads.join(ch_host))
//                 MAP2REF_SWF_SYM_PE(MAP2REF_SWF_CORAL_PE.out.non_mapped_reads.join(ch_sym))
//                     EXTRACT_BACTERIA_PE(MAP2REF_SWF_SYM_PE.out.non_mapped_reads)

//                 FASTP_REPORT_PE(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))
//                 fastp_parse_pe_ch = FASTP_REPORT_PE.out.report_json.concat(
//                     EXTRACT_BACTERIA_PE.out.report_json,
//                     MAP2REF_SWF_CORAL_PE.out.report_json,
//                     MAP2REF_SWF_SYM_PE.out.report_json)

//                     names_pe = fastp_parse_pe_ch.collect{it[0]}.toList()   
//                     paths_pe = fastp_parse_pe_ch.collect{it[1]}.toList()
//                     names_pe.merge(paths_pe).set { fastp_parse_ch_pe }
//                     FASTP_PLOT_PE(fastp_parse_ch_pe)

//             MAP2REF_SWF_CORAL(seq_reads_ont_ch.join(ch_host))
//                 MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads.join(ch_sym))        
//                     EXTRACT_BACTERIA_ONT(MAP2REF_SWF_SYM.out.non_mapped_reads)
//                 FASTP_REPORT_ONT(seq_reads_ont_ch)
//                 fastp_parse_ont_ch = FASTP_REPORT_ONT.out.report_json.concat(
//                     MAP2REF_SWF_CORAL.out.report_json,
//                     MAP2REF_SWF_SYM.out.report_json,
//                     EXTRACT_BACTERIA_ONT.out.report_json)

//                     names_ont = fastp_parse_ont_ch.collect{it[0]}.toList()   
//                     paths_ont = fastp_parse_ont_ch.collect{it[1]}.toList()
//                     names_ont.merge(paths_ont).set { fastp_parse_ch_ont }
//                     FASTP_PLOT_ONT(fastp_parse_ch_ont)
// }