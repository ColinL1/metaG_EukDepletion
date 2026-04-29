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
include { EXTRACT_BACTERIA; EXTRACT_BACTERIA_ONT; EXTRACT_BACTERIA_FA; EXTRACT_BACTERIA_PE } from '../subworkflows/extract_bacteria.nf'
include { TRIM_PE } from '../subworkflows/trim_illumina.nf'

params.nodes = "/share/databases/kaiju/nodes.dmp"
params.kaiju_db = "/share/databases/kaiju/refseq/kaiju_db_refseq.fmi"
params.barplot_script = "/home/colinl/metaG/Git/metaG_EukDepletion/R-scripts/Bar_plot_kaiju.r"
// params.parse_script = "/home/colinl/metaG/Git/metaG_EukDepletion/parse_json_cl.py"
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
include { fastp_report as fastp_report1; fastp_report as fastp_report2 } from '../modules/fastp_stats.nf'
// include { FASTP_PARSE } from '../modules/FASTP_PARSE.nf'
include { FASTP_PLOT as fastp_plot_pe; FASTP_PLOT as fastp_plot_ont } from '../modules/FASTP_PARSE.nf'
// include { fasta2fastq } from '../modules/ONT_prep.nf'

workflow TEST_MAPPING_2_X {
    take: 
        input_file

    main:
        Channel.fromPath(input_file)
            .splitCsv(header:true, quote: '\"')
            .map {row -> tuple(row.Sample_name, row.Sample_name, tuple((row.R1), (row.R2)), (row.Seqencing_type), (row.Entry))}
            .filter { it[3] == 'Illumina' }
            .filter { it[4] == 'TEST_MAPPING_2' }
            .map { tuple(it[0], it[1], it[2], it[3]) }
            .set{ seq_reads_pe_ch }

        Channel.fromPath(input_file)
            .splitCsv(header:true, quote: '\"')
            .map {row -> tuple(row.Sample_name, row.Sample_name, (row.R), (row.Seqencing_type), (row.Entry))}
            .filter { it[3] == 'ONT' }
            .filter { it[4] == 'TEST_MAPPING_2' }
            .map { tuple(it[0], it[1], it[2], it[3]) }
            .set{ seq_reads_ont_ch }

        Channel.fromPath(input_file)
            .splitCsv(header:true, quote: '\"')
            .map {row -> tuple(row.Sample_name, (row.reference_host_path))}
            .map {tuple (it[0], it[1].split('/')[-1], it[1])}
            .set{ch_host}
            
        Channel.fromPath(input_file)
            .splitCsv(header:true, quote: '\"')
            .map {row -> tuple(row.Sample_name, (row.reference_Sym_path))}
            .map {tuple (it[0], it[1].split('/')[-1], it[1])}
            .set{ch_sym}

            kaiju_nodes = file(params.nodes)
            kaiju_db = file(params.kaiju_db)

            TRIM_PE(seq_reads_pe_ch)
                fastp_report1(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))

            MAP2REF_SWF_CORAL_PE(TRIM_PE.out.trimmed_reads.join(ch_host))
            
            MAP2REF_SWF_SYM_PE(TRIM_PE.out.trimmed_reads.join(ch_sym))        
            
            EXTRACT_BACTERIA_PE(TRIM_PE.out.trimmed_reads)

            fastp_parse_ch1 = fastp_report1.out.report_json.concat(
                    MAP2REF_SWF_CORAL_PE.out.report_json,
                    MAP2REF_SWF_SYM_PE.out.report_json,
                    EXTRACT_BACTERIA_PE.out.report_json)

                    names_pe = fastp_parse_ch1.collect{it[0]}.toList()   
                    paths_pe = fastp_parse_ch1.collect{it[1]}.toList()
                    names_pe.merge(paths_pe).set { fastp_parse_ch_pe }

            MAP2REF_SWF_CORAL(seq_reads_ont_ch.join(ch_host))
            
            MAP2REF_SWF_SYM(seq_reads_ont_ch.join(ch_sym))        
            
            EXTRACT_BACTERIA_ONT(seq_reads_ont_ch)

                fastp_report2(seq_reads_ont_ch)
                fastp_parse_ch2 = fastp_report2.out.report_json.concat(
                    MAP2REF_SWF_CORAL.out.report_json,
                    MAP2REF_SWF_SYM.out.report_json,
                    EXTRACT_BACTERIA_ONT.out.report_json)

                    names_ont = fastp_parse_ch2.collect{it[0]}.toList()   
                    paths_ont = fastp_parse_ch2.collect{it[1]}.toList()
                    names_ont.merge(paths_ont).set { fastp_parse_ch_ont }

                    fastp_plot_pe(fastp_parse_ch_pe)
                    fastp_plot_ont(fastp_parse_ch_ont)
}
