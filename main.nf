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

// Channel.fromFilePairs(params.input).set { seq_reads_pe_ch }
// Channel.fromPath(params.input).map {tuple( it.name.split('.fastq.gz')[0], it )}.set { seq_reads_ont_ch }

/*
========================================================================================
    Set params for kaiju database and minimap2 index references files
========================================================================================
*/
// import kaiju subworkflow 
include { MAP2REF_SWF as MAP2REF_SWF_CORAL; MAP2REF_SWF as MAP2REF_SWF_SYM } from './subworkflows/map2ref_extract.nf'
include { MAP2REF_SWF_PE as MAP2REF_SWF_CORAL_PE; MAP2REF_SWF_PE as MAP2REF_SWF_SYM_PE } from './subworkflows/map2ref_extract.nf'
include { EXTRACT_BACTERIA_ONT; EXTRACT_BACTERIA_PE } from './subworkflows/extract_bacteria.nf'// EXTRACT_BACTERIA; EXTRACT_BACTERIA_ONT; EXTRACT_BACTERIA_FA; EXTRACT_BACTERIA_PE
include { TRIM_PE } from './subworkflows/trim_illumina.nf'

params.nodes = "/share/databases/kaiju/nodes.dmp"
params.kaiju_db = "/share/databases/kaiju/refseq/kaiju_db_refseq.fmi"
params.barplot_script = "/home/colinl/metaG/Git/metaG_EukDepletion/R-scripts/Bar_plot_kaiju.r"
// params.parse_script = "/home/colinl/metaG/Git/metaG_EukDepletion/parse_json_cl.py"
params.parse_script = "/home/colinl/metaG/Git/metaG_EukDepletion/parse_json_cl_3f.py"
params.sample_metadata_sheet = "/home/colinl/metaG/Git/metaG_EukDepletion/sample_metadata_sheet.csv"

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
// include { FASTP_REPORT as fastp_report1; FASTP_REPORT as fastp_report2 } from './modules/fastp_stats.nf'
include { FASTP_PARSE } from './modules/fastp_parse.nf'
include { FASTP_PLOT as FASTP_PLOT_PE; FASTP_PLOT as FASTP_PLOT_ONT } from './modules/fastp_parse.nf'
include { fasta2fastq } from './modules/ONT_prep.nf'

// TODO import module for assembly variant (currently unused.) 
include { MEGAHIT_SE; MEGAHIT_PE } from './modules/assembly_megahit.nf'

// TODO check for importance of this for ONT workflow
// Channel.fromPath(params.input+'/*.fastq.gz').map {tuple( it.name.split('.fastq.gz')[0], it )}.set { seq_reads_ch }
// // Channel.fromPath(params.contigs).map {tuple( (it.name.split('.fasta')[0]), it )}.set { seq_contigs_ch }
// Channel.fromPath(params.input+'/*.fasta').map {tuple( (it.name.split('.fasta')[0]).split('_')[0..3].join("_"), it )}.set { seq_contigs_ch } // version that cleans the names (specific to my use case)

// include { MAP2REF_BWA_PE  } from './subworkflows/extract_bacteria.nf'
// include { MAP2REF_BWA_PE as MAP2REF_BWA_PE_CORAL; MAP2REF_BWA_PE as MAP2REF_BWA_PE_SYM } from './subworkflows/map2ref_extract.nf'
// include { MAP2REF_S_SWF_PE as MAP2REF_S_SWF_CORAL_PE; MAP2REF_S_SWF_PE as MAP2REF_S_SWF_SYM_PE } from './subworkflows/map2ref_extract.nf'

// workflow test_ch {

//     // TODO: find a way to use the filter function
//         Channel.fromPath(params.input)
//         .splitCsv(header:true, quote: '\"')
//         // Sample_name,R1,R2,R,Seqencing_type,Entry,reference_host_path,reference_Sym_path
//         .map {row -> tuple(row.Sample_name, row.Sample_name, tuple((row.R1), (row.R2)), (row.Seqencing_type))}
//         .filter { it[3] == 'Illumina' }
//         .set{ seq_reads_pe_ch }

//         Channel.fromPath(params.input)
//         .splitCsv(header:true, quote: '\"')
//         // Sample_name,R1,R2,R,Seqencing_type,Entry,reference_host_path,reference_Sym_path
//         .map {row -> tuple(row.Sample_name, (row.reference_host_path))}
//         .map {tuple (it[0], it[1].split('/')[-1], it[1])}
//         .set{ch_host}
        
//         Channel.fromPath(params.input)
//         .splitCsv(header:true, quote: '\"')
//         // Sample_name,R1,R2,R,Seqencing_type,Entry,reference_host_path,reference_Sym_path
//         .map {row -> tuple(row.Sample_name, (row.reference_Sym_path))}
//         .map {tuple (it[0], it[1].split('/')[-1], it[1])}
//         .set{ch_sym}

//         Channel.fromPath(params.input)
//         .splitCsv(header:true, quote: '\"')
//         // Sample_name,R1,R2,R,Seqencing_type,Entry,reference_host_path,reference_Sym_path
//         .map {row -> tuple(row.Sample_name, row.Sample_name, (row.R), (row.Seqencing_type))}
//         // .filter( it[3] == 'ONT' )
//         .filter { it[3] == 'ONT' }
//         .set{ seq_reads_ont_ch }
// }

// TODO: add switch for kaiju vs kraken bracken
// log.info """\

//     metaG kingdom taxonomy - NF   PIPELINE
//     =======================================
//     kaiju_db   : ${params.kaiju_db}
//     kraken2    : ${params.kraken_db}
//     reads      : ${params.input}
//     outdir     : ${params.outdir}

//     """


include { NCBI_DOWNLOAD } from './modules/NCBI_Datasets.nf' 
include { MINIMAP2_INDEX } from './modules/minimap2_index.nf' 

workflow test_dataset_creation {
    Channel.fromPath(params.input)
        .splitCsv(header:true, quote: '\"')
        .map {row -> tuple(row.Sample_name, row.reference_host_path, row.Seqencing_type, row.mapping_arg)}
        .unique{it[1]}
        .map {tuple (it[0], it[0], it[1], it[2], it[3])}
        .set{ch_host}

        if (ch_host)
            ch_host.view()
            ch_host.filter{ it[4] == "reference"}.set{ download_ncbi_ch }
            if (download_ncbi_ch)

        if (!ch_host)
            print("not_existing")
    // NCBI_DOWNLOAD()
    // MINIMAP2_INDEX()
}

workflow test_ch_input_GOOD { //THIS works as intended
        

//        def file = new File(params.input)
        // def rows = file.readLines().tail()*.split(',')
        // int total = rows.size()
        // Set type = rows.collect { it[5] + ' ' + it[4] }
        
  //      file.withReader { r ->
      //  def rows = []
    //    def reader = new CSVReaderHeaderAware(r)
   //     while ((next = reader.readMap())) rows << next
     //   Set type = rows.collect { it.Seqencing_type + ' ' + it.Entry }
   // }
        
     //   println(type)
        
    //     Channel.fromPath(params.input)
    //         .splitCsv(header:true, quote: '\"')
    //         .map {row -> tuple((row.Entry))}
    //         .filter { it[0] == 'TEST_MAPPING_2' }
    //         .ifEmpty (entry_point = false)

    // if ( entry_point == true) {
    //     println("something else") 
    //  } else {
    //     print("test working")
    //         // Channel.fromPath(params.input)
    //         //     .splitCsv(header:true, quote: '\"')
    //         //     .map{row -> tuple(row.Sample_name, row.Sample_name, tuple((row.R1), (row.R2)), (row.Seqencing_type), (row.Entry))}
    //         //     .filter{ it[3] == 'Illumina' }
    //         //     .filter{ it[4] == 'TEST_MAPPING_1' }
    //         //     .map{ tuple(it[0], it[1], it[2], it[3]) }
    //         //     .set{ seq_reads_pe_ch }

    //         // Channel.fromPath(params.input)
    //         //     .splitCsv(header:true, quote: '\"')
    //         //     .map {row -> tuple(row.Sample_name, row.Sample_name, (row.R), (row.Seqencing_type), (row.Entry))}
    //         //     .filter { it[3] == 'ONT' }
    //         //     .filter { it[4] == 'TEST_MAPPING_1' }
    //         //     .map { tuple(it[0], it[1], it[2], it[3]) }
    //         //     .set{ seq_reads_ont_ch }
    //     Channel.fromPath(params.input)
    //         .splitCsv(header:true, quote: '\"')
    //         .map {row -> tuple((row.Seqencing_type))}
    //         .filter { it[0] == 'Illumina' }
    //         .ifEmpty (pe_illumina_entry = false)

    //     Channel.fromPath(params.input)
    //         .splitCsv(header:true, quote: '\"')
    //         .map {row -> tuple((row.Seqencing_type))}
    //         .filter { it[0] == 'ONT' }
    //         .ifEmpty (ont_entry = false)
            
    //     // if(seq_reads_ont_ch && seq_reads_pe_ch) {
    //     //     println("both exist")
    //     // } else 
    //     if(pe_illumina_entry == true) {
    //         println("ONT not present ==> do illumina only ")
    //         // tst1.view()
    //         print(pe_illumina_entry)

    //     }
    //     if(ont_entry == true) {
    //         println("PE illumina not present ==>  do ONT only ")
    //         // tst.view()
    //         print(ont_entry)

    //     } 
    //     else { println("do both?")
    //     }
    //     }
        //     println(!(seq_reads_ont_ch && seq_reads_pe_ch))
        // //     print("no input!!")
            
        //     println( (seq_reads_ont_ch && !seq_reads_pe_ch))
        // //     print("RUN ONLY ONT")
            
        //     println((!seq_reads_ont_ch))
        //     println((!seq_reads_pe_ch))
        // //     print("RUN ONLY PE")
            
        //     println((seq_reads_ont_ch && seq_reads_pe_ch))
        //     print("RUN both")
        // seq_reads_ont_ch.view()
        // seq_reads_pe_ch.view()
}



// empty_ch = Channel.empty()
        Channel.fromPath(params.input)
            .splitCsv(header:true, quote: '\"')
            .map {row -> tuple(row.Sample_name, row.Sample_name, tuple((row.R1), (row.R2)), (row.Seqencing_type), (row.Entry))}
            .filter { it[3] == 'Illumina' }
            .filter { it[4] == 'TEST_MAPPING_1' }
            .map { tuple(it[0], it[1], it[2], it[3]) }
            // .ifEmpty('empty')
            .set{ seq_reads_pe_ch }

        Channel.fromPath(params.input)
            .splitCsv(header:true, quote: '\"')
            .map {row -> tuple(row.Sample_name, row.Sample_name, (row.R), (row.Seqencing_type), (row.Entry))}
            .filter { it[3] == 'ONT' }
            .filter { it[4] == 'TEST_MAPPING_1' }
            .map { tuple(it[0], it[1], it[2], it[3]) }
            .ifEmpty("empty")
            .set{ seq_reads_ont_ch }

        Channel.fromPath(params.input)
            .splitCsv(header:true, quote: '\"')
            .map {row -> tuple(row.Sample_name, (row.reference_host_path))}
            .map {tuple (it[0], it[1].split('/')[-1].split('.mmi')[0], it[1])}
            .set{ch_host}

        Channel.fromPath(params.input)
            .splitCsv(header:true, quote: '\"')
            .map {row -> tuple(row.Sample_name, (row.reference_Sym_path))}
            .map {tuple (it[0], it[1].split('/')[-1].split('.mmi')[0], it[1])}
            .set{ch_sym}

// import mapping subworkflow 
include { MAPPING_CONSECUTIVE_PE; MAPPING_PARALLEL_PE  } from './subworkflows/mapping.nf'
include { MAPPING_CONSECUTIVE; MAPPING_PARALLEL} from './subworkflows/mapping.nf'
include { MAPPING_CONSECUTIVE as MAPPING_CONSECUTIVE_CONTIGS_PE ; MAPPING_PARALLEL as MAPPING_PARALLEL_CONTIGS_PE } from './subworkflows/mapping.nf'
include { MAPPING_CONSECUTIVE as MAPPING_CONSECUTIVE_CONTIGS_ONT ; MAPPING_PARALLEL as MAPPING_PARALLEL_CONTIGS_ONT } from './subworkflows/mapping.nf'

workflow {
    //map PE reads in parallel and hierarchically.
    TRIM_PE(seq_reads_pe_ch)
        MAPPING_CONSECUTIVE_PE(TRIM_PE.out.trimmed_reads, ch_host, ch_sym)
        MAPPING_PARALLEL_PE(TRIM_PE.out.trimmed_reads, ch_host, ch_sym)
        
    //map PE contigs in parallel and hierarchically.
    MEGAHIT_PE(TRIM_PE.out.trimmed_reads)
        MAPPING_CONSECUTIVE_CONTIGS_PE(MEGAHIT_PE.out, ch_host, ch_sym)
        MAPPING_PARALLEL_CONTIGS_PE(MEGAHIT_PE.out, ch_host, ch_sym)

    //map ONT contigs in parallel and hierarchically.
    MEGAHIT_SE(seq_reads_ont_ch)
                MAPPING_CONSECUTIVE_CONTIGS_ONT(MEGAHIT_SE.out, ch_host, ch_sym)
                MAPPING_PARALLEL_CONTIGS_ONT(MEGAHIT_SE.out, ch_host, ch_sym)
    
    //map ONT reads in parallel and hierarchically.
    MAPPING_CONSECUTIVE(seq_reads_ont_ch, ch_host, ch_sym)
    MAPPING_PARALLEL(seq_reads_ont_ch, ch_host, ch_sym)
}

include {FASTP_REPORT} from './modules/fastp_stats.nf'
//TODO :  TO BE optimised. (adding database mv to shared memory)
workflow KRBR_PE{
    TRIM_PE(seq_reads_pe_ch)
        FASTP_REPORT(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))
            KRAKEN_BRACKEN_PE(TRIM_PE.out.trimmed_reads)
            //TODO add qc
}
workflow KRBR_SE{
        FASTP_REPORT(seq_reads_ont_ch)
            KRAKEN_BRACKEN(seq_reads_ont_ch, FASTP_REPORT.out.report_json)
            //TODO add qc
}

workflow TestWorkflow {
    // Print the output of seq_reads_pe_ch
    seq_reads_pe_ch.subscribe { sample ->
        println("seq_reads_pe_ch: $sample")
    }

    // // Print the output of seq_reads_ont_ch
    // seq_reads_ont_ch.subscribe { sample ->
    //     println("seq_reads_ont_ch: $sample")
    // }

    // Print the output of ch_host
    ch_host.subscribe { host ->
        println("ch_host: $host")
    }

    // Print the output of ch_sym
    ch_sym.subscribe { sym ->
        println("ch_sym: $sym")
    }
}

// // Run the TestWorkflow
// TestWorkflow.run()

        Channel.fromPath(params.input)
            .splitCsv(header:true, quote: '\"')
            .map {row -> tuple(row.Sample_name, (row.reference_host_path_cont))}
            .map {tuple (it[0], it[1].split('/')[-1].split('.mmi')[0], it[1])}
            .set{ch_host_contigs}

        Channel.fromPath(params.input)
            .splitCsv(header:true, quote: '\"')
            .map {row -> tuple(row.Sample_name, (row.reference_Sym_path_cont))}
            .map {tuple (it[0], it[1].split('/')[-1].split('.mmi')[0], it[1])}
            .set{ch_sym_contigs}


workflow PE_CONSECUTIVE {
    //map PE reads in parallel and hierarchically.
    TRIM_PE(seq_reads_pe_ch)
        MAPPING_CONSECUTIVE_PE(TRIM_PE.out.trimmed_reads, ch_host, ch_sym, seq_reads_pe_ch)
        // MAPPING_PARALLEL_PE(TRIM_PE.out.trimmed_reads, ch_host, ch_sym)
        
    //map PE contigs in parallel and hierarchically.
    MEGAHIT_PE(TRIM_PE.out.trimmed_reads)
        MAPPING_CONSECUTIVE_CONTIGS_PE(MEGAHIT_PE.out, ch_host_contigs, ch_sym_contigs)
        // MAPPING_PARALLEL_CONTIGS_PE(MEGAHIT_PE.out, ch_host, ch_sym)
    // //map ONT contigs in parallel and hierarchically.
    // MEGAHIT_SE(seq_reads_ont_ch)
    //             MAPPING_CONSECUTIVE_CONTIGS_ONT(MEGAHIT_SE.out, ch_host, ch_sym)
    //             MAPPING_PARALLEL_CONTIGS_ONT(MEGAHIT_SE.out, ch_host, ch_sym)
    
    // //map ONT reads in parallel and hierarchically.
    // MAPPING_CONSECUTIVE(seq_reads_ont_ch, ch_host, ch_sym)
    // MAPPING_PARALLEL(seq_reads_ont_ch, ch_host, ch_sym)
}



// workflow TEST_MAPPING_1 {
//         params.seq = "illumina"
//         params.type = "reads"

//         kaiju_nodes = file(params.nodes)
//         kaiju_db = file(params.kaiju_db)

//         TRIM_PE(seq_reads_pe_ch)
//             FASTP_REPORT(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))

//         MAP2REF_SWF_CORAL_PE(TRIM_PE.out.trimmed_reads, ref_host_ch)
//             MAP2REF_SWF_SYM_PE(MAP2REF_SWF_CORAL_PE.out.non_mapped_reads, ref_sym_ch)        
//                 EXTRACT_BACTERIA_PE(MAP2REF_SWF_SYM_PE.out.non_mapped_reads)

//                 fastp_parse_ch = FASTP_REPORT.out.report_json.concat(
//                 EXTRACT_BACTERIA_PE.out.report_json,
//                 MAP2REF_SWF_CORAL_PE.out.report_json,
//                 MAP2REF_SWF_SYM_PE.out.report_json)
//                     // fastp_parse_ch.collect{it[1]}.toList().view()
//                 names = fastp_parse_ch.collect{it[0]}.toList()   
//                 paths = fastp_parse_ch.collect{it[1]}.toList()
//                 names.merge(paths).set { fastp_parse_ch_2 }
//                     FASTP_PLOT(fastp_parse_ch_2)
// }

// workflow TEST_MAPPING_2 {
//         params.seq = "illumina"
//         params.type = "reads"

//         kaiju_nodes = file(params.nodes)
//         kaiju_db = file(params.kaiju_db)

//         TRIM_PE(seq_reads_pe_ch)
//             FASTP_REPORT(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))

//         MAP2REF_SWF_CORAL_PE(TRIM_PE.out.trimmed_reads, ref_host_ch)
        
//         MAP2REF_SWF_SYM_PE(TRIM_PE.out.trimmed_reads, ref_sym_ch)        
        
//         EXTRACT_BACTERIA_PE(TRIM_PE.out.trimmed_reads)

//                 fastp_parse_ch = FASTP_REPORT.out.report_json.concat(
//                 EXTRACT_BACTERIA_PE.out.report_json,
//                 MAP2REF_SWF_CORAL_PE.out.report_json,
//                 MAP2REF_SWF_SYM_PE.out.report_json)
//                     // fastp_parse_ch.collect{it[1]}.toList().view()
//                 names = fastp_parse_ch.collect{it[0]}.toList()   
//                 paths = fastp_parse_ch.collect{it[1]}.toList()
//                 names.merge(paths).set { fastp_parse_ch_2 }
//                     FASTP_PLOT(fastp_parse_ch_2)

// }


// //TODO : WORKING AS IS. TO BE optimised. 
// workflow KAIJU_PE_ILLUMINA {
//     // main:
//         kaiju_nodes = file(params.nodes)
//         kaiju_db = file(params.kaiju_db)
// // TODO check integration with fastqc + multiqc and optional make contigs assembly version. and fastp
//     // seq_reads_pe_ch.view()
//         TRIM_PE(seq_reads_pe_ch)
//             // FASTP_REPORT(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))
//             FASTP_REPORT(TRIM_PE.out.trimmed_reads)
//         // EXTRACT_BACTERIA_PE(TRIM_PE.out.trimmed_reads)
//             MAP2REF_SWF_CORAL_PE(TRIM_PE.out.trimmed_reads, ref_host_ch)
//                 MAP2REF_SWF_SYM_PE(MAP2REF_SWF_CORAL_PE.out.non_mapped_reads, ref_sym_ch)        
//         EXTRACT_BACTERIA_PE(MAP2REF_SWF_SYM_PE.out.non_mapped_reads)

//                 fastp_parse_ch = FASTP_REPORT.out.report_json.concat(
//                 EXTRACT_BACTERIA_PE.out.report_json,
//                 MAP2REF_SWF_CORAL_PE.out.report_json,
//                 MAP2REF_SWF_SYM_PE.out.report_json)
//                     // fastp_parse_ch.collect{it[1]}.toList().view()
//                 names = fastp_parse_ch.collect{it[0]}.toList()   
//                 paths = fastp_parse_ch.collect{it[1]}.toList()
//                 names.merge(paths).set { fastp_parse_ch_2 }
//                     FASTP_PLOT(fastp_parse_ch_2)
//     // emit:
//     //     EXTRACT_BACTERIA_PE.out.bacteria_reads
//     //     MAP2REF_SWF_CORAL_PE.out.mapped_reads
//     //     MAP2REF_SWF_SYM_PE.out.mapped_reads
//     //     EXTRACT_BACTERIA_PE.out.report_json
//     //     MAP2REF_SWF_CORAL_PE.out.report_json
//     //     MAP2REF_SWF_SYM_PE.out.report_json
//     //     FASTP_PARSE.out.report_csv
// }

// //TODO : WORKING AS IS. TO BE optimised. 
// workflow EXTRACT_CORAL_SYMB {
//     // main:
//         kaiju_nodes = file(params.nodes)
//         kaiju_db = file(params.kaiju_db)
// // TODO check integration with fastqc + multiqc and optional make contigs assembly version. and fastp
//     // seq_reads_pe_ch.view()
//         TRIM_PE(seq_reads_pe_ch)
//             // FASTP_REPORT(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))
//             FASTP_REPORT(TRIM_PE.out.trimmed_reads)
//         // EXTRACT_BACTERIA_PE(TRIM_PE.out.trimmed_reads)
//             MAP2REF_S_SWF_CORAL_PE(TRIM_PE.out.trimmed_reads, ref_host_ch)
//                 MAP2REF_S_SWF_SYM_PE(MAP2REF_S_SWF_CORAL_PE.out.non_mapped_reads, ref_sym_ch)        

//                 fastp_parse_ch = FASTP_REPORT.out.report_json.concat(
//                 // EXTRACT_BACTERIA_PE.out.report_json,
//                 MAP2REF_S_SWF_CORAL_PE.out.report_json,
//                 MAP2REF_S_SWF_SYM_PE.out.report_json)
//                     // fastp_parse_ch.collect{it[1]}.toList().view()
//                 names = fastp_parse_ch.collect{it[0]}.toList()   
//                 paths = fastp_parse_ch.collect{it[1]}.toList()
//                 names.merge(paths).set { fastp_parse_ch_2 }
//                     FASTP_PLOT(fastp_parse_ch_2)
// }

// // //TODO : WORKING AS IS. TO BE optimised. 
// workflow EXTRACT_CORAL_SYMB_BWA {
//     // main:
//         kaiju_nodes = file(params.nodes)
//         kaiju_db = file(params.kaiju_db)
// // TODO check integration with fastqc + multiqc and optional make contigs assembly version. and fastp
//     // seq_reads_pe_ch.view()
//         TRIM_PE(seq_reads_pe_ch)
//             // FASTP_REPORT(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))
//             FASTP_REPORT(TRIM_PE.out.trimmed_reads)
//         // EXTRACT_BACTERIA_PE(TRIM_PE.out.trimmed_reads)
//             MAP2REF_BWA_PE_CORAL(TRIM_PE.out.trimmed_reads, ref_host_ch)
//                 MAP2REF_BWA_PE_SYM(MAP2REF_BWA_PE_CORAL.out.non_mapped_reads, ref_sym_ch)        

//                 fastp_parse_ch = FASTP_REPORT.out.report_json.concat(
//                 // EXTRACT_BACTERIA_PE.out.report_json,
//                 MAP2REF_BWA_PE_CORAL.out.report_json,
//                 MAP2REF_BWA_PE_SYM.out.report_json)
//                     // fastp_parse_ch.collect{it[1]}.toList().view()
//                 names = fastp_parse_ch.collect{it[0]}.toList()   
//                 paths = fastp_parse_ch.collect{it[1]}.toList()
//                 names.merge(paths).set { fastp_parse_ch_2 }
//                     FASTP_PLOT(fastp_parse_ch_2)
// }

//TODO : WORKING AS IS. TO BE optimised. 
workflow KAIJU_CONTIGS_PE {
// TODO check integration with fastqc + multiqc and optional make contigs assembly version. and fastp
    // seq_reads_pe_ch.view()

        TRIM_PE(seq_reads_pe_ch)
            // FASTP_REPORT(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))
            FASTP_REPORT(TRIM_PE.out.trimmed_reads)
        MEGAHIT_PE(TRIM_PE.out.trimmed_reads)
        EXTRACT_BACTERIA_FA(MEGAHIT_PE.out.contigs_fa)
            MAP2REF_SWF_CORAL(EXTRACT_BACTERIA_FA.out.non_bacteria_reads, ref_host_ch)
                MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        
                
                fastp_parse_ch = FASTP_REPORT.out.report_json.concat(
                EXTRACT_BACTERIA_FA.out.report_json,
                MAP2REF_SWF_CORAL.out.report_json,
                MAP2REF_SWF_SYM.out.report_json)
                    // fastp_parse_ch.collect{it[1]}.toList().view()
                    // FASTP_PARSE(fastp_parse_ch.collect{it[1]}.toList())
                names = fastp_parse_ch.collect{it[0]}.toList() 
                paths = fastp_parse_ch.collect{it[1]}.toList()
                names.merge(paths).set { fastp_parse_ch_2 }
                    FASTP_PLOT(fastp_parse_ch_2)
}

workflow ONT_ANNA_SYM_CORAL {

// TODO check integration with fastqc + multiqc and optional make contigs assembly version. and fastp
    // seq_reads_pe_ch.view()
        // TRIM_PE(seq_reads_pe_ch)
            // FASTP_REPORT(seq_reads_pe_ch.concat(TRIM_PE.out.trimmed_reads))
            FASTP_REPORT(seq_reads_ont_ch)
        // MEGAHIT_PE(TRIM_PE.out.trimmed_reads)
        // EXTRACT_BACTERIA_FA(seq_reads_ch)
            MAP2REF_SWF_CORAL(seq_reads_ont_ch, ref_host_ch)
                MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        
                fastp_parse_ch = FASTP_REPORT.out.report_json.concat(
                // EXTRACT_BACTERIA_FA.out.report_json,
                MAP2REF_SWF_CORAL.out.report_json,
                MAP2REF_SWF_SYM.out.report_json)
                    // fastp_parse_ch.collect{it[1]}.toList().view()
                    // FASTP_PARSE(fastp_parse_ch.collect{it[1]}.toList())
                names = fastp_parse_ch.collect{it[0]}.toList()   
                paths = fastp_parse_ch.collect{it[1]}.toList()
                names.merge(paths).set { fastp_parse_ch_2 }
                    FASTP_PLOT(fastp_parse_ch_2)
}

// workflow {
//     // TODO: add R or python plots (long term add)
//     KAIJU_PE_ILLUMINA()
//     KRBR_PE_ILLUMINA()

// }
//TODO: cleanup and add same as above for ONT reads.
// workflow {
// //     fasta2fastq(seq_contigs_ch)
// //         FASTP_REPORT(seq_reads_pe_ch.out.converted_fasta)
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
//                 FASTP_PARSE(fastp_parse_ch.collect{it[1]}.toList())
// }


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

// workflow ONT_VF {

//     TRIM_NANOPORE(seq_reads_ch)
//         EXTRACT_BACTERIA(TRIM_NANOPORE.out.trimmed_reads)
//             MAP2REF_SWF_CORAL(EXTRACT_BACTERIA.out.non_bacteria_reads, ref_host_ch)
//                 MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        

//             fastp_parse_ch = TRIM_NANOPORE.out.report_json.concat(
//             EXTRACT_BACTERIA.out.report_json,
//             MAP2REF_SWF_CORAL.out.report_json,
//             MAP2REF_SWF_SYM.out.report_json)
//                 // fastp_parse_ch.collect{it[1]}.toList().view()
//                 FASTP_PARSE(fastp_parse_ch.collect{it[1]}.toList())
// }

// workflow CONTIGS {

//     fastp_fasta_report(seq_contigs_ch)
//     EXTRACT_BACTERIA_FA(seq_contigs_ch)
//         MAP2REF_SWF_CORAL(EXTRACT_BACTERIA_FA.out.non_bacteria_reads, ref_host_ch)
//             MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        

//             fastp_parse_ch = fastp_fasta_report.out.report_json.concat(
//             EXTRACT_BACTERIA_FA.out.report_json,
//             MAP2REF_SWF_CORAL.out.report_json,
//             MAP2REF_SWF_SYM.out.report_json)
//                 // fastp_parse_ch.collect{it[1]}.toList().view()
//                 FASTP_PARSE(fastp_parse_ch.collect{it[1]}.toList())
// }

// workflow  {
    
//     fasta2fastq(seq_contigs_ch)
//         FASTP_REPORT(seq_reads_pe_ch.out.converted_fasta)

//     TRIM_NANOPORE(seq_reads_ch)
//         EXTRACT_BACTERIA(TRIM_NANOPORE.out.trimmed_reads.concat(fasta2fastq.out.converted_fasta))
//             MAP2REF_SWF_CORAL(EXTRACT_BACTERIA.out.non_bacteria_reads, ref_host_ch)
//                 MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        

//             fastp_parse_ch = TRIM_NANOPORE.out.report_json.concat(
//             FASTP_REPORT.out.report_json,
//             EXTRACT_BACTERIA.out.report_json,
//             MAP2REF_SWF_CORAL.out.report_json,
//             MAP2REF_SWF_SYM.out.report_json)
//                 // fastp_parse_ch.collect{it[1]}.toList().view()
//                 FASTP_PARSE(fastp_parse_ch.collect{it[1]}.toList())
// }

// workflow NON_TRIMMED_assembly  {
    
//     MEGAHIT_SE(seq_reads_ch)
//     fasta2fastq(seq_contigs_ch.concat(MEGAHIT_SE.out.contigs_fa))
//         FASTP_REPORT(fasta2fastq.out.converted_fasta)

//         EXTRACT_BACTERIA(fasta2fastq.out.converted_fasta)
//             MAP2REF_SWF_CORAL(EXTRACT_BACTERIA.out.non_bacteria_reads, ref_host_ch)
//                 MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        

//             fastp_parse_ch = FASTP_REPORT.out.report_json.concat(
//             EXTRACT_BACTERIA.out.report_json,
//             MAP2REF_SWF_CORAL.out.report_json,
//             MAP2REF_SWF_SYM.out.report_json)
//                 // fastp_parse_ch.collect{it[1]}.toList().view()
//                 FASTP_PARSE(fastp_parse_ch.collect{it[1]}.toList())
// }

// workflow NON_TRIMMED  {
    
//     fasta2fastq(seq_contigs_ch)
//         FASTP_REPORT(fasta2fastq.out.converted_fasta.concat(seq_reads_ch))

//         EXTRACT_BACTERIA(fasta2fastq.out.converted_fasta.concat(seq_reads_ch))
//             MAP2REF_SWF_CORAL(EXTRACT_BACTERIA.out.non_bacteria_reads, ref_host_ch)
//                 MAP2REF_SWF_SYM(MAP2REF_SWF_CORAL.out.non_mapped_reads, ref_sym_ch)        

//             fastp_parse_ch = FASTP_REPORT.out.report_json.concat(
//             EXTRACT_BACTERIA.out.report_json,
//             MAP2REF_SWF_CORAL.out.report_json,
//             MAP2REF_SWF_SYM.out.report_json)
//                 // fastp_parse_ch.collect{it[1]}.toList().view()
//                 FASTP_PARSE(fastp_parse_ch.collect{it[1]}.toList())
// }


// workflow PE_ILLUMINA  {
//     // seq_reads_pe_ch.view()
//     FASTP_REPORT(seq_reads_pe_ch)
//     EXTRACT_BACTERIA_PE(seq_reads_pe_ch)
//         MAP2REF_SWF_CORAL_PE(EXTRACT_BACTERIA_PE.out.non_bacteria_reads, ref_host_ch)
//             MAP2REF_SWF_SYM_PE(MAP2REF_SWF_CORAL_PE.out.non_mapped_reads, ref_sym_ch)        

//             fastp_parse_ch = FASTP_REPORT.out.report_json.concat(
//             EXTRACT_BACTERIA_PE.out.report_json,
//             MAP2REF_SWF_CORAL_PE.out.report_json,
//             MAP2REF_SWF_SYM_PE.out.report_json)
//                 fastp_parse_ch.collect{it[1]}.toList().view()
//                 // FASTP_PARSE(fastp_parse_ch.collect{it[1]}.toList())
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


workflow PARTIAL_INPUT{
//     countSeqReadsONT = seq_reads_ont_ch.count().into { count ->
//     if (count == 0) {
//         println ("The channel is empty")
//         // return true // Or do something indicating the channel is empty
//     } else {
//         println ("The channel has $count elements")
//         // return false // Or do something indicating the channel is not empty
//     }
// }
// Check the output of seq_reads_ont_ch.view { it[3] }
// seq_reads_ont_ch.view { it[3] }.collect().subscribe { value ->
//     println "Values in the channel: $value"
//     print("$value" == '[empty]')
// }   
isONTPresent = seq_reads_ont_ch.view { it[3] }.any { it == 'empty' }
// print(isONTPresent)
// print(seq_reads_ont_ch.any { it == 'empty' })
    if (isONTPresent)
    MAPPING_CONSECUTIVE_PE(seq_reads_pe_ch, ch_host, ch_sym)
    else if( seq_reads_pe_ch == 'empty' )
    MAPPING_CONSECUTIVE(seq_reads_ont_ch, ch_host, ch_sym)
    else
    seq_reads_ont_ch.view()
    // // seq_reads_ont_ch.count()
    // print(isChannelEmpty)
    // seq_reads_ont_ch.view()
    // print(seq_reads_ont_ch.view().size())
    // // print(seq_reads_ont_ch == 'empty')
    // print(exist_ont)
    // MAPPING_CONSECUTIVE_PE(seq_reads_pe_ch, ch_host, ch_sym)
    // MAPPING_CONSECUTIVE(seq_reads_ont_ch, ch_host, ch_sym)
    // MAPPING_PARALLEL(seq_reads_pe_ch, seq_reads_ont_ch, ch_host, ch_sym)
}
