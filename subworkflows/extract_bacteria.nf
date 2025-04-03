#!/usr/bin/env nextflow

/*
========================================================================================
    Run kaiju and extract bacteria matching reads
========================================================================================
*/

params.nodes = "/share/databases/kaiju/nodes.dmp"
params.kaiju_db = "/share/databases/kaiju/refseq/kaiju_db_refseq.fmi"

/*
========================================================================================
    Include Modules
========================================================================================
*/

include { KAIJU_SE; KAIJU_PE; } from '../modules/kaiju_multi.nf' //kaiju_multi; kaiju_multi_pe
include { SEQKIT_SPLIT_BAC; SEQKIT_SPLIT_BAC_PE } from '../modules/split_seqkit.nf'
include { FASTP_REPORT; fastp_fasta_report } from '../modules/fastp_stats.nf'

/*
========================================================================================
    Workflow EXTRACT_BACTERIA
========================================================================================
*/



workflow EXTRACT_BACTERIA_ONT {
    take: 
        file
    main:
        
        KAIJU_SE(file)
            KAIJU_SE.out.report_kaiju_out.join(file)
            .map { tuple (it[0], it[1], it[2], it[5], it[6])}
            .set { kaiju_out_ch }
                SEQKIT_SPLIT_BAC(kaiju_out_ch)
                FASTP_REPORT(SEQKIT_SPLIT_BAC.out.concat())

    emit:
        bacteria_reads = SEQKIT_SPLIT_BAC.out.bacteria_reads
        non_bacteria_reads = SEQKIT_SPLIT_BAC.out.non_bacteria_reads
        report_json = FASTP_REPORT.out.report_json
}

workflow EXTRACT_BACTERIA_PE {
    take: 
        file
    main:

        KAIJU_PE(file)

        KAIJU_PE.out.report_kaiju_out.join(file)
            .map { tuple (it[0], it[1], it[2], it[5], it[6])}
            .set { kaiju_out_ch }

            SEQKIT_SPLIT_BAC_PE(kaiju_out_ch)
                FASTP_REPORT(SEQKIT_SPLIT_BAC_PE.out.concat())

    emit:
        bacteria_reads = SEQKIT_SPLIT_BAC_PE.out.bacteria_reads
        non_bacteria_reads = SEQKIT_SPLIT_BAC_PE.out.non_bacteria_reads
        report_json = FASTP_REPORT.out.report_json
}

// // // individual tester #TODO:remove once complete
// workflow  {
//     Channel.fromFilePairs("/home/colinl/metaG/Git/metaG_EukDepletion/input/tests/cor_test_ill/*_{1,2}.fq.gz"). set {input_fq}
//     // // input_fq.view()
//         // names = input_fq.collect{it[0]}.toList()
//         // path1 = input_fq.collect{it[1][0]}.toList()
//         // path2 = input_fq.collect{it[1][1]}.toList()
//         // output = input_fq.collect{it[0]+ "_out"}.toList()
//         //     // input_fq.map{tuple (it[0], it[1][1], it[1][0], it[0]+ "_out")}
//         //     // .filter("*_1.fq.gz")'${1}_2.fq.gz
//         //     // .map{ tuple ( it == "*_1.fq.gz", it == "*_2.fq.gz")}
//         //     // .collect()
//         //     // .view()
//         // names
//         //     .merge(path1)
//         //     .merge(path2)
//         //     .merge(output)
//         //     .set { kaiju_multi_input_ch }
//         //     kaiju_multi_input_ch.view()
//     // kaiju_multi_input_ch.view()
//     // EXTRACT_BACTERIA_FA(input_fq)
//     // EXTRACT_BACTERIA_FA.out.bacteria_reads.view()
//     EXTRACT_BACTERIA_PE(input_fq)
// }

// workflow EXTRACT_BACTERIA_PE_MULTI {
//     take: 
//         file
//     main:
//         names = file.collect{it[0]}.toList()
//         path1 = file.collect{it[1][0]}.toList()
//         path2 = file.collect{it[1][1]}.toList()
//         output = file.collect{it[0]+ "_out"}.toList()

//         names
//             .merge(path1)
//             .merge(path2)
//             .merge(output)
//             .set { kaiju_multi_input_ch }

//         kaiju_multi_pe(kaiju_multi_input_ch)
//             kaiju_multi_pe.out.kaiju_out.flatten().map{ tuple (it.name.replaceAll(/_out/, ""), it)}.set{ kaiju_multi_out_clean_ch }
//             SEQKIT_SPLIT_BAC_PE(kaiju_multi_out_clean_ch.join(file))
//                 FASTP_REPORT(SEQKIT_SPLIT_BAC_PE.out.concat())

//     emit:
//         bacteria_reads = SEQKIT_SPLIT_BAC_PE.out.bacteria_reads
//         non_bacteria_reads = SEQKIT_SPLIT_BAC_PE.out.non_bacteria_reads
//         report_json = FASTP_REPORT.out.report_json
// }
// workflow CH_TEST {
    
//     Channel.fromPath("/home/colinl/metaG/Git/metaG_EukDepletion/input.csv")
//     .splitCsv(header:true, quote: '\"')
//     .map {row -> tuple(row.Sample_name, row.Sample_name, (row.R), (row.Seqencing_type), (row.Entry))}
//     .filter { it[3] == 'ONT' }
//     .filter { it[4] == 'TEST_MAPPING_2' }
//     .map { tuple(it[0], it[1], it[2], it[3]) }
//     .set{ seq_reads_ont_ch }
        
//         KAIJU_SE(seq_reads_ont_ch)
//             KAIJU_SE.out.report_kaiju_out.join(seq_reads_ont_ch)
//             .map { tuple (it[0], it[1], it[2], it[5], it[6])}
//             .set { kaiju_out_ch }

//         KAIJU_SE.out.report_kaiju_out.join(seq_reads_ont_ch)
//         .set { kaiju_out_ch_2 }
// // 
//         kaiju_out_ch.view()
//         // kaiju_out_ch_2.view()
//                 // SEQKIT_SPLIT_BAC(kaiju_out_ch)
//                 // FASTP_REPORT(SEQKIT_SPLIT_BAC.out.concat())

// }
// workflow EXTRACT_BACTERIA {
//     take: 
//         file
//     main:
//         names = file.collect{it[0]}.toList()
//             paths = file.collect{it[1]}.toList()
//             output = file.collect{it[0]+ "_out"}.toList()
            
//             names
//                 .merge(paths)
//                 .merge(output)
//                 .set { kaiju_multi_input_ch }
        
//         kaiju_multi(kaiju_multi_input_ch)
//             kaiju_multi.out.kaiju_out.flatten().map{ tuple (it.name.replaceAll(/_out/, ""), it)}.set{ kaiju_multi_out_clean_ch }
//                 SEQKIT_SPLIT_BAC(kaiju_multi_out_clean_ch.join(file))
//                 FASTP_REPORT(SEQKIT_SPLIT_BAC.out.concat())

//     emit:
//         bacteria_reads = SEQKIT_SPLIT_BAC.out.bacteria_reads
//         non_bacteria_reads = SEQKIT_SPLIT_BAC.out.non_bacteria_reads
//         report_json = FASTP_REPORT.out.report_json
// }


// workflow EXTRACT_BACTERIA_FA {
//     take: 
//         file
//     main:
//         names = file.collect{it[0]}.toList()
//             paths = file.collect{it[1]}.toList()
//             output = file.collect{it[0]+ "_out"}.toList()
            
//             names
//                 .merge(paths)
//                 .merge(output)
//                 .set { kaiju_multi_input_ch }
        
//         kaiju_multi(kaiju_multi_input_ch)
//             kaiju_multi.out.kaiju_out.flatten().map{ tuple (it.name.replaceAll(/_out/, ""), it)}.set{ kaiju_multi_out_clean_ch }
//                 SEQKIT_SPLIT_BAC(kaiju_multi_out_clean_ch.join(file))
//                     fastp_fasta_report(SEQKIT_SPLIT_BAC.out.concat())
    
//     emit:
//         bacteria_reads = SEQKIT_SPLIT_BAC.out.bacteria_reads
//         non_bacteria_reads = SEQKIT_SPLIT_BAC.out.non_bacteria_reads
//         report_json = fastp_fasta_report.out.report_json
// }
