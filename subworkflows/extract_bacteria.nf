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

include { kaiju_pe; kaiju_multi; kaiju_multi_pe } from '../modules/kaiju_multi.nf'
include { split_bac; split_bac_pe } from '../modules/split_seqkit.nf'
include { fastp_report; fastp_fasta_report } from '../modules/fastp_stats.nf'

/*
========================================================================================
    Workflow EXTRACT_BACTERIA
========================================================================================
*/

workflow EXTRACT_BACTERIA {
    take: 
        file
    main:
        names = file.collect{it[0]}.toList()
            paths = file.collect{it[1]}.toList()
            output = file.collect{it[0]+ "_out"}.toList()
            
            names
                .merge(paths)
                .merge(output)
                .set { kaiju_multi_input_ch }
        
        kaiju_multi(kaiju_multi_input_ch)
            kaiju_multi.out.kaiju_out.flatten().map{ tuple (it.name.replaceAll(/_out/, ""), it)}.set{ kaiju_multi_out_clean_ch }
                split_bac(kaiju_multi_out_clean_ch.join(file))
                fastp_report(split_bac.out.concat())

    emit:
        bacteria_reads = split_bac.out.bacteria_reads
        non_bacteria_reads = split_bac.out.non_bacteria_reads
        report_json = fastp_report.out.report_json
}

workflow EXTRACT_BACTERIA_PE {
    take: 
        file
    main:
        // names = file.collect{it[0]}.toList()
        // path1 = file.collect{it[1][0]}.toList()
        // path2 = file.collect{it[1][1]}.toList()
        // output = file.collect{it[0]+ "_out"}.toList()

        // names
        //     .merge(path1)
        //     .merge(path2)
        //     .merge(output)
        //     .set { kaiju_multi_input_ch }
        kaiju_pe(file)
            // kaiju_multi_pe.out.kaiju_out.flatten().map{ tuple (it.name.replaceAll(/_out/, ""), it)}.set{ kaiju_multi_out_clean_ch }
            split_bac_pe(kaiju_pe.out.report_kaiju_out.join(file))
                fastp_report(split_bac_pe.out.concat())

    emit:
        bacteria_reads = split_bac_pe.out.bacteria_reads
        non_bacteria_reads = split_bac_pe.out.non_bacteria_reads
        report_json = fastp_report.out.report_json
}

workflow EXTRACT_BACTERIA_FA {
    take: 
        file
    main:
        names = file.collect{it[0]}.toList()
            paths = file.collect{it[1]}.toList()
            output = file.collect{it[0]+ "_out"}.toList()
            
            names
                .merge(paths)
                .merge(output)
                .set { kaiju_multi_input_ch }
        
        kaiju_multi(kaiju_multi_input_ch)
            kaiju_multi.out.kaiju_out.flatten().map{ tuple (it.name.replaceAll(/_out/, ""), it)}.set{ kaiju_multi_out_clean_ch }
                split_bac(kaiju_multi_out_clean_ch.join(file))
                    fastp_fasta_report(split_bac.out.concat())
    
    emit:
        bacteria_reads = split_bac.out.bacteria_reads
        non_bacteria_reads = split_bac.out.non_bacteria_reads
        report_json = fastp_fasta_report.out.report_json
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
//             split_bac_pe(kaiju_multi_out_clean_ch.join(file))
//                 fastp_report(split_bac_pe.out.concat())

//     emit:
//         bacteria_reads = split_bac_pe.out.bacteria_reads
//         non_bacteria_reads = split_bac_pe.out.non_bacteria_reads
//         report_json = fastp_report.out.report_json
// }