#!/usr/bin/env nextflow

// Set params for reads and power 
params.reads = "$baseDir/input/*.fastq.gz"
params.contigs = "$baseDir/input/*.fasta"
params.outdir = "$baseDir/results_ont/"
// params.cpusHigh = "40"
// params.cpusVHigh = "100"
// params.cpusMin = "4"
// params.memMax = '500 GB'

//temporary fix to multiply the reference channel to match the number of samples
params.n_samples = 120

// // params for mmeseqs2
// params.diamond_db_nr = "/home/colinl/database/diamond/nr.dmnd"
// params.meganizer_db = "/home/colinl/database/megan/megan-map-Feb2022.db"
// params.py_scripts = "/home/colinl/metaG/Git/metaG_methods/blastX/py_script"

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

include { ont_trim } from './modules/ONT_prep.nf'
include { kaiju } from './modules/kaiju_multi.nf'
include { split_bac } from './modules/split_seqkit.nf'

include { map2ref; map2ref as map2ref1 } from './modules/minimap2.nf' 
include { split_bam; split_bam as split_bam1 } from './modules/split_bam.nf'

include { fastp_report } from './modules/fastp_stats.nf'
include { fastp_parse } from './modules/fastp_parse.nf'

workflow test {
    main:
        // kaiju(seq_reads_ch)
        // split_bac_ch = kaiju.out.report_kaiju_out.join(seq_reads_ch)
        //     split_bac_ch.view()
        seq_contigs_ch.view()
}

workflow trimm_nanopore {
    main:
        ont_trim(seq_reads_ch)
            fastp_report(ont_trim.out.trimmed_fastq)
    
    emit:
        trimmed_reads = ont_trim.out.trimmed_fastq
        report_json = fastp_report.out.report_json
}

workflow extract_bacteria {
    main:
        trimm_nanopore()
        kaiju(trimm_nanopore.out.trimmed_reads)
        split_bac_ch = kaiju.out.report_kaiju_out.join(trimm_nanopore.out.trimmed_reads)
            split_bac(split_bac_ch)
            fastp_ch = split_bac.out.bacteria_reads.concat(split_bac.out.non_bacteria_reads)
            fastp_report(fastp_ch)
    
    emit:
        bacteria_reads = split_bac.out.bacteria_reads
        non_bacteria_reads = split_bac.out.non_bacteria_reads
        report_json = fastp_report.out.report_json
        
        trimmed_reads = trimm_nanopore.out.trimmed_reads
        report_json_og = trimm_nanopore.out.report_json
}

workflow extract_corals {
    main:
    // multiply reference file tuple by the number of samples present (temporary solution)
        Channel.fromPath(params.ref_coral) 
            .map {tuple( it.name.split('.mmi')[0], it )}
            .flatMap {it * params.n_samples }
            .collate(2)
            .set { ref_coral_ch }
        
        extract_bacteria()
        map2ref(extract_bacteria.out.non_bacteria_reads, ref_coral_ch)
            split_bam(map2ref.out.mapp_file, ref_coral_ch)
                fastp_report(split_bam.out.mapped_reads.concat(split_bam.out.unmapped_reads))

    emit:
        coral_reads = split_bam.out.mapped_reads
        non_coral_reads = split_bam.out.unmapped_reads
        report_json = fastp_report.out.report_json

        bacteria_reads = extract_bacteria.out.bacteria_reads
        non_bacteria_reads = extract_bacteria.out.non_bacteria_reads
        report_json_bc = extract_bacteria.out.report_json
        trimmed_reads = extract_bacteria.out.trimmed_reads
        report_json_og = extract_bacteria.out.report_json
}

workflow extract_sym {
    main:
    // multiply reference file tuple by the number of samples present (temporary solution)
        Channel.fromPath(params.ref_sym)
            .map {tuple( it.name.split('.mmi')[0], it )}
            .flatMap {it * params.n_samples }
            .collate(2)
            .set { ref_sym_ch }
        
        extract_corals()
        map2ref(extract_corals.out.non_coral_reads, ref_sym_ch)
            split_bam(map2ref.out.mapp_file, ref_sym_ch)
                fastp_report(split_bam.out.mapped_reads.concat(split_bam.out.unmapped_reads))

    emit:
        symbiodiniaceae_reads = split_bam.out.mapped_reads
        non_symbiodiniaceae_reads = split_bam.out.unmapped_reads
        report_json = fastp_report.out.report_json

        coral_reads = extract_corals.out.coral_reads
        non_coral_reads = extract_corals.out.non_coral_reads
        report_json_co = extract_corals.out.report_json
        bacteria_reads = extract_corals.out.bacteria_reads
        non_bacteria_reads = extract_corals.out.non_bacteria_reads
        report_json_bc = extract_corals.out.report_json
        trimmed_reads = extract_corals.out.trimmed_reads
        report_json_og = extract_corals.out.report_json
}

workflow reports_parse {
    main:
        fastp_parse()
    
    emit:
        report_csv = fastp_parse.out.report_csv
}


workflow ONT_pipeline {
    main:
        mode = 'ONT'
        extract_sym()
        
        fastp_parse(extract_sym.out.report_json.concat(
            extract_sym.out.report_json_co,
            extract_sym.out.report_json_bc,
            extract_sym.out.report_json_og))
}

workflow {
    main:
        extract_bacteria(seq_reads_ch)
        extract_corals(extract_bacteria.out.non_bacteria_reads)
        extract_sym(extract_corals.out.non_coral_reads)
        
        reports_fastp(seq_reads_ch
                .concat(extract_sym.out.symbiodiniaceae_reads,
                extract_sym.out.non_symbiodiniaceae_reads,
                extract_corals.out.coral_reads,
                extract_corals.out.non_coral_reads,
                extract_bacteria.out.bacteria_reads,
                extract_bacteria.out.non_bacteria_reads))

}

workflow fasta_pipeline{
    
    main:
        mode = 'illumina_contigs'

        extract_bacteria(seq_contigs_ch)
        extract_corals(extract_bacteria.out.non_bacteria_reads)
        extract_sym(extract_corals.out.non_coral_reads)
        
        reports_fastp(seq_reads_ch
                .concat(extract_sym.out.symbiodiniaceae_reads,
                extract_sym.out.non_symbiodiniaceae_reads,
                extract_corals.out.coral_reads,
                extract_corals.out.non_coral_reads,
                extract_bacteria.out.bacteria_reads,
                extract_bacteria.out.non_bacteria_reads))

}

// workflow {
//     main:
//         kaiju(seq_reads_ch)
//         split_bac_ch = kaiju.out.report_kaiju_out.join(seq_reads_ch)
//             split_bac(split_bac_ch)

//         Channel.fromPath(params.ref_coral)
//             .map {tuple( it.name.split('.mmi')[0], it )}
//             .flatMap {it * params.n_samples }
//             .collate(2)
//             .set { ref_coral_ch }


//         map2ref(split_bac.out.non_bacteria_reads, ref_coral_ch)
//             split_bam(map2ref.out.mapp_file, ref_coral_ch)

//         Channel.fromPath(params.ref_sym)
//             .map {tuple( it.name.split('.mmi')[0], it )}
//             .flatMap {it * params.n_samples }
//             .collate(2)
//             .set { ref_sym_ch }

//         map2ref1(split_bam.out.unmapped_reads, ref_sym_ch)
//             split_bam1(map2ref1.out.mapp_file, ref_sym_ch)

//                 fastp_report(split_bam1.out.mapped_reads.concat(split_bam1.out.unmapped_reads, split_bam.out.unmapped_reads, split_bam.out.mapped_reads, split_bac.out.bacteria_reads, split_bac.out.non_bacteria_reads, seq_reads_ch))

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

    // def subject = 'Pipeline execution summary'
    // def recipient = 'luigi.colin@uni-konstanz.de'

    // // ['mailx', '-s', subject, recipient].execute() << """

    // Pipeline execution summary
    // ---------------------------
    // Completed at: ${workflow.complete}
    // Duration    : ${workflow.duration}
    // Success     : ${workflow.success}
    // workDir     : ${workflow.workDir}
    // exit status : ${workflow.exitStatus}
    // Error report: ${workflow.errorReport ?: '-'}
    // """
	log.info ( workflow.success ? "\nDone! --> $params.outdir\n" : "Oops .. something went wrong" )

}