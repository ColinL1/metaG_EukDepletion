#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Set params for reads and power 
params.reads = "$baseDir/input/*.fastq.gz"
params.outdir = "$baseDir/results_full_1/"
params.cpusHigh = "40"
params.cpusVHigh = "100"
params.cpusMin = "4"
params.memMax = '500 GB'

//temporary fix to multiply the reference channel to match the number of samples
params.n_samples = 108

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
    cpus       : ${params.cpusMin} - ${params.cpusHigh}
    memory     : ${params.memMax}

    """

Channel.fromPath(params.reads).map {tuple( it.name.split('.fastq.gz')[0], it )}.set { seq_file_ch }

include { kaiju } from './modules/kaiju_multi.nf'
include { split_bac } from './modules/split_seqkit.nf'

include { map2ref; map2ref as map2ref1 } from './modules/minimap2.nf' 
include { split_bam; split_bam as split_bam1 } from './modules/split_bam.nf'

include { fastp_report } from './modules/fastp_stats.nf'

workflow {
    main:
        kaiju(seq_file_ch)
        split_bac_ch = kaiju.out.report_kaiju_out.join(seq_file_ch)
            split_bac(split_bac_ch)

        Channel.fromPath(params.ref_coral)
            .map {tuple( it.name.split('.mmi')[0], it )}
            .flatMap {it * params.n_samples }
            .collate(2)
            .set { ref_coral_ch }


        map2ref(split_bac.out.non_bacteria_reads, ref_coral_ch)
            split_bam(map2ref.out.mapp_file, ref_coral_ch)

        Channel.fromPath(params.ref_sym)
            .map {tuple( it.name.split('.mmi')[0], it )}
            .flatMap {it * params.n_samples }
            .collate(2)
            .set { ref_sym_ch }

        map2ref1(split_bam.out.unmapped_reads, ref_sym_ch)
            split_bam1(map2ref1.out.mapp_file, ref_sym_ch)

                    fastp_report(split_bam1.out.mapped_reads.concat(split_bam1.out.unmapped_reads, split_bam.out.unmapped_reads, split_bam.out.mapped_reads, split_bac.out.bacteria_reads, split_bac.out.non_bacteria_reads, seq_file_ch))

}

workflow test {
    main:
        kaiju(seq_file_ch)
        split_bac_ch = kaiju.out.report_kaiju_out.join(seq_file_ch)
            split_bac_ch.view()
}

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

    def subject = 'Pipeline execution summary - Diamond BlastX'
    def recipient = 'luigi.colin@uni-konstanz.de'

    ['mailx', '-s', subject, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
	log.info ( workflow.success ? "\nDone! --> $params.outdir\n" : "Oops .. something went wrong" )

}