process megahit_pe {
    tag "${sample}"
    publishDir "$params.outdir/illumina_contigs", mode: 'symlink'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path ("${sample}/final.contigs.fa"), emit: contigs_fa


    script:
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} --out-dir ${sample} --k-min 27 --k-max 127 --k-step 10 --num-cpu-threads ${task.cpus}
    """
}

process megahit_se {
    tag "${sample}"
    publishDir "$params.outdir/ONT_contigs", mode: 'symlink'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path ("${sample}/final.contigs.fa"), emit: contigs_fa

    script:
    """
    megahit -r ${reads}  --out-dir ${sample} --k-min 27 --k-max 127 --k-step 10 --num-cpu-threads ${task.cpus}
    """
}

workflow {
    Channel.fromFilePairs("/home/colinl/metaG/Git/metaG_EukDepletion/input/tests/cor_test_ill/*_{1,2}.fq.gz"). set {input_fq}
    megahit_pe(input_fq)
    megahit_pe.out.contigs_fa.view()
}