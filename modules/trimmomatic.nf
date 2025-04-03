#!/usr/bin/env nextflow
// Trimmomatic illumina reads trimming

process TRIMMOMATIC {
    tag "${sample}"
    // cpus "${params.cpusHigh}"
    // memory "${params.memMax}"
    publishDir "$params.outdir/Trimmed_fastq/", mode: 'symlink'

    input: 
    tuple val(sample), val(base_name), path(reads), val (seq_type)

    output:
    tuple val(sample), val("${base_name}.TRIM"), path("${base_name}_TRIM_paried_{1,2}.fq.gz"), val (seq_type), emit: trimmed_fq_out
    tuple val(sample), val("${base_name}.TRIM"), path("${base_name}_TRIM_unparied_{1,2}.fq.gz"), val (seq_type), emit: singleton_fq_out
    tuple val(sample), val(base_name), path("${base_name}_trimmomatic.log"), val (seq_type), emit: log_trim

    script:
    """
    trimmomatic PE -threads ${task.cpus} ${reads[0]} ${reads[1]} ${base_name}_TRIM_paried_1.fq.gz ${base_name}_TRIM_unparied_1.fq.gz ${base_name}_TRIM_paried_2.fq.gz ${base_name}_TRIM_unparied_2.fq.gz ILLUMINACLIP:/home/colinl/metaG/Git/metaG_EukDepletion/modules/adapters/Truseq_V.3_edited.fa:2:30:10:8:TRUE HEADCROP:5 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 2>> ${base_name}_trimmomatic.log 1>>file.out
    """
    stub:
    """
    touch ${base_name}_TRIM_paried_1.fq.gz
    touch ${base_name}_TRIM_paried_2.fq.gz
    touch ${base_name}_TRIM_unparied_1.fq.gz
    touch ${base_name}_TRIM_unparied_2.fq.gz
    touch ${base_name}_trimmomatic.log
    """
}

workflow {
    Channel.fromFilePairs("/home/colinl/metaG/Git/metaG_EukDepletion/concat_fastqs/*_{1,2}.fq.gz")
        | map { id, reads ->
            (species, replicate, method, buffer, other, unspecified) = id.tokenize("_")
            meta = [
                id:id,
                // single_end:'PE',
                species:species,
                replicate:replicate,
                method:method,
                buffer:buffer,
                other:other,
                unspecified:unspecified,
            ]
            [meta, reads]
    }   //| view()
    | TRIMMOMATIC()
}