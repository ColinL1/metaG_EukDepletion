#!/usr/bin/env nextflow

process clean_names {
    tag "${sample}"
    cpus "${params.cpusMin}"
    publishDir "$params.outdir/fixed_names/", mode: 'symlink'

    input: 
    tuple val(sample), path(reads)

    output:
    path ("${sample}_clean_name.fastq.gz"), emit: clean_names_reads

    script:
    """
    bioawk -v FILE="${sample}" -c fastx '{ print \$name , FILE}' < ${reads} > ${sample}_new_keyfile.txt
    seqkit replace -k ${sample}_new_keyfile.txt -p '^(\\S+)' -r '{kv}_\$1' ${reads} -j ${task.cpus} -o ${sample}_clean_name-temp.fastq >/dev/null 2>/dev/null
    seqtk seq -l 0 ${sample}_clean_name-temp.fastq > ${sample}_clean_name.fastq
    rm ${sample}_clean_name-temp.fastq ${sample}_new_keyfile.txt
    pigz -p ${task.cpus} ${sample}_clean_name.fastq
    """
}

process clean_names_fasta {
    tag "${sample}"
    cpus "${params.cpusMin}"
    publishDir "$params.outdir/fixed_names/", mode: 'symlink'

    input: 
    tuple val(sample), path(reads)

    output:
    path ("${sample}_clean_name.fastq.gz"), emit: clean_names_reads

    script:
    """
    bioawk -v FILE="${sample}" -c fastx '{ print \$name , FILE}' < ${reads} > ${sample}_new_keyfile.txt
    seqkit replace -k ${sample}_new_keyfile.txt -p '^(\\S+)' -r '{kv}_\$1' ${reads} -j ${task.cpus} -o ${sample}_clean_name-temp.fastq >/dev/null 2>/dev/null
    seqtk seq -l 0 ${sample}_clean_name-temp.fastq > ${sample}_clean_name.fastq
    seqtk seq -a ${sample}_clean_name.fastq > ${sample}_clean_name.fa
    rm ${sample}_clean_name-temp.fastq ${sample}_new_keyfile.txt
    pigz -p ${task.cpus} ${sample}_clean_name.fa
    """
}

process medaka {
    tag "${sample}"
    cpus "${params.cpusMin}"
    publishDir "$params.outdir/medaka/", mode: 'symlink'
    
    input: 
    tuple val (sample), val(sample_assembly), val(contigs), val(reads)

    output:
    tuple val("${sample_assembly}_medaka"), path ("${sample_assembly}_medaka/consensus.fasta"), emit: medaka_consensus
    

    script:
    """
    medaka_consensus -i ${reads} -d ${contigs} -o ${sample_assembly}_medaka -t ${task.cpus} -m r941_min_sup_g507
    """
}

process unicycler {
    tag "${sample}"
    cpus "${params.cpusHigh}"
    errorStrategy "ignore"
    publishDir "$params.outdir/unicycler/", mode: 'symlink'
    
    input: 
    tuple val(sample), path(reads)

    output:
    tuple val ("${sample}"), val("${sample}_unicycler"), path ("${sample}_unicycler/assembly.fasta"), emit: unicycler_assembly
    // tuple val("${sample}"), path ("${sample}_medaka/consensus.fasta"), emit: medaka_consensus
    

    script:
    """
    unicycler -l ${reads} -o ${sample}_unicycler -t ${task.cpus} --keep 3
    """
}

process canu {
    tag "${sample}"
    cpus "${params.cpusHigh}"
    errorStrategy "ignore"
    publishDir "$params.outdir/canu/", mode: 'symlink'
    
    input: 
    tuple val(sample), path(reads)

    output:
    tuple val("${sample}_canu"), path ("${sample}_canu/${sample}.contigs.fasta"), emit: canu_assembly
    

    script:
    """
    canu -p ${sample} -d ${sample}_canu genomeSize=4.8m correctedErrorRate=0.16 maxInputCoverage=50 stopOnLowCoverage=5 maxThreads=${task.cpus} -nanopore ${reads} 
    """
}

process flye {
    tag "${sample}"
    cpus "${params.cpusHigh}"
    errorStrategy "ignore"
    publishDir "$params.outdir/flye/", mode: 'symlink'
    
    input: 
    tuple val(sample), path(reads)

    output:
    tuple val("${sample}_flye"), path ("${sample}_flye/assembly.fasta"), emit: flye_assembly
    

    script:
    """
    flye --nano-hq ${reads} --out-dir ${sample}_flye --genome-size 4.8m --threads ${task.cpus} --meta
    """
}

process megahit {
    tag "${sample}"
    cpus "${params.cpusHigh}"
    publishDir "$params.outdir/megahit/", mode: 'symlink'

    input:
    tuple val(sample), path(reads) 

    output:
    tuple val("${sample}_megahit"), path ("${sample}_megahit/final.contigs.fa"), emit: megahit_assembly

    script:
    """
    megahit -r ${reads} --out-dir ${sample}_megahit --presets meta-sensitive --num-cpu-threads ${task.cpus}
    """
}

process consensus_diamond_daa {
    tag "${sample}"
    cpus "${params.cpusMedium}"
    publishDir "$params.outdir/diamond_megan/", mode: 'symlink'
    
    input: 
    tuple val(sample), val(contigs)

    output:
    tuple val("${sample}"), val("${sample}.daa"), emit: blastx_meganizer_results
    

    script:
    """
    diamond blastx -q ${contigs} -d ${params.diamond_db_nr} -o ${sample}.daa -F 15 -f 100 --range-culling --top 10 -p ${task.cpus}
    """
}


process diamond_txt {
    tag "${sample}"
    cpus "${params.cpusMedium}"
    publishDir "$params.outdir/diamond_txt/", mode: 'symlink'
    
    input: 
    tuple val(sample), path(reads)

    output:
    tuple val("${sample}"), path("${sample}_blastx_out.txt"), emit: blastx_results
    

    script:
    """
    diamond blastx -d ${params.diamond_db_nr} -q ${reads} -o ${sample}_blastx_out.txt -p ${task.cpus}  --evalue 0.00001 -f 6 qseqid sseqid pident length evalue bitscore staxids stitle qcovhsp 
    """
}

process diamond_daa {
    tag "${sample}"
    cpus "${params.cpusMedium}"
    publishDir "$params.outdir/diamond_daa/", mode: 'symlink'
    
    input: 
    tuple val(sample), path(reads)

    output:
    tuple val("${sample}"), path ("${sample}_diamond_out.daa"), emit: blastx_meganizer_results
    

    script:
    """
    diamond blastx -d ${params.diamond_db_nr} -q ${reads} -o ${sample}_diamond_out.daa -p ${task.cpus} -f 100
    """
}

process meganizer {
    tag "${sample}"
    cpus "${params.cpusMin}"
    maxRetries 5
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    // publishDir "$params.outdir/diamond_daa/meganizer", mode: 'symlink'
    
    input: 
    tuple val(sample), path(blastx_daa)

    // output:
    // tuple val("${sample}"), path ("${sample}_diamond_out.daa"), emit: meganizer_results

    script:
    """
    daa-meganizer -i ${blastx_daa} -mdb ${params.meganizer_db} -lg -t ${task.cpus}
    """
}


// process add_kingdom {
//     tag "${sample}"
//     cpus "${params.cpusMin}"
//     publishDir "$params.outdir/diamond_txt/Kingdom", mode: 'symlink'
    
//     input: 
//     tuple val(sample), val(blastx)

//     output:
//     path ("${sample}_blastx_kingdom.txt"), emit: kingdoms_blastx
    

//     script:
//     """
//     python ${params.py_scripts}/add_kingdoms_blastx.py  -b  ${blastx}  -o ${sample}_blastx_kingdom.txt
//     """
// }