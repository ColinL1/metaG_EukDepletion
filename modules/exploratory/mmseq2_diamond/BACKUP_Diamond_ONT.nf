#!/usr/bin/env nextflow

process clean_names {
    tag "${sample}"
    cpus "${params.cpusMin}"
    publishDir "$params.outdir/fixed_names/", mode: 'symlink'

    input: 
    tuple val(sample), val(reads)

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
    tuple val(sample), val(reads)

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

process unicycler {
    tag "${sample}"
    cpus "${params.cpusVHigh}"
    errorStrategy "ignore"
    publishDir "$params.outdir/unicycler/", mode: 'symlink'
    
    input: 
    tuple val(sample), val(reads)

    output:
    tuple val("${sample}"), path ("${sample}/assembly.fasta"), emit: unicycler_assembly
    tuple val("${sample}"), path ("${sample}_medaka/consensus.fasta"), emit: medaka_consensus
    

    script:
    """
    unicycler -l ${reads} -o ${sample} -t ${task.cpus} --keep 3
    medaka_consensus -i ${reads} -d ${sample}/assembly.fasta -o ${sample}_medaka -t ${task.cpus} -m r941_min_sup_g507
    """
}

process canu {
    tag "${sample}"
    cpus "${params.cpusVHigh}"
    errorStrategy "ignore"
    publishDir "$params.outdir/canu/", mode: 'symlink'
    
    input: 
    tuple val(sample), val(reads)

    output:
    tuple val("${sample}"), path ("${sample}_canu/${sample}.contigs.fasta"), emit: canu_assembly
    tuple val("${sample}"), path ("${sample}_medaka/consensus.fasta"), emit: medaka_consensus
    

    script:
    """
    canu -p ${sample} -d ${sample}_canu genomeSize=4.8m correctedErrorRate=0.16 maxInputCoverage=50 stopOnLowCoverage=5 maxThreads=${task.cpus} -nanopore ${reads} 
    medaka_consensus -i ${reads} -d ${sample}_canu/${sample}.contigs.fasta -o ${sample}_medaka -t ${task.cpus} -m r941_min_sup_g507
    """
}

process flye {
    tag "${sample}"
    cpus "${params.cpusVHigh}"
    errorStrategy "ignore"
    publishDir "$params.outdir/flye/", mode: 'symlink'
    
    input: 
    tuple val(sample), val(reads)

    output:
    tuple val("${sample}"), path ("${sample}_flye/assembly.fasta"), emit: flye_assembly
    tuple val("${sample}"), path ("${sample}_medaka/consensus.fasta"), emit: medaka_consensus
    

    script:
    """
    flye --nano-hq ${reads} --out-dir ${sample}_flye --genome-size 4.8m --threads ${task.cpus} --meta
    medaka_consensus -i ${reads} -d ${sample}_flye/assembly.fasta -o ${sample}_medaka -t ${task.cpus} -m r941_min_sup_g507
    """
}

process megahit {
    tag "${sample}"
    cpus "${params.cpusVHigh}"
    publishDir "$params.outdir/megahit/", mode: 'symlink'

    input:
    tuple val(sample), path(reads) 

    output:
    tuple val("${sample}"), path ("${sample}_megahit/final.contigs.fa"), emit: megahit_assembly
    tuple val("${sample}"), path ("${sample}_medaka/consensus.fasta"), emit: medaka_consensus

    script:
    """
    megahit -r ${reads} --out-dir ${sample}_megahit --presets meta-sensitive --num-cpu-threads ${task.cpus}
    medaka_consensus -i ${reads} -d ${sample}_megahit/final.contigs.fa -o ${sample}_medaka -t ${task.cpus} -m r941_min_sup_g507
    """
}

process consensus_diamond_meganizer {
    tag "${sample}"
    cpus "${params.cpusHigh}"
    publishDir "$params.outdir/diamond/", mode: 'symlink'
    
    input: 
    tuple val(sample), val(contigs)

    output:
    tuple val("${sample}"), path ("${sample}_diamond_out.daa"), emit: blastx_meganizer_results
    

    script:
    """
    diamond blastx -q ${contigs} -d ${params.diamond_db_nr} -o ${sample}.daa -F 15 -f 100 --range-culling --top 10 -p ${task.cpus}
    """
}


process diamond_blastx {
    tag "${sample}"
    cpus "${params.cpusHigh}"
    publishDir "$params.outdir/diamond_blastx/", mode: 'symlink'
    
    input: 
    tuple val(sample), val(reads)

    output:
    tuple val("${sample}"), path ("${sample}_blastx_out.txt"), emit: blastx_results
    

    script:
    """
    diamond blastx -d ${params.diamond_db_nr} -q ${reads} -o ${sample}_blastx_out.txt -p ${task.cpus}  --evalue 0.00001 -f 6 qseqid sseqid pident length evalue bitscore staxids stitle qcovhsp 
    """
}

process diamond_meganizer {
    tag "${sample}"
    cpus "${params.cpusHigh}"
    publishDir "$params.outdir/diamond_meganizer/", mode: 'symlink'
    
    input: 
    tuple val(sample), val(reads)

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
    publishDir "$params.outdir/diamond_meganizer/meganizer", mode: 'symlink'
    
    input: 
    tuple val(sample), val(reads)

    // output:
    // tuple val("${sample}"), path ("${sample}_diamond_out.daa"), emit: meganizer_results

    script:
    """
    daa-meganizer -i ${reads} -mdb ${params.meganizer_db} -lg -t ${task.cpus}
    """
}


// process add_kingdom {
//     tag "${sample}"
//     cpus "${params.cpusMin}"
//     publishDir "$params.outdir/diamond_blastx/Kingdom", mode: 'symlink'
    
//     input: 
//     tuple val(sample), val(blastx)

//     output:
//     path ("${sample}_blastx_kingdom.txt"), emit: kingdoms_blastx
    

//     script:
//     """
//     python ${params.py_scripts}/add_kingdoms_blastx.py  -b  ${blastx}  -o ${sample}_blastx_kingdom.txt
//     """
// }