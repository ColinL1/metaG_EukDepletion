#!/usr/bin/env nextflow
// use NCBI dataset to downloa genome reference files for mapping

process NCBI_DOWNLOAD {
    tag "${sample}"

    publishDir "$params.outdir/reference_genomes/", mode: 'symlink'

    input: 
    tuple val(sample), val(base_name), val(taxon), val (args), 

    output:
    tuple val(sample), val(taxon), path("${taxon}_ncbi_*genome.fna"), emit: reference_genomes

    script:
    if(args == "none")
    """
    datasets download genome taxon ${taxon} --filename ${taxon}.zip
    unzip ${taxon}.zip
    cat ncbi_dataset/data/*/*.fna > ${taxon}_ncbi_genome.fna
    rm ${taxon}.zip
    rm -rf ncbi_dataset
    """
    else if (args == 'reference')
    """
    datasets download genome taxon ${taxon} --filename ${taxon}.zip --reference
    unzip ${taxon}.zip
    cat ncbi_dataset/data/*/*.fna > ${taxon}_ncbi_reference_genome.fna
    rm ${taxon}.zip
    rm -rf ncbi_dataset
    """
    stub:
    """
    touch ${taxon}_ncbi_genome.fna
    """
}