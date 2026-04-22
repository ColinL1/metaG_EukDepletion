#!/usr/bin/env nextflow

process sort_contigs_kingdom {
    tag "${sample}"
    cpus "${params.cpusMin}"
    publishDir "$params.outdir/diamond_blastx/Kingdom", mode: 'symlink'
    
    input: 
    tuple val(sample), val(fasta)
    tuple val(blast_name), val(blastx)

    output:
    path ("${sample}"), emit: sorted_fasta
    

    script:
    """
    python ${params.py_scripts}/sort_contigs_kingdom.py -c ${fasta} -b ${sample}_kingdom.txt -o ${sample}
    """
}


// process sort_contigs_genus {
//     tag "${sample}"
//     cpus "${params.cpusMin}"
//     publishDir "$params.outdir/diamond_blastx/Kingdom", mode: 'symlink'
    
//     input: 
//     tuple val(sample), val(fasta)
//     //tuple val(blast_name), val(blastx)

//     output:
//     path ("${sample}_blastx_kingdom.txt"), emit: kingdoms_blastx
    

//     script:
//     """
//     grep ">" ${fasta} | awk '{print \$1}' | sed 's/>//g' > ${fasta}_ids.txt
    
//     grep -f ${fasta}_ids.txt $SAMPLE_NAME.kingdom > $SAMPLE_NAME_1.kingdom

    
    
//     python ${params.py_scripts}/sort_contigs_kingdom.py -c ${fasta} -b ${sample}_kingdom.txt -o ${sample}
//     """
// }