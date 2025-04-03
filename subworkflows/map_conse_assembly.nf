#!/usr/bin/env nextflow


//TODO: check if possible to inherit conditional directly from queue channel  
process MINIMAP2_MAP {
    tag "${meta.id}"
    label "process_high"
    publishDir "${
        meta.species == 'H2' ? "$baseDir/mapping/Aiptasia/${meta.id}/bam/" :
        meta.species == 'F003' ? "$baseDir/mapping/Aiptasia/${meta.id}/bam/" :
        meta.species == 'Porites' ? "$baseDir/mapping/Corals/${meta.id}/bam/" :
        meta.species == 'Acropora' ? "$baseDir/mapping/Corals/${meta.id}/bam/" :
        meta.species == 'Pocillopora' ? "$baseDir/mapping/Corals/${meta.id}/bam/" :
        meta.species == 'unknown' ? 'unknown_samples' :
        'other_samples'}", mode: 'symlink'

    // publishDir "$baseDir/mapping/${meta.species}/bam/${meta.id}/", mode: 'symlink'
    // publishDir "$baseDir/mapping/all_aiptasiidae/${meta.id}/reads/${meta.species}/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)
    // tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: mapp_file
    // tuple val(meta), path("*.unmapped_{1,2}.fq.gz"), emit: unmapped_reads
    // tuple val(meta), path("*.mapped_{1,2}.fq.gz"), emit: mapped_reads

    script:
    if( "${meta.species}" == 'H2' )
        """
        minimap2 -t ${task.cpus} -ax asm5 ${params.ref_file_aip} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.aiptasiidae.bam
        """
    else if( "${meta.species}" == 'F003' )
        """
        minimap2 -t ${task.cpus} -ax asm5 ${params.ref_file_aip} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.aiptasiidae.bam
        """
    else if( "${meta.species}" == 'Porites' )
        """
	    minimap2 -t ${task.cpus} -ax asm5 ${params.ref_file_scl} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.scleractina.bam
        
        """
    else if ( "${meta.species}" == 'Acropora' )
        """
	    minimap2 -t ${task.cpus} -ax asm5 ${params.ref_file_scl} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.scleractina.bam
        """
    else if ( "${meta.species}" == 'Pocillopora' )
        """
	    minimap2 -t ${task.cpus} -ax asm5 ${params.ref_file_scl} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.scleractina.bam
        """
    else
    error "Invalid alignment mode: ${meta.species}"

    stub:
    if( "${meta.species}" == 'H2' )
        """
        touch ${meta.id}.${meta.species}-map-H2.bam
        touch ${meta.id}
        """
    else if( "${meta.species}" == 'F003' )
        """
        touch ${meta.id}.${meta.species}-map-F003.bam
        touch ${meta.id}
        """
    else if( "${meta.species}" == 'Porites' )
        """
        touch ${meta.id}.${meta.species}-map-Porites.bam
        touch ${meta.id}
        """
    else if ( "${meta.species}" == 'Acropora' )
        """
        touch ${meta.id}.${meta.species}-map-Acropora.bam
        touch ${meta.id}
        """
    else if ( "${meta.species}" == 'Pocillopora' )
        """
        touch ${meta.id}.${meta.species}-map-Pocillopora.bam
        touch ${meta.id}
        """
    else
    error "Invalid alignment mode: ${meta.species}"

}

params.ref_file_aip = '/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_minimap2_contigs/all_aiptasiidae.mmi'
params.ref_file_scl = '/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_minimap2_contigs/all_scleractinia.mmi'

process SPLIT_BAM_1 {
    tag "${meta.id}"
    label "med_mem"
    // publishDir "$baseDir/mapping/all_scleractina/${meta.id}/bam/${meta.species}/", mode: 'symlink'
    // publishDir "$baseDir/mapping/${meta.species}/reads/${meta.id}/", mode: 'symlink'
    publishDir "${
    meta.species == 'H2' ? "$baseDir/mapping/Aiptasia/${meta.id}/reads/" :
    meta.species == 'F003' ? "$baseDir/mapping/Aiptasia/${meta.id}/reads/" :
    meta.species == 'Porites' ? "$baseDir/mapping/Corals/${meta.id}/reads/" :
    meta.species == 'Acropora' ? "$baseDir/mapping/Corals/${meta.id}/reads/" :
    meta.species == 'Pocillopora' ? "$baseDir/mapping/Corals/${meta.id}/reads/" :
    meta.species == 'unknown' ? 'unknown_samples' :
    'other_samples'}", mode: 'symlink'

    input:
    tuple val(meta), path(bam)
    // tuple val(meta2), path(index)

    output:
    // tuple val(meta), path("*.bam"), emit: mapp_file
    tuple val(meta), path("*.unmapped.fa.gz"), emit: unmapped_reads
    tuple val(meta), path("*.mapped.fa.gz"), emit: mapped_reads

    script:
    """
    samtools view -F 0x900 -F 0x4 ${bam} --threads ${task.cpus} -o ${meta.id}.scleractina.mapped.fa
    samtools view -f 0x4 -F 0x900 ${bam} --threads ${task.cpus} -o ${meta.id}.scleractina.unmapped.fa
    pigz -p ${task.cpus} ${meta.id}.scleractina.mapped.fa ${meta.id}.scleractina.unmapped.fa
    """
    stub:
    """
    touch ${meta.id}.unmapped.fa.gz
    touch ${meta.id}.mapped.fa.gz
    """
}


process MINIMAP2_MAP_2 {
    tag "${meta.id}"
    label "process_high"
    publishDir "$baseDir/mapping/Symbiodiniaceae/${meta.id}/bam/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)
    // tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: mapp_file
    // tuple val(meta), path("*.unmapped_{1,2}.fq.gz"), emit: unmapped_reads
    // tuple val(meta), path("*.mapped_{1,2}.fq.gz"), emit: mapped_reads

    script:
	"""
	minimap2 -t ${task.cpus} -ax asm5 ${params.ref_file_2} ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.symbiodiniaceae.bam
	"""
    stub:
    """
    touch ${meta.id}.bam
    touch ${meta.id}
    """
}

params.ref_file_2 = '/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_minimap2_contigs/all_symbiodiniaceae.mmi'

process SPLIT_BAM_2 {
    tag "${meta.id}"
    label "med_mem"
    publishDir "$baseDir/mapping/Symbiodiniaceae/${meta.id}/reads/", mode: 'symlink'

    input:
    tuple val(meta), path(bam)
    // tuple val(meta2), path(index)

    output:
    // tuple val(meta), path("*.bam"), emit: mapp_file
    tuple val(meta), path("*.unmapped.fa.gz"), emit: unmapped_reads
    tuple val(meta), path("*.mapped.fa.gz"), emit: mapped_reads

    script:
    """
    samtools view -F 0x900 -F 0x4 ${bam} --threads ${task.cpus} -o ${meta.id}.symbiodiniaceae.mapped.fa
    samtools view -f 0x4 -F 0x900 ${bam} --threads ${task.cpus} -o ${meta.id}.symbiodiniaceae.unmapped.fa
    pigz -p ${task.cpus} ${meta.id}.symbiodiniaceae.mapped.fa ${meta.id}.symbiodiniaceae.unmapped.fa
    """
    stub:
    """
    touch ${meta.id}.unmapped.fa.gz
    touch ${meta.id}.mapped.fa.gz
    """
}


include {KAIJU_SE} from '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/kaiju.nf'

params.nodes = '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/bacteria_kaiju/nodes.dmp'
params.names = '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/bacteria_kaiju/names.dmp'
params.kaiju_db = '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/bacteria_kaiju/kaiju_db_refseq.fmi'


params.input = "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/assembly/assembly/*/*/final.contigs.fa"

Channel
    .fromPath(params.input)
    .map { reads ->
        (species, replicate, method, buffer, other, unspecified) = reads.getParent().getName().tokenize("_")
        meta = [
            id:reads.getParent().getName().minus(~/_TRIM_paried/),
            // single_end:'PE',
            species:species,
            replicate:replicate,
            method:method,
            buffer:buffer,
            other:other,
            unspecified:unspecified,
        ]
        [meta, reads]
    }
    .map { meta, reads ->
        def newmap = [
            species: meta.species == "Po" ? "Porites" :
                    meta.species == "Por" ? "Porites" :
                    meta.species == "Ac" ? "Acropora" :
                    meta.species == "Acro" ? "Acropora" :
                    meta.species == "F003" ? "F003" :
                    meta.species == "F3" ? "F003" :
                    meta.species == "H2" ? "H2" :
                    meta.species == "Poci" ? "Pocillopora" :
                    meta.species == "Pr" ? "Pocillopora" :
                    "Unknown"
        ]
        [meta + newmap, reads]
    }
    .set { samples_coral }

workflow {
    // samples_coral.view()
    MINIMAP2_MAP(samples_coral) 
        SPLIT_BAM_1(MINIMAP2_MAP.out.mapp_file)
            MINIMAP2_MAP_2(SPLIT_BAM_1.out.unmapped_reads) // was run against mapped reads. not good. 
                SPLIT_BAM_2(MINIMAP2_MAP_2.out.mapp_file)
                    KAIJU_SE(SPLIT_BAM_2.out.unmapped_reads)
}



// workflow {
//     Channel.fromPath(params.input)
//         | map { reads ->
//             (species, replicate, method, buffer, other, unspecified) = reads.getParent().getName().tokenize("_")
//             meta = [
//                 id:reads.getParent().getName().minus(~/_TRIM_paried/),
//                 // single_end:'PE',
//                 species:species,
//                 replicate:replicate,
//                 method:method,
//                 buffer:buffer,
//                 other:other,
//                 unspecified:unspecified,
//             ]
//             [meta, reads]
//     }
//     | view()
// }