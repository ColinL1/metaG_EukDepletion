#!/usr/bin/env nextflow

process BOWTIE2_MAP {
    tag "${meta.id}"
    label "process_high"
    publishDir "${
            meta.species == 'H2' ? "$baseDir/mapping/Aiptasia/${meta.id}/" :
            meta.species == 'F003' ? "$baseDir/mapping/Aiptasia/${meta.id}/" :
            meta.species == 'Porites' ? "$baseDir/mapping/Corals/${meta.id}/" :
            meta.species == 'Acropora' ? "$baseDir/mapping/Corals/${meta.id}/" :
            meta.species == 'Pocillopora' ? "$baseDir/mapping/Corals/${meta.id}/" :
            meta.species == 'unknown' ? 'unknown_species' :
            'other_samples'}", mode: 'symlink'

    input:
    tuple val(meta), path(reads)
    // tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: mapp_file
    tuple val(meta), path("*.unmapped_{1,2}.fq.gz"), emit: unmapped_reads
    tuple val(meta), path("*.mapped_{1,2}.fq.gz"), emit: mapped_reads

    script:
    if( "${meta.species}" == 'H2' )
        """
        INDEX=`find -L ${params.ref_file_aip} -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
        [ -z "\$INDEX" ] && INDEX=`find -L ${params.ref_file_aip} -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
        [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

        bowtie2 -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} --un-conc-gz ${meta.id}.scleractina.unmapped.fq.gz --al-conc-gz ${meta.id}.scleractina.mapped.fq.gz -p ${task.cpus} | samtools view -S -b > ${meta.id}.bam
        ls *.1.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.1.gz/_1.fq.gz/'g) ; done
        ls *.2.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.2.gz/_2.fq.gz/'g) ; done
        """
    else if( "${meta.species}" == 'F003' )
        """
        INDEX=`find -L ${params.ref_file_aip} -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
        [ -z "\$INDEX" ] && INDEX=`find -L ${params.ref_file_aip} -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
        [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

        bowtie2 -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} --un-conc-gz ${meta.id}.scleractina.unmapped.fq.gz --al-conc-gz ${meta.id}.scleractina.mapped.fq.gz -p ${task.cpus} | samtools view -S -b > ${meta.id}.bam
        ls *.1.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.1.gz/_1.fq.gz/'g) ; done
        ls *.2.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.2.gz/_2.fq.gz/'g) ; done
        """
    else if( "${meta.species}" == 'Porites' )
        """
        INDEX=`find -L ${params.ref_file_scl} -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
        [ -z "\$INDEX" ] && INDEX=`find -L ${params.ref_file_scl} -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
        [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

        bowtie2 -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} --un-conc-gz ${meta.id}.scleractina.unmapped.fq.gz --al-conc-gz ${meta.id}.scleractina.mapped.fq.gz -p ${task.cpus} | samtools view -S -b > ${meta.id}.bam
        ls *.1.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.1.gz/_1.fq.gz/'g) ; done
        ls *.2.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.2.gz/_2.fq.gz/'g) ; done        
        """
    else if ( "${meta.species}" == 'Acropora' )
        """
        INDEX=`find -L ${params.ref_file_scl} -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
        [ -z "\$INDEX" ] && INDEX=`find -L ${params.ref_file_scl} -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
        [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

        bowtie2 -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} --un-conc-gz ${meta.id}.scleractina.unmapped.fq.gz --al-conc-gz ${meta.id}.scleractina.mapped.fq.gz -p ${task.cpus} | samtools view -S -b > ${meta.id}.bam
        ls *.1.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.1.gz/_1.fq.gz/'g) ; done
        ls *.2.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.2.gz/_2.fq.gz/'g) ; done
        """
    else if ( "${meta.species}" == 'Pocillopora' )
        """
        INDEX=`find -L ${params.ref_file_scl} -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
        [ -z "\$INDEX" ] && INDEX=`find -L ${params.ref_file_scl} -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
        [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

        bowtie2 -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} --un-conc-gz ${meta.id}.scleractina.unmapped.fq.gz --al-conc-gz ${meta.id}.scleractina.mapped.fq.gz -p ${task.cpus} | samtools view -S -b > ${meta.id}.bam
        ls *.1.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.1.gz/_1.fq.gz/'g) ; done
        ls *.2.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.2.gz/_2.fq.gz/'g) ; done
        """
    else
    error "Invalid alignment mode: ${meta.species}"

    stub:
    if( "${meta.species}" == 'H2' )
        """
        touch ${meta.id}-map-H2.bam
        touch ${meta.id}-map-H2.unmapped_1.fq.gz ${meta.id}-map-H2.unmapped_2.fq.gz
        touch ${meta.id}-map-H2.mapped_1.fq.gz ${meta.id}-map-H2.mapped_2.fq.gz
        """
    else if( "${meta.species}" == 'F003' )
        """
        touch ${meta.id}-map-F003.bam
        touch ${meta.id}-map-F003.unmapped_1.fq.gz ${meta.id}-map-F003.unmapped_2.fq.gz
        touch ${meta.id}-map-F003.mapped_1.fq.gz ${meta.id}-map-F003.mapped_2.fq.gz
        """
    else if( "${meta.species}" == 'Porites' )
        """
        touch ${meta.id}-map-Porites.bam
        touch ${meta.id}-map-Porites.unmapped_1.fq.gz ${meta.id}-map-Porites.unmapped_2.fq.gz
        touch ${meta.id}-map-Porites.mapped_1.fq.gz ${meta.id}-map-Porites.mapped_2.fq.gz
        """
    else if ( "${meta.species}" == 'Acropora' )
        """
        touch ${meta.id}-map-Acropora.bam
        touch ${meta.id}-map-Acropora.unmapped_1.fq.gz ${meta.id}-map-Acropora.unmapped_2.fq.gz
        touch ${meta.id}-map-Acropora.mapped_1.fq.gz ${meta.id}-map-Acropora.mapped_2.fq.gz
        """
    else if ( "${meta.species}" == 'Pocillopora' )
        """
        touch ${meta.id}-map-Pocillopora.bam
        touch ${meta.id}-map-Pocillopora.unmapped_1.fq.gz ${meta.id}-map-Pocillopora.unmapped_2.fq.gz
        touch ${meta.id}-map-Pocillopora.mapped_1.fq.gz ${meta.id}-map-Pocillopora.mapped_2.fq.gz
        """
    else
    error "Invalid alignment mode: ${meta.species}"

}

params.ref_file_aip = '/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_bowtie2/all_aiptasiidae*'
params.ref_file_scl = '/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_bowtie2/all_scleractina*'

process BOWTIE2_MAP_2 {
    tag "${meta.id}"
    label "process_high"
    publishDir "$baseDir/mapping/symbiodiniaceae/${meta.id}/", mode: 'symlink'

    input:
    tuple val(meta), path(reads)
    // tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: mapp_file
    tuple val(meta), path("*.unmapped_{1,2}.fq.gz"), emit: unmapped_reads
    tuple val(meta), path("*.mapped_{1,2}.fq.gz"), emit: mapped_reads

    script:
    """
    INDEX=`find -L ${params.ref_file_2} -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ${params.ref_file_2} -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

    bowtie2 -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} --un-conc-gz ${meta.id}.symbiodiniaceae.unmapped.fq.gz --al-conc-gz ${meta.id}.symbiodiniaceae.mapped.fq.gz -p ${task.cpus} | samtools view -S -b > ${meta.id}.bam
    ls *.1.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.1.gz/_1.fq.gz/'g) ; done
    ls *.2.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.2.gz/_2.fq.gz/'g) ; done
    """
    stub:
    """
    touch ${meta.id}.bam
    touch ${meta.id}.unmapped_1.fq.gz ${meta.id}.unmapped_2.fq.gz
    touch ${meta.id}.mapped_1.fq.gz ${meta.id}.mapped_2.fq.gz
    """
}

params.ref_file_2 = '/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_bowtie2/all_symbiodiniaceae*'

include {KAIJU_PE} from '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/kaiju.nf'

params.nodes = '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/bacteria_kaiju/nodes.dmp'
params.names = '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/bacteria_kaiju/names.dmp'
params.kaiju_db = '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/bacteria_kaiju/kaiju_db_refseq.fmi'

params.input = "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/input/*_{1,2}.fq.gz"

Channel
    .fromFilePairs(params.input)
    .map { id, reads ->
        (species, replicate, method, buffer, other, unspecified) = id.tokenize("_")
        meta = [
            id:id,
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
            BOWTIE2_MAP(samples_coral)
                BOWTIE2_MAP_2(BOWTIE2_MAP.out.unmapped_reads)
                    KAIJU_PE(BOWTIE2_MAP_2.out.unmapped_reads)
}
