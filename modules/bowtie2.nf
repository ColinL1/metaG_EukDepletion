#!/usr/bin/env nextflow

process BOWTIE2_MAP_HOST {
    tag "${meta.id}"
    label "process_high"
    conda "bioconda::bowtie2 bioconda::samtools"
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'
    publishDir path: { "${params.outdir}/reads_based/mapping/${meta.species}/${meta.id}/" }, mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: mapp_file
    tuple val(meta), path("*.unmapped_{1,2}.fq.gz"), emit: unmapped_reads
    tuple val(meta), path("*.mapped_{1,2}.fq.gz"), emit: mapped_reads

    script:
    if( "${meta.species}" == 'H2' )
        """
        INDEX=`find -L ${params.ref_bt2_aip} -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
        [ -z "\$INDEX" ] && INDEX=`find -L ${params.ref_bt2_aip} -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
        [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

        bowtie2 --mm -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} --un-conc-gz ${meta.id}.host.unmapped.fq.gz --al-conc-gz ${meta.id}.host.mapped.fq.gz -p ${task.cpus} | samtools view -S -b > ${meta.id}.host.bam
        ls *.1.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.1.gz/_1.fq.gz/'g) ; done
        ls *.2.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.2.gz/_2.fq.gz/'g) ; done
        """
    else if( "${meta.species}" == 'F003' )
        """
        INDEX=`find -L ${params.ref_bt2_aip} -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
        [ -z "\$INDEX" ] && INDEX=`find -L ${params.ref_bt2_aip} -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
        [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

        bowtie2 --mm -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} --un-conc-gz ${meta.id}.host.unmapped.fq.gz --al-conc-gz ${meta.id}.host.mapped.fq.gz -p ${task.cpus} | samtools view -S -b > ${meta.id}.host.bam
        ls *.1.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.1.gz/_1.fq.gz/'g) ; done
        ls *.2.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.2.gz/_2.fq.gz/'g) ; done
        """
    else if( "${meta.species}" == 'Porites' )
        """
        INDEX=`find -L ${params.ref_bt2_scl} -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
        [ -z "\$INDEX" ] && INDEX=`find -L ${params.ref_bt2_scl} -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
        [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

        bowtie2 --mm -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} --un-conc-gz ${meta.id}.host.unmapped.fq.gz --al-conc-gz ${meta.id}.host.mapped.fq.gz -p ${task.cpus} | samtools view -S -b > ${meta.id}.host.bam
        ls *.1.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.1.gz/_1.fq.gz/'g) ; done
        ls *.2.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.2.gz/_2.fq.gz/'g) ; done        
        """
    else if ( "${meta.species}" == 'Acropora' )
        """
        INDEX=`find -L ${params.ref_bt2_scl} -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
        [ -z "\$INDEX" ] && INDEX=`find -L ${params.ref_bt2_scl} -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
        [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

        bowtie2 --mm -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} --un-conc-gz ${meta.id}.host.unmapped.fq.gz --al-conc-gz ${meta.id}.host.mapped.fq.gz -p ${task.cpus} | samtools view -S -b > ${meta.id}.host.bam
        ls *.1.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.1.gz/_1.fq.gz/'g) ; done
        ls *.2.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.2.gz/_2.fq.gz/'g) ; done
        """
    else if ( "${meta.species}" == 'Pocillopora' )
        """
        INDEX=`find -L ${params.ref_bt2_scl} -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
        [ -z "\$INDEX" ] && INDEX=`find -L ${params.ref_bt2_scl} -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
        [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

        bowtie2 --mm -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} --un-conc-gz ${meta.id}.host.unmapped.fq.gz --al-conc-gz ${meta.id}.host.mapped.fq.gz -p ${task.cpus} | samtools view -S -b > ${meta.id}.host.bam
        ls *.1.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.1.gz/_1.fq.gz/'g) ; done
        ls *.2.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.2.gz/_2.fq.gz/'g) ; done
        """
    else
        error "Invalid alignment mode: ${meta.species}"
    stub:
    if( "${meta.species}" == 'H2' )
        """
        touch ${meta.id}.host.bam
        touch ${meta.id}.host.unmapped_1.fq.gz ${meta.id}.host.unmapped_2.fq.gz
        touch ${meta.id}.host.mapped_1.fq.gz ${meta.id}.host.mapped_2.fq.gz
        """
    else if( "${meta.species}" == 'F003' )
        """
        touch ${meta.id}.host.bam
        touch ${meta.id}.host.unmapped_1.fq.gz ${meta.id}.host.unmapped_2.fq.gz
        touch ${meta.id}.host.mapped_1.fq.gz ${meta.id}.host.mapped_2.fq.gz
        """
    else if( "${meta.species}" == 'Porites' )
        """
        touch ${meta.id}.host.bam
        touch ${meta.id}.host.unmapped_1.fq.gz ${meta.id}.host.unmapped_2.fq.gz
        touch ${meta.id}.host.mapped_1.fq.gz ${meta.id}.host.mapped_2.fq.gz
        """
    else if ( "${meta.species}" == 'Acropora' )
        """
        touch ${meta.id}.host.bam
        touch ${meta.id}.host.unmapped_1.fq.gz ${meta.id}.host.unmapped_2.fq.gz
        touch ${meta.id}.host.mapped_1.fq.gz ${meta.id}.host.mapped_2.fq.gz
        """
    else if ( "${meta.species}" == 'Pocillopora' )
        """
        touch ${meta.id}.host.bam
        touch ${meta.id}.host.unmapped_1.fq.gz ${meta.id}.host.unmapped_2.fq.gz
        touch ${meta.id}.host.mapped_1.fq.gz ${meta.id}.host.mapped_2.fq.gz
        """
    else
        error "Invalid alignment mode: ${meta.species}"
}

process BOWTIE2_MAP_SYM {
    tag "${meta.id}"
    label "process_high"
    conda "bioconda::bowtie2 bioconda::samtools"
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'
    publishDir path: { "${params.outdir}/reads_based/mapping/symbiodiniaceae/${meta.id}/" }, mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: mapp_file
    tuple val(meta), path("*.unmapped_{1,2}.fq.gz"), emit: unmapped_reads
    tuple val(meta), path("*.mapped_{1,2}.fq.gz"), emit: mapped_reads

    script:
    """
    INDEX=`find -L ${params.ref_bt2_sym} -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ${params.ref_bt2_sym} -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

    bowtie2 --mm -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} --un-conc-gz ${meta.id}.symbiodiniaceae.unmapped.fq.gz --al-conc-gz ${meta.id}.symbiodiniaceae.mapped.fq.gz -p ${task.cpus} | samtools view -S -b > ${meta.id}.symbiodiniaceae.bam
    ls *.1.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.1.gz/_1.fq.gz/'g) ; done
    ls *.2.gz | while read line ; do mv \$line \$(echo \$line | sed s'/.fq.2.gz/_2.fq.gz/'g) ; done
    """
    stub:
    """
    touch ${meta.id}.symbiodiniaceae.bam
    touch ${meta.id}.symbiodiniaceae.unmapped_1.fq.gz ${meta.id}.symbiodiniaceae.unmapped_2.fq.gz
    touch ${meta.id}.symbiodiniaceae.mapped_1.fq.gz ${meta.id}.symbiodiniaceae.mapped_2.fq.gz
    """
}
