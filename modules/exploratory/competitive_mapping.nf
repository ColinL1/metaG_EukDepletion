#!/usr/bin/env nextflow
// Competitive mapping against a combined reference. Not used by main.nf.

process BOWTIE2_MAP_COMPETITIVE {
    tag "${meta.id}"
    label "process_high"
    conda "bioconda::bowtie2 bioconda::samtools"
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.bam"),                        emit: bam
    tuple val(meta), path("${meta.id}.primary.bam"),                emit: bam_primary
    tuple val(meta), path("${meta.id}.competitive_counts.tsv"),     emit: counts
    tuple val(meta), path("*.unmapped_{1,2}.fq.gz"),                emit: unmapped_reads

    script:
    """
    # Find combined coral+sponge index (supports both .bt2 and .bt2l)
    INDEX=`find -L ${params.competitive_index} -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ${params.competitive_index} -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

    # Map to combined reference, retain all alignments in BAM
    bowtie2 --mm \\
        -x \$INDEX \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --un-conc-gz ${meta.id}.unmapped.fq.gz \\
        -p ${task.cpus} \\
        | samtools view -bS -F 4 - \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.bam

    # Rename unmapped files to expected pattern
    mv ${meta.id}.unmapped.fq.1.gz ${meta.id}.unmapped_1.fq.gz
    mv ${meta.id}.unmapped.fq.2.gz ${meta.id}.unmapped_2.fq.gz

    samtools view -b -F 2308 ${meta.id}.bam \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.primary.bam
    samtools index ${meta.id}.primary.bam

    echo -e "sample\\tgroup\\tread_count" > ${meta.id}.competitive_counts.tsv
    samtools idxstats ${meta.id}.primary.bam \\
        | awk -v sample="${meta.id}" '
            /^CORAL_/  { coral  += \$3 }
            /^Scleractinia_/  { coral += \$3 }
            /^Porifera_/ { sponge += \$3 }
            /^Symbiodiniaceae_/ { symbiodiniaceae += \$3 }
            END {
                print "${meta.id}\tcoral\t"  (coral  ? coral  : 0)
                print "${meta.id}\tsponge\t" (sponge ? sponge : 0)
                print "${meta.id}\tsymbiodiniaceae\t" (symbiodiniaceae ? symbiodiniaceae : 0)
            }
        ' >> ${meta.id}.competitive_counts.tsv
    """

    stub:
    """
    touch ${meta.id}.bam
    touch ${meta.id}.primary.bam
    touch ${meta.id}.competitive_counts.tsv
    touch ${meta.id}.unmapped_1.fq.gz ${meta.id}.unmapped_2.fq.gz
    """
}

process COMPETITIVE_MAP_ONT {
    tag "${meta.id}"
    label "process_high"
    conda "bioconda::minimap2 bioconda::samtools"
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:b411340b52d82a9c276d87c7a3dcffc880be762f-0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: mapp_file
    tuple val(meta), path("*.competitive_counts.tsv"), emit: counts

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont ${params.competitive_index_ont} \\
    ${reads} --split-prefix=tmp | samtools view -S -b > ${meta.id}.bam

    samtools view -b -F 2308 ${meta.id}.bam \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.primary.bam
    samtools index ${meta.id}.primary.bam

    echo -e "sample\\tgroup\\tread_count" > ${meta.id}.competitive_counts.tsv
    samtools idxstats ${meta.id}.primary.bam \\
        | awk -v sample="${meta.id}" '
            /^CORAL_/  { coral  += \$3 }
            /^Scleractinia_/  { coral += \$3 }
            /^Porifera_/ { sponge += \$3 }
            /^Symbiodiniaceae_/ { symbiodiniaceae += \$3 }
            END {
                print "${meta.id}\\tcoral\\t"  (coral  ? coral  : 0)
                print "${meta.id}\\tsponge\\t" (sponge ? sponge : 0)
                print "${meta.id}\\tsymbiodiniaceae\\t" (symbiodiniaceae ? symbiodiniaceae : 0)
            }
        ' >> ${meta.id}.competitive_counts.tsv

    """
    stub:
    """
    touch ${meta.id}.bam
    touch ${meta.id}.competitive_counts.tsv
    touch ${meta.id}.primary.bam
    """
}
process COMPETITIVE_MAP_MINPE {
    tag "${meta.id}"
    label "process_high"
    conda "bioconda::minimap2 bioconda::samtools"
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:b411340b52d82a9c276d87c7a3dcffc880be762f-0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: mapp_file
    tuple val(meta), path("*.competitive_counts.tsv"), emit: counts

    script:
    """
    minimap2 -t ${task.cpus} -ax sr ${params.competitive_index} \\
    ${reads[0]} ${reads[1]} --split-prefix=tmp | samtools view -S -b > ${meta.id}.bam

    samtools view -b -F 2308 ${meta.id}.bam \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.primary.bam
    samtools index ${meta.id}.primary.bam

    echo -e "sample\\tgroup\\tread_count" > ${meta.id}.competitive_counts.tsv
    samtools idxstats ${meta.id}.primary.bam \\
        | awk -v sample="${meta.id}" '
            /^CORAL_/  { coral  += \$3 }
            /^Scleractinia_/  { coral += \$3 }
            /^Porifera_/ { sponge += \$3 }
            /^Symbiodiniaceae_/ { symbiodiniaceae += \$3 }
            END {
                print "${meta.id}\\tcoral\\t"  (coral  ? coral  : 0)
                print "${meta.id}\\tsponge\\t" (sponge ? sponge : 0)
                print "${meta.id}\\tsymbiodiniaceae\\t" (symbiodiniaceae ? symbiodiniaceae : 0)
            }
        ' >> ${meta.id}.competitive_counts.tsv

    """
    stub:
    """
    touch ${meta.id}.bam
    touch ${meta.id}.competitive_counts.tsv
    touch ${meta.id}.primary.bam
    """
}
