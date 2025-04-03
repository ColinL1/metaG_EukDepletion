#!/usr/bin/env nextflow


//TODO: check if possible to inherit conditional directly from queue channel  

params.ref_file_aip = '/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_minimap2_contigs/all_aiptasiidae.mmi'
params.ref_file_scl = '/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_minimap2_contigs/all_scleractinia.mmi'

params.ref_file_2 = '/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_minimap2_contigs/all_symbiodiniaceae.mmi'

include { MINIMAP2_MAP_HOST ; MINIMAP2_MAP_SYM } from "$baseDir/modules/minimap2.nf"
include { SPLIT_BAM_HOST ; SPLIT_BAM_SYM } from "$baseDir/modules/samtools_split.nf"
include {KAIJU_SE} from "$baseDir/modules/kaiju.nf"
include { MEGAHIT_PE } from "$baseDir/modules/assembly_megahit.nf"

params.nodes = '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/bacteria_kaiju/nodes.dmp'
params.names = '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/bacteria_kaiju/names.dmp'
params.kaiju_db = '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/bacteria_kaiju/kaiju_db_refseq.fmi'

params.input = "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/assembly/assembly/*/*/final.contigs.fa"

// Channel
//     .fromPath(params.input)
//     .map { reads ->
//         (species, replicate, method, buffer, other, unspecified) = reads.getParent().getName().tokenize("_")
//         meta = [
//             id:reads.getParent().getName().minus(~/_TRIM_paried/),
//             // single_end:'PE',
//             species:species,
//             replicate:replicate,
//             method:method,
//             buffer:buffer,
//             other:other,
//             unspecified:unspecified,
//         ]
//         [meta, reads]
//     }
//     .map { meta, reads ->
//         def newmap = [
//             species: meta.species == "Po" ? "Porites" :
//                     meta.species == "Por" ? "Porites" :
//                     meta.species == "Ac" ? "Acropora" :
//                     meta.species == "Acro" ? "Acropora" :
//                     meta.species == "F003" ? "F003" :
//                     meta.species == "F3" ? "F003" :
//                     meta.species == "H2" ? "H2" :
//                     meta.species == "Poci" ? "Pocillopora" :
//                     meta.species == "Pr" ? "Pocillopora" :
//                     "Unknown"
//         ]
//         [meta + newmap, reads]
//     }
//     .set { samples_coral }

workflow MAP_CONSE_ASSEMBLY {
    take: 
        reads
    main:
    MEGAHIT_PE(reads)
        MINIMAP2_MAP_HOST(MEGAHIT_PE.out.contigs_fa) 
            SPLIT_BAM_HOST(MINIMAP2_MAP_HOST.out.mapp_file)
                MINIMAP2_MAP_SYM(SPLIT_BAM_HOST.out.unmapped_reads) // was run against mapped reads. not good. 
                    SPLIT_BAM_SYM(MINIMAP2_MAP_SYM.out.mapp_file)
                        KAIJU_SE(SPLIT_BAM_SYM.out.unmapped_reads)
    emit:
    kaiju_report = KAIJU_SE.out.report_kaiju_out

    bowtie2_host_bam = MINIMAP2_MAP_HOST.out.mapp_file
    bowtie2_host_map_fq = SPLIT_BAM_HOST.out.mapped_reads
    bowtie2_host_unmap_fq = SPLIT_BAM_HOST.out.unmapped_reads
    
    bowtie2_sym_bam = MINIMAP2_MAP_SYM.out.mapp_file
    bowtie2_sym_map_fq = SPLIT_BAM_SYM.out.mapped_reads
    bowtie2_sym_unmap_fq = SPLIT_BAM_SYM.out.unmapped_reads
}
