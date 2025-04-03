#!/usr/bin/env nextflow

params.ref_file_aip = '/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_bowtie2/all_aiptasiidae*'
params.ref_file_scl = '/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_bowtie2/all_scleractina*'

params.ref_file_2 = '/home/colinl/metaG/Git/metaG_EukDepletion/references/ref_genomes_bowtie2/all_symbiodiniaceae*'

include { BOWTIE2_MAP_HOST ; BOWTIE2_MAP_SYM } from "$baseDir/modules/bowtie2.nf"
include {KAIJU_PE} from '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/kaiju.nf'

params.nodes = '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/bacteria_kaiju/nodes.dmp'
params.names = '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/bacteria_kaiju/names.dmp'
params.kaiju_db = '/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/bacteria_kaiju/kaiju_db_refseq.fmi'

params.input_fq = "/home/colinl/metaG/Git/metaG_EukDepletion/manual_piplines/reads/input/*_{1,2}.fq.gz"

// Channel
//     .fromFilePairs(params.input_fq)
//     .map { id, reads ->
//         (species, replicate, method, buffer, other, unspecified) = id.tokenize("_")
//         meta = [
//             id:id,
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


workflow MAP_CONSE_PE {
    take: 
        reads
    main:
        BOWTIE2_MAP_HOST(reads)
            BOWTIE2_MAP_SYM(BOWTIE2_MAP_HOST.out.unmapped_reads)
                KAIJU_PE(BOWTIE2_MAP_SYM.out.unmapped_reads)
    emit:
        kaiju_report = KAIJU_PE.out.report_kaiju_out
        
        bowtie2_host_bam = BOWTIE2_MAP_HOST.out.mapp_file
        bowtie2_host_map_fq = BOWTIE2_MAP_HOST.out.mapped_reads
        bowtie2_host_unmap_fq = BOWTIE2_MAP_HOST.out.unmapped_reads

        bowtie2_sym_bam = BOWTIE2_MAP_SYM.out.mapp_file
        bowtie2_sym_map_fq = BOWTIE2_MAP_SYM.out.mapped_reads
        bowtie2_sym_unmap_fq = BOWTIE2_MAP_SYM.out.unmapped_reads
}
