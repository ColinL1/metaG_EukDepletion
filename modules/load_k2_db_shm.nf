#!/usr/bin/env nextflow
// load the kr2 db into shared memory, then remove it after running k2. NOT WORKING 

process load_in_shm {
    tag ""
    label "min_mem"

    input: 
    path(params.kraken_db)

    output:
    path(""), emit: alpha_diversity_report

    script:
    """
        FILE=/dev/shm/*.k2d
    if [ -f "$FILE" ]; then
        echo "db in memory"
    else 
        cp ${params.kraken_db}/* /dev/shm/.
    fi
    """
}


process clean_shm {
    tag ""
    label "min_mem"

    input: 
    tuple val(sample), path(bracken_report)

    output:
    path(""), emit: alpha_diversity_report

    script:
    """
    rm -rf /dev/shm/* 
    """
}