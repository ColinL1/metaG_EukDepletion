#!/usr/bin/env nextflow
// Load the Kraken2 db into shared memory before K2 runs, then clean up after.

process load_in_shm {
    tag "load_k2_db"
    label "min_mem"

    // val trigger: collected signal from all mapping steps; gates this process
    // so the DB is only loaded once all bowtie2/minimap2 work is done.
    input:
    val trigger

    output:
    val true, emit: db_loaded

    script:
    """
    if [ -f /dev/shm/${params.kraken_dbName}/hash.k2d ]; then
        echo "DB already in shared memory"
    else
        cp -rv ${params.kraken_db} /dev/shm/${params.kraken_dbName}
        echo "DB copied to /dev/shm/${params.kraken_dbName}"
    fi
    """
    stub:
    """
    sleep 3 # Simulate time taken to copy DB; adjust as needed
    echo "[stub] Skipping DB copy to /dev/shm/${params.kraken_dbName}"
    """
}

process clean_shm {
    tag "clean_k2_db"
    label "min_mem"

    // Receives a collected list of all K2 outputs so this only runs after
    // every K2 job (across both sub-workflows) has completed.
    input:
    val(all_k2_done)

    script:
    """
    rm -rf /dev/shm/${params.kraken_dbName}
    echo "Shared memory cleaned"
    """
    stub:
    """
    sleep 3 # Simulate time taken to clean up; adjust as needed
    echo "[stub] Skipping /dev/shm/${params.kraken_dbName} removal"
    """
}