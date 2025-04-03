#!/usr/bin/env nextflow

process FASTP_PARSE {
    // tag "${sample_name}"
    // cpus "${params.cpusMin}"
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    publishDir "$params.outdir/fastp/", mode: 'symlink'
    label 'min_mem'


    input: 
    // tuple val(sample_name), path(reads)
    tuple val(reports), path(reports)
    // tuple val(sample), val(base_name), path ("${base_name}_fastp.html"), val (seq_type), emit: report_html


    output:
    path ("report.csv"), emit: report_csv

    script:
    """
    python /home/colinl/metaG/Git/metaG_EukDepletion/parse_json_cl.py -j . -o report.csv
    """
    stub:
    """
    touch report.csv
    """    
}

process FASTP_PLOT { //TODO: fix error   -- Check script './subworkflows/../modules/fastp_parse.nf' at line: 31 or see '.nextflow.log' file for more details
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    publishDir "$params.outdir/out/${seq_type[0]}", mode: 'symlink'
    label 'min_mem'

    input: 
    // tuple val(reports), path(reports)
    tuple val(sample), path (reports), val (seq_type)

    // path (sample_metadata_sheet)
    // path(parse_script)
    
    output:
    path ("report.csv"), emit: report_csv
    path ("*.png"), emit: barplot_png
    path ("*.svg"), emit: barplot_svg

    script:
    """
    find . -name "*_fastp.json" -exec rename 's/_fastp//' {} \\;
    find . -name "*_TRIM*" -exec rename 's/_TRIM/.TRIM/' {} \\;

    find . -name "*.all_scleractina.mapped.json" -exec rename 's/.all_scleractina.mapped./.scleractina./' {} \\;
    find . -name "*.all_aiptasiidae.mapped.json" -exec rename 's/.all_aiptasiidae.mapped./.aiptasiidae./' {} \\;

    find . -name "*.all_scleractina.unmapped.all_symbiodiniaceae.mapped.json" -exec rename 's/.all_scleractina.unmapped.all_symbiodiniaceae.mapped./.symbiodiniaceae./' {} \\;
    find . -name "*.all_aiptasiidae.unmapped.all_symbiodiniaceae.mapped.json" -exec rename 's/.all_aiptasiidae.unmapped.all_symbiodiniaceae.mapped./.symbiodiniaceae./' {} \\;

    find . -name "*.all_symbiodiniaceae.mapped.json" -exec rename 's/.all_symbiodiniaceae.mapped./.symbiodiniaceae./' {} \\;

    find . -name "*.all_scleractina.unmapped.all_symbiodiniaceae.unmapped.bacteria.json" -exec rename 's/.all_scleractina.unmapped.all_symbiodiniaceae.unmapped.bacteria./.bacteria./' {} \\;
    find . -name "*.all_aiptasiidae.unmapped.all_symbiodiniaceae.unmapped.bacteria.json" -exec rename 's/.all_aiptasiidae.unmapped.all_symbiodiniaceae.unmapped.bacteria./.bacteria./' {} \\;
    
    find . -name "*.all_scleractina.unmapped.all_symbiodiniaceae.unmapped.non-bacteria.json" -exec rename 's/.all_scleractina.unmapped.all_symbiodiniaceae.unmapped.non-bacteria./.other./' {} \\;
    find . -name "*.all_aiptasiidae.unmapped.all_symbiodiniaceae.unmapped.non-bacteria.json" -exec rename 's/.all_aiptasiidae.unmapped.all_symbiodiniaceae.unmapped.non-bacteria./.other./' {} \\;
    
    rm *.all_scleractina.unmapped.all_symbiodiniaceae.unmapped.json *.all_scleractina.unmapped.json *.all_aiptasiidae.unmapped.json *.all_symbiodiniaceae.unmapped.json *.non-bacteria.json
    
    python ${params.parse_script} -j . -o report.csv
    
    Rscript ${params.barplot_script} -f report.csv -m ${params.sample_metadata_sheet} -s ${seq_type[0]} -r reads
    
    """
    stub:
    """
    touch report.csv
    touch stub_photos_empty.png
    touch stub_photos_empty.svg
    echo "$seq_type" > checkme_please.txt
    """

}

// Add unrecognized words to the custom dictionary
    // cspell:ignore fastp scleractina aiptasiidae symbiodiniaceae
    // ls -1 *_fastp.json | while read line ;do mv \$line \$(echo \$line | sed s'/_fastp//'g) ; done
    // ls -1 *_TRIM* | while read line ;do mv \$line \$(echo \$line | sed s'/_TRIM/.TRIM/'g) ; done
    // ls -1 *.all_scleractina.mapped.json | while read line ;do mv \$line \$(echo \$line | sed s'/.all_scleractina.mapped./.scleractina./'g) ; done
    // ls -1 *.all_aiptasiidae.mapped.json | while read line ;do mv \$line \$(echo \$line | sed s'/.all_aiptasiidae.mapped./.aiptasiidae./'g) ; done
    // ls -1 *.all_scleractina.unmapped.all_symbiodiniaceae.mapped.json | while read line ;do mv \$line \$(echo \$line | sed s'/.all_scleractina.unmapped.all_symbiodiniaceae.mapped./.symbiodiniaceae./'g) ; done
    // ls -1 *.all_aiptasiidae.unmapped.all_symbiodiniaceae.mapped.json | while read line ;do mv \$line \$(echo \$line | sed s'/.all_aiptasiidae.unmapped.all_symbiodiniaceae.mapped./.symbiodiniaceae./'g) ; done
    // ls -1 *.all_symbiodiniaceae.mapped.json | while read line ;do mv \$line \$(echo \$line | sed s'/.all_symbiodiniaceae.mapped./.symbiodiniaceae./'g) ; done
    // ls -1 *.all_scleractina.unmapped.all_symbiodiniaceae.unmapped.bacteria.json | while read line ;do mv \$line \$(echo \$line | sed s'/.all_scleractina.unmapped.all_symbiodiniaceae.unmapped.bacteria./.bacteria./'g) ; done
    // ls *.all_scleractina.unmapped.all_symbiodiniaceae.unmapped.non-bacteria.json| while read line ; do mv \$line \$(echo \$line | sed s'/.all_scleractina.unmapped.all_symbiodiniaceae.unmapped.non-bacteria./.other./'g); done
    // rm *.all_scleractina.unmapped.all_symbiodiniaceae.unmapped.json *.all_scleractina.unmapped.json *.all_symbiodiniaceae.unmapped.json *.non-bacteria.json
    
// mamba env create --file /home/colinl/metaG/Git/metaG_EukDepletion/environment.yaml

process FASTP_PLOT_CONTIGS { //TODO: fix error   -- Check script './subworkflows/../modules/fastp_parse.nf' at line: 31 or see '.nextflow.log' file for more details
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    publishDir "$params.outdir/out/${seq_type[0]}", mode: 'symlink'
    label 'min_mem'

    input: 
    // tuple val(reports), path(reports)
    tuple val(sample), path (reports), val (seq_type)

    // path (sample_metadata_sheet)
    // path(parse_script)
    
    output:
    path ("report.csv"), emit: report_csv
    path ("*.png"), emit: barplot_png
    path ("*.svg"), emit: barplot_svg

    script:
    """
    find . -name "*_fastp.json" -exec rename 's/_fastp//' {} \\;
    find . -name "*_TRIM*" -exec rename 's/_TRIM/.TRIM/' {} \\;
    find . -name "*.TRIM.contigs.json" -exec rename 's/.TRIM.contigs./.TRIM./' {} \\;

    find . -name "*.contigs.corals.mapped.json" -exec rename 's/.contigs.corals.mapped./.scleractina./' {} \\;
    find . -name "*.contigs.aiptasiidae.mapped.json" -exec rename 's/.contigs.aiptasiidae.mapped./.aiptasiidae./' {} \\;

    find . -name "*.contigs.corals.unmapped.all_symbiodiniaceae_v2.mapped.json" -exec rename 's/.contigs.corals.unmapped.all_symbiodiniaceae_v2.mapped./.symbiodiniaceae./' {} \\;
    find . -name "*.contigs.aiptasiidae.unmapped.all_symbiodiniaceae_v2.mapped.json" -exec rename 's/.contigs.aiptasiidae.unmapped.all_symbiodiniaceae_v2.mapped./.symbiodiniaceae./' {} \\;

    find . -name "*.contigs.corals.unmapped.all_symbiodiniaceae_v2.unmapped.bacteria.json" -exec rename 's/.contigs.corals.unmapped.all_symbiodiniaceae_v2.unmapped.bacteria./.bacteria./' {} \\;
    find . -name "*.contigs.aiptasiidae.unmapped.all_symbiodiniaceae_v2.unmapped.bacteria.json" -exec rename 's/.contigs.aiptasiidae.unmapped.all_symbiodiniaceae_v2.unmapped.bacteria./.bacteria./' {} \\;
    
    find . -name "*.contigs.corals.unmapped.all_symbiodiniaceae_v2.unmapped.non-bacteria.json" -exec rename 's/.contigs.corals.unmapped.all_symbiodiniaceae_v2.unmapped.non-bacteria./.other./' {} \\;
    find . -name "*.contigs.aiptasiidae.unmapped.all_symbiodiniaceae_v2.unmapped.non-bacteria.json" -exec rename 's/.contigs.aiptasiidae.unmapped.all_symbiodiniaceae_v2.unmapped.non-bacteria./.other./' {} \\;
    
    rm *.TRIM.contigs.corals.unmapped.all_symbiodiniaceae_v2.unmapped.json *.TRIM.contigs.corals.unmapped.json 
    rm *.TRIM.contigs.aiptasiidae.unmapped.all_symbiodiniaceae_v2.unmapped.json *.TRIM.contigs.aiptasiidae.unmapped.json 
    
    python ${params.parse_script} -j . -o report.csv
    
    Rscript ${params.barplot_script} -f report.csv -m ${params.sample_metadata_sheet} -s ${seq_type[0]} -r contigs
    
    """
    stub:
    """
    touch report.csv
    touch stub_photos_empty.png
    touch stub_photos_empty.svg
    echo "$seq_type" > checkme_please.txt
    """

}