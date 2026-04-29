process BRACKEN_PE {
	tag "$meta.id"
	label 'min_mem'
	publishDir "${params.outdir}/bracken/${meta.id}/", mode: 'symlink'
	conda "bioconda::bracken"
	maxRetries params.bracken_threshold
	errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }

	input:
	tuple val(meta), path(kraken_report)

	output:
	tuple val(meta), path("${meta.id}.bracken.tsv"), emit: bracken
	tuple val(meta), path("${meta.id}.bracken.report"), emit: report

	script:
	def threshold = params.bracken_threshold - (task.attempt - 1)
	"""
	bracken \\
		-d /dev/shm/${params.kraken_dbName} \\
		-i ${kraken_report} \\
		-o ${meta.id}.bracken.tsv \\
		-w ${meta.id}.bracken.report \\
        -t ${threshold} \\
		-r ${params.bracken_read_len_pe} \\
		-l ${params.bracken_tax_level}
	"""
	stub:
	"""
	touch ${meta.id}.bracken.tsv
	touch ${meta.id}.bracken.report
	"""
}


process BRACKEN_ONT {
	tag "$meta.id"
	label 'min_mem'
	publishDir "${params.outdir}/bracken/${meta.id}/", mode: 'symlink'
	conda "bioconda::bracken"
	// Retry up to 3 times (attempt 1 = threshold 5, 2 = threshold 3, 3 = threshold 0)
    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.bracken_threshold

	input:
	tuple val(meta), path(kraken_report)

	output:
	tuple val(meta), path("${meta.id}.bracken.tsv"), emit: bracken
	tuple val(meta), path("${meta.id}.bracken.report"), emit: report

	script:
	def threshold = params.bracken_threshold - (task.attempt - 1)
	"""
	bracken \\
		-d /dev/shm/${params.kraken_dbName} \\
		-i ${kraken_report} \\
		-o ${meta.id}.bracken.tsv \\
		-w ${meta.id}.bracken.report \\
        -t ${threshold} \\
		-r ${params.bracken_read_len_ont} \\
		-l ${params.bracken_tax_level}
	"""
	stub:
	"""
	touch ${meta.id}.bracken.tsv
	touch ${meta.id}.bracken.report
	"""
}
