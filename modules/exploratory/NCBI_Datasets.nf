// mock example of what this module should eventually look like, not implemented yet
// Idea currently parked for future development
process NCBI_DOWNLOAD {
	tag "$accession"

	input:
	val accession

	output:
	path "${accession}.zip"

	script:
	"""
	datasets download genome accession '${accession}' \
		--filename '${accession}.zip'
	"""
}
