// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GUNZIP {
    tag '$archive'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) },
        enabled: options.publish

    conda (params.conda ? "bioconda::samtools=1.12" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.12--h9aed4be_1"
    } else {
        container "quay.io/biocontainers/samtools:1.12--h9aed4be_1"
    }

    input:
    path archive

    output:
    path "$gunzip",       emit: gunzip
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    gunzip       = archive.toString() - '.gz'
    """
    gunzip -f $options.args $archive

    echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//' > ${software}.version.txt
    """
}
