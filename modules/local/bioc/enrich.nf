// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process BIOC_ENRICH {
    tag "$meta.id"
    label 'process_high'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::bioconductor-clusterprofiler=3.18.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-clusterprofiler:3.18.1--r40hdfd78af_0"
    } else {
        container "quay.io/biocontainers/bioconductor-clusterprofiler:3.18.1--r40hdfd78af_0"
    }

    input:
    tuple val(meta), val(bin_size), path(diff)

    output:
    tuple val(meta), val(bin_size), path("enrichment/*"), emit: enrichment
    path "*.version.txt"                                , emit: version

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    install_packages.r clusterProfiler pathview biomaRt optparse
    enrich.r -s ${params.species} $options.args

    # *.version.txt
    """
}
