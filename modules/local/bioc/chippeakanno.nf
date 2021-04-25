// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process BIOC_CHIPPEAKANNO {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::bioconductor-chippeakanno=3.24.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-chippeakanno:3.24.1--r40hdfd78af_0"
    } else {
        container "quay.io/biocontainers/bioconductor-chippeakanno:3.24.1--r40hdfd78af_0"
    }

    input:
    tuple val(meta), val(bin_size), path(diff)

    output:
    tuple val(meta), val(bin_size), path("*.anno.csv"), emit: anno
    path "*.version.txt"                              , emit: version

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    annopeaks.r -s ${params.species}
    """
}
