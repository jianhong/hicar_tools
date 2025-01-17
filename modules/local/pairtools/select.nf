// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process PAIRTOOLS_SELECT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) },
        enabled: options.publish

    conda (params.conda ? "bioconda::pairtools=0.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pairtools:0.3.0--py37hb9c2fc3_5"
    } else {
        container "quay.io/biocontainers/pairtools:0.3.0--py37hb9c2fc3_5"
    }

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.${options.args2}.pairs.gz"), emit: sel
    tuple val(meta), path("*.${options.args3}.pairs.gz"), emit: rest
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    pairtools select \\
        "$options.args" \\
        -o ${prefix}.${options.args2}.pairs.gz \\
        --output-rest ${prefix}.${options.args3}.pairs.gz \\
        ${input}

    echo \$(pairtools --version 2>&1) | sed 's/pairtools.*version //' > ${software}.version.txt
    """
}
