// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process PAIRTOOLS_PARSE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) },
        enabled: options.publish

    conda (params.conda ? "bioconda::pairtools=0.3.0 bioconda::samtools=1.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pairtools:0.3.0--py37hb9c2fc3_5"
    } else {
        container "quay.io/biocontainers/pairtools:0.3.0--py37hb9c2fc3_5"
    }

    input:
    tuple val(meta), path(bam)
    path(chromsizes)

    output:
    tuple val(meta), path("*.pairsam.gz"), emit: pairsam
    tuple val(meta), path("*.pairsam.stat"), emit: stat
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    pairtools parse -c ${chromsizes} \\
        --assembly ${params.species} \\
        $options.args \\
        --output-stats ${prefix}.raw.pairsam.stat \\
        -o ${prefix}.raw.pairsam.gz \\
        ${prefix}.bam

    echo \$(pairtools --version 2>&1) | sed 's/pairtools.*version //' > ${software}.version.txt
    """
}
