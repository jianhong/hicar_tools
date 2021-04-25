// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process DIFFHICAR {
    tag "$bin_size"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:bin_size) }

    conda (params.enable_conda ? "bioconda::bioconductor-deseq2=1.30.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-deseq2:1.30.1--r40h399db7b_0"
    } else {
        container "quay.io/biocontainers/bioconductor-deseq2:1.30.1--r40h399db7b_0"
    }

    input:
    tuple val(meta), val(bin_size), path(bedpe)

    output:
    tuple val(meta), val(bin_size), path("diffhic/*"), emit: diff
    path "*.version.txt"                             , emit: version

    script:
    """
    echo '${metadata}' > designtab.txt
    diffhicar.r -d 'designtab.txt' \\
        $options.args
    ## must output the packaes version as *.version.txt
    """
}
