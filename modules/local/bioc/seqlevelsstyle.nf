// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process SEQLEVELS_STYLE {
    tag "$bed"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') },
        enabled: options.publish

    conda (params.enable_conda ? "bioconda::bioconductor-genomeinfodb=1.26.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-genomeinfodb:1.26.4--r40hdfd78af_0"
    } else {
        container "quay.io/biocontainers/bioconductor-genomeinfodb:1.26.4--r40hdfd78af_0"
    }

    input:
    path bed

    output:
    stdout emit: seqlevels_style
    path "*.version.txt"          , emit: version

    script:
    def software = "GenomeInfoDb"

    """
    seqlevelsstyle.r $bed > r.log.txt 2>&1
    cat tmp.txt
    """
}
