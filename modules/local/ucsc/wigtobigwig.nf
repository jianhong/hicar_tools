// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process UCSC_WIGTOBIGWIG {
    tag '$wig'
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') },
        enabled: options.publish

    conda (params.conda ? "bioconda::ucsc-wigtobigwig=377" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ucsc-wigtobigwig:377--h0b8a92a_2"
    } else {
        container "quay.io/biocontainers/ucsc-wigtobigwig:377--h0b8a92a_2"
    }

    input:
    path wig
    path chromsizes

    output:
    path "*.bw"                   , emit: bw
    path "*.version.txt"          , emit: version

    script:
    def software = 'ucsc-wigtobigwig'

    """
    wigToBigWig \\
        $options.args \\
        $wig \\
        $chromsizes \\
        ${wig.getSimpleName()}.bw

    echo \$(wigToBigWig 2>&1) | sed 's/wigToBigWig v //; s/ - Convert.*\$//' > ${software}.version.txt
    """
}
