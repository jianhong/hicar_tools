// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process UCSC_BIGWIGAVERAGEOVERBED {
    tag "$bin_size"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) },
        enabled: options.publish

    conda (params.conda ? "bioconda::ucsc-bigwigaverageoverbed=377" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ucsc-bigwigaverageoverbed:377--h0b8a92a_2"
    } else {
        container "quay.io/biocontainers/ucsc-bigwigaverageoverbed:377--h0b8a92a_2"
    }

    input:
    tuple val(bin_size), path(cut)
    path mappability

    output:
    tuple val(bin_size), path("*.map.tab")           , emit: tab
    path "*.version.txt"                             , emit: version

    script:
    def software = "ucsc-bigWigAverageOverBed"
    """
    # there is a bug that bigWigAverageOverBed can not handle ensembl seqlevels style.
    bigWigAverageOverBed $mappability $cut ${cut.getSimpleName()}.map.tab

    echo \$(bigWigAverageOverBed 2>&1) | sed 's/bigWigAverageOverBed v//; s/ - Compute.*\$//' > ${software}.version.txt
    """
}
