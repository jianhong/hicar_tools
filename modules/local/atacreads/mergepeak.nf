// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process MERGE_PEAK {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) },
        enabled: options.publish

    conda (params.conda ? "bioconda::bedtools=2.30.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_1"
    } else {
        container "quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_1"
    }

    input:
    path peak

    output:
    path "ATAC_merged_peak.bed"    , emit: peak
    path  "*.version.txt"          , emit: version

    script:
    def software = "bedtools"
    """
    cat *.narrowPeak | cut -f1-3 | sort -k1,1 -k2,2n |
    bedtools merge $options.args \\
        -i stdin > ATAC_merged_peak.bed

    echo \$(bedtools --version) | sed -e "s/bedtools v//g" > ${software}.version.txt
    """
}
