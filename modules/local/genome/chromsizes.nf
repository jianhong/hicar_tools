// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process CHROMSIZES {
    tag "$fasta"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) },
        enabled: options.publish

    conda (params.conda ? "bioconda::samtools=1.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"
    } else {
        container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    }

    input:
    path fasta

    output:
    path '*.sizes'      , emit: sizes
    path '*.fai'        , emit: fai
    path "*.version.txt", emit: version

    script:
    def software = 'samtools'
    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
    """
}
