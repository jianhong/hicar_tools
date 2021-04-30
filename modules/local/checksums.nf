// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CHECKSUMS {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) },
        enabled: options.publish

    conda (params.conda ? "bioconda::sed=4.7.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/sed:4.7.0"
    } else {
        container "quay.io/biocontainers/sed:4.7.0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    path "md5.*.txt", emit: md5

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    touch md5.${prefix}.txt
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
    gunzip -c ${prefix}_1.fastq.gz > ${prefix}_1.fastq
    ${params.md5sum} ${prefix}_1.fastq >>md5.${prefix}.txt
    gunzip -c ${prefix}_2.fastq.gz > ${prefix}_2.fastq
    ${params.md5sum} ${prefix}_2.fastq >>md5.${prefix}.txt
    if [ "${meta.md5_1}" != "null" ] && [ "${meta.md5_1}" != "" ]; then
        md5=(\$(${params.md5sum} ${prefix}_1.fastq.gz))
        if [ "\$md5" != "${meta.md5_1}" ]
        then
            echo "${meta.id} has checksum ${meta.md5_1}, but we got checksum \$md5!"
            exit 128
        fi
    fi
    if [ "${meta.md5_2}" != "null" ] && [ "${meta.md5_2}" != "" ]; then
        md5=(\$(${params.md5sum} ${prefix}_2.fastq.gz))
        if [ "\$md5" != "${meta.md5_2}" ]
        then
            echo "${meta.id} has checksum ${meta.md5_2}, but we got checksum \$md5!"
            exit 128
        fi
    fi
    """
}
