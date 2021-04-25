// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process SHIFTREADS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) },
        enabled: options.publish


    input:
    tuple val(meta), path(pair)

    output:
    tuple val(meta), path("*.bed.gz"), emit: bed
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    zcat < $pair | \\
    awk 'BEGIN {{OFS="\t"}} ;  /^[^#]/ {{ {{ if (\$7 == "+") {{\$5 = \$5 + 4}} else if (\$7 == "-") {{\$5 = \$5 - 5}}  print \$4, \$5, \$5+1, "*", "*", \$7}} }}' | \\
    sort -k1,1 -k2,2n | uniq  | gzip -nc > ${prefix}.R2.ATAC.bed.gz

    echo \$(awk --version 2>&1) | sed 's/^.*awk //; s/Using.*\$//' > ${software}.version.txt
    """
}
