/*
 * Createing Stats for mapping results
 */
include { initOptions } from './functions'
params.options = [:]
options        = initOptions(params.options)

include { SAMTOOLS_SORT                } from '../../nf-core/software/samtools/sort/main'                     addParams(options: options.samtools_sort)
include { SAMTOOLS_INDEX               } from '../../nf-core/software/samtools/index/main'                    addParams(options: options.samtools_index)
include { SAMTOOLS_STATS               } from '../../nf-core/software/samtools/stats/main'                    addParams(options: options.samtools_stats)
include { SAMTOOLS_IDXSTATS            } from '../../nf-core/software/samtools/idxstats/main'                 addParams(options: options.samtools_idxstats)
include { SAMTOOLS_FLAGSTAT            } from '../../nf-core/software/samtools/flagstat/main'                 addParams(options: options.samtools_flagstat)

workflow BAM_STAT {
    take:
    bam          // channel: [ val(meta), path(bam) ]

    main:
    ch_version = SAMTOOLS_SORT(bam).version
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    ch_bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0])
    SAMTOOLS_STATS(ch_bam_bai)
    SAMTOOLS_FLAGSTAT(ch_bam_bai)
    SAMTOOLS_IDXSTATS(ch_bam_bai)

    emit:
    stats    = SAMTOOLS_STATS.out.stats           // channel: [ val(meta), [ stats ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat     // channel: [ val(meta), [ flagstat ] ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats     // channel: [ val(meta), [ idxstats ] ]
    version  = ch_version                         // channel: [ path(version) ]
}
