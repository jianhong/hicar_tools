/*
 * pair the proper mapped pairs
 */
include { initOptions } from './functions'
params.options = [:]
options        = initOptions(params.options)

include { PAIRTOOLS_DEDUP    } from '../pairtools/dedup'    addParams(options: options.paritools_dedup)
include { PAIRTOOLS_FLIP     } from '../pairtools/flip'     addParams(options: options.pairtools_flip)
include { PAIRTOOLS_PARSE    } from '../pairtools/parse'    addParams(options: options.pairtools_parse)
include { PAIRTOOLS_RESTRICT } from '../pairtools/restrict' addParams(options: options.pairtools_restrict)
include { PAIRTOOLS_SELECT   } from '../pairtools/select'   addParams(options: options.pairtools_select)
include { PAIRTOOLS_SELECT
       as PAIRTOOLS_SELECT2  } from '../pairtools/select'   addParams(options: options.pairtools_select_long)
include { PAIRTOOLS_SORT     } from '../pairtools/sort'     addParams(options: options.pairtools_sort)
include { PAIRIX             } from '../pairix/pairix'      addParams(options: options.pairix)
include { READS_STAT         } from '../reads_stat'         addParams(options: options.reads_stat)
include { READS_SUMMARY      } from '../reads_summary'      addParams(options: options.reads_summary)
include { PAIRSQC            } from '../pairix/pairsqc'     addParams(options: options.pairsqc)
include { PAIRSPLOT          } from '../pairix/pairsplot'   addParams(options: options.pairsplot)

workflow PAIRTOOLS_PAIRE {
  take:
  ch_bam      // channel: [ val(meta), [bam] ]
  chromsizes  // channel: [ path(chromsizes) ]
  frag        // channel: [ path(fragment) ]

  main:
  //raw pairs, output raw.pairsam
  PAIRTOOLS_PARSE(ch_bam, chromsizes)
  // select valid pairs, output sorted.pairs
  PAIRTOOLS_FLIP(PAIRTOOLS_PARSE.out.pairsam, chromsizes)
  PAIRTOOLS_SELECT(PAIRTOOLS_FLIP.out.flip)
  PAIRTOOLS_SORT(PAIRTOOLS_SELECT.out.sel)
  // remove duplicate pairs, output dedup.pairs
  PAIRTOOLS_DEDUP(PAIRTOOLS_SORT.out.sorted)
  // remove same fragment pairs, output samefrag.pairs, valid.pairs <- like HiC pairs
  PAIRTOOLS_RESTRICT(PAIRTOOLS_DEDUP.out.pairs, frag)
  PAIRTOOLS_SELECT2(PAIRTOOLS_RESTRICT.out.restrict)
  // make index for valid.pairs
  PAIRIX(PAIRTOOLS_SELECT2.out.rest)
  //reads information
  PAIRTOOLS_PARSE.out.stat
                 .map{meta, stat -> [meta.id, meta, stat]}
                 .join(PAIRTOOLS_DEDUP.out.stat.map{meta, stat -> [meta.id, stat]})
                 .map{id, meta, raw, dedup -> [meta, raw, dedup ]}
                 .set{ reads_stat }
  READS_STAT(reads_stat)
  READS_SUMMARY(READS_STAT.out.stat.map{it[1]}.collect())
  PAIRSQC(PAIRIX.out.index, chromsizes)
  PAIRSPLOT(PAIRSQC.out.qc)

  emit:
  pair = PAIRIX.out.index               // channel: [ val(meta), [valid.pair.gz], [valid.pair.gz.px] ]
  stat = READS_SUMMARY.out.summary      // channel: [ val(meta), [summary] ]
  raw = PAIRTOOLS_PARSE.out.pairsam     // channel: [ val(meta), [pairsam] ]
  version = PAIRTOOLS_PARSE.out.version // channel: [ path(version) ]
}
