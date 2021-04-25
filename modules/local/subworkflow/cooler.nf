/*
 * pair the proper mapped pairs
 */
include { initOptions } from './functions'
params.options = [:]
options        = initOptions(params.options)

include { COOLER_CLOAD   } from '../cooler/cload'   addParams(options: options.cooler_cload)
include { COOLER_MERGE   } from '../cooler/merge'   addParams(options: options.cooler_merge)
include { COOLER_ZOOMIFY } from '../cooler/zoomify' addParams(options: options.cooler_zoomify)
include { COOLER_DUMP    } from '../cooler/dump'    addParams(options: options.cooler_dump_per_group)
include { COOLER_DUMP
    as COOLER_DUMP_SAMPLE} from '../cooler/dump'    addParams(options: options.cooler_dump_per_sample)

workflow COOLER {
  take:
  valid_pairs  // channel: [ val(meta), val(bin), [pairs], [pairs.px] ]
  chromsizes   // channel: [ path(chromsizes) ]

  main:
  // HiC-like contact matrix
  ch_version = COOLER_CLOAD(valid_pairs, chromsizes).version
  // Merge contacts
  COOLER_CLOAD.out.cool
              .map{
                meta, bin, cool ->
                  [meta.group, bin, cool]
              }
              .groupTuple(by:[0, 1])
              .map{group, bin, cool -> [[id:group, bin:bin], cool]}
              .set{ch_cooler}
  COOLER_MERGE(ch_cooler)
  // create mcooler file for visualization
  COOLER_ZOOMIFY(COOLER_MERGE.out.cool)
  // dump long.intra.bedpe for each group for MAPS to call peaks
  COOLER_DUMP(COOLER_MERGE.out.cool)
  // dump long.intra.bedpe for each sample
  COOLER_DUMP_SAMPLE(COOLER_CLOAD.out.cool.map{ meta, bin, cool -> [[id:meta.id, group:meta.group, bin:bin], cool]})

  emit:
  mcool       = COOLER_ZOOMIFY.out.mcool     // channel: [ val(meta), [mcool] ]
  bedpe       = COOLER_DUMP.out.bedpe        // channel: [ val(meta), [bedpe] ]
  samplebedpe = COOLER_DUMP_SAMPLE.out.bedpe // channel: [ val(meta), [bedpe] ]
  version     = ch_version                   // channel: [ path(version) ]
}
