/*
 * Createing Genomic Features Files
 */
include { initOptions } from './functions'
params.options = [:]
options        = initOptions(params.options)

include { MAPS_MAPS             } from '../maps/maps'         addParams(options: options.maps_maps)
include { MAPS_CALLPEAK         } from '../maps/callpeak'     addParams(options: options.maps_callpeak)
include { MAPS_REFORMAT         } from '../maps/reformat'     addParams(options: options.maps_reformat)

workflow MAPS_PEAK {
  take:
  reads        // channel: [ meta, bin_size, path(macs2), path(long_bedpe), path(short_bed), path(background) ]

  main:
  //create parameter table
  //input=val(meta), val(bin_size), path(macs2), path(long_bedpe), path(short_bed), path(background)
  //maps from bedpe
  ch_version = MAPS_MAPS(reads).version
  //regression and peak calling
  peak = MAPS_CALLPEAK(MAPS_MAPS.out.maps).peak
  ch_version = ch_version.mix(MAPS_CALLPEAK.out.version)
  //peak formatting
  MAPS_REFORMAT(peak)
  ch_version = ch_version.mix(MAPS_REFORMAT.out.version)

  emit:
  peak         = MAPS_REFORMAT.out.bedpe      // channel: [ path(bedpe) ]
  version      = ch_version                   // channel: [ path(version) ]
}
