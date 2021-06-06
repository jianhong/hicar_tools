#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/hicar
========================================================================================
 nf-core/hicar Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/hicar
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
def json_schema = "$baseDir/nextflow_schema.json"
def unexpectedParams = []
if (params.validate_params) {
    unexpectedParams = NfcoreSchema.validateParameters(params, json_schema, log)
}
// Show help message
if (params.help) {
    def command = "nextflow run jianhong/hicar --input 'input.csv' -profile docker"
    log.info Headers.nf_core(workflow, params.monochrome_logs)
    log.info Schema.params_help(json_schema, command)
    exit 0
}
////////////////////////////////////////////////////

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Configurable variables
params.fasta      = Checks.get_genome_attribute(params, 'fasta')
params.bwa_index  = Checks.get_genome_attribute(params, 'bwa')
params.gtf        = Checks.get_genome_attribute(params, 'gtf')
params.gene_bed   = Checks.get_genome_attribute(params, 'bed12')
params.macs_gsize = Checks.get_genome_attribute(params, 'macs_gsize')
params.deep_gsize = Checks.get_genome_attribute(params, 'deep_gsize')
params.blacklist  = Checks.get_genome_attribute(params, 'blacklist')
anno_readme       = Checks.get_genome_attribute(params, 'readme')
params.species    = Checks.get_genome_attribute(params, 'species')
params.autosomal  = Checks.get_genome_attribute(params, 'autosomal')
params.mappability= Checks.get_genome_attribute(params, 'mappability')

def RE_cutsite = [
    "mboi": "^GATC",
    "dpnii": "^GATC",
    "bglii": "^GATCT",
    "hindiii": "^AGCTT",
    "cviqi": "^TAC"]
if (!params.enzyme.toLowerCase() in RE_cutsite){
    exit 1, "Not supported yet!"
}
params.restriction_sites = RE_cutsite[params.enzyme.toLowerCase()]

if (!params.autosomal){
    exit 1, "Not supported yet!"
}
if (!params.macs_gsize){
    exit 1, "Not supported yet!"
}

if(params.enable_conda || params.conda){
    params.enable_conda = true
    params.conda = true
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

// Header log info
def summary = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

// Check the hostnames against configured profiles
checkHostname()

/*
 * Create a channel for input read files
 */
ch_input = Channel.fromPath("${params.input}").splitCsv(header: true, sep:",")

// index.Rmd; TODO
//ch_index_docs = file("$projectDir/docs/index.Rmd", checkIfExists: true)

// software collection
ch_software_versions = Channel.empty()

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
def getParam(modules, module) {
    return modules[module]?:[:]
}
def getSubWorkFlowParam(modules, mods) {
    def Map options = [:]
    mods.each{
      val ->
        options[val] = modules[val]?:[:]
    }
    return options
}

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { CHECKSUMS              } from './modules/local/checksums'
include { PREPARE_GENOME         } from './modules/local/subworkflow/preparegenome' addParams(options: getSubWorkFlowParam(modules, ['gunzip', 'gtf2bed', 'chromsizes', 'genomefilter', 'bwa_index', 'gffread', 'digest_genome']))
include { BAM_STAT               } from './modules/local/subworkflow/bam_stats' addParams(options: getSubWorkFlowParam(modules, ['samtools_sort', 'samtools_index', 'samtools_stats', 'samtools_flagstat', 'samtools_idxstats']))
include { PAIRTOOLS_PAIRE        } from './modules/local/subworkflow/pairtools' addParams(options: getSubWorkFlowParam(modules, ['paritools_dedup', 'pairtools_flip', 'pairtools_parse', 'pairtools_restrict', 'pairtools_select', 'pairtools_select_long', 'pairtools_sort', 'pairix', 'reads_stat', 'reads_summary', 'pairsqc', 'pairsplot']))
include { COOLER                 } from './modules/local/subworkflow/cooler' addParams(options: getSubWorkFlowParam(modules, ['cooler_cload', 'cooler_merge', 'cooler_zoomify', 'cooler_dump_per_group', 'cooler_dump_per_sample']))
include { ATAC_PEAK              } from './modules/local/subworkflow/callatacpeak' addParams(options: getSubWorkFlowParam(modules, ['pairtools_select', 'pairtools_select_short', 'merge_reads', 'shift_reads', 'macs2_atac', 'dump_reads_per_group', 'dump_reads_per_sample', 'merge_peak']))
include { MAPS_MULTIENZYME       } from './modules/local/subworkflow/multienzyme'   addParams(options: getSubWorkFlowParam(modules, ['maps_cut', 'maps_fend', 'genmap_mappability', 'ucsc_wigtobigwig', 'maps_mapability', 'maps_merge', 'maps_feature']))
include { MAPS_PEAK              } from './modules/local/subworkflow/maps_peak' addParams(options: getSubWorkFlowParam(modules, ['maps_maps', 'maps_callpeak', 'maps_reformat']))
include { DIFFHICAR              } from './modules/local/bioc/diffhicar' addParams(options: getParam(modules, 'diffhicar'))
include { BIOC_CHIPPEAKANNO      } from './modules/local/bioc/chippeakanno' addParams(options: getParam(modules, 'chippeakanno'))
include { BIOC_ENRICH            } from './modules/local/bioc/enrich' addParams(options: getParam(modules, 'enrichment'))
include { GET_SOFTWARE_VERSIONS  } from './modules/local/get_software_version'

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

include { FASTQC   } from './modules/nf-core/software/fastqc/main'
include { CUTADAPT } from './modules/nf-core/software/cutadapt/main' addParams(options: getParam(modules, 'cutadapt'))
include { BWA_MEM  } from './modules/nf-core/software/bwa/mem/main'  addParams(options: getParam(modules, 'bwa_mem'))

def rtitle = ''
def rfilename = ''
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    rtitle = "--title \"${workflow.runName}\""
    rfilename = "--filename " + workflow.runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report"
}
include { MULTIQC  } from './modules/nf-core/software/multiqc/main'  addParams(options: [args:"$rtitle $rfilename"])

ch_reads = ch_input.map{
  row ->
      fastq1 = file(row.remove("fastq_1"), checkIfExists: true)
      fastq2 = file(row.remove("fastq_2"), checkIfExists: true)
      meta = row
      meta.id = row.group + "_REP" + row.replicate
      [meta, [fastq1, fastq2]]
}
//ch_reads.view()
cool_bin = Channel.fromList(params.cool_bin.tokenize('_'))

/*
 * HiCAR workflow
 */
workflow{
  // check the input fastq files are correct
  CHECKSUMS(ch_reads)

  PREPARE_GENOME()
  ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.version.ifEmpty(null))

  // Read QC
  if(!params.skip_fastqc){
    FASTQC(ch_reads)
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.ifEmpty(null))
  }

  // trimming
  CUTADAPT(ch_reads)
  ch_software_versions = ch_software_versions.mix(CUTADAPT.out.version.ifEmpty(null))

  // mapping
  BWA_MEM(
    CUTADAPT.out.reads,
    PREPARE_GENOME.out.bwa_index
  )
  ch_software_versions = ch_software_versions.mix(BWA_MEM.out.version.ifEmpty(null))
  // mapping stats
  BAM_STAT(BWA_MEM.out.bam)
  ch_software_versions = ch_software_versions.mix(BAM_STAT.out.version.ifEmpty(null))

  // filter reads, output pair (like hic pair), raw (pair), and stats
  PAIRTOOLS_PAIRE(
    BWA_MEM.out.bam,
    PREPARE_GENOME.out.chrom_sizes,
    PREPARE_GENOME.out.digest_genome
  )
  ch_software_versions = ch_software_versions.mix(PAIRTOOLS_PAIRE.out.version.ifEmpty(null))

  // combine bin_size and create cooler file, and dump long_bedpe
  cool_bin.combine(PAIRTOOLS_PAIRE.out.pair)
          .map{bin, meta, pair, px -> [meta, bin, pair, px]}
          .set{cool_input}
  COOLER(
    cool_input,
    PREPARE_GENOME.out.chrom_sizes
  )
  ch_software_versions = ch_software_versions.mix(COOLER.out.version.ifEmpty(null))

  // calling ATAC peaks, output ATAC narrowPeak and reads in peak
  ATAC_PEAK(
    PAIRTOOLS_PAIRE.out.raw
  )
  ch_software_versions = ch_software_versions.mix(ATAC_PEAK.out.version.ifEmpty(null))

  // calling distal peaks: [ meta, bin_size, path(macs2), path(long_bedpe), path(short_bed), path(background) ]
  background = MAPS_MULTIENZYME(PREPARE_GENOME.out.fasta, cool_bin, PREPARE_GENOME.out.chrom_sizes).bin_feature
  ch_software_versions = ch_software_versions.mix(MAPS_MULTIENZYME.out.version.ifEmpty(null))
  reads_peak   = ATAC_PEAK.out.reads
                        .map{ meta, reads ->
                                [meta.id, reads]} // here id is group
                        .combine(ATAC_PEAK.out.mergedpeak)// group, reads, peaks
                        .cross(COOLER.out.bedpe.map{[it[0].id, it[0].bin, it[1]]})// group, bin, bedpe
                        .map{ short_bed, long_bedpe -> //[bin_size, group, macs2, long_bedpe, short_bed]
                                [long_bedpe[1], short_bed[0], short_bed[2], long_bedpe[2], short_bed[1]]}
  background.cross(reads_peak)
              .map{ background, reads -> //[group, bin_size, macs2, long_bedpe, short_bed, background]
                    [[id:reads[1]], background[0], reads[2], reads[3], reads[4], background[1]]}
              .set{ maps_input }
  //maps_input.view()
  MAPS_PEAK(maps_input)
  ch_software_versions = ch_software_versions.mix(MAPS_PEAK.out.version.ifEmpty(null))

  // Differential analysis
  if(!params.skip_diff_analysis){
    MAPS_PEAK.out.peak //[]
             .map{meta, bin_size, peak -> [bin_size, peak]}
             .groupTuple()
             .cross(COOLER.out.samplebedpe.map{[it[0].bin, it[1]]}.groupTuple())
             .map{ peak, long_bedpe -> [peak[0], peak[1].flatten(), long_bedpe[1].flatten()] }//bin_size, meta, peak, long_bedpe
             .groupTuple()
             .map{[it[0], it[1].flatten().unique(), it[2].flatten()]}
             .set{ch_diffhicar}
    //ch_diffhicar.view()
    DIFFHICAR(ch_diffhicar)
    ch_software_versions = ch_software_versions.mix(DIFFHICAR.out.version.ifEmpty(null))
    //annotation
    if(!params.skip_peak_annotation){
        BIOC_CHIPPEAKANNO(DIFFHICAR.out.diff, PREPARE_GENOME.out.gtf)
        ch_software_versions = ch_software_versions.mix(BIOC_CHIPPEAKANNO.out.version.ifEmpty(null))
        BIOC_ENRICH(BIOC_CHIPPEAKANNO.out.anno)
        ch_software_versions = ch_software_versions.mix(BIOC_ENRICH.out.version.ifEmpty(null))
    }
  }

  //annotation


  // Parse software version numbers
  ch_software_versions.flatten()
                      .map{[it.getSimpleName(), it]}
                      .groupTuple()
                      .map{it[1][0]}
                      .collect()
                      .set{ch_version}
  GET_SOFTWARE_VERSIONS(ch_version)

  if(!params.skip_multiqc){
      MULTIQC(
        GET_SOFTWARE_VERSIONS.out.ch_software_versions_yaml
                             .concat(ch_multiqc_config)
                             .concat(ch_multiqc_custom_config.collect().ifEmpty([]))
                             .concat(FASTQC.out.zip.map{it[1]}.collect().ifEmpty([]))
                             .concat(CUTADAPT.out.log.map{it[1]}.collect().ifEmpty([]))
                             .concat(BAM_STAT.out.stats.map{it[1]}.collect().ifEmpty([]))
                             .concat(BAM_STAT.out.flagstat.map{it[1]}.collect().ifEmpty([]))
                             .concat(BAM_STAT.out.idxstats.map{it[1]}.collect().ifEmpty([]))
                             .concat(ATAC_PEAK.out.xls.map{it[1]}.collect().ifEmpty([]))
                             .concat(PAIRTOOLS_PAIRE.out.stat.collect().ifEmpty([]))
                             .collect()
      )
  }
}

/*
 * Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {
    def multiqc_report = []
    Completion.email(workflow, params, summary, run_name, projectDir, multiqc_report, log)
    Completion.summary(workflow, params, log)
}

workflow.onError {
    // Print unexpected parameters
    for (p in unexpectedParams) {
        log.warn "Unexpected parameter: ${p}"
    }
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
