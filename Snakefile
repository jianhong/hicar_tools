shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")


import collections

configfile: "config.yaml"

BWA_INDEX      = config['BWA_INDEX']
chromsizes     = config['chromsizes']
genome         = config['genome']
frag_path      = config['frag_path']
genome_version = config['species']
tmp            = config['tmp']
group_text     = config['grouptext']
cool_bin       = config['cool_bin']   ## two bin resolution

smooth_window  = 150
shiftsize      = -75
pval_thresh    = 0.01

FILES      = json.load(open(config['SAMPLES_JSON']))
SAMPLES    = sorted(FILES.keys())
maps_background = '/datacommons/ydiaolab/genome_ref/maps_feature/hg38_cviqi/hg38/CviQI/F_GC_M_CviQI_10Kb_el.txt' ## updatelater

## prepare the key sample
def load_group(grouptxt):
    '''
    group per sample key
    '''
    groupid = collections.defaultdict(list)
    with open(grouptxt,'r') as f:
        for line in f:
            line = line.strip()
            if not line =='':
                group, value = line.split("\t")
                groupid[group].append(value)
    return groupid

groupid = load_group(group_text)


TARGETS = []
cooler = expand("coolers-{genome}/{sample}.{cool_bin}.cool", sample = SAMPLES, genome = genome, cool_bin = cool_bin)
read_summary = expand("reads_sumamry_info-{genome}/{sample}.read_summary", sample = SAMPLES, genome = genome)

# merged_cooler = expand("merged_cooler-{genome}/{sample}.merged.{cool_bin}.cool", \
    # sample = groupid.keys(), genome = genome, cool_bin = cool_bin)
# merged_R2 = expand("merged_R2_reads-{genome}/{sample}.merged.ATAC.bed.gz", \
#     sample = groupid.keys(), genome = genome)

merged_peaks = expand("merged_R2_reads-{genome}/{sample}.merged.ATAC.bed.gz",  sample = groupid.keys(), genome = genome)
peaks  = expand("macs2_peak_{sample}/{sample}_{genome}_peaks.narrowPeak" , sample = SAMPLES, genome = genome) ## no need for each sample

TARGETS.extend(expand("macs2_peak_{sample}/{sample}_{genome}_peaks.narrowPeak", genome = genome, sample = groupid.keys()))
TARGETS.extend(expand("longreads/{sample}-{cool_bin}-{genome}/{sample}.chr17.long.intra.bedpe", genome = genome, sample = groupid.keys(), cool_bin =  cool_bin))
TARGETS.extend(expand("short_reads/{sample}-{genome}/{sample}.chr1.shrt.vip.bed",  genome = genome, sample = groupid.keys()))

# TARGETS.extend(expand("reads_folder/{sample}-{cool_bin}-{genome}/{sample}.chr17.long.intra.bedpe",genome = genome, sample = groupid.keys(), cool_bin =  cool_bin))
TARGETS.extend(expand("maps_result/{sample}-{cool_bin}-{genome}/MAPS_output/{sample}/{sample}.maps.qc", genome = genome, sample = groupid.keys(), cool_bin =  10000))
# mcool = expand("merged_cooler-{genome}/merged.{cool_bin}.mcool", genome = genome, cool_bin = cool_bin)

TARGETS.extend(read_summary)
# TARGETS.extend(merged_cooler)
# TARGETS.extend(merged_R2)
# TARGETS.extend(merged_peaks)



localrules: all, read_info ,create_sym_folder

rule all:
    input: TARGETS
        
rule fastq_preprocessing:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output: 
        "00_CviQI_trimed_fq/{sample}_r1.fq.gz", 
        "00_CviQI_trimed_fq/{sample}_r2.fq.gz"
    threads: 12
    message: "keep read1 start with CviQI cutting site "
    log:
         "00_log/{sample}.cutadapt"
    shell:
        """
        cutadapt -Z -j {threads} -e 0 --no-indels \
        --action none  --discard-untrimmed \
        -g ^TAC \
        -o {output[0]} -p {output[1]}  {input[0]} {input[1]}  2> {log} 
        """

rule bwa_mem_mapping:
    input:
        r1 ="00_CviQI_trimed_fq/{sample}_r1.fq.gz",
        r2 = "00_CviQI_trimed_fq/{sample}_r2.fq.gz"
    output: "01_bam/{sample}.bam"
    threads: 24
    message: "bwa {input}: {threads} threads"
    log:
        "00_log/{sample}.bwa"
    shell:
        """
        module load BWA
        module load samtools
        bwa mem  -SP -t {threads} {BWA_INDEX} {input} | samtools view -bS - > {output}  2> {log}
        """



rule prase_pairs_from_BAM_files: ## no flip to makesure the R1 R2 position for the peak calling
    input:  "01_bam/{sample}.bam"
    output: "pairs-{genome}/{sample}.raw.pairsam.gz", "pairs-{genome}/{sample}.raw.pairsam.stat"
    message: "prase bam {input} "
    threads: 2
    shell:
        """
        module load samtools
        pairtools parse -c {chromsizes}  \
        --assembly {genome} --min-mapq 10 \
        --max-molecule-size 2000 --max-inter-align-gap 20 \
        --walks-policy mask  --no-flip --drop-seq --drop-sam  \
        --output-stats {output[1]} -o {output[0]}  {input} 
        """

rule select_valid_pairs:
    input:  "pairs-{genome}/{sample}.raw.pairsam.gz"
    output: "pairs-{genome}/{sample}.sorted.pairs.gz"
    message: "flip and sort {input} "
    threads: 8
    shell:
        """
         pairtools flip -c {chromsizes} {input} | \
         pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")' | \
         pairtools sort  --nproc 8  --memory 15G  --tmpdir {tmp} -o {output}
        """

rule remove_duplicate_pairs:
    input:  "pairs-{genome}/{sample}.sorted.pairs.gz"
    output: "filtered-{genome}/{sample}.dedup.pairs.gz" ,"filtered-{genome}/{sample}.dedup.pairs.stat"
    message: "dedup to filted {input} "
    threads: 5
    shell:
        """
         pairtools dedup --max-mismatch 1 --method max -o {output[0]} {input} \
         --output-stats  {output[1]}
        """



rule remove_same_fragment_pairs:
    input:  "filtered-{genome}/{sample}.dedup.pairs.gz"
    output: valid = "filtered-{genome}/{sample}.valid.pairs.gz", same_f = "filtered-{genome}/{sample}.samefrag.pairs.gz"
    message: "selected pairsam {input} "
    threads: 5
    shell:
        """
        pairtools restrict -f {frag_path} {input} |  \
        pairtools select '(COLS[-6]==COLS[-3]) and (chrom1==chrom2)' \
        --output-rest {output[valid]} -o {output[same_f]} 
        """

rule make_index:
    input:  "filtered-{genome}/{sample}.valid.pairs.gz"
    output: "filtered-{genome}/{sample}.valid.pairs.gz.px2"
    message: "dedup to filted {input} "
    threads: 5
    shell:
        """
         pairix -p pairs  {input}
        """
        
rule HiC_contact_matrices_Cooler:
    input:  "filtered-{genome}/{sample}.valid.pairs.gz", "filtered-{genome}/{sample}.valid.pairs.gz.px2"
    output: "coolers-{genome}/{sample}.{cool_bin}.cool"
    message: "cooler {input} "
    # params: res = {cool_bin} 
    threads: 10
    shell:
        """
        cooler cload pairix --assembly hg38 --nproc {threads} \
        --max-split 2 {chromsizes}:{wildcards.cool_bin} {input[0]} {output}
        """


rule read_info:  
    input:  "pairs-{genome}/{sample}.raw.pairsam.stat", "filtered-{genome}/{sample}.dedup.pairs.stat"
    output: "reads_sumamry_info-{genome}/{sample}.read_summary"
    threads: 1
    script:
        "script/read_summary.R"



rule extract_R2_ATAC_reads:
    input:  "pairs-{genome}/{sample}.raw.pairsam.gz"
    output: "peaks-{genome}/{sample}.longRange_Trans.pairs.gz", "peaks-{genome}/{sample}.short.pairs.gz"
    message: "flip to filted {input} "
    threads: 8
    shell:
        """
         pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")' {input} | \
         pairtools select '(chrom1==chrom2) and (abs(pos1 - pos2) < 1e4)'  -o {output[1]}  --output-rest {output[0]} 
        """

rule ATAC_reads_Tn5_shifting_duplicate_remove:  
    input:  "peaks-{genome}/{sample}.longRange_Trans.pairs.gz"
    output: "peaks-{genome}/{sample}.R2.ATAC.bed.gz"
    threads: 1
    shell:
        """
        export TMPDIR={tmp}
        zcat {input} | \
        awk ' BEGIN {{OFS="\\t"}} ;  /^[^#]/ {{ {{ if ($7 == "+") {{$5 = $5 + 4}} else if ($7 == "-") {{$5 = $5 - 5}}  print $4, $5, $5+1, "*", "*", $7}} }} ' |\
        sort -k1,1 -k2,2n | uniq  | gzip -nc > {output}
        """


def get_atac_input(wildcards):
    return expand("peaks-{genome}/{sample}.R2.ATAC.bed.gz", genome = wildcards.genome, sample = groupid[wildcards.sample])


rule merged_ATAC_peaks:  
    input:  get_atac_input ## sample group
    output: "merged_R2_reads-{genome}/{sample}.merged.ATAC.bed.gz"
    threads: 1
    shell:
        """
        zcat {input} | \
        sort -k1,1 -k2,2n |  gzip -nc > {output}
        """



# rule merged_ATAC_peaks:  
#     input:  expand("peaks-{genome}/{sample}.R2.ATAC.bed.gz", genome = genome, sample = lambda wildcards: groupid[wildcards.group]) ## sample group
#     output: "merged_R2_reads-{genome}/{group}.merged.ATAC.bed.gz"
#     threads: 1
#     shell:
#         """
#         zcat {input} | \
#         sort -k1,1 -k2,2n |  gzip -nc > {output}
#         """

def get_cooler_input(wildcards):  ## group functin to define the input sample
    return expand("coolers-{genome}/{sample}.{cool_bin}.cool", cool_bin = wildcards.cool_bin , 
        genome = wildcards.genome, sample = groupid[wildcards.sample] )

rule merge_cooler:
    input: get_cooler_input
    output: "merged_cooler-{genome}/{sample}.merged.{cool_bin}.cool"
    threads:1
    shell:
        """
        cooler merge {output} {input}
        """

rule mcooler:
    input:  expand("merged_cooler-{genome}/merged.{cool_bin}.cool", genome = genome, cool_bin = cool_bin)
    output: expand("merged_cooler-{genome}/merged.{cool_bin}.mcool", genome = genome, cool_bin = cool_bin)
    threads:16
    params: res = "5000,10000,50000,100000,250000,500000,1000000,2500000,5000000,10000000"
    shell:
        """
        cooler zoomify -n 16 -o {output} -r {params.res} {input}
        """

rule ATAC_macs2_peaks:
    input:  "merged_R2_reads-{genome}/{sample}.merged.ATAC.bed.gz"
    output: "macs2_peak_{sample}/{sample}_{genome}_peaks.narrowPeak"
    threads:1
    params: name = "{sample}_{genome}"
    shell:
        """
        macs2 callpeak -t {input} -f BED -n {params.name}  -g {genome_version} --qval {pval_thresh} \
        --shift {shiftsize} --extsize {smooth_window} --nomodel -B --SPMR --keep-dup all --call-summits \
        --outdir macs2_peak_{wildcards.sample} 
        """

rule dump_ATAC_reads:
    input: "merged_R2_reads-{genome}/{sample}.merged.ATAC.bed.gz"
    output: "short_reads/{sample}-{genome}/{sample}.chr1.shrt.vip.bed"
    threads: 1
    params: "short_reads/{sample}-{genome}"
    shell:
        """
        mkdir -p {params}
        zcat {input[0]} | awk -v setname={wildcards.sample} -v outdir={params}  -F $"\\t"  'BEGIN {{OFS=FS}};{{print $1,$2,$3 > outdir"/"setname"."$1".shrt.vip.bed" }}' 
        """


rule dump_long_reads:
    input:  "merged_cooler-{genome}/{sample}.merged.{cool_bin}.cool"
    output: "longreads/{sample}-{cool_bin}-{genome}/{sample}.chr17.long.intra.bedpe"
    threads: 1
    params: "longreads/{sample}-{cool_bin}-{genome}"
    shell:
        """
        mkdir  -p {params}
        cooler dump -t pixels -H --join {input} | \
        awk -v setname={wildcards.sample} -v outdir={params}   -F $"\\t" '{{if($1 == $4) {{print > outdir"/"setname"."$1".long.intra.bedpe"}} }}'
        """

rule create_sym_folder:
    input:  "longreads/{sample}-{cool_bin}-{genome}/{sample}.chr17.long.intra.bedpe",
    "short_reads/{sample}-{genome}/{sample}.chr1.shrt.vip.bed"
    output: directory("reads_folder/{sample}-{cool_bin}-{genome}")  ## symbol maybe treated as directory
    threads: 1
    params: "reads_folder/{sample}-{cool_bin}-{genome}"
    shell:
        """
        mkdir  -p {params}
        ln -sr longreads/{wildcards.sample}-{wildcards.cool_bin}-{wildcards.genome}/{wildcards.sample}.*.long.intra.bedpe {params}
        ln -sr short_reads/{wildcards.sample}-{wildcards.genome}/{wildcards.sample}.*.shrt.vip.bed {params}
        """

rule run_maps:
    input: reads_folder = "reads_folder/{sample}-{cool_bin}-{genome}",
     mac2_peak = "macs2_peak_{sample}/{sample}_{genome}_peaks.narrowPeak"
    output: "maps_result/{sample}-{cool_bin}-{genome}/MAPS_output/{sample}/{sample}.maps.qc"
    threads: 1
    params: outdir = "maps_result/{sample}-{cool_bin}-{genome}",  readsdir = "reads_folder/{sample}-{cool_bin}-{genome}"
    shell:
        """
        mkdir -p maps_result/{wildcards.sample}-{wildcards.cool_bin}-{wildcards.genome}/MAPS_output/{wildcards.sample}
        bash maps.sh {wildcards.sample} {input.mac2_peak} {params.outdir} {params.readsdir} {wildcards.genome} {wildcards.cool_bin} {maps_background}
        """


