#!/usr/bin/env Rscript

#######################################################################
#######################################################################
## Created on Nov. 10, 2020 call DiffBind 3.0
## Copyright (c) 2020 Jianhong Ou (jianhong.ou@gmail.com)
#######################################################################
#######################################################################
c("optparse", "ChIPpeakAnno", "clusterProfiler", "pathview", "biomaRt")
# set libPath to pwd
pwd <- getwd()
pwd <- file.path(pwd, "lib")
dir.create(pwd)
.libPaths(c(pwd, .libPaths()))

library(ChIPpeakAnno)
library(clusterProfiler)
library(pathview)
library(biomaRt)
library(optparse)

option_list <- list(make_option(c("-s", "--species"), type="character", default=NULL, help="species", metavar="string"),
                    make_option(c("-n", "--ucscname"), type="character", default=NULL, help="ucscname", metavar="string"),
                    make_option(c("-c", "--cores"), type="integer", default=1, help="Number of cores", metavar="integer"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$species)){
  print_help(opt_parser)
  stop("Please provide species.", call.=FALSE)
}

attr <- c("hsapiens_homolog_associated_gene_name")
scientificName <- c("GRCh37"="Homo sapiens",
                    "GRCh38"="Homo sapiens",
                    "GRCm38"="Mus musculus",
                    "TAIR10"="Arabidopsis thaliana",
                    "EB2"="Bacillus subtilis 168",
                    "UMD3.1"="Bos taurus",
                    "WBcel235"="Caenorhabditis elegans",
                    "CanFam3.1"="Canis familiaris",
                    "GRCz10"="Danio rerio",
                    "BDGP6"="Drosophila melanogaster",
                    "EquCab2"="Equus caballus",
                    "EB1"="Escherichia coli str. K12a",
                    "Galgal4"="Gallus gallus",
                    "Gm01"="Glycine max",
                    "Mmul_1"="Macaca mulatta",
                    "IRGSP-1.0"="Oryza sativa japonica",
                    "CHIMP2.1.4"="Pan troglodytes",
                    "Rnor_6.0"="Rattus norvegicus",
                    "R64-1-1"="Saccharomyces cerevisiae",
                    "EF2"="Schizosaccharomyces pombe",
                    "Sbi1"="Sorghum bicolor",
                    "Sscrofa10.2"="Sus scrofa",
                    "AGPv3"="Zea mays",
                    "hg38"="Homo sapiens",
                    "hg19"="Homo sapiens",
                    "mm10"="Mus musculus",
                    "bosTau8"="Bos taurus",
                    "ce10"="Caenorhabditis elegans",
                    "canFam3"="Canis familiaris",
                    "danRer10"="Danio rerio",
                    "dm6"="Drosophila melanogaster",
                    "equCab2"="Equus caballus",
                    "galGal4"="Gallus gallus",
                    "panTro4"="Pan troglodytes",
                    "rn6"="Rattus norvegicus",
                    "sacCer3"="Saccharomyces cerevisiae",
                    "susScr3"="Sus scrofa")
if(opt$ucscname %in% names(scientificName)){
  scientificName <- scientificName[opt$ucscname]
}else{
  if(opt$species %in% names(scientificName)){
    scientificName <- scientificName[opt$species]
  }else{
    stop("Not a valid genome for enrichment analysis.")
  }
}

scientificName2martTab <- function(.ele){
  .ele <- tolower(strsplit(.ele)[[1]])
  paste0(substring(.ele[1], 1, 1), .ele[2])
}
org <- egOrgMap(scientificName)
lib <- .libPaths()
if(file.access(lib[1], mode=2)!=0){
  pwd <- getwd()
  pwd <- file.path(pwd, "lib")
  dir.create(pwd)
  .libPaths(c(pwd, lib))
}
BiocManager::install(org, update = FALSE, ask = FALSE)
library(org, character.only = TRUE)
organism <- ChIPpeakAnno:::.findKEGGRESTOrganismName(org)
org <- get(org)

shortStrs <- function(strs, len=60){
  if(length(strs)==0) return(strs)
  strs <- as.character(strs)
  shortStr <- function(str, len=60){
    stopifnot(length(str)==1)
    stopifnot(is.character(str))
    if(nchar(str)<=len) return(str)
    strs <- strsplit(str, " ")[[1]]
    nc <- nchar(strs)
    nclast <- nc[length(nc)] + 3
    paste0(substring(str, first = 1, last = len-nclast), "...", strs[length(strs)])
  }
  strs <- sapply(strs, shortStr, len=len)
  make.unique(strs)
}

files <- dir(".", "DiffBind.res..*.all.csv", recursive = TRUE, full.names = TRUE)
gmt <- "ftp.broadinstitute.org://pub/gsea/gene_sets/c2.all.v7.2.symbols.gmt"
for(file in files){
  data <- read.csv(file, row.names = 1)
  if(length(data$gene)!=nrow(data)){
    data$gene <- data$symbol
  }
  data.s <- data[data$FDR<0.05, , drop=FALSE]
  if(nrow(data.s)>1){
    gene.df <- bitr(data.s$gene,
                    fromType=ifelse(grepl("^ENS", as.character(data.s$gene)[1]),
                                    "ENSEMBL", "SYMBOL"),
                    toType = c("ENTREZID", "SYMBOL"), OrgDb = org)
    
    ego <- sapply(c("BP", "MF", "CC"), function(.onto){
      enrichGO(  gene = gene.df$ENTREZID,
                 OrgDb = org,
                 ont = .onto,
                 readable = TRUE
      )
    })
    pff <- file.path("enrichment", basename(dirname(file)))
    dir.create(pff, recursive = TRUE)
    
    null <- mapply(ego, names(ego), FUN=function(.ele, .name){
      write.csv(.ele, file.path(pff, paste0("GO.", .name, ".enrichment.for.fdr0.05.csv")))
      
      .ele <- as.data.frame(.ele)
      if(nrow(.ele)>1){
        .ele$qvalue <- -log10(.ele$p.adjust)
        plotdata <- .ele[!is.na(.ele$qvalue), c("Description", "qvalue", "Count")]
        if(nrow(plotdata)>20) plotdata <- plotdata[1:20, ]
        plotdata$Description <- shortStrs(plotdata$Description)
        ggplot(plotdata, aes(x=reorder(Description, -qvalue), y=qvalue, fill=Count, label=Count)) +
          scale_fill_gradient2(low = muted("blue"), high = muted("red"), oob = scales::squish) +
          geom_bar(stat="identity") + scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
          geom_text(vjust=-.1) +
          xlab("") + ylab("-log10(p-value)") +
          theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
        ggsave(file.path(pff, paste("GO.", .name, ".enrichment.for.fdr0.05.top.pdf", sep = ".")), width = 6, height = 6)
      }
    })
    
    kk <- enrichKEGG(gene = gene.df$ENTREZID, organism = organism)
    kk <- as.data.frame(kk)
    eid <- strsplit(kk$geneID, "\\/")
    symbol <- lapply(eid, function(.ele) gene.df[match(.ele, gene.df$ENTREZID), "SYMBOL"])
    symbol <- sapply(symbol, paste, collapse="/")
    kk$geneSYM <- symbol
    write.csv(kk, file.path(pff, "KEGGenrichment.for.fdr0.05.csv"))
    
    ds <- data$log2FoldChange
    gene.df <- bitr(data$gene,
                    fromType=ifelse(grepl("^ENS", as.character(data.s$gene)[1]),
                                    "ENSEMBL", "SYMBOL"),
                    toType = c("ENTREZID", "SYMBOL"), OrgDb = org)
    names(ds) <- gene.df[match(data$gene, gene.df$ENSEMBL), "ENTREZID"] 
    ds <- ds[!is.na(names(ds))]
    ds <- ds[!is.na(ds)]
    
    p <- file.path(pff, "pathview")
    dir.create(p)
    for (pid in kk$ID[-seq.int(45)]) {
      tryCatch(pv.out <- pathview(gene.data = ds, pathway.id = pid, 
                                  species=organism, kegg.dir = p,
                                  kegg.native=TRUE),
               error=function(.e) message(.e))
    }
    pngs <- dir(".", "pathview.png")
    file.rename(pngs, file.path(p, pngs))
    
    rnk <- file.path(pff, sub(".csv", ".preranked.rnk", file))
    if(scientificName!="Homo sapiens"){
      mart <- useMart("ensembl", paste0(scientificName2martTab(scientificName), "_gene_ensembl"))
      bm <- getBM(values = unique(data$ensembl_id), 
                  attributes = c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"),
                  filters = "ensembl_gene_id",
                  mart = mart)
      data$hsapiens_homolog_associated_gene_name <- 
        bm[match(data$ensembl_id, bm$ensembl_gene_id), 
           "hsapiens_homolog_associated_gene_name"]
      data.rnk <- data[data$stat, c("hsapiens_homolog_associated_gene_name", "stat")]
    }else{
      data.rnk <- data[data$stat, c("gene", "stat")]
    }
    colnames(data.rnk)[1] <- c("hsapiens_gene_name")
    data.rnk <- data[data$stat, c("hsapiens_gene_name", "stat")]
    data.rnk <- data.rnk[!is.na(data.rnk[, 1]), ]
    data.rnk <- data.rnk[data.rnk[, 1]!="", ]
    write.table(data.rnk, file = rnk, quote=FALSE, row.names = FALSE, sep="\t")
    rpt_label <- "c2.all.v7.2"
    outfolder <- file.path(pff, "enrichment")
    
    cmd <- paste("gsea-cli.sh GSEAPreranked -gmx", gmt, "-norm meandiv -nperm 1000 -rnk",
                 rnk, "-scoring_scheme weighted -rpt_label", rpt_label, 
                 "-create_svgs true -make_sets true -plot_top_x 100 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out",
                 outfolder)
    tryCatch(system(cmd), error=function(e){message(e)})
  }
}











