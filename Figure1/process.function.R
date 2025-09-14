library(pacman)
p_load(BiocIO,tools,hbmcbioR,JASPAR2020,readxl,ggtext,htmlwidgets,bootnet,tidygraph,ggraph,igraph,
       motifStack,ggtreeExtra,ggnewscale,ggtree,ShortRead,GenomicAlignments,Rsamtools,BiocParallel,
       DESeq2,enrichplot,DOSE,qpdf,Cairo,clusterProfiler,org.At.tair.db,entropy,BSgenome.Athaliana.TAIR.TAIR9,
       GenomicFeatures,AnnotationDbi,karyoploteR,regioneR,RColorBrewer,
       doParallel,parallel,iterators,foreach,ggpubr,gridExtra,openxlsx,DiffBind,SummarizedExperiment,
       Biobase,MatrixGenerics,rtracklayer,fjComm,GenomicRanges,TFBSTools,universalmotif,motifmatchr,
       Biostrings,GenomeInfoDb,XVector,IRanges,S4Vectors,stats4,BiocGenerics,ggseqlogo,matrixStats,
       reshape2,magrittr,lubridate,forcats,stringr,dplyr,purrr,readr,tidyr,tibble,tidyverse,pacman,
       plotly,ggplot2,glue,ggthemes,grid,tools,stats,graphics
)


mergePeakFiles <- function(path, pattern = "narrowPeak$") {
  peakfiles <- list.files(path = path, pattern = pattern, full.names = TRUE)
  if (length(peakfiles) == 0) {
    stop("No peak files found in the given path.")
  }
  g <- rtracklayer::import(peakfiles[1])
  if (length(peakfiles) > 1) {
    for (i in 2:length(peakfiles)) {
      g2 <- rtracklayer::import(peakfiles[i])
      a <- findOverlaps(g, g2)
      g2 <- g2[-unique(a@to)]
      g <- c(g, g2)
    }
  }

  
  SELEX_6mer_ic<-function(seqs,gapMin=0,gapMax=10)
{
  seqs=fjComm::rmdup(as.data.frame(seqs))[[1]]
  cnts_3mer=fjComm::kmerCntBit_rc(strings = seqs,k=3,diffLen = T,collapse = T,asDf = T,all_possible_k = T,pseudo = 5,rc_combine = T,rmdup = F,rc_k_uniq = F) %>% mutate(freq=counts/sum(counts))
  cnts_3mer=cnts_3mer$freq%>% set_names(cnts_3mer$kmer)
  cnts_g6mer=fjComm::gkmerCntBit_rc(strings = seqs,gapNo = 1,k = 3,gapMins = gapMin, gapMaxs = gapMax, pseudo = 5,diffLen = T,posInfo = F,all_possible_k = T,rmdup = F,melt_result = T,rc_combine = T,rc_k_uniq = F) %>% mutate(freq=counts/sum(counts),kmer1=str_sub(kmer,1,3),kmer2=str_sub(kmer,4,-1))
  cnts_g6mer %<>% mutate(exp_freq=cnts_3mer[kmer1]*cnts_3mer[kmer2]/(gapMax-gapMin+1))
  with(cnts_g6mer,freq*log2(freq/exp_freq)) %>% sum ## KL divergence
}
