callpeak_memechip <- function(dap,input,outdir,pixHalo){
  mclapply(1:nrow(bamlist),function(x){
    name <- glue("{dap$TF[x]}_{dap$well[x]}_{dap$batch[x]}")
    cmd <- glue("macs2 callpeak -t {dap$bamfile[x]} -c {input} -n {name} -f BAMPE --outdir {outdir} -g 1.3e8")
    system(cmd)
  },mc.cores = 30)
  
  files <- list.files(outdir,"bed$") %>% paste0(outdir,.)
  cl=makeCluster(30)
  registerDoParallel(cl)
  foreach(x=seq(length(files)))%dopar%{
    library(tools)
    library(glue)
    library(BiocIO)
    library(rtracklayer)
    library(dplyr)
    library(Biostrings)
    library(BSgenome.Athaliana.TAIR.TAIR9)
    library(IRanges)
    library(rtracklayer)
    df <- import.bed(files[x])
    df <- df[-(findOverlaps(df,pixHalo)@from %>% unique())]
    df <- resize(df,200,fix="center")
    df <- df[grep("^Chr", df@seqnames),]
    if (length(df)<=600) {
      df <- df
    }else{
      df <- df[order(df$score,decreasing = TRUE)][1:600]    }
    seqlevels(df, pruning.mode="coarse") <- seqlevels(BSgenome.Athaliana.TAIR.TAIR9::Athaliana)
    seq <- getSeq(BSgenome.Athaliana.TAIR.TAIR9::Athaliana,df)
    name <- file_path_sans_ext(files[x]) %>% basename()# %>% str_extract(pattern = ".*WRKY[0-9]+")
    names(seq) <- df$name
    cmd <- glue("mkdir {ourdir}{name}")
    system(cmd)
    path = glue("{ourdir}{name}.fa")
    writeXStringSet(seq,filepath = path,append = TRUE,format = "fasta")
    cmd <- glue("/wrk/chenhao/meme/bin/meme-chip {ourdir}{name}.fa \\
              -dna -oc {ourdir}{name}/ -meme-mod anr \\
              -streme-nmotifs 0 -meme-nmotifs 3 -order 2 -spamo-skip -fimo-skip")
    system(cmd)
  }
  stopImplicitCluster()
  stopCluster(cl)
}


rename_duplicates <- function(strings) {
  counts <- table(strings)
  seen <- setNames(rep(1, length(counts)), names(counts))
  new_strings <- character(length(strings))
  
  for (i in seq_along(strings)) {
    current_string <- strings[i]
    if (counts[current_string] > 1) {
      seen_count <- seen[current_string]
      new_strings[i] <- paste(current_string, seen_count, sep = "_")
      seen[current_string] <- seen_count + 1
    } else {
      new_strings[i] <- paste(current_string, 1, sep = "_")
    }
  }
  
  return(new_strings)
}
#计算k-mer 并带gap 的富集####
gkmer_enrich<-function(file,len=8,gapMax,pseudo)
{
  w63=file
  w63k=fjComm::gkmerCntBit_rc(w63,gapNo = 1,gapMins = 0,gapMaxs = gapMax,k = len/2,diffLen = T,posInfo = F,all_possible_k = T,pseudo = pseudo,rc_combine = T)
  w63s=universalmotif::shuffle_sequences(w63 %>% DNAStringSet()) %>% as.character()
  w63sk=fjComm::gkmerCntBit_rc(w63s,gapNo = 1,gapMins = 0,gapMaxs = gapMax,k = len/2,diffLen = T,posInfo = F,all_possible_k = T,pseudo = pseudo,rc_combine = T)
  w63e=w63k %>% mutate(counts=counts/w63sk$counts)
  w63e
}
cnt_kmer_rc<-function(file,len=8,pseudo=200)
{
  w63= file %>% readLines()
  w63k=fjComm::kmerCntBit(w63,len,diffLen = T, collapse = T, asDf = T, all_possible_k = T, pseudo = pseudo)
  w63s=universalmotif::shuffle_sequences(w63 %>% DNAStringSet()) %>% as.character()
  w63sk=fjComm::kmerCntBit(w63s,len,diffLen = T, collapse = T, asDf = T, all_possible_k = T, pseudo = pseudo)
  w63sk <-  w63sk[which(w63sk$kmer %in% w63k$kmer),]
  w63e=w63k %>% mutate(counts=(counts-w63sk$counts)/w63sk$counts) #减背景再除背景
  w63e$counts=(w63e$counts + w63e[match(w63e$kmer,w63e$kmer %>% revComp()),]$counts)/2
  w63e
}

  #去除gr对象中和自身重叠的gr对象，不包括自己
remove_gr_grangse <- function(b,minoverlap=50){
  overlaps <- findOverlaps(b,minoverlap = minoverlap, drop.self = TRUE, drop.redundant = TRUE)
  keep <- setdiff(seq_along(b), unique(subjectHits(overlaps)))
  gr_no_overlaps <- b[keep]
  gr_no_overlaps
}

frip <- function(reads, peaks, singleEnd=T) {
  reads <- BamFile(reads)
  overlaps <- summarizeOverlaps(
    peaks,
    reads,
    mode="IntersectionNotEmpty",
    ignore.strand=T,
    singleEnd=singleEnd,
    count.mapped.reads=T
  )
  readsInPeaks <- sum(assay(overlaps))
  if (class(reads) %in% c("BamViews", "BamFile")) {
    allReads <- colData(overlaps)$mapped
  } else {
    allReads <- colData(overlaps)$records
  }
  result <- readsInPeaks/allReads
  
  return(result)
}
#
enrich_in_DNAstringSet <- function(PWMatrix,p_seq,n_seq,background="subject",match="matches",
                                     pseudo=20,p.cutoff=1e-05,...){
  match <- motifmatchr::matchMotifs(pwms = PWMatrix, p_seq,
                                    out = match,bg =background,p.cutoff = p.cutoff)
  match <- match@assays@data@listData$motifMatches %>% colSums()
  df <- motifmatchr::matchMotifs(pwms =PWMatrix,subject = n_seq,
                                 out =match,p.cutoff=p.cutoff,bg=background)
  match_s <- df@assays@data@listData$motifMatches %>% colSums()
  result <- (match+pseudo )/(match_s+pseudo)
  result  
}

motif.enrich.in.TSS <- function(motifs_pwm,seqs,width=4000){
  #'@param motifs_pwm description PWMatrixList object
  #'@param  seqs description DNAstringset object
  Txdb <- TxDb.Athaliana.BioMart.plantsmart51
  seqlevels(Txdb) <- c(paste0("Chr", 1:5),"Mt","Pt")
  seqlevels(Txdb,pruning.mode="coarse")<- seqlevels(Txdb)[str_detect(seqlevels(Txdb),"[cC][hH][rR][0-9]+")]
  promoters_=promoters(Txdb,upstream = width,downstream = width) %>% GenomicRanges::trim() %>% {.[width(.)==2*width]}
  genome=readDNAStringSet("/wrk/data/genome/yz_genome_data/aragenome/ara.fa")
  seqs=genome[promoters_]
  seqnum=length(seqs)
  motifs_rev=lapply(motifs_pwm, function(motif){
    # motif@profileMatrix %>% .[4:1,] %>% set_rownames(qw("A C G T")) %>%  universalmotif::create_motif() %>% {.@name=motif@name;.}
    motif@profileMatrix %>% .[,ncol(.):1] %>% set_rownames(qw("A C G T")) %>%  universalmotif::create_motif() %>% {.@name=motif@name;.}
  })
  motifs_rev_pwm=map(motifs_rev, ~convert_motifs(.x,"TFBSTools-PWMatrix")) %>% {do.call(PWMatrixList,.)}
  motif_hits=motifmatchr::matchMotifs(motifs_pwm, seqs %>% as.character(),
                                      out = "positions",bg ="subject",p.cutoff = 1e-4, 
                                      ranges=GenomicRanges::GRanges(seqnames="Chr1",ranges=IRanges(start = rep(1,seqnum),width = rep(2*width,seqnum) ))  )
  motif_hits_rev=motifmatchr::matchMotifs(motifs_rev_pwm, seqs %>% as.character(), 
                                          out = "positions",bg ="subject",p.cutoff = 1e-4, 
                                          ranges=GenomicRanges::GRanges(seqnames="Chr1",ranges=IRanges(start = rep(1,seqnum),width = rep(2*width,seqnum) ))  )
  
  results=lapply(seq_along(motifs_pwm),function(ind){
    (motif_hits[[ind]]@ranges@start) %>% table() %>% as.data.frame() %>% mutate(motif=names(motifs_pwm)[ind])   #name 是null  
  })
  result_df=do.call(rbind,results) %>% set_colnames(qw("spacing freq motif"))
  result_df %<>% pivot_wider(names_from = spacing,values_from = freq)
  result_df %<>% column_to_rownames("motif")
  
  results_rev=lapply(seq_along(motifs_pwm),function(ind){
    (motif_hits_rev[[ind]]@ranges@start) %>% table() %>% as.data.frame() %>% mutate(motif=names(motifs_pwm)[ind])
  })
  result_rev_df=do.call(rbind,results_rev) %>% set_colnames(qw("spacing freq motif"))
  result_rev_df %<>% pivot_wider(names_from = spacing,values_from = freq)
  result_rev_df %<>% column_to_rownames("motif")
  
  
  
  trans_dist=result_df %>% as.matrix() %>% {.[is.na(.)]=0;.} %>% sweep(1,rowSums(.),"/") %>%
    apply(1,function(x)zoo::rollmean(x,200)) %>% t
  trans_dist_rev=result_rev_df %>% as.matrix() %>% {.[is.na(.)]=0;.} %>% sweep(1,rowSums(.),"/") %>%
    apply(1,function(x)zoo::rollmean(x,200)) %>% t
  
  trans_dist=trans_dist-trans_dist_rev
  col_draw=colorspace::sequential_hcl(200,palette = "PuBuGn",rev = T)#
  ComplexHeatmap::Heatmap(trans_dist,col = col_draw, cluster_rows = F, show_row_dend = F, cluster_columns =F,use_raster = T,
                          row_names_gp = gpar(fontsize=7),
                          column_labels = c("-4kb",rep("",(length(trans_dist[1,])/2)-1 %>% round()),
                                            "TSS",rep("",(length(trans_dist[1,])/2)-1 %>% round()),"4kb")    )
  
}
