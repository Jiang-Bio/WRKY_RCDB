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




remove_reverse_complementary <- function(dna_seqs) {
  seqs <- as.character(dna_seqs)
  seqs_seen <- character()
  seqs_to_remove <- character()
  reverse_complement <- function(dna) {
    rev_comp <- reverseComplement(DNAString(dna))
    return(as.character(rev_comp))
  }
  for (seq in seqs) {
    rev_comp_seq <- reverse_complement(seq)
    
    if (rev_comp_seq %in% seqs_seen) {
      seqs_to_remove <- c(seqs_to_remove, seq)
    } else {
      seqs_seen <- c(seqs_seen, seq, rev_comp_seq)
    }
  }
  
  unique_seqs <- setdiff(seqs, seqs_to_remove)
  return(DNAStringSet(unique_seqs))
}

reverse_if_motif_found <- function(dna_seqs, motif = "TTGAC") {
  result <- DNAStringSet()
  for (i in 1:length(dna_seqs)) {
    seq <- dna_seqs[i]
    seq_str <- as.character(seq)

    motif_pos <- gregexpr(motif, seq_str, fixed = TRUE)[[1]]
    if (motif_pos[[1]]==c(-1)) {
      result <- c(result, seq)
    }else if(length(motif_pos) == 1 || length(motif_pos) == 2|| length(motif_pos) == 3){
      result <- c(result, reverseComplement(seq))
    }else{
      result <- c(result, seq)
    }
  }
  return(result)
}

classify_motif_pattern <- function(dna_seqs, motifs = c("GTCAA", "GTCTT")) {
  motifs_dna <- DNAStringSet(motifs)
  motifs_rc <- as.character(reverseComplement(motifs_dna))
  motif_len <- width(motifs_dna[1])  
  
  result <- character(length(dna_seqs))
  
  for (i in seq_along(dna_seqs)) {
    seq <- as.character(dna_seqs[i])
    
    all_matches <- data.frame(start = integer(), dir = character(), motif = character())
    

    for (motif in motifs) {
      pos <- unlist(gregexpr(motif, seq, fixed = TRUE))
      if (pos[1] != -1) {
        all_matches <- rbind(all_matches, 
                             data.frame(start = pos, dir = "F", motif = motif))
      }
    }
    
    for (motif in motifs_rc) {
      pos <- unlist(gregexpr(motif, seq, fixed = TRUE))
      if (pos[1] != -1) {
        all_matches <- rbind(all_matches, 
                             data.frame(start = pos, dir = "R", motif = motif))
      }
    }
    
    if (nrow(all_matches) == 0) {
      result[i] <- "no motif"
      next
    }
    
    if (nrow(all_matches) > 1) {
      all_matches <- all_matches[order(all_matches$start), ]
      filtered <- all_matches[1, , drop = FALSE]
      for (j in 2:nrow(all_matches)) {
        if (all_matches$start[j] - filtered$start[nrow(filtered)] >= motif_len) {
          filtered <- rbind(filtered, all_matches[j, ])
        }
      }
    } else {
      filtered <- all_matches
    }
    
    n_hits <- nrow(filtered)
    
    if (n_hits == 1) {
      result[i] <- "monomer"
    } else if (n_hits == 2) {
      gap <- filtered$start[2] - filtered$start[1] - motif_len
      dir1 <- filtered$dir[1]
      dir2 <- filtered$dir[2]
      
      if (dir1 == "F" && dir2 == "F") {
        result[i] <- paste0("DR", gap)
      } else if (dir1 == "F" && dir2 == "R") {
        result[i] <- paste0("ER", gap)
      } else if (dir1 == "R" && dir2 == "F") {
        result[i] <- paste0("IR", gap)
      } else if (dir1 == "R" && dir2 == "R") {
        result[i] <- paste0("DR", gap)  # 反向+反向也认为是DR
      } else {
        result[i] <- "unknown"
      }
    } else {
      result[i] <- "multiple motifs"
    }
  }
  
  return(result)
}
sort_motif_classes <- function(classes) {
  extract_number <- function(x) {
    as.numeric(sub("^[A-Z]+", "", x))
  }
  types <- c("monomer", "DR", "ER", "IR","multiple motifs","no motif")
  type_order <- match(sub("[0-9]+$", "", classes), types)
  number_order <- ifelse(classes == "monomer", 0, extract_number(classes))
  
  sorted <- order(type_order, number_order)
  return(sorted)
}



plot_track <- function(zoom_region,bwfile,y_max=40,track_color="#80B1D3",peak_color="#FB8072",sitename="TF",
                       show_gene=TRUE,show_peak=NULL,genome=NULL,ylables=c("10","40"),tick.pos=c(10,40),
                       label_margin=0.08,show_tatget_gene=NULL){
  #'@param zoom_region Grange object
  #'@param show_peak Grange object
  #'@param show_tatget_gene gene ID
  if (is.null(genome)) {
    genome <- readDNAStringSet("/wrk/data/genome/yz_genome_data/aragenome/ara.fa")}
  
  genome <- GRanges(seqnames = names(genome),
                    ranges = IRanges(start = rep(1, length(genome)), width = width(genome)),
                    strand = rep("*", length(genome)))
  
  #plot边距
  pp <- getDefaultPlotParams(plot.type=1)
  pp$leftmargin <- 0.15
  pp$topmargin <- 15
  pp$bottommargin <- 15
  pp$ideogramheight <- 5
  pp$data1inmargin <- 10
  if (show_gene) {
    txdbbb <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
    seqlevels(txdbbb, pruning.mode="coarse") <- seqlevels(BSgenome.Athaliana.TAIR.TAIR9::Athaliana)

    kp <- plotKaryotype(zoom = zoom_region,genome = genome,plot.type =4,plot.params = pp)
    genes.data <- makeGenesDataFromTxDb(txdbbb,
                                        karyoplot=kp,
                                        plot.transcripts = TRUE, 
                                        plot.transcripts.structure = TRUE)
    genes.data <- mergeTranscripts(genes.data)
    if (!is.null(show_tatget_gene)) {
      message(genes.data$genes$gene_id)
      num <- which(genes.data$genes$gene_id==show_tatget_gene)
      genes.data$genes <- genes.data$genes[num]
      genes.data$transcripts <- genes.data$transcripts[num]
      genes.data$coding.exons <- genes.data$coding.exons[num]
      genes.data$non.coding.exons <- genes.data$non.coding.exons[num]
    }
    
    kpPlotGenes(kp, data=genes.data, r0=0, r1=0.05,add.strand.marks = TRUE,mark.height = 0.5,
                introns.col = "gray",mark.width = 1,mark.distance = 0.5,marks.col = "gray")
  }else {
    kp <- plotKaryotype(zoom = g,genome = genome,plot.type =4,plot.params = pp)
  }
  kpAddBaseNumbers(kp, tick.dist = 1000,tick.len = 2.5 ,minor.tick.len = 1,minor.tick.dist = 500,
                   add.units = TRUE,units = "kb" ,cex=0.5, digits = 2)
  
  if (!is.null(show_peak)) {
    peakinfo <-  show_peak 
    kpPlotRegions(kp, peakinfo, col=peak_color, r0=0.15, r1=0.18)
    kpAddLabels(kp, labels = glue("{sitename} \n binding site"), r0=0.15, r1=0.18, cex = 0.8)
  }
  
  for(i in seq_len(length(bwfile))) {
    bigwig.file <-bwfile[i]
    at <- autotrack(i, length(bwfile), r0=0.3, r1=1,margin = 0.4)
    kp <- kpPlotBigWig(kp, data=bigwig.file, ymax=y_max,r0=at$r0, r1=at$r1,col = track_color)
    computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    kpAxis(kp, r0=at$r0, r1=at$r1,cex=0.5, ymin=0, ymax=y_max,tick.pos = tick.pos,
           labels = ylables )
    kpAddLabels(kp, labels = names(bwfile)[i], r0=at$r0, r1=at$r1,cex=0.5,srt=0,
                pos=3, offset = 10,label.margin = label_margin)
  }
}


TF_logo <- function(pwmlist,outfile,weigh=NULL){
  #'@ pwmlist  PWMatrixList
  motif_pics=BiocParallel::bplapply(1:length(pwmlist),function(x){
    t_pwmlist <- universalmotif::convert_motifs(pwmlist[[x]])
    motif_pic=suppressMessages(
      universalmotif::view_motifs(t_pwmlist,show.positions = F,tryRC = T,
                                  min.mean.ic = 0,fit.to.height = 1,names.pos = "right",
                                  normalise.scores = TRUE,min.position.ic = 0,
                                  text.size = 4,min.height = 0,min.overlap = 20)+
        ylab("")+
        scale_y_continuous() +
        ggplot2::theme(axis.ticks.y =element_blank(),
                       axis.line.y = element_blank(),
                       axis.text.y = element_blank())  )
    ggplotGrob(motif_pic)
  },BPPARAM = BiocParallel::MulticoreParam(workers = 30))
  
  if (length(pwmlist)>100) {
    blank_plot <- ggplot() +
      theme_void()+
      xlim(0, 6) +
      ylim(0, 100)#)#+xlim(0,2)
    anno_plots <- c()
    anno_texts <- c()
    for (i in 1:100) {
      cmd <- glue("annotation_custom(motif_pics[[{i}]],
                    xmin =0.5,
                    xmax={0.5+ncol(pwmlist[[i]]@profileMatrix)*0.1},
                    ymin={99-i},
                    ymax={99-(i-2.5)})")
      cmd2 <- glue("annotate('text', 
                    x = 0.4, 
                    y = {100.5-(i)}, 
                    label = '{pwmlist[[i]]@name}',size = 2)")
      anno_plots <- append(anno_plots,cmd)
      anno_texts <- append(anno_texts,cmd2)
    }
    for (i in 101:length(pwmlist)) {
      cmd <- glue("annotation_custom(motif_pics[[{i}]],xmin =3.5,xmax={3.5+
                ncol(pwmlist[[i]]@profileMatrix)*0.1},ymin={99-i+100},ymax={99-(i-2.5)+100})")
      cmd2 <- glue("annotate('text', x = 3.4, y = {100.5-(i)+100}, label = '{pwmlist[[i]]@name}',size = 2)")
      anno_plots <- append(anno_plots,cmd)
      anno_texts <- append(anno_texts,cmd2)
    }
    anno_plots %<>% paste0(collapse = "+")
    anno_texts %<>% paste0(collapse = "+")
    plot_text=paste0("blank_plot+ ",anno_plots,"+",anno_texts)
    plot_1 <- eval(parse(text = plot_text))
    ggsave(filename = outfile,plot = plot_1,width = 5,height = 20)
  }else{
    if (is.null(weigh)) {
      blank_plot <- ggplot() +
        theme_void()+
        xlim(0, 3) +
        ylim(0, 100)
    }else{
      blank_plot <- ggplot() +
        theme_void()+
        xlim(0, weigh) +
        ylim(0, 100)
    }
    
    anno_plots <- c()
    anno_texts <- c()
    for (i in 1:length(pwmlist)) {
      cmd <- glue("annotation_custom(motif_pics[[{i}]],xmin =0.5,xmax={0.5+
                length(pwmlist[[i]]@profileMatrix[1,])*0.1},ymin={99-i},ymax={99-(i-2.5)})")
      cmd2 <- glue("annotate('text', x = 0.4, y = {100.5-(i)}, label = '{pwmlist[[i]]@name}',size = 1)")
      anno_plots <- append(anno_plots,cmd)
      anno_texts <- append(anno_texts,cmd2)
    }
    anno_plots %<>% paste0(collapse = "+")
    anno_texts %<>% paste0(collapse = "+")
    plot_text=paste0("blank_plot+ ",anno_plots,"+",anno_texts)
    plot_1 <- eval(parse(text = plot_text))
    ggsave(filename = outfile,plot = plot_1,width = 5,height = 20)
  }
  
}
