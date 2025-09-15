  df_1 <- read_xlsx("~/WRKY_For_analysis.xlsx",sheet = "ampDAP")[1:94,]
  df_2 <- read_xlsx("~/WRKY_For_analysis.xlsx",sheet = "M_ampDAP")
  df_S <- read_xlsx("~/WRKY_For_analysis.xlsx",sheet = "SELEX")
  df_M <- read_xlsx("~/WRKY_For_analysis.xlsx",sheet = "M_SELEX")
  
  name <- intersect(df_1$TF,df_2$TF)
  name <- intersect(df_M$TF,df_S$TF)
  name <- intersect(intersect(df_1$TF,df_2$TF), intersect(df_M$TF,df_S$TF))
  
  df_1 <- df_1[which(sapply(1:nrow(df_1), function(x){ df_1$TF[x]%in%name && df_1$QC[x] > 0.045})),]
  df_2 <- df_2[which(sapply(1:nrow(df_2), function(x){ df_2$TF[x]%in%name && df_2$QC[x] > 0.045})),]
  
  df_M <- df_M[which(sapply(1:nrow(df_M), function(x){ df_M$TF[x]%in%name && df_M$KL_divs_6mer[x] > 0.018})),]
  df_S <- df_S[which(sapply(1:nrow(df_S), function(x){ df_S$TF[x]%in%name && df_S$KL_divs_6mer[x] > 0.018})),]
  name <- intersect(intersect(df_1$TF,df_2$TF), intersect(df_M$TF,df_S$TF))
  
  
  
  wrky <- "WRKY70"
  kmer_num <- 2000
  peak_num <- 1000
  #kmer 、enrichment
  for (i in 1:length(name)) {
    file <- df_S[df_S$TF==name[i],]
    seqs <- getSeq_fqfachrFile(file$clean_reads[which.max(file$KL_divs_6mer)])
    seqlens <- 101L #seqs %>% nchar %>% table() %>% which.max() %>% names %>% as.integer()
    seqs <- fjComm::length_adjust(seqs, seqlens)
    Selex <- gkmer_enrich(file = seqs,len = 8,gapMax = 8,pseudo = 200)
    
    file <- df_M[df_M$TF==name[i],]
    seqs <- getSeq_fqfachrFile(file$clean_reads[which.max(file$KL_divs_6mer)])
    seqlens <- 101L #seqs %>% nchar %>% table() %>% which.max() %>% names %>% as.integer()
    seqs <- fjComm::length_adjust(seqs, seqlens)
    M_Selex <- gkmer_enrich(file = seqs,len = 8,gapMax = 8,pseudo = 200)
    
    M_Selex <- M_Selex[order(M_Selex$counts,decreasing = T),]
    Selex <-     Selex[order(Selex$counts,decreasing = T),]
    save(Selex,M_Selex,file = glue("~/ROC_curve/{name[i]}_kmer_enrich_score.Rdata"))
  }
  
  
  # save(Selex,M_Selex,file = glue("~/2024_WRKY/5_CALLPEAK/ROC_curve/{wrky}_kmer_enrich_score.Rdata"))
  # load(glue("~/ROC_curve/{wrky}_kmer_enrich_score.Rdata"))
  
  # Selex$counts <- log2((Selex$counts)/max(log2(Selex$counts)))
  # M_Selex$counts <- log2(M_Selex$counts)/max(log2(M_Selex$counts))
  # Selex <-  Selex[!duplicated(Selex$counts), ][1:kmer_num,]
  
  ### @@@
  # Selex <-  Selex[!duplicated(Selex$counts), ][1:kmer_num,]
  
  
  kmer_file <- list.files("~/ROC_curve/","_kmer_enrich_score.Rdata")
  
  for (i in 1:length(kmer_file)) {
    load(paste0("~/ROC_curve/",kmer_file[i],recycle0 = ""))
    
    Selex <-  Selex[1:kmer_num,]
    Selex$kemr_g <-  sapply(1:nrow(Selex), function(x){
      if(Selex$gap[x] %>% gsub("n","",.) %>% as.numeric() >0){
        paste0(str_split(string = Selex$kmer[x],"")[[1]][1:4] %>% paste0(collapse = ""),
               rep("N",Selex$gap[x] %>% gsub("n","",.) %>% as.numeric()) %>% paste0(collapse = ""),
               str_split(string = Selex$kmer[x],"")[[1]][5:8] %>% paste0(collapse = ""),collapse = "")
      }else{  Selex$kmer[x]}  })
    
    ### @@@
    # M_Selex <-  M_Selex[!duplicated(M_Selex$counts), ][1:kmer_num,]
    M_Selex <-  M_Selex[1:kmer_num,]
    M_Selex$kemr_g <-  sapply(1:nrow(M_Selex), function(x){
      if(M_Selex$gap[x] %>% gsub("n","",.) %>% as.numeric() >0){
        paste0(str_split(string =  M_Selex$kmer[x],"")[[1]][1:4] %>% paste0(collapse = ""),
               rep("N",M_Selex$gap[x] %>% gsub("n","",.) %>% as.numeric()) %>% paste0(collapse = ""),
               str_split(string = M_Selex$kmer[x],"")[[1]][5:8] %>% paste0(collapse = ""),collapse = "")
      }else{   M_Selex$kmer[x]}  })
    
    #peak sequence
    pixHalo <- readRDS("~/pixhalo.rds")
    genome <- getSeq(BSgenome.Athaliana.TAIR.TAIR9)[1:5]
    genome <- GRanges(seqnames = names(genome),
                      ranges = IRanges(start = rep(1, length(genome)), width = width(genome)),
                      strand = rep("*", length(genome)))
    wrky <- str_extract(kmer_file[i],"WRKY[0-9]+")
    peak_1 <- df_1[df_1$TF==wrky,] %>% arrange(desc(QC)) %>% .$callpeak %>% .[1] %>% rtracklayer::import()
    peak_2 <- df_2[df_2$TF==wrky,] %>% arrange(desc(QC)) %>% .$callpeak %>% .[1] %>% rtracklayer::import()
    
    ### @@@
    for (ii in 1:2) {
      used_peak_file=ii
      df=list(peak_1,peak_2)[[used_peak_file]]
      
      df <- df[-findOverlaps(df,pixHalo,minoverlap = 0.1)@from %>% unique()]
      df <- df[grep("^Chr", df@seqnames),]
      negSet <- setdiff(genome,df)
      df <- df[order(df$signalValue,decreasing = TRUE)][1:peak_num]
      df %<>% resize(200,"center")
      seq <- getSeq(BSgenome.Athaliana.TAIR.TAIR9,df)
      n_seq <- getSeq(BSgenome.Athaliana.TAIR.TAIR9,negSet)#%>% shuffle_sequences()
      
      negSet__ <- lapply(seq_along(genome), function(chr){
        data.frame(
          seqnames = genome[chr]@seqnames,
          start = seq(1, genome[chr]@ranges@width, 200),
          strand = "*",
          end = c(seq(1, genome[chr]@ranges@width, 200)[-1] - 1, genome[chr]@ranges@width))
      }) %>% do.call(rbind, .)
      
      negSet_gr <- makeGRangesFromDataFrame(negSet__, ignore.strand = T, seqinfo = NULL, seqnames.field = "seqnames",
                                            start.field = "start", end.field = "end", strand.field = "strand")
      lap <- findOverlaps(negSet_gr, df, minoverlap = 1)
      negSet_gr <- negSet_gr[-(lap@from) %>% unique()] #%>% getSeq(BSgenome.Athaliana.TAIR.TAIR9,.)
      
      
      GenomicRanges::mcols(negSet_gr)$score <- 0
      # SELEX pred
      skmer <- DNAStringSet(Selex$kemr_g)
      s <- lapply(1:length(skmer), function(x){
        #获取匹配到同一个bin的所有kmer的富集值的和
        ### @@@
        hits=vcountPattern(pattern = skmer[[x]],subject =seq ,fixed = F ,max.mismatch = 0)
        Selex$counts[x]*hits
      })
      s <- do.call(cbind,s) %>% apply(1,sum)
      library(doParallel)
      library(foreach)
      cl=makeCluster(20)
      registerDoParallel(cl)
      GR <- foreach(x=seq(length(skmer)),.combine = function(y,z){
        score <- if (is.null(y)) z else y + z
      })%dopar%{
        library(pacman)
        pacman::p_load(Biostrings,GenomicRanges,IRanges,BiocGenerics,magrittr)
        df <- vmatchPattern(pattern =  skmer[[x]],subject = n_seq,fixed = F,max.mismatch = 0) %>% as.data.frame()
        df$group_name <- as.character(negSet@seqnames)[df$group]
        df$start <- df$start+negSet@ranges@start[df$group]
        df$end <- df$end+negSet@ranges@start[df$group]
        GR <- makeGRangesFromDataFrame(df, ignore.strand = T, seqinfo = NULL, seqnames.field = "group_name",
                                       start.field = "start", end.field = "end") #%>% resize( width(.) + 200, fix = "center")# %>% remove_gr_grangse(minoverlap = 1)
        z <- negSet_gr
        a <- findOverlaps(z,GR) %>% as.data.frame()
        a <- a[!duplicated(a$queryHits ), ]
        z[a$queryHits]$score <- z[a$queryHits]$score + Selex$counts[x]
        z$score
      }
      stopImplicitCluster()
      stopCluster(cl)
      n <- GR#negSet_gr$score
      
      data <- c(s,n) %>%as.data.frame() %>% `colnames<-`("predictions")
      data$labels <- c(rep(1,length(s)),rep(0,length(negSet_gr$score)))
      data$predictions[is.na(data$predictions)] =0
      # data$predictions <- log2(data$predictions)
      data$predictions[data$predictions=="-Inf"]=0
      S= data
      
      
      # mSELEX pred
      skmer <- DNAStringSet(M_Selex$kemr_g)
      s <- lapply(1:length(skmer), function(x){
        #获取匹配到同一个bin的所有kmer的富集值的和
        ### @@@
        hits=vcountPattern(pattern = skmer[[x]],subject =seq ,fixed = F ,max.mismatch = 0)
        M_Selex$counts[x]*hits
      })
      s <- do.call(cbind,s) %>% apply(1,sum)
      cl=makeCluster(20)
      registerDoParallel(cl)
      GR <- foreach(x=seq(length(skmer)),
                    .combine = function(y,z){
                      score <- if (is.null(y)) z else y + z
                    })%dopar%{
                      library(pacman)
                      pacman::p_load(Biostrings,GenomicRanges,IRanges,BiocGenerics,magrittr)
                      df <- vmatchPattern(pattern =  skmer[[x]],subject = n_seq,fixed = F,max.mismatch = 0) %>% as.data.frame()
                      df$group_name <- as.character(negSet@seqnames)[df$group]
                      df$start <- df$start+negSet@ranges@start[df$group]
                      df$end <- df$end+negSet@ranges@start[df$group]
                      GR <- makeGRangesFromDataFrame(df, ignore.strand = T, seqinfo = NULL, seqnames.field = "group_name",
                                                     start.field = "start", end.field = "end") #%>% resize( width(.) + 200, fix = "center")# %>% remove_gr_grangse(minoverlap = 1)
                      z <- negSet_gr
                      a <- findOverlaps(z,GR) %>% as.data.frame()
                      a <- a[!duplicated(a$queryHits ), ]
                      z[a$queryHits]$score <- z[a$queryHits]$score + M_Selex$counts[x]
                      z$score
                    }
      stopImplicitCluster()
      stopCluster(cl)
      n <-GR# negSet_gr$score
      
      data <- c(s,n) %>%as.data.frame() %>% `colnames<-`("predictions")
      data$labels <- c(rep(1,length(s)),rep(0,length(negSet_gr$score)))
      data$predictions[is.na(data$predictions)] =0
      # data$predictions <- log2(data$predictions)
      data$predictions[data$predictions=="-Inf"]=0
      M= data
      name <- str_extract(kmer_file[i],"WRKY[0-9]+")
      save(S,M,file = glue("~/ROC_curve/{ii}_{name}_S_M_score.Rdata"))
      
    }
    
    
  }
  
  name <- list.files("/ROC_curve/",pattern = "_S_M_score.Rdata") %>% str_extract("WRKY[0-9]+") %>% unique()
  mclapply( 1:length(name),FUN = function(i){
    wrky <- name[i]
    file <- list.files("/ROC_curve/",pattern = glue("{wrky}.S_M_score.Rdata")) %>% paste0("~/ROC_curve/",.)
    load(file[1])

    Q3 <- quantile(S$predictions, 0.99); upper_bound <- Q3 *2  ;S$predictions[S$predictions>upper_bound] <- mean(S$predictions)
    Q3 <- quantile(M$predictions, 0.99); upper_bound <- Q3 *2  ;M$predictions[M$predictions>upper_bound] <- mean(M$predictions)
    pr_data <-  list(S = pr_curve_fj(S$predictions, S$labels),
                     M = pr_curve_fj(M$predictions, M$labels))
    Df_All <- map2(pr_data, names(pr_data), function(pr, prName) {
      pr[, 1:2] %>% as_tibble() %>% mutate(group = prName) } ) %>% do.call(rbind, .)
    
    p=ggplot(Df_All,aes(x = Recall, y = Precision,color=group))  +
      geom_line(linewidth=1) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F) +
      gg_theme_Publication() +
      ggplot2::theme(axis.ticks = element_blank()) +
      scale_color_manual(values = brewer.pal(3,"Set1")) +
      scale_y_continuous(breaks = c(0,1)) +
      scale_x_continuous(breaks = c(0, 1))+
      ggplot2::theme(axis.title.x = element_text(size = 15),
                     axis.title.y = element_text(size = 15),
                     axis.text.x =  element_text(size = 15),
                     axis.text.y =  element_text(size = 15))
    
    load(file[2])
    Q3 <- quantile(S$predictions, 0.99); upper_bound <- Q3 *10  ;S$predictions[S$predictions>upper_bound] <- mean(S$predictions)
    Q3 <- quantile(M$predictions, 0.99); upper_bound <- Q3 *10  ;M$predictions[M$predictions>upper_bound] <- mean(M$predictions)
    pr_data <-  list(S = pr_curve_fj(S$predictions, S$labels),
                     M = pr_curve_fj(M$predictions, M$labels))
    Df_All <- map2(pr_data, names(pr_data), function(pr, prName) {
      pr[, 1:2] %>% as_tibble() %>% mutate(group = prName) } ) %>% do.call(rbind, .)
    
    q=ggplot(Df_All,aes(x = Recall, y = Precision,color=group))  +
      geom_line(linewidth=1) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F) +
      gg_theme_Publication() +
      ggplot2::theme(axis.ticks = element_blank()) +
      scale_color_manual(values = brewer.pal(3,"Set1")) +
      scale_y_continuous(breaks = c(0,1)) +
      scale_x_continuous(breaks = c(0, 1))+
      ggplot2::theme(axis.title.x = element_text(size = 15),
                     axis.title.y = element_text(size = 15),
                     axis.text.x =  element_text(size = 15),
                     axis.text.y =  element_text(size = 15))
    plot <- plot_grid(p,q,nrow = 1)
    ggsave(filename = glue("/ROC_curve/{wrky}_PR_curve.pdf"),plot = plot,width = 8,height = 4)
  },mc.cores = 2)
  
  
  name <- c("WRKY4","WRKY14","WRKY15","WRKY28",
            "WRKY29","WRKY35","WRKY71")
  auc <- mclapply(1:length(name),FUN = function(x){
    wrky <- name[x]
    file <- list.files("/ROC_curve/",pattern = glue("{wrky}.S_M_score.Rdata")) %>% paste0("~/ROC_curve/",.)
    load(file[1])
    Q3 <- quantile(S$predictions, 0.99); upper_bound <- Q3 *10  ;S$predictions[S$predictions>upper_bound] <- mean(S$predictions)
    Q3 <- quantile(M$predictions, 0.99); upper_bound <- Q3 *10  ;M$predictions[M$predictions>upper_bound] <- mean(M$predictions)
    pr_data <-  list(S = pr_curve_fj(S$predictions, S$labels),
                     M = pr_curve_fj(M$predictions, M$labels))
    amp_s <- integrate(f = approxfun(pr_data$S$Recall, pr_data$S$Precision), lower = min(pr_data$S$Recall), upper = max(pr_data$S$Recall), subdivisions = 2000)$value
    amp_m <- integrate(f = approxfun(pr_data$M$Recall, pr_data$M$Precision), lower = min(pr_data$M$Recall), upper = max(pr_data$M$Recall), subdivisions = 2000)$value
    
    load(file[2])
    Q3 <- quantile(S$predictions, 0.99); upper_bound <- Q3 *10  ;S$predictions[S$predictions>upper_bound] <- mean(S$predictions)
    Q3 <- quantile(M$predictions, 0.99); upper_bound <- Q3 *10  ;M$predictions[M$predictions>upper_bound] <- mean(M$predictions)
    pr_data <-  list(S = pr_curve_fj(S$predictions, S$labels),
                     M = pr_curve_fj(M$predictions, M$labels))
    mamp_s <- integrate(f = approxfun(pr_data$S$Recall, pr_data$S$Precision), lower = min(pr_data$S$Recall), upper = max(pr_data$S$Recall), subdivisions = 2000)$value
    mamp_m <- integrate(f = approxfun(pr_data$M$Recall, pr_data$M$Precision), lower = min(pr_data$M$Recall), upper = max(pr_data$M$Recall), subdivisions = 2000)$value
    
    area <- c(amp_s,amp_m,mamp_s,mamp_m)
    area
  },mc.cores = 2)
  
  auc <- auc %>% largeListToDf() %>% as.matrix()
  rownames(auc) <- c("amp_s","amp_m","mamp_s","mamp_m")
  colnames(auc) <- name
  save(auc,file = "2024_WRKY/5_CALLPEAK/ROC_curve/oneoff/AUC.Rdata")
  
  df <- auc[1:2,] %>% `ROWNAMES<-`(c("SELEX","MSELEX")) %>% as.matrix() %>% melt()
  p1 <- ggplot(df,aes(x=Var2,y=value,color=Var1,group=Var1))+
    geom_line(linewidth=1)+
    labs(title = "ampdap")+xlab("")+ylab("auc")+
    scale_y_continuous(limits  = c(0,0.3)) + 
    scale_color_manual(values = brewer.pal(3,"Set1")[1:2] %>% rev())+
    ggplot2::theme(axis.title.x = element_text(size = 15),
                   axis.title.y = element_text(size = 15),
                   axis.text.x =  element_text(size = 10,angle = 90),
                   axis.text.y =  element_text(size = 10))
  p1
  df <- auc[3:4,] %>% `ROWNAMES<-`(c("SELEX","MSELEX")) %>% as.matrix() %>% melt()
  p2 <- ggplot(df,aes(x=Var2,y=value,color=Var1,group=Var1))+
    geom_line(linewidth=1)+
    labs(title = "mampdap")+xlab("")+ylab("auc")+
    scale_y_continuous(limits  = c(0,0.3)) + 
    scale_color_manual(values = brewer.pal(3,"Set1")[1:2] %>% rev())+
    ggplot2::theme(axis.title.x = element_text(size = 15),
                   axis.title.y = element_text(size = 15),
                   axis.text.x =  element_text(size = 10,angle = 90),
                   axis.text.y =  element_text(size = 10))
  p2
  plot <- plot_grid(p1,p2)
  plot
  ggsave(filename = "/ROC_curve/AUC_area.pdf",plot,width = 5,height = 2)
  auc[3,]<auc[4,]
  auc[1,]>auc[2,]
  
  
  dat <- auc %>% melt()
  dat$paired <- dat$Var2
  dat$Var1 <- as.character(dat$Var1)
  dat$Var2 <- "amp"
  dat$Var2[grep("mamp",dat$Var1)] <- "mamp"
  dat$Var1[grep("_s",dat$Var1)] <- "aSELEX"
  dat$Var1[grep("_m",dat$Var1)] <- "mSELEX"
  plot_auc <- ggplot(dat, aes(x = Var1, y = value, color = Var1)) +
    geom_boxplot(outlier.size = 2) +
    geom_line(aes(group = paired), color = "grey80", size = 0.5) +
    geom_point(size = 2) +
    scale_y_continuous(limits = c(0,0.3))+
    scale_color_manual(values =  c("#00AFBB", "#E7B800"))+
    stat_compare_means(comparisons = list(c("aSELEX", "mSELEX")), 
                       paired = T, size = 2.5,label = "p.format",tip.length = 0)+ 
    facet_wrap(~ Var2, scales = 'free_y', nrow = 1)+ylab("AUC")+xlab("")
  
  plot_auc
