
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

# -----------------------------------------------------------
#fig1A####
# -----------------------------------------------------------
       
       df <- read.table(
         "~/2024_WRKY/fig_data.summary/subfamily.txt",
         header = TRUE
       )
       
       
       data <- data.frame(matrix(nrow = 35, ncol = 6))
       colnames(data) <- c("subfamily", "1", "2", "3", "4", "group")
       data[, "subfamily"] <- rep(unique(df$Subfamily), 5)
       data[, "group"] <- rep(c("selex", "meselex", "dap", "amp", "medap"), each = 7)
       
       # Category labels
       chars <- c("1", "2", "3", "4")
       # Legend:
       # 1 = this study
       # 2 = previous study
       # 3 = shared
       # 4 = not available
       
       
       fill_counts <- function(df_col, offset) {
         for (i in seq_along(unique(df$Subfamily))) {
           idx <- which(df$Subfamily == unique(df$Subfamily)[i])
           a <- table(df_col[idx])
           # Ensure all categories exist
           for (char in chars) {
             if (!(char %in% names(a))) a[char] <- 0
           }
           data[i + offset, 2:5] <<- a[order(names(a))]
         }
       }
       
       
       fill_counts(df$SELEX, 0)
       fill_counts(df$Me.SELEX, 7)
       fill_counts(df$DAP, 14)
       fill_counts(df$ampDAP, 21)
       fill_counts(df$MethylampDAP, 28)
       
       
       data$R <- rowSums(data[, 2:5])
       
       data <- data %>%
         mutate(
           x = rep(seq(2, 14, 2), 5),
           y = as.numeric(as.character(gl(5, 7, labels = c(2, 4, 6, 8, 10))))
         )
       
       # Swap columns "2" and "3" (for visualization order)
       data[, c(3, 4)] <- data[, c(4, 3)]
       
       
       df1 <- data.frame(
         x = seq(2, 14, 2),
         y = 11.5,
         label = unique(df$Subfamily)
       )
       
       df2 <- data.frame(
         x = seq(2, 14, 4),
         y = 0,
         label = c("this study", "shared", "previous study", "not available")
       )
       
       df3 <- data.frame(
         x = 16,
         y = c(2, 4, 6, 8, 10),
         label = c("SELEX", "MESELEX", "DAP", "ampDAP", "MethylampDAP")
       )
       
       
       p <- ggplot() +
         geom_scatterpie(
           data = data,
           aes(x, y, group = group, r = 0.7),
           cols = c("1", "2", "3", "4")
         ) +
         coord_equal() +
         theme_void() +
         theme(legend.position = "none") +
         scale_fill_manual(values = c("#DF7380", "#EFB466", "#90C895", "#DADADA")) +
         geom_text(data = df3, aes(x = x, y = y, label = label)) +
         geom_label(data = df1, aes(x = x - 0.1, y = y, label = label)) +
         geom_label(
           data = df2, aes(x = x, y = y, label = label),
           fill = c("#DF7380", "#EFB466", "#90C895", "#DADADA")
         )
       
       # Save to PDF
       ggsave(filename = outfile, plot = p,width = 10, height = 8)



# -----------------------------------------------------------
#fig1B####
# -----------------------------------------------------------


## Compare DAP, ampDAP, and m_ampDAP peaks with public datasets

       
       library(rtracklayer)
       library(ChIPpeakAnno)
       library(eulerr)

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
  
  return(g)
}

public_dap <- mergePeakFiles("2024_WRKY/5_CALLPEAK/Public_DAP_ampDAP/ampDAP/callpeak/")
dap <- mergePeakFiles("2024_WRKY/5_CALLPEAK/DAP_callpeak_minus_pixHalo/temp/")
public_ampdap <- mergePeakFiles("2024_WRKY/5_CALLPEAK/Public_DAP_ampDAP/ampDAP/callpeak/")
ampdap <- mergePeakFiles("2024_WRKY/5_CALLPEAK/AMPdap_callpeak_minus_pixHalo/temp/")

 pixHalo <- readRDS("~/2024_WRKY/5_CALLPEAK/pixhalo.rds")
 dap <- dap[-findOverlaps(dap,pixHalo)@from %>% unique()]
       a <- findOverlapsOfPeaks(unique(dap), unique(public_dap))
       dat <- c(
         "this_study" = a$peaklist$unique.dap. %>% length(),
         "previous_study" = a$peaklist$unique.public_dap. %>% length(),
         "this_study&previous_study" = a$peaklist$`unique.dap.///unique.public_dap.` %>% length()
       )
       p1 <- plot(
         euler(dat),
         quantities = TRUE,
         fills = list(fill = c("#FF6E72", "#F7BC66", "#61CA9D")),
         edges = list(col = "black", lty = 1, lwd = 2),
         main = "DAP"
       )
# -----------------------------------------------------------
#fig1C####
# -----------------------------------------------------------
## mSELEX vs SELEX induced divergence analysis
       library(doParallel)
       library(dplyr)
       library(readr)
       library(universalmotif)
       library(fjComm)
       library(magrittr)
       library(stringr)

SELEX_6mer_ic<-function(seqs,gapMin=0,gapMax=10)
{
  seqs=fjComm::rmdup(as.data.frame(seqs))[[1]]
  cnts_3mer=fjComm::kmerCntBit_rc(strings = seqs,k=3,diffLen = T,collapse = T,asDf = T,all_possible_k = T,pseudo = 5,rc_combine = T,rmdup = F,rc_k_uniq = F) %>% mutate(freq=counts/sum(counts))
  cnts_3mer=cnts_3mer$freq%>% set_names(cnts_3mer$kmer)
  cnts_g6mer=fjComm::gkmerCntBit_rc(strings = seqs,gapNo = 1,k = 3,gapMins = gapMin, gapMaxs = gapMax, pseudo = 5,diffLen = T,posInfo = F,all_possible_k = T,rmdup = F,melt_result = T,rc_combine = T,rc_k_uniq = F) %>% mutate(freq=counts/sum(counts),kmer1=str_sub(kmer,1,3),kmer2=str_sub(kmer,4,-1))
  cnts_g6mer %<>% mutate(exp_freq=cnts_3mer[kmer1]*cnts_3mer[kmer2]/(gapMax-gapMin+1))
  with(cnts_g6mer,freq*log2(freq/exp_freq)) %>% sum ## KL divergence
}

KL_divs=fjComm::mclapply(reads_files,function(file){seqs=readLines(file);SELEX_6mer_ic(seqs)},mc.cores = 20) %>% unlist()

       KL_dist_cutoff <- 0.018
       selex_files <- selex_files[selex_files$KL_divs_6mer > KL_dist_cutoff, ]
       m_selex_files <- m_selex_files[m_selex_files$KL_divs_6mer > KL_dist_cutoff, ]
       name <- intersect(m_selex_files$TF, selex_files$TF) %>% setdiff(too_weak)

       m_df <- m_selex_files[m_selex_files$TF %in% name, ] %>%
         group_by(TF) %>%
         arrange(desc(KL_divs_6mer)) %>%
         slice_head(n = 1) %>%
         ungroup()
       m_df$TF <- paste0(m_df$TF, "_m")
       
       s_df <- selex_files[selex_files$TF %in% name, ] %>%
         group_by(TF) %>%
         arrange(desc(KL_divs_6mer)) %>%
         slice_head(n = 1) %>%
         ungroup()
       s_df$TF <- paste0(s_df$TF, "_s")
       
       merge_df <- bind_rows(s_df, m_df)
       
       # -------------------------------
       # Parallel GK-mer enrichment calculation
       # -------------------------------
       cl <- makeCluster(20)
       registerDoParallel(cl)
       
       result <- foreach(i = seq(nrow(merge_df))) %dopar% {
         library(readr)
         library(universalmotif)
         library(fjComm)
         gk <- fjComm::gkmer_enrich(merge_df$clean_reads[i] %>% read_lines(), len = 8)
       }
       
       stopImplicitCluster()
       stopCluster(cl)
       
       # Combine counts
       thrid <- lapply(result, function(x) x$counts) %>% do.call(cbind, .)
       colnames(thrid) <- merge_df$TF
       
       # -------------------------------
       # Pearson correlation heatmap
       # -------------------------------
       s_cor <- cor(thrid[, 1:(nrow(merge_df)/2)] %>% as.matrix())
       m_cor <- cor(thrid[, (nrow(merge_df)/2+1):(nrow(merge_df))] %>% as.matrix())
       
       heatmap <- m_cor - s_cor
       heatmap %<>%
         set_rownames(rownames(heatmap) %>% str_replace("_m", "")) %>%
         set_colnames(colnames(heatmap) %>% str_replace("_m", ""))
       
       plot <- ggheat(heatmap, clustering = "both") +
         scale_fill_gradientn(
           colours = colorspace::diverge_hcl(100, palette = "Tropic"),
           limits = c(-0.7, 0.7),
           oob = scales::squish,
           name = "Pearson's r\n(mSELEX-SELEX)"
         ) +
         ylab("") + xlab("") +
         gg_theme_Publication(7) +
         gg_axis_x_label_angle(45)
       
       
       # -------------------------------
       # Spearman correlation heatmap
       # -------------------------------
       s_cor <- cor(thrid[, 1:(nrow(merge_df)/2)] %>% as.matrix(), method = "spearman")
       m_cor <- cor(thrid[, (nrow(merge_df)/2+1):(nrow(merge_df))] %>% as.matrix(), method = "spearman")
       
       heatmap <- m_cor - s_cor
       heatmap %<>%
         set_rownames(rownames(heatmap) %>% str_replace("_m", "")) %>%
         set_colnames(colnames(heatmap) %>% str_replace("_m", ""))
       
       plot <- ggheat(heatmap, clustering = "both") +
         scale_fill_gradientn(
           colours = colorspace::diverge_hcl(100, palette = "Tropic"),
           limits = c(-0.4, 0.4),
           oob = scales::squish,
           name = "Spearman's rho\n(mSELEX-SELEX)"
         ) +
         ylab("") + xlab("") +
         gg_theme_Publication(7) +
         gg_axis_x_label_angle(45)
       




# -----------------------------------------------------------
#fig1D####
# -----------------------------------------------------------

  KL_dist_cutoff=0 #0.018
  too_weak=c("WRKY11","WRKY46","WRKY23")
  # TF1="WRKY33";TF2="WRKY35"
  TF1="WRKY71";TF2="WRKY65"
  
  load("~/2024_WRKY/Me-selex/KL_divs_6mer_QC_1.Rdata")
  m_selex_files <- selex_files
  load("~/2024_WRKY/Me-selex/KL_divs_6mer_QC_2.Rdata")
  selex_files <- selex_files[selex_files$KL_divs_6mer>KL_dist_cutoff,]
  m_selex_files <- m_selex_files[m_selex_files$KL_divs_6mer>KL_dist_cutoff,]
  
  name <- intersect(m_selex_files$TF,selex_files$TF) %>% setdiff(too_weak) #
  m_df <- m_selex_files[m_selex_files$TF%in%name &m_selex_files$KL_divs_6mer>KL_dist_cutoff,] %>% group_by(TF) %>% arrange(desc(KL_divs_6mer)) %>% dplyr::slice_head(n = 1) %>% ungroup()
  m_df$TF <- paste0(m_df$TF,"_m")
  s_df <- selex_files[selex_files$TF%in%name &selex_files$KL_divs_6mer>KL_dist_cutoff,]%>% group_by(TF) %>% arrange(desc(KL_divs_6mer)) %>% dplyr::slice_head(n = 1) %>% ungroup()
  s_df$TF <- paste0(s_df$TF,"_s")
  merge_df <- bind_rows(s_df,m_df)
  
  
  wrky33 <- merge_df %>% dplyr::filter(TF %in% c(paste0(c(TF1,TF2),c("_s")),paste0(c(TF1,TF2),c("_m"))))
  
  dis <- fjComm::mclapply(1:4, function(i){fjComm::gkmer_enrich(wrky33$clean_reads[i] %>% read_lines) },mc.cores = 4)
  dis_merge <- lapply(dis, function(x) x$counts) %>% do.call(cbind, .) %>% `colnames<-`(wrky33$TF)
  dis_merge <- as.data.frame(dis_merge)
  dis_merge$kmer <- dis[[1]]$kmer
  dis_merge$gap <- dis[[1]]$gap
  dis_merge <- dis_merge[order(dis_merge[[1]]),]
  dis_merge <- dis_merge %>% mutate(is_duplicate = (lag(.[[1]]) == .[[1]] & lag(.[[2]]) == .[[2]] & lag(.[[3]]) == .[[3]] & lag(.[[4]]) == .[[4]])) %>% filter(is.na(is_duplicate) | !is_duplicate) %>%
    select(-is_duplicate) 
  # dis_merge[,1] <- dis_merge[,1]/1193242
  # dis_merge[,2] <- dis_merge[,2]/773170
  # dis_merge[,3] <- dis_merge[,3]/131897
  # dis_merge[,4] <- dis_merge[,4]/578626
  
  # dis_merge$fold_s <- (dis_merge$WRKY33_s/dis_merge$WRKY35_s) %>% log2()
  # dis_merge$fold_m <- (dis_merge$WRKY33_m/dis_merge$WRKY35_m) %>% log2()
  # plot_df <- dis_merge[,c(3:4,7,8,5)]
  # plot_df$max <- apply(plot_df[,1:2], 1,function(x)(max(x)))
  #
  # plot_df <- plot_df[,3:6]
  # plot_df$color <- "A"
  # plot_df$color[plot_df$fold_s>1]="B"
  # plot_df$color[plot_df$fold_s+1<0]="C"
  # a <- plot_df[which(plot_df$color=="A") %>% sample(20000),]
  # plot_df <- plot_df[-c(which(plot_df$color=="A") ),]
  # plot_df <- rbind(plot_df,a)
  
  plot_s=ggplot(dis_merge,aes(y=(dis_merge[[2]]), x=(dis_merge[[1]])))+geom_hex(bins=70)+
    scale_fill_gradientn(colours = colorspace::sequential_hcl(50,palette = "Mako"),name="Point density",values = c(0,0.01,1))+
    scale_x_continuous(expand = c(0,0.1))+scale_y_continuous(expand = c(0,0.1))+xlab(paste0("8-mer enrichment (",colnames(dis_merge)[1]%>% str_replace("_.",""),")"))+ylab(paste0("8-mer enrichment (",colnames(dis_merge)[2]%>% str_replace("_.",""),""))+gg_theme_Publication(7)
  #
  plot_m=ggplot(dis_merge,aes(y=(dis_merge[[4]]), x=(dis_merge[[3]])))+geom_hex(bins=70)+
    scale_fill_gradientn(colours = colorspace::sequential_hcl(50,palette = "Mako"),name="Point density",values = c(0,0.01,1))+
    scale_x_continuous(expand = c(0,0.1))+scale_y_continuous(expand = c(0,0.1))+xlab(paste0("8-mer enrichment (",colnames(dis_merge)[3]%>% str_replace("_.",""),")"))+ylab(paste0("8-mer enrichment (",colnames(dis_merge)[2]%>% str_replace("_.",""),""))+gg_theme_Publication(7)
  
  # gg_save_pdf(plot_s,7,5,filename = "plot/xyplot_35_33_selex.pdf")
  # gg_save_pdf(plot_m,7,5,filename = "plot/xyplot_35_33_mselex.pdf")
  gg_save_pdf(plot_s,7,5,filename = "plot/xyplot_71_65_selex.pdf")
  gg_save_pdf(plot_m,7,5,filename = "plot/xyplot_71_65_mselex.pdf")

