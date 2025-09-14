
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
         "~/subfamily.txt",
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
