# -----------------------------------------------------------
#fig3A####
# -----------------------------------------------------------
## -------------------------------
## Motif alignment, clustering, and circular tree visualization with stars
## -------------------------------

library(universalmotif)
library(motifStack)
library(ggtree)
library(ggnewscale)
library(ggplot2)
library(fjComm)
library(stringr)
library(dplyr)
library(ade4)

# Load PWM list
pwmlist <- read_jaspar("2024_WRKY/8_motif_net/ALL_mono_P_S_T.pwm")

# Align motifs
pwmlist_alig <- DNAmotifAlignment(
  threshold = 0,
  pfms = pwmlist %>% convert_motifs(class = "motifStack-pfm"),
  revcomp = rep(FALSE, length(pwmlist))
)

# Trim PWM matrices (remove flanking columns)
for (i in seq_along(pwmlist_alig)) {
  pwmlist_alig[[i]]@mat <- pwmlist_alig[[i]]@mat[, -c(1:2, 16:17)]
}

# Rename motifs for clarity
c <- pwmlist_alig
for (i in seq_along(c)) {
  c[[i]]@name <- c[[i]]@name %>%
    gsub("P", "_1", .) %>%
    gsub("S", "_2", .) %>%
    gsub("T", "_3", .) %>%
    gsub("F", "_4", .)
  names(c)[i] <- c[[i]]@name
}

# Construct motif correlation matrix for clustering
p <- data.frame()
for (i in seq_along(c)) {
  p1 <- t(rbind(c[[i]]@mat[1, ], c[[i]]@mat[2, ], c[[i]]@mat[3, ], c[[i]]@mat[4, ]))
  p <- if (i == 1) p1 else cbind(p, p1)
}
colnames(p) <- names(c)
p <- cor(p, method = "spearman")

# Hierarchical clustering
a <- dist(x = p, method = "maximum")
hc <- hclust(a)

# Convert to phylog object for circular tree plotting
phylog <- ade4::hclust2phylog(hc, add.tools = TRUE)
p <- ggtree(phylog, layout = "circular", size = 0.6) +
  geom_text(aes(label = node), size = 1.5) +
  geom_tiplab(hjust = -0.5, size = 4.2, align = TRUE, offset = 0) +
  ggtree::rotate(220) # Rotate tree

# Prepare star plot data for each PWM
dt <- lapply(c, function(pwm) {
  df <- data.frame()
  for (i in seq_len(ncol(pwm@mat))) {
    id <- pwm@name
    group <- which.max(pwm@mat[, i]) %>% names()
    size <- pwm@mat[which.max(pwm@mat[, i]), i]
    strain <- i
    df <- rbind(df, c(id, group, size, strain))
  }
  colnames(df) <- c("id", "group", "size", "strain")
  df$size <- as.numeric(df$size)
  df$strain <- as.numeric(df$strain)
  df
}) %>% do.call(rbind, .)
row.names(dt) <- NULL

# Separate motif names into two groups
p1 <- names(c)[str_detect(names(c), "m")]
p2 <- names(c)[!str_detect(names(c), "m")]

# Generate geom_strip strings
text1 <- paste0("geom_strip(\"", p1, "\", \"", p1, "\", offset = 0.06, barsize = 20, extend = 0.6, color = \"#00BFC4\", offset.text = 3)") %>% paste0(collapse = "+")
text2 <- paste0("geom_strip(\"", p2, "\", \"", p2, "\", offset = 0.06, barsize = 20, extend = 0.6, color = \"#F8766D\", offset.text = 3)") %>% paste0(collapse = "+")

# Generate class strips based on clustering
tree <- cutree(hc, h = 0.37)
colors <- c("#FFFFB3","#E5C494","#B3DE69","#85CBBF","#F0766D","#80A9CB","#B6B2D2","#F4C5DD","#FF9D00","#F35EAB","#6A3D9A")
class_strips <- lapply(seq_along(colors), function(i){
  class_nodes <- names(tree[tree == i])
  paste0("geom_strip(\"", class_nodes, "\", \"", class_nodes, "\", offset = 0.15, barsize = 10, extend = 0.6, color = \"", colors[i], "\", angle = 10, hjust = \"center\", fontsize = 4, offset.text = 2)") %>% paste0(collapse = "+")
}) %>% paste0(collapse = "+")

# Combine all strips
plot_text <- paste0("p + ", text1, "+", text2, "+", class_strips)
plot <- eval(parse(text = plot_text))

# Add stars to circular tree
p1 <- plot + geom_fruit(
  data = dt,
  geom = geom_star,
  mapping = aes(x = strain, y = id, fill = group, size = size),
  starstroke = 0,
  starshape = 13,
  offset = 0.24,
  pwidth = 0.30
) +
  scale_size_continuous(range = c(0, 4), limits = c(min(dt$size), max(dt$size)), breaks = c(1, 2, 3)) +
  scale_fill_manual(values = c("#109648", "#255C99", "#F7B32B", "#D62839"),
                    guide = guide_legend(keywidth = 0.4, keyheight = 0.4, order = 4,
                                         override.aes = list(starstroke = 0.3))) +
  geom_tiplab(hjust = -0.5, size = 4.2, align = TRUE, offset = 0) +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.key.size = unit(10, "lines")
  )

# Save plot
ggsave(filename = "2024_WRKY/8_motif_net/plot/fig3a.pdf", plot = p1, width = 20, height = 20)


# -----------------------------------------------------------
#fig3B####
# -----------------------------------------------------------

# =============================

library(openxlsx)
library(magrittr)
library(stringr)
library(foreach)
library(doParallel)
library(Biostrings)
library(motifmatchr)
library(fjComm)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(reshape2)

# =============================
prepare_motifs <- function(pwm_file, n_motifs = 11){
  pwm <- read_jaspar(pwm_file) %>%
    convert_motifs("TFBSTools-PWMatrix") %>%
    do.call(PWMatrixList, .)
  pwm <- pwm[1:n_motifs]
  names(pwm@listData) <- sapply(seq_len(n_motifs), function(x){pwm@listData[[x]]@name})
  return(pwm)
}

prepare_selex_files <- function(excel_file){
  selex <- read.xlsx(excel_file, sheet = "SELEX")
  m_selex <- read.xlsx(excel_file, sheet = "M_SELEX")
  
  selex$TF <- paste0(selex$TF, "_s")
  m_selex$TF <- paste0(m_selex$TF, "_m")
  
  bind_rows(selex, m_selex)
}

# =============================

compute_motif_enrichment <- function(selex_files, pwm, pseudo = 1000, ncores = 30){
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  rich_list <- foreach(i = seq_len(nrow(selex_files))) %dopar% {
    library(magrittr)
    library(Biostrings)
    library(motifmatchr)
    library(fjComm)
    
    seqs <- fjComm::getSeq_fqfachrFile(selex_files$for_motif_disc[i])
    seqlen <- names(which.max(table(nchar(seqs)))) %>% as.integer()
    seqs <- fjComm::length_adjust(seqs, seqlen)
    seqs <- c(seqs, revComp(seqs)) %>% DNAStringSet()
    
    a <- matchMotifs(pwm, seqs, out = "matches", p.cutoff = 1e-5)
    a_counts <- colSums(a@assays@data@listData$motifMatches)
    
    b <- matchMotifs(pwm, shuffle_sequences(seqs), out = "matches", p.cutoff = 1e-5)
    b_counts <- colSums(b@assays@data@listData$motifMatches)
    
    counts <- (a_counts + pseudo) / (b_counts + pseudo)
    names(counts) <- names(pwm)
    counts
  }
  
  stopImplicitCluster()
  stopCluster(cl)
  
  names(rich_list) <- paste0(selex_files$TF, "_", selex_files$repeats)
  return(rich_list)
}


# =============================
process_enrichment_matrix <- function(rich_list, pwm_names){
  df <- largeListToDf(rich_list) %>% t()
  colnames(df) <- pwm_names
  
  df <- t(apply(df, 1, function(x){ log2(x + 1) / max(log2(x + 1)) }))
  df
}


# =============================
plot_heatmap_bar <- function(df, group1_idx, group2_idx, output_file){
  df1 <- df[group1_idx, ]
  df2 <- df[group2_idx, ]
  df2 <- df2[order(rownames(df2)), ]
  
  SELEX <- colMeans(df1) - 2 * colMeans(df1)
  MESELEX <- colMeans(df2)
  
  # 柱状图
  bar_df <- rbind(SELEX, MESELEX) %>% melt()
  bar_plot <- ggplot(bar_df, aes(x = Var2, y = value, fill = Var1)) +
    geom_col(position = "identity") +
    coord_flip() +
    scale_y_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1, 1)) +
    scale_fill_manual(values = c("SELEX" = "#F8766D", "MESELEX" = "#00BFC4")) +
    theme_bw()
  

  heatmap_plot <- function(mat){
    df_melt <- melt(as.matrix(mat))
    ggplot(df_melt, aes(x = Var2, y = Var1, fill = value)) +
      geom_tile() +
      scale_fill_gradientn(colours = brewer.pal(9, "PuBu"), limits = c(0, 1), name = "log2x/log2max") +
      xlab("") + ylab("") +
      theme_bw() + theme(panel.grid = element_blank())
  }
  
  p1 <- heatmap_plot(df1)
  p2 <- heatmap_plot(df2)
  

  plot <- plot_grid(p1, p2, bar_plot, ncol = 3)
  ggsave(filename = output_file, plot = plot, width = 18, height = 5)
}

# =============================

pwm <- prepare_motifs("2024_WRKY/8_motif_net/class11_pwm.txt")
selex_files <- prepare_selex_files("2024_WRKY/WRKY_For_analysis.xlsx")


rich_in_101N <- compute_motif_enrichment(selex_files, pwm, pseudo = 1000, ncores = 30)
save(rich_in_101N, file = "2024_WRKY/8_motif_net/class11_in_all_sample.Rdata")

df <- process_enrichment_matrix(rich_in_101N, names(pwm))

plot_heatmap_bar(df, group1_idx = 1:99, group2_idx = 100:139, 
                 output_file = "2024_WRKY/8_motif_net/plot/class11_in_all_sample_1.pdf")













# -----------------------------------------------------------
#fig3D####
# -----------------------------------------------------------




# =============================
# 1. Load genome and prepare motif
# =============================
library(Biostrings)
library(TFBSTools)
library(motifmatchr)
library(GenomicRanges)
library(rtracklayer)
library(reshape2)
library(ggplot2)
library(cowplot)
library(purrr)
library(readxl)
library(dplyr)

# Load genome sequences (first 5 chromosomes for demo)
genomeDna <- readDNAStringSet("~/aragenome/ara.fa")[1:5]

# Parameters
p.cutoff <- 1e-04
disp_rng <- 500

# Create PRE motif
WRKY_PRE <- create_motif("TACTGCGCTTAGT") %>% convert_motifs("TFBSTools-PWMatrix")

# =============================
# 2. Load DAP-seq files and filter
# =============================
mDAPfiles <- read_excel("~/2024_WRKY/WRKY_For_analysis.xlsx", sheet = 5, col_names = TRUE) %>%
  filter(`repeat` == 1)

aDAPfiles <- read_excel("~/2024_WRKY/WRKY_For_analysis.xlsx", sheet = 4, col_names = TRUE) %>%
  filter(`repeat` == 1) %>%
  .[str_detect(.$TF, "WRKY"), ]

# =============================
# 3. Match motif in genome
# =============================
matches_ <- matchMotifs(WRKY_PRE, genomeDna, out = "positions", p.cutoff = p.cutoff, bg = "even")

# Convert matches to GRanges for each chromosome
matches_ <- lapply(seq_along(matches_), function(li) {
  Match <- map2(
    c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"), 1:5,
    ~GRanges(
      seqnames = .x,
      ranges = ranges(matches_[[li]][[.y]]),
      strand = matches_[[li]][[.y]]@elementMetadata$strand
    )
  ) %>% suppressWarnings() %>%
    do.call(c, .) %>%
    `seqlengths<-`(seqlengths(genomeDna)) %>%
    GenomicRanges::trim()
})
names(matches_) <- names(matches_)

# =============================
# 4. Plot PRE motif logo
# =============================
plot_PRE <- getSeq(genomeDna, matches_[[1]]) %>%
  create_motif() %>%
  {ggseqlogo(.@motif) + gg_theme_Publication()}

gg_save_pdf(plot_PRE, 6, 3, "~/work/WRKY_paper/Enrich_PRE_DAP/",
            paste0("motif_PRE_p", p.cutoff))

# =============================
# 5. Convert bigwig to coverage
# =============================
aBwCvg <- map(aDAPfiles$bwfile, ~rtracklayer::import(.x) %>% coverage(., weight = .$score)) %>%
  `names<-`(aDAPfiles$TF)

mBwCvg <- map(mDAPfiles$bwfile, ~rtracklayer::import(.x) %>% coverage(., weight = .$score)) %>%
  `names<-`(mDAPfiles$TF)

# Get TFs present in both datasets
TFs <- intersect(names(mBwCvg), names(aBwCvg))

# =============================
# 6. Generate coverage plots for each TF
# =============================
plotList <- lapply(TFs, function(TF) {
  matchesExt <- matches_[[1]] %>% resize(2 * disp_rng + 1, "center") %>% .[width(.) == 2 * disp_rng + 1]
  
  # Normalize coverage to sum = 1
  Df <- data.frame(
    mAmp = mBwCvg[[TF]][matchesExt] %>% as.matrix() %>% colSums() %>% {./sum(.)},
    amp  = aBwCvg[[TF]][matchesExt] %>% as.matrix() %>% colSums() %>% {./sum(.)},
    Pos  = -disp_rng:disp_rng
  ) %>% melt("Pos")
  
  ggplot(Df, aes(Pos, value, color = variable)) +
    geom_line() +
    # geom_smooth(se = FALSE, method = "loess", span = 0.005, linewidth = 0.5) +
    gg_theme_Publication() +
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      plot.background = element_blank(),
      panel.background = element_blank()
    ) +
    labs(title = TF)
})

# Combine plots in grid
plot <- do.call(cowplot::plot_grid, plotList)
gg_save_pdf(plot, 6 * 5, 4 * 4, "~/work/WRKY_paper/Enrich_PRE_DAP/",
            paste0("PRE_p", p.cutoff))

# =============================
# 7. Coverage plots without y-axis
# =============================
plotList_nony <- lapply(TFs, function(TF) {
  matchesExt <- matches_[[1]] %>% resize(2 * disp_rng + 1, "center") %>% .[width(.) == 2 * disp_rng + 1]
  
  Df <- data.frame(
    mAmp = mBwCvg[[TF]][matchesExt] %>% as.matrix() %>% colSums() %>% {./sum(.)},
    amp  = aBwCvg[[TF]][matchesExt] %>% as.matrix() %>% colSums() %>% {./sum(.)},
    Pos  = -disp_rng:disp_rng
  ) %>% melt("Pos")
  
  ggplot(Df, aes(Pos, value, color = variable)) +
    geom_line() +
    scale_colour_manual(values = c("#5475A5", "#B92F74")) +
    # geom_smooth(se = FALSE, method = "loess", span = 0.005, linewidth = 0.5) +
    gg_theme_Publication() +
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      plot.background = element_blank(),
      panel.background = element_blank(),
      axis.text.y = element_blank()
    ) +
    labs(title = TF)
})

plot_nony <- do.call(cowplot::plot_grid, plotList_nony)
gg_save_pdf(plot_nony, 24, 16, "~/work/WRKY_paper/Enrich_PRE_DAP/",
            paste0("PRE_p", p.cutoff, "_nony"))





# -----------------------------------------------------------
#fig3F####
# -----------------------------------------------------------



# =========================================
# 1. Load required libraries and set parameters
# =========================================
library(Biostrings)
library(TFBSTools)
library(universalmotif)
library(motifmatchr)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(magrittr)
library(zoo)
library(colorspace)
library(ComplexHeatmap)

width <- 5000
# Load Arabidopsis TxDb
Txdb <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
seqlevels(Txdb) <- c(paste0("Chr", 1:5), "Mt", "Pt")
seqlevels(Txdb, pruning.mode="coarse") <- seqlevels(Txdb)[str_detect(seqlevels(Txdb), "[cC][hH][rR][0-9]+")]

# Load genome sequence
genome <- readDNAStringSet("~/aragenome/ara.fa")

# =========================================
# 2. Extract promoter sequences (±5kb)
# =========================================
promoters_ <- promoters(Txdb, upstream = width, downstream = width) %>%
  GenomicRanges::trim() %>%
  {.[width(.) == 2 * width]}  # only keep exact width

seqs <- genome[promoters_]
seqnum <- length(seqs)

# =========================================
# 3. Load motifs and clean PWM
# =========================================
class_motif <- read_jaspar("~/2024_WRKY/8_motif_net/class11_pwm.txt") %>%
  convert_motifs("TFBSTools-PWMatrix")

# Remove columns 11:13 from profileMatrix
for (i in seq_along(class_motif)) {
  class_motif[[i]]@profileMatrix <- class_motif[[i]]@profileMatrix[, -c(11:13)]
}

motifs_pwm <- do.call(PWMatrixList, class_motif)

# =========================================
# 4. Define motif enrichment function around TSS
# =========================================
motif.enrich.in.TSS <- function(motifs_pwm, width = 5000) {

  # Create reversed motifs for potential reverse complement analysis
  motifs_rev <- lapply(motifs_pwm, function(motif) {
    motif@profileMatrix[, ncol(motif@profileMatrix):1] %>%
      set_rownames(c("A", "C", "G", "T")) %>%
      universalmotif::create_motif() %>%
      {.@name = motif@name; .}
  })
  motifs_rev_pwm <- map(motifs_rev, ~convert_motifs(.x, "TFBSTools-PWMatrix")) %>%
    {do.call(PWMatrixList, .)}

  # Match motifs to promoter sequences
  motif_hits <- motifmatchr::matchMotifs(
    motifs_pwm,
    seqs %>% as.character(),
    out = "positions",
    bg = "subject",
    p.cutoff = p.cutoff,
    ranges = GRanges(
      seqnames = "Chr1",
      ranges = IRanges(start = rep(1, seqnum), width = rep(2 * width, seqnum))
    )
  )

  # Summarize motif hits per position
  results <- lapply(seq_along(motifs_pwm), function(ind) {
    motif_hits[[ind]]@ranges@start %>%
      table() %>%
      as.data.frame() %>%
      mutate(motif = motifs_pwm[[ind]]@name)
  })

  result_df <- do.call(rbind, results) %>%
    set_colnames(c("spacing", "freq", "motif")) %>%
    pivot_wider(names_from = spacing, values_from = freq) %>%
    column_to_rownames("motif")

  # Normalize and smooth using rolling mean
  trans_dist <- result_df %>%
    as.matrix() %>%
    {.[is.na(.)] <- 0; .} %>%
    sweep(1, rowSums(.), "/") %>%
    apply(1, function(x) zoo::rollmean(x, 20)) %>%
    t

  # Color palette for heatmap
  col_draw <- colorspace::sequential_hcl(200, palette = "PuBuGn", rev = TRUE)

  # Plot heatmap
  ComplexHeatmap::Heatmap(
    trans_dist,
    col = col_draw,
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    cluster_columns = FALSE,
    use_raster = TRUE,
    row_names_gp = gpar(fontsize = 12)
  )
}

# =========================================
# 5. Run motif enrichment and save to PDF
# =========================================
pdf(file = "2024_WRKY/8_motif_net/motif_around_in_TSS_1e-4.pdf", width = 5, height = 5)
motif.enrich.in.TSS(motifs_pwm, width = width)
dev.off()



# -----------------------------------------------------------
#fig3G RNA Analysis####
# -----------------------------------------------------------
{
  # ------------------------
  # 1. Load RNA-seq count data
  # ------------------------
  df1 <- read.table(
    "/wrk/manana/work/Nana4_3_2__qcfs_RNAseq/pre_processed/result/Nana4_3_2__qcfs_RNAseq.count",
    header = TRUE
  )[, 1:14]
  
  df2 <- read.table(
    "/wrk/manana/work/Nana4_2__Cotyledon_Root_RNASeq/pre_processed/result/Nana4_2__Cotyledon_Root_RNASeq.count",
    header = TRUE
  )[, 9:10]
  
  # Combine datasets and extract relevant columns
  df <- cbind(df1, df2)
  count_table <- df[, 7:16]
  colnames(count_table) <- c(
    "Cotyledon","Cotyledon","Flower","Flower","Silique","Silique",
    "Stem","Stem","Root","Root"
  ) %>% rename_duplicates()
  
  # Calculate TPM
  rnaseq_tpm <- t(t(expr1) / colSums(expr1)) * 1e6
  rnaseq_tpm <- cbind(df[, 1:6], rnaseq_tpm)
  rnaseq_tpm$Geneid <- gsub(".Araport11.447", replacement = "", rnaseq_tpm$Geneid)

  # ------------------------
  # 2. Parallel motif enrichment for top 10% genes per tissue
  # ------------------------
  library(parallel)
  library(parallelly)
  library(foreach)
  library(doParallel)
  
  cl <- makeCluster(20)
  registerDoParallel(cl)
  
  result <- foreach(x = seq(1:10)) %dopar% {
    # Load libraries inside foreach
    library(SummarizedExperiment)
    library(universalmotif)
    library(dplyr)
    library(Biostrings)
    library(motifmatchr)
    library(fjComm)
    
    pseudo <- 1000
    
    # Identify top 10% up- and down-regulated genes
    up_gene <- rnaseq_tpm[order(rnaseq_tpm[, x + 6], decreasing = TRUE), 1][1:round(nrow(rnaseq_tpm) * 0.10)]
    down_gene <- rnaseq_tpm[order(rnaseq_tpm[, x + 6], decreasing = FALSE), 1][1:round(nrow(rnaseq_tpm) * 0.10)]
    
    # Extract promoters
    library(TxDb.Athaliana.BioMart.plantsmart51)
    promoter <- promoters(TxDb.Athaliana.BioMart.plantsmart51, upstream = 300, downstream = 300) %>%
      GenomicRanges::trim() %>%
      {.[width(.) == 600]}
    promoter$tx_name <- gsub("\\.[0-9]+", "", promoter$tx_name)
    
    # Get sequences for up-regulated genes
    up_gene <- promoter[which(promoter$tx_name %in% up_gene)] %>% {.[!duplicated(.$tx_name)]}
    seqlevels(up_gene, pruning.mode = "coarse") <- seqlevels(BSgenome.Athaliana.TAIR.TAIR9::Athaliana)
    up_seq <- getSeq(BSgenome.Athaliana.TAIR.TAIR9::Athaliana, up_gene)
    
    # Get sequences for down-regulated genes
    down_gene <- promoter[which(promoter$tx_name %in% down_gene)] %>% {.[!duplicated(.$tx_name)]}
    seqlevels(down_gene, pruning.mode = "coarse") <- seqlevels(BSgenome.Athaliana.TAIR.TAIR9::Athaliana)
    down_seq <- getSeq(BSgenome.Athaliana.TAIR.TAIR9::Athaliana, down_gene)
    
    # Motif matching for up-regulated genes
    df <- motifmatchr::matchMotifs(pwms = pwm, subject = up_seq, out = "matches", p.cutoff = 5e-05, bg = "subject")
    match <- df@assays@data@listData$motifMatches %>% colSums()
    up_data <- match %>% as.data.frame() %>% `colnames<-`("value")
    up_data$group <- colnames(rnaseq_tpm)[x + 6]
    up_data$motif <- names(pwm.l)
    
    # Motif matching for down-regulated genes
    df <- motifmatchr::matchMotifs(pwms = pwm, subject = down_seq, out = "matches", p.cutoff = 5e-05, bg = "subject")
    match <- df@assays@data@listData$motifMatches %>% colSums()
    down_data <- match %>% as.data.frame() %>% `colnames<-`("value")
    down_data$group <- colnames(rnaseq_tpm)[x + 6]
    down_data$motif <- names(pwm.l)
    
    # Calculate ratio up/down
    up_data$value <- up_data$value / down_data$value
    up_data
  }
  
  stopImplicitCluster()
  stopCluster(cl)
  
  # ------------------------
  # 3. Combine results and organize factors
  # ------------------------
  a <- bind_rows(result)
  a <- mutate(a,
              motif = fct_relevel(motif, "class11", "class10", "class9", "class8", "class7",
                                  "class6", "class5", "class4","class3","class2","class1"))
  a$group <- a$group %>% gsub(".[0-9]+_[0-9]+_", "", .) %>%
    fct_relevel("Flower_1", "Flower_2", "Cotyledon_1", "Cotyledon_2",
                "Silique_1", "Silique_2", "Root_1", "Root_2", "Stem_1", "Stem_2")
  
  # ------------------------
  # 4. Plot heatmap
  # ------------------------
  col_draw <- colorspace::diverge_hcl(200, palette = "Blue-Red3", rev = FALSE)
  p <- ggplot(a, aes(x = group, y = as.factor(motif), fill = log2(value))) +
    geom_tile() +
    coord_fixed(ratio = 1) +
    scale_fill_gradientn(colors = col_draw, name = "log2 Enrichment", limits = c(-1, 1)) +
    theme(
      axis.text.x = element_text(angle = 30, size = 10, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 15)
    ) +
    labs(title = "10%")
  
  ggsave("~/2024_WRKY/8_motif_net/diff_tissus_RNAseq_top15per/RNA_3.pdf", plot = p, width = 6, height = 6)
}


# -----------------------------------------------------------
#fig3G ATI Analysis####
# -----------------------------------------------------------


          
{
  # ------------------------
  # 1. Load ATI data
  # ------------------------
  file <- read.table("2024_WRKY/8_motif_net/diff_tissues_ATI.txt", header = TRUE)
  file$tissue <- rename_duplicates(file$tissue)
  
  # ------------------------
  # 2. Parallel motif enrichment
  # ------------------------
  cl <- makeCluster(21)
  registerDoParallel(cl)
  
  result <- foreach(x = seq(nrow(file))) %dopar% {
    library(SummarizedExperiment)
    library(universalmotif)
    library(dplyr)
    library(Biostrings)
    library(motifmatchr)
    library(fjComm)
    
    pseudo <- 1000
    
    # Read sequences from file
    seq <- readLines(file$clean_reads_30n[x]) %>% DNAStringSet()
    
    # Motif matching
    df <- motifmatchr::matchMotifs(pwms = pwm, subject = seq, out = "matches", p.cutoff = 1e-05, bg = "subject")
    match <- df@assays@data@listData$motifMatches %>% colSums()
    
    # Motif matching with shuffled sequences
    df_s <- motifmatchr::matchMotifs(pwms = pwm, subject = seq %>% universalmotif::shuffle_sequences(),
                                     out = "matches", p.cutoff = 1e-05, bg = "subject")
    match_s <- df_s@assays@data@listData$motifMatches %>% colSums()
    
    # Calculate enrichment ratio
    heat_data <- (match + pseudo) / (match_s + pseudo) %>% as.data.frame() %>% `colnames<-`("value")
    heat_data$group <- file$tissue[x]
    heat_data$motif <- names(pwm.l)
    heat_data
  }
  
  stopImplicitCluster()
  stopCluster(cl)
  
  # ------------------------
  # 3. Save results
  # ------------------------
  save(result, file = "2024_WRKY/8_motif_net/diff_tissue_ATI.Rdata")
}


# -----------------------------------------------------------
#fig3G ATAC Analysis####
# -----------------------------------------------------------

 load("~/2024_WRKY/8_motif_net/diff_tissue_ATAC.Rdata")
  
  result <- bind_rows(result) 
  shuff_result <- bind_rows(shuff_result)
  shuff_result$group %<>%  paste0("Shuffle_",.)
  result <- rbind(shuff_result,result)
  result <- mutate(result,motif = fct_relevel(motif, 
                                              "class1", "class2", "class3", 
                                              "class4", "class5", "class6", 
                                              "class7", "class8","class9","class10","class11"))
  result <- mutate(result,group = fct_relevel(group,
                                              "flower1","flower2","flower3",
                                              "fruit1", "fruit2", "fruit3",
                                              "leaf1",  "leaf2",  "leaf3",
                                              "root1",  "root2",  "root3",
                                              "stem1",  "stem2",  "stem3",
                                              "Shuffle_flower1","Shuffle_flower2","Shuffle_flower3",
                                              "Shuffle_fruit1", "Shuffle_fruit2", "Shuffle_fruit3",
                                              "Shuffle_leaf1",  "Shuffle_leaf2",  "Shuffle_leaf3",
                                              "Shuffle_root1",  "Shuffle_root2",  "Shuffle_root3",
                                              "Shuffle_stem1",  "Shuffle_stem2",  "Shuffle_stem3"))
  
  plot1 <- ggplot(data = result,mapping = aes(x=group,y = as.factor(motif),fill=log2(value)))+
    geom_tile()+ 
    coord_fixed(ratio = 1)+
    scale_fill_gradientn(colors =  colorspace::diverge_hcl(200,palette = "Blue-Red3",rev = F),name = "log2 Enrichment",limits=c(-1,1))+
    theme(
      axis.text.x = element_text(angle = 30,size = 10,hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title =element_text(size=15)
    )    +labs(title = "ATAC pseudo=1000")
  plot1
  ggsave(plot = plot1,"~/2024_WRKY/8_motif_net/plot/diff_tissue_ATAC.pdf",width = 8,height = 8)

          
          
