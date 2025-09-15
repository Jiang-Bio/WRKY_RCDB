############################################################
# Script: WRKY35 motif analysis and coverage plotting
#   1. Load WRKY35 motifs from JASPAR
#   2. Import genome fasta and peak files
#   3. Filter peaks by removing background overlaps
#   4. Match motifs across the genome
#   5. Extract coverage signals and plot
############################################################

# -------------------------
# Load required packages
# -------------------------
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(BiocGenerics)
library(BiocIO)
library(dplyr)
library(purrr)
library(stringr)
library(reshape2)
library(ggplot2)
library(cowplot)

# 读取 motif
WRKY35 <- read_jaspar("~/WRKY35_s1_m1m2.txt")
names(WRKY35) <- map(WRKY35, ~ .x@name)
WRKY35 <- convert_motifs(WRKY35, "TFBSTools-PWMatrix") %>%
  do.call(PWMatrixList, .)

# 读取基因组序列
genomeDna <- readDNAStringSet("/wrk/data/genome/yz_genome_data/aragenome/ara.fa")[1:5]
p.cutoff  <- 1e-04
disp_rng  <- 500

# Peak 文件
peakFiles <- c(
  DAP    = "~/WRKY35_B11_Nana3_3_3__DAP_WRKY_WG_peaks.narrowPeak",
  AmpDAP = "~/WRKY35_B10_Nana3_3__AmpDAP_WRKY_WG_50ul_5ul_peaks.narrowPeak",
  mAmpDAP = "~/WRKY35_C5_peaks.narrowPeak"
)

# BigWig 文件
bwFiles <- c(
  DAP    = "~/bwfiles/B11_23_WRKY35_WG_5ul_D2.bs10.bw",
  AmpDAP = "~/bwfiles/B10_22_WRKY35_WG_5ul_A1.bs10.bw",
  mAmpDAP = "~/bwfiles/C5_29_WRKY35_mA1.bs10.bw"
)

# 处理 peaks
peakGR <- lapply(seq_along(peakFiles), function(peak) {
  dap_peak <- rtracklayer::import(peakFiles[peak])
  seqlevels(dap_peak, pruning.mode = "coarse") <- seqlevels(genomeDna)
  seqlengths(dap_peak) <- seqlengths(genomeDna)

  ## 移除背景
  if (str_detect(names(peakFiles)[peak], "^DAP")) {
    pixHalo_d <- BiocIO::import("~/callpeak/G3_75_Dap3_4_input1_peaks.narrowPeak")
    pixHalo_d_ov <- findOverlaps(pixHalo_d, dap_peak)
    overlap_width <- pintersect(
      pixHalo_d[pixHalo_d_ov@from],
      dap_peak[pixHalo_d_ov@to]
    ) %>% width()
    overlaps <- pixHalo_d_ov@to[
      (overlap_width / width(pixHalo_d[pixHalo_d_ov@from])) > 0.5 &
      (overlap_width / width(dap_peak[pixHalo_d_ov@to])) > 0.5
    ]
    dap_peak <- dap_peak[-overlaps]

  } else if (str_detect(names(peakFiles)[peak], "^AmpDAP")) {
    pixHalo_a <- BiocIO::import("~/callpeak/input_F4_Nana3_2__WRKY_WG_50ul_5ul_DAPAmpDAP_peaks.narrowPeak")
    pixHalo_a_ov <- findOverlaps(pixHalo_a, dap_peak)
    overlap_width <- pintersect(
      pixHalo_a[pixHalo_a_ov@from],
      dap_peak[pixHalo_a_ov@to]
    ) %>% width()
    overlaps <- pixHalo_a_ov@to[
      (overlap_width / width(pixHalo_a[pixHalo_a_ov@from])) > 0.5 &
      (overlap_width / width(dap_peak[pixHalo_a_ov@to])) > 0.5
    ]
    dap_peak <- dap_peak[-overlaps]

  } else if (str_detect(names(peakFiles)[peak], "^mAmpDAP")) {
    pixHalo_m <- BiocIO::import("~/callpeak_with_input/Nana3_8_A4_pIXHalo_peaks.narrowPeak")
    pixHalo_m_ov <- findOverlaps(pixHalo_m, dap_peak)
    overlap_width <- pintersect(
      pixHalo_m[pixHalo_m_ov@from],
      dap_peak[pixHalo_m_ov@to]
    ) %>% width()
    overlaps <- pixHalo_m_ov@to[
      (overlap_width / width(pixHalo_m[pixHalo_m_ov@from])) > 0.5 &
      (overlap_width / width(dap_peak[pixHalo_m_ov@to])) > 0.5
    ]
    dap_peak <- dap_peak[-overlaps]
  }

  dap_peak <- dap_peak[order(dap_peak$qValue, decreasing = TRUE)] # %>% head(600)
}) %>% `names<-`(names(peakFiles))

# 合并所有 peaks
allPeaks <- BiocGenerics::Reduce("c", peakGR)

# Motif 匹配
WRKY35Match <- matchMotifs(
  WRKY35,
  genomeDna,
  out = "positions",
  p.cutoff = p.cutoff,
  bg = "even"
)

WRKY35Match <- lapply(seq_along(WRKY35Match), function(li) {
  Match <- map2(c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"), 1:5,
    ~ GRanges(
      seqnames = .x,
      ranges   = ranges(WRKY35Match[[li]][[.y]]),
      strand   = WRKY35Match[[li]][[.y]]@elementMetadata$strand
    )
  ) %>%
    { suppressWarnings(do.call(c, .)) } %>%
    `seqlengths<-`(seqlengths(genomeDna)) %>%
    GenomicRanges::trim()
}) %>% `names<-`(names(WRKY35Match))

# 绘图
plotList <- map(c("s1", "m1", "m2"), function(motif_) {
  matches_ <- WRKY35Match[[motif_]]

  # 取各文库与 motif 交集区间，为比较文库间差异需统一区间
  matches_ <- matches_[findOverlaps(matches_, allPeaks)@from %>% unique()]
  matches_ <- resize(matches_, width = 2 * disp_rng + 1, fix = "center")
  matches_ <- matches_[GenomicRanges::trim(matches_) %>% width() %>% {. == 2 * disp_rng + 1}]

  plotDf <- data.frame(
    # DAP = genome_cvgs$DAP[matches_] %>% as.matrix() %>% {./peakLen["DAP"]} %>% colSums(),
    AmpDAP  = genome_cvgs$AmpDAP[matches_] %>% as.matrix() %>% { ./length(matches_) } %>% colSums(),
    mAmpDAP = genome_cvgs$mAmpDAP[matches_] %>% as.matrix() %>% { ./length(matches_) } %>% colSums(),
    Pos     = -disp_rng:disp_rng
  )

  plotDf <- melt(plotDf, "Pos") %>%
    dplyr::rename(Lib = variable, Cov = value)

  ggplot(plotDf) +
    geom_line(aes(x = Pos, y = Cov, color = Lib)) +
    gg_theme_Publication() +
    theme(
      axis.line    = element_blank(),
      panel.border = element_rect(fill = NA, color = "black")
    ) +
    labs(title = motif_)
})

# 合并绘图
p <- cowplot::plot_grid(plotlist = plotList, nrow = 1)
