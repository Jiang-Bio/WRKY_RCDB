# =============================================
# figure 5A
# =============================================







# =============================================
# figure 5B
# =============================================

# Load required packages
pacman::p_load(universalmotif, fjComm, ggplot2, cowplot, reshape2, dplyr, Biostrings, GenomicRanges, rtracklayer)

# ------------------------
# 1. Load motif and convert to PWMatrix
# ------------------------
WRKY35 <- read_jaspar("/wrk/chenhao/work/WRKY_paper/genome_coverage_for_wrky35/WRKY35_s1_m12.txt")
names(WRKY35) <- map(WRKY35, ~.x@name)
WRKY35 <- convert_motifs(WRKY35, "TFBSTools-PWMatrix") %>% do.call(PWMatrixList, .)

# ------------------------
# 2. Load genome and set parameters
# ------------------------
genomeDna <- readDNAStringSet("/wrk/data/genome/yz_genome_data/aragenome/ara.fa")[1:5]
p.cutoff <- 1e-4
disp_rng <- 500  # window around motif for plotting

# ------------------------
# 3. Define peak files and bigWig files
# ------------------------
peakFiles <- c(
  DAP = "~/2024_WRKY/5_CALLPEAK/DAP_callpeak_minus_pixHalo/callpeak/WRKY35_B11_Nana3_3_3__DAP_WRKY_WG_peaks.narrowPeak", 
  AmpDAP = "~/2024_WRKY/5_CALLPEAK/AMPdap_callpeak_minus_pixHalo/callpeak/WRKY35_B10_Nana3_3__AmpDAP_WRKY_WG_50ul_5ul_peaks.narrowPeak", 
  mAmpDAP = "~/2024_WRKY/6_Me-DAP/callpeak_with_input/WRKY35_C5_peaks.narrowPeak"
)

bwFiles <- c(
  DAP = "/wrk/manana/work/Nana3_3_3__DAP_WRKY_WG/pre_processed/bwfiles/B11_23_WRKY35_WG_5ul_D2.bs10.bw",
  AmpDAP = "/wrk/manana/work/Nana3_3__AmpDAP_WRKY_WG_50ul_5ul/pre_processed/bwfiles/B10_22_WRKY35_WG_5ul_A1.bs10.bw",
  mAmpDAP = "/wrk/manana/work/Nana3_8__WRKY_methyl_ampDAP/pre_processed/bwfiles/C5_29_WRKY35_mA1.bs10.bw"
)

# ------------------------
# 4. Import peaks and remove "pixHalo" background
# ------------------------
peakGR <- lapply(seq_along(peakFiles), function(peak) {
  dap_peak <- rtracklayer::import(peakFiles[peak])
  seqlevels(dap_peak, pruning.mode = "coarse") <- seqlevels(genomeDna)
  seqlengths(dap_peak) <- seqlengths(genomeDna)

  # Remove background peaks based on library type
  lib_name <- names(peakFiles)[peak]
  if (str_detect(lib_name, "^DAP")) {
    pixHalo <- BiocIO::import("~/2024_WRKY/5_CALLPEAK/DAP_callpeak_minus_pixHalo/callpeak/G3_75_Dap3_4_input1_peaks.narrowPeak")
  } else if (str_detect(lib_name, "^AmpDAP")) {
    pixHalo <- BiocIO::import("~/2024_WRKY/5_CALLPEAK/AMPdap_callpeak_minus_pixHalo/callpeak/input_F4_Nana3_2__WRKY_WG_50ul_5ul_DAPAmpDAP_peaks.narrowPeak")
  } else if (str_detect(lib_name, "^mAmpDAP")) {
    pixHalo <- BiocIO::import("~/2024_WRKY/6_Me-DAP/callpeak_with_input/Nana3_8_A4_pIXHalo_peaks.narrowPeak")
  }

  pixHalo_ov <- findOverlaps(pixHalo, dap_peak)
  overlap_width <- pintersect(pixHalo[pixHalo_ov@from], dap_peak[pixHalo_ov@to]) %>% width()
  overlaps <- pixHalo_ov@to[
    (overlap_width / width(pixHalo[pixHalo_ov@from])) > 0.5 &
      (overlap_width / width(dap_peak[pixHalo_ov@to])) > 0.5
  ]
  dap_peak <- dap_peak[-overlaps]

  # Sort by qValue
  dap_peak[order(dap_peak$qValue, decreasing = TRUE)]
}) %>% `names<-`(names(peakFiles))

# ------------------------
# 5. Map motifs to genome
# ------------------------
WRKY35Match <- matchMotifs(WRKY35, genomeDna, out = "positions", p.cutoff = p.cutoff, bg = "even")

WRKY35Match <- lapply(seq_along(WRKY35Match), function(li) {
  map2(c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"), 1:5, 
       ~GRanges(seqnames = .x, ranges = ranges(WRKY35Match[[li]][[.y]]),
                strand = WRKY35Match[[li]][[.y]]@elementMetadata$strand)) %>%
    {suppressWarnings(do.call(c, .))} %>%
    `seqlengths<-`(seqlengths(genomeDna)) %>%
    GenomicRanges::trim()
}) %>% `names<-`(names(WRKY35Match))

# ------------------------
# 6. Generate coverage plots
# ------------------------
plotList <- map(c("s1", "m1", "m2"), function(motif_) {
  matches_ <- WRKY35Match[[motif_]]

  # Intersect with allPeaks for uniform comparison
  matches_ <- matches_[findOverlaps(matches_, allPeaks)@from %>% unique()]
  matches_ <- resize(matches_, width = 2 * disp_rng + 1, fix = "center")
  matches_ <- matches_[GenomicRanges::trim(matches_) %>% width() == 2 * disp_rng + 1]

  plotDf <- data.frame(
    AmpDAP = genome_cvgs$AmpDAP[matches_] %>% as.matrix() %>% {./length(matches_)} %>% colSums(),
    mAmpDAP = genome_cvgs$mAmpDAP[matches_] %>% as.matrix() %>% {./length(matches_)} %>% colSums(),
    Pos = -disp_rng:disp_rng
  ) %>% melt("Pos") %>% dplyr::rename(Lib = variable, Cov = value)

  ggplot(plotDf) + 
    geom_line(aes(x = Pos, y = Cov, color = Lib)) +
    gg_theme_Publication() +
    theme(axis.line = element_blank(),
          panel.border = element_rect(fill = NA, color = "black")) +
    labs(title = motif_)
})

# Combine plots into a single figure
p <- cowplot::plot_grid(plotlist = plotList, nrow = 1)
gg_save_pdf(p, 18, 3, "/wrk/chenhao/work/WRKY_paper/genome_coverage_for_wrky35/",
            paste0("WRKY35_p", p.cutoff))





