 WRKY35 <- read_jaspar("~/WRKY35_s1_m1m2.txt")
  names(WRKY35) <- map(WRKY35, ~.x@name)
  WRKY35 <- convert_motifs(WRKY35, "TFBSTools-PWMatrix") %>% do.call(PWMatrixList, .)
  
  genomeDna <- readDNAStringSet("/wrk/data/genome/yz_genome_data/aragenome/ara.fa")[1:5]
  p.cutoff <- 1e-04
  disp_rng <- 500 
  
  peakFiles <- c(DAP = "~/WRKY35_B11_Nana3_3_3__DAP_WRKY_WG_peaks.narrowPeak", 
                 AmpDAP = "~/WRKY35_B10_Nana3_3__AmpDAP_WRKY_WG_50ul_5ul_peaks.narrowPeak", 
                 mAmpDAP = "~/WRKY35_C5_peaks.narrowPeak")
  allPeaks <- BiocGenerics::Reduce("c", peakGR) ## 合并各文库区间
  ## 用于展示区间coverage信息
  bwFiles <- c(DAP = "~/bwfiles/B11_23_WRKY35_WG_5ul_D2.bs10.bw",
               AmpDAP = "~/bwfiles/B10_22_WRKY35_WG_5ul_A1.bs10.bw",
               mAmpDAP = "~/bwfiles/C5_29_WRKY35_mA1.bs10.bw")
  
  peakGR <- lapply(seq_along(peakFiles), function(peak) {
    dap_peak <- rtracklayer::import(peakFiles[peak])
    seqlevels(dap_peak, pruning.mode = "coarse") <- seqlevels(genomeDna)
    seqlengths(dap_peak) <- seqlengths(genomeDna)
    ## remove background
    if (str_detect(names(peakFiles)[peak], "^DAP")) {
      pixHalo_d <-  BiocIO::import("~/callpeak/G3_75_Dap3_4_input1_peaks.narrowPeak")
      pixHalo_d_ov <- findOverlaps(pixHalo_d, dap_peak)
      overlap_width <- pintersect(pixHalo_d[pixHalo_d_ov@from], dap_peak[pixHalo_d_ov@to]) %>% width()
      overlaps <- pixHalo_d_ov@to[(overlap_width / width(pixHalo_d[pixHalo_d_ov@from])) > 0.5 & (overlap_width / width(dap_peak[pixHalo_d_ov@to])) > 0.5]
      dap_peak <- dap_peak[-overlaps]
    } else if (str_detect(names(peakFiles)[peak], "^AmpDAP")) {
      pixHalo_a <-  BiocIO::import("~/callpeak/input_F4_Nana3_2__WRKY_WG_50ul_5ul_DAPAmpDAP_peaks.narrowPeak")
      pixHalo_a_ov <- findOverlaps(pixHalo_a, dap_peak)
      overlap_width <- pintersect(pixHalo_a[pixHalo_a_ov@from], dap_peak[pixHalo_a_ov@to]) %>% width()
      overlaps <- pixHalo_a_ov@to[(overlap_width / width(pixHalo_a[pixHalo_a_ov@from])) > 0.5 & (overlap_width / width(dap_peak[pixHalo_a_ov@to])) > 0.5]
      dap_peak <- dap_peak[-overlaps]
    } else if (str_detect(names(peakFiles)[peak], "^mAmpDAP")) {
      pixHalo_m <-  BiocIO::import("~/callpeak_with_input/Nana3_8_A4_pIXHalo_peaks.narrowPeak")
      pixHalo_m_ov <- findOverlaps(pixHalo_m, dap_peak)
      overlap_width <- pintersect(pixHalo_m[pixHalo_m_ov@from], dap_peak[pixHalo_m_ov@to]) %>% width()
      overlaps <- pixHalo_m_ov@to[(overlap_width / width(pixHalo_m[pixHalo_m_ov@from])) > 0.5 & (overlap_width / width(dap_peak[pixHalo_m_ov@to])) > 0.5]
      dap_peak <- dap_peak[-overlaps]
    }
    dap_peak <- dap_peak[order(dap_peak$qValue, decreasing = T)] # %>% head(600)
  } ) %>% `names<-`(names(peakFiles))
  
  WRKY35Match <- matchMotifs(WRKY35, genomeDna, out = "positions", p.cutoff = p.cutoff, bg = "even")
  WRKY35Match <- lapply(seq_along(WRKY35Match), function(li) {
    Match <- map2(c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"), 1:5, 
                  ~GRanges(seqnames = .x, ranges = ranges(WRKY35Match[[li]][[.y]]), strand = WRKY35Match[[li]][[.y]]@elementMetadata$strand)) %>% 
      {suppressWarnings(do.call(c, .))} %>% `seqlengths<-`(seqlengths(genomeDna)) %>% GenomicRanges::trim()
  } ) %>% `names<-`(names(WRKY35Match))
  
  
  plotList <- map(c("s1", "m1", "m2"), function(motif_) {
    matches_ <- WRKY35Match[[motif_]]
    
    # 取各文库（总）与motif交集区间，为比较文库间差异，需统一各文库区间
    matches_ <- matches_[findOverlaps(matches_, allPeaks)@from %>% unique()] %>% resize(width = 2 * disp_rng + 1, fix = "center")
    matches_ <- matches_[GenomicRanges::trim(matches_) %>% width() %>% {. == 2 * disp_rng + 1}]
    
    plotDf <- data.frame( # DAP = genome_cvgs$DAP[matches_] %>% as.matrix() %>% {./peakLen["DAP"]} %>% colSums(), 
      AmpDAP = genome_cvgs$AmpDAP[matches_] %>% as.matrix() %>% {./length(matches_)} %>% colSums(), 
      mAmpDAP = genome_cvgs$mAmpDAP[matches_] %>% as.matrix() %>% {./length(matches_)} %>% colSums(), Pos = -disp_rng:disp_rng)
    plotDf <- melt(plotDf, "Pos") %>% dplyr::rename(Lib = variable, Cov = value)
    ggplot(plotDf) + 
      geom_line(aes(x = Pos, y = Cov, color = Lib)) + 
      gg_theme_Publication() + 
      theme(axis.line = element_blank(), panel.border = element_rect(fill = NA, color = "black")) + 
      labs(title = motif_)
  })
  
  p <- cowplot::plot_grid(plotlist = plotList, nrow = 1)



