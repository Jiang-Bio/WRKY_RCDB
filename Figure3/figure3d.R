 
  genomeDna <- readDNAStringSet("~/aragenome/ara.fa")[1:5]
  p.cutoff <- 1e-04
  disp_rng <- 500
  WRKY_PRE <- create_motif("TACTGCGCTTAGT") %>% convert_motifs("TFBSTools-PWMatrix")
  
  mDAPfiles=readxl::read_excel("~/WRKY_For_analysis.xlsx",sheet = 5,col_names = T) %>% dplyr::filter(`repeat`==1)
  aDAPfiles=readxl::read_excel("~/WRKY_For_analysis.xlsx",sheet = 4,col_names = T) %>% 
    dplyr::filter(`repeat` == 1) %>% .[str_detect(.$TF, "WRKY"),]
  
  matches_<- matchMotifs(WRKY_PRE, genomeDna, out = "positions", p.cutoff = p.cutoff, bg = "even")
  matches_ <- lapply(seq_along(matches_), function(li) {
    Match <- map2(c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"), 1:5,
                  ~GRanges(seqnames = .x, ranges = ranges(matches_[[li]][[.y]]), strand = matches_[[li]][[.y]]@elementMetadata$strand)) %>%
      {suppressWarnings(do.call(c, .))} %>% #resize(width = 1000, fix = "center") %>%
      `seqlengths<-`(seqlengths(genomeDna)) %>% GenomicRanges::trim()
  } ) %>% `names<-`(names(matches_ ))
  
  plot_PRE <- getSeq(genomeDna, matches_[[1]]) %>% create_motif() %>% {ggseqlogo(.@motif) + gg_theme_Publication()}

  
  ## 将 bw 改为 coverage
  aBwCvg <- map(aDAPfiles$bwfile, ~rtracklayer::import(.x) %>% coverage(., weight = .$score)) %>% 
    `names<-`(aDAPfiles$TF)
  mBwCvg <- map(mDAPfiles$bwfile, ~rtracklayer::import(.x) %>% coverage(., weight = .$score)) %>% 
    `names<-`(mDAPfiles$TF)
  
  TFs <- intersect(names(mBwCvg), names(aBwCvg))

  ## 无 y 轴坐标显示
  plotList_nony <- lapply(TFs, function(TF) {
    matchesExt <- matches_[[1]] %>% resize(2 * disp_rng + 1, "center") %>% .[width(.) == 2 * disp_rng + 1]
    Df <- data.frame(mAmp = mBwCvg[[TF]][matchesExt] %>% as.matrix() %>% colSums() %>% {./sum(.)}, 
                     amp = aBwCvg[[TF]][matchesExt] %>% as.matrix() %>% colSums() %>% {./sum(.)}, 
                     Pos = -disp_rng:disp_rng) %>% melt("Pos")
    plot <- ggplot(Df, aes(Pos, value, color = variable)) + 
      geom_line() + 
      scale_colour_manual(values = c("#5475A5", "#B92F74")) + 
      # geom_smooth(se = F, method = "loess", span = 0.005, linewidth = .5) +
      gg_theme_Publication() + 
      theme(axis.line = element_blank(), panel.border = element_rect(fill = NA, color = "black"), 
            plot.background = element_blank(), panel.background = element_blank(), 
            axis.text.y = element_blank()) + 
      labs(title = TF)
    plot
  })
