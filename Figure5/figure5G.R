  genomeDna <- readDNAStringSet("/wrk/data/genome/yz_genome_data/aragenome/ara.fa")[1:5]
  p.cutoff <- 1e-05
  disp_rng=1000
  
  KL_dist_cutoff=0.018
  merge_df="placeholder" %>% {
    load("~/KL_divs_6mer_QC_1.Rdata")
    load("~/KL_divs_6mer_QC_2.Rdata")
    selex_files <- selex_files[selex_files$KL_divs_6mer>KL_dist_cutoff,]
    m_selex_files <- m_selex_files[m_selex_files$KL_divs_6mer>KL_dist_cutoff,]
    name <- intersect(m_selex_files$TF,selex_files$TF)
    m_df <- m_selex_files[m_selex_files$TF%in%name &m_selex_files$KL_divs_6mer>KL_dist_cutoff,] %>% group_by(TF) %>% arrange(desc(KL_divs_6mer)) %>% dplyr::slice_head(n = 1) %>% ungroup()
    m_df$TF <- paste0(m_df$TF,"_m")
    s_df <- selex_files[selex_files$TF%in%name &selex_files$KL_divs_6mer>KL_dist_cutoff,]%>% group_by(TF) %>% arrange(desc(KL_divs_6mer)) %>% dplyr::slice_head(n = 1) %>% ungroup()
    s_df$TF <- paste0(s_df$TF,"_s")
    merge_df <- bind_rows(s_df,m_df); merge_df
  }
  
  # seed: AAGTNAAC
  TF1="WRKY35"

  selexfile=merge_df %>% dplyr::filter(TF=="{TF1}_s" %>% glue()) %>% .$clean_reads
  seqs_s=read_lines(selexfile) %>% fjComm::length_adjust(101L)
  pfm_s=fjComm::pfm_from_seed(seqs_or_file = seqs_s,seed1 = "AAGT",gapLen = 1,seed2 = "AAC",two_strands = T,flankLen = 3,all_start_with_specified_gap = T,rmdup_fg = T)

  pfm_s %<>% apply(2, function(x)x/sum(x))
  
  pfm_list=map(1:ncol(pfm_s),function(x){pfm_s[,x]=0.25;pfm_s}) %>% set_names(paste0("pfm",1:14))
  pfm_list= lapply(pfm_list, function(x) create_motif(x %>% set_rownames(c("A","C","G","T"))) %>% convert_motifs("TFBSTools-PWMatrix")) %>% do.call(PWMatrixList, .)
  
  matches_<- matchMotifs(pfm_list, genomeDna, out = "positions", p.cutoff = p.cutoff, bg = "even")
  matches_ <- lapply(seq_along(matches_), function(li) {
    Match <- map2(c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"), 1:5,
                  ~GRanges(seqnames = .x, ranges = ranges(matches_[[li]][[.y]]), strand = matches_[[li]][[.y]]@elementMetadata$strand)) %>%
      {suppressWarnings(do.call(c, .))} %>% #resize(width = 1000, fix = "center") %>%
      `seqlengths<-`(seqlengths(genomeDna)) %>% GenomicRanges::trim()
  } ) %>% `names<-`(names(matches_ ))
  
  
  mDAPfiles=readxl::read_excel("~/WRKY_For_analysis.xlsx",sheet = 5,col_names = T) %>% dplyr::filter(TF==TF1 & `repeat`==1)
  aDAPfiles=readxl::read_excel("~/WRKY_For_analysis.xlsx",sheet = 4,col_names = T) %>% dplyr::filter(TF==TF1 & `repeat`==1)
  
  mPeak=rtracklayer::import(mDAPfiles$callpeak)
  aPeak=rtracklayer::import(aDAPfiles$callpeak)
  allPeaks=c(mPeak,aPeak)
  
  mCvg=rtracklayer::import(mDAPfiles$bwfile,as="Rle")
  aCvg=rtracklayer::import(aDAPfiles$bwfile,as="Rle")
  

  ind=8
  matches_curr=matches_[[ind]]#[overlapsAny(matches_[[ind]], allPeaks)]
  match_seqs=genomeDna[matches_curr]
  letter_curr=match_seqs %>% as.matrix() %>% .[,ind]
  
  ## resize rngs to calc cvg ##
  matches_curr %<>% resize(2*disp_rng+1,"center")
  filter_=GenomicRanges::trim(matches_curr) %>% width() %>% {.==2*disp_rng+1}
  matches_curr= matches_curr[filter_]
  letter_curr=letter_curr[filter_]
  plotdf=aCvg[matches_curr] %>% as.matrix() %>% set_rownames(letter_curr) %>% set_colnames(-disp_rng:disp_rng) %>% melt() %>% set_colnames(qw("letter pos signal"))
  plotdf %<>% group_by(letter,pos) %>% summarise(signal=sum(signal)/n()) %>% ungroup() %>% mutate(letter=factor(letter,levels=qw("A C G T")))
  plot_ampDAP <- ggplot(plotdf)+
    geom_line(aes(pos,signal,color=letter))+ #%>% plotly::ggplotly()
    scale_colour_manual(values = c("#109648","#265C99","#FFA300","#D62839"))+
    gg_theme_Publication()+
    theme(axis.line = element_blank(), panel.border = element_rect(fill = NA, color = "black"), 
          plot.background = element_blank(), panel.background = element_blank())

  ## mAmpDAP
  plotdf_methylampDAP=mCvg[matches_curr] %>% as.matrix() %>% set_rownames(letter_curr) %>% set_colnames(-disp_rng:disp_rng) %>% melt() %>% set_colnames(qw("letter pos signal"))
  plotdf_m %<>% group_by(letter,pos) %>% summarise(signal=sum(signal)/n()) %>% ungroup() %>% mutate(letter=factor(letter,levels=qw("A C G T")))
  plot_m <- ggplot(plotdf_m)+
    geom_line(aes(pos,signal,color=letter))+ #%>% plotly::ggplotly()
    scale_colour_manual(values = c("#109648","#265C99","#FFA300","#D62839"))+
    gg_theme_Publication()+
    theme(axis.line = element_blank(), panel.border = element_rect(fill = NA, color = "black"), 
          plot.background = element_blank(), panel.background = element_blank())

  
  
