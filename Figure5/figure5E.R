



    KL_dist_cutoff=0.012
    wb <- loadWorkbook("~/WRKY_For_analysis.xlsx")
    selex_files  <- read.xlsx(wb, sheet = "SELEX") %>% as.tibble()
    m_selex_files<- read.xlsx(wb, sheet = "M_SELEX") %>% as.tibble()
    
    selex_files <- selex_files[selex_files$KL_divs_6mer>KL_dist_cutoff,]
    m_selex_files <- m_selex_files[m_selex_files$KL_divs_6mer>KL_dist_cutoff,]
    
    name <- intersect(m_selex_files$TF,selex_files$TF) #%>% setdiff(too_weak) #除去信号太弱的

    m_df <- m_selex_files[m_selex_files$TF%in%name ,] %>% group_by(TF) %>% arrange(desc(KL_divs_6mer)) %>% dplyr::slice_head(n = 1) %>% ungroup()
    m_df$TF <- paste0(m_df$TF,"_m")
   
    s_df <- selex_files[selex_files$TF%in%name ,]%>% group_by(TF) %>% arrange(desc(KL_divs_6mer)) %>% dplyr::slice_head(n = 1) %>% ungroup()
    s_df$TF <- paste0(s_df$TF,"_s")
    merge_df <- bind_rows(s_df,m_df)
    
    
    TF1="WRKY29"
    list <- mclapply(c(1:22), function(x){
      TF1 <- merge_df$TF[x] %>% gsub("_s","",.)
      
      wrky33 <- merge_df %>% dplyr::filter(TF %in% c(paste0(c(TF1),c("_s")),paste0(c(TF1),c("_m")))) %>% arrange(desc(TF))
      
      seqs_s=read_lines(wrky33$clean_reads[1]) %>% fjComm::length_adjust(101L)
      pfm_s=fjComm::pfm_from_seed(seqs_or_file = seqs_s,seed1 = "AAGT",gapLen = 0,seed2 = "CAAC",two_strands = T,flankLen = 3,all_start_with_specified_gap = T,rmdup_fg = T)
      pfm_s %<>% apply(2, function(x)x/sum(x))
      seqs_m=read_lines(wrky33$clean_reads[2]) %>% fjComm::length_adjust(101L)
      pfm_m=fjComm::pfm_from_seed(seqs_or_file = seqs_m,seed1 = "AAGT",gapLen = 0,seed2 = "CAAC",two_strands = T,flankLen = 3,all_start_with_specified_gap = T,rmdup_fg = T)
      pfm_m %<>% apply(2, function(x)x/sum(x))
      pseudo=0.1
      C_strand_plus=log2( (pfm_m[2,]+pseudo)/   (pfm_s[2,]+pseudo) )
      C_strand_minus=log2((pfm_m[3,]+pseudo)/  (pfm_s[3,]+pseudo) )
      C_effects= c(C_strand_plus,C_strand_minus) %>% as.data.frame() %>% t() %>% set_rownames(TF1) %>% melt()
      abs_max <- c(abs(pfm_m[2,]-pfm_s[2,]),abs(pfm_m[3,]-pfm_s[3,])) %>% as.data.frame() %>% t() %>% set_rownames(c("size"))
      df <- cbind(C_effects,t(abs_max))
      
      library(cowplot)
      blank_=fjComm::createDummy() %>% ggplotGrob()
      pfmlogo=ggplotGrob(ggseqlogo_lab(pfm_s)+theme_void())
      prow1=plot_grid(pfmlogo,blank_,nrow = 1,rel_widths = c(0.67,0.33))
      pfmlogo=ggplotGrob(ggseqlogo_lab(pfm_m)+theme_void())
      prow2=plot_grid(pfmlogo,blank_,nrow = 1,rel_widths = c(0.67,0.33))
      plot_grid(prow1,prow2,ncol = 1)
      
      plot1 <- plot_grid(prow1,prow2,blank_,ncol = 1,rel_heights = c(0.60,0.35,0.05),greedy = F)
      ggsave(filename = glue("~/2024_WRKY/11_diff_bind_all_WRKY_amp_dap_methy/plot/LOGO_{TF1}_new_version.pdf"),plot1,width = 3.26,height = 1.34)
      df
    },mc.cores = 22)
    
    
    df <- do.call(rbind,list)
    ord <- df$Var1[df$Var2==11][df$value[df$Var2==11] %>% order()]
    df$Var1 <- factor(df$Var1, levels = c(ord %>% as.character()))
    df <- df[which(df$Var1 %in% name),]
    
    df <- df %>%
      mutate(
        group = case_when(
          abs(value) < 0.2 ~ "medium",
          value >0.2  ~ "increase",
          TRUE ~ "decrease"
        )
      )
    
    df$Var1[df$Var2==11& df$group=="medium"] %>% as.character()
    plot <- ggplot(df, aes(x = Var2, y = Var1)) +
      geom_tile(color = "white", fill = "white") +
      geom_point(
        aes(size = size + 2, fill = value),  
        alpha = 0.99,
        shape = 21,         
        color = "black",       
        stroke = 0.5    
      ) +
      scale_fill_gradientn(   
        colors = diverge_hcl(n = 100, palette = "Blue-Red3"),
        name = "Log2 fold change\n(mSELEX / SELEX)", 
        limits = c(-2.5, 2.5)
      ) +
      coord_equal() +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme_minimal() +labs(subtitle = "KL_dist_cutoff=0.012  abs(value) < 0.3 ~ medium")+
      theme(
        axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.ticks  = element_blank(),
        axis.text.y  = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 10),
        axis.title  = element_blank(),
        legend.title  = element_blank(),
        legend.key  = element_blank(),
        legend.text  = element_text(color = "black", size = 9),
        legend.spacing.x  = unit(0.1, 'cm'),
        legend.key.width  = unit(0.5, 'cm'),
        legend.key.height  = unit(0.5, 'cm'),
        legend.background  = element_blank()
      )


  
