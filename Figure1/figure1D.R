KL_dist_cutoff=0 #0.018
  TF1="WRKY71";TF2="WRKY65"
  
  load("~/2024_WRKY/Me-selex/KL_divs_6mer_QC_1.Rdata")
  m_selex_files <- selex_files
  load("~/2024_WRKY/Me-selex/KL_divs_6mer_QC_2.Rdata")
  selex_files <- selex_files[selex_files$KL_divs_6mer>KL_dist_cutoff,]
  m_selex_files <- m_selex_files[m_selex_files$KL_divs_6mer>KL_dist_cutoff,]
  
  name <- intersect(m_selex_files$TF,selex_files$TF) %>% setdiff(too_weak) 
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

  plot_s=ggplot(dis_merge,aes(y=(dis_merge[[2]]), x=(dis_merge[[1]])))+geom_hex(bins=70)+
    scale_fill_gradientn(colours = colorspace::sequential_hcl(50,palette = "Mako"),name="Point density",values = c(0,0.01,1))+
    scale_x_continuous(expand = c(0,0.1))+scale_y_continuous(expand = c(0,0.1))+xlab(paste0("8-mer enrichment (",colnames(dis_merge)[1]%>% str_replace("_.",""),")"))+ylab(paste0("8-mer enrichment (",colnames(dis_merge)[2]%>% str_replace("_.",""),""))+gg_theme_Publication(7)
  #
  plot_m=ggplot(dis_merge,aes(y=(dis_merge[[4]]), x=(dis_merge[[3]])))+geom_hex(bins=70)+
    scale_fill_gradientn(colours = colorspace::sequential_hcl(50,palette = "Mako"),name="Point density",values = c(0,0.01,1))+
    scale_x_continuous(expand = c(0,0.1))+scale_y_continuous(expand = c(0,0.1))+xlab(paste0("8-mer enrichment (",colnames(dis_merge)[3]%>% str_replace("_.",""),")"))+ylab(paste0("8-mer enrichment (",colnames(dis_merge)[2]%>% str_replace("_.",""),""))+gg_theme_Publication(7)
