m_selex_files= read_xlsx("2024_WRKY/WRKY_For_analysis.xlsx",sheet = "M_SELEX")
reads_files=m_selex_files$clean_reads

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
reads_num=fjComm::mclapply(reads_files,function(file){seqs=readLines(file);seqs=fjComm::rmdup(as.data.frame(seqs))[[1]];length(seqs)},mc.cores = 10) %>% unlist()
m_selex_files %<>% mutate(KL_divs_6mer=KL_divs, total_reads=reads_num)

p_selex_QC_stat=ggplot(m_selex_files)+geom_histogram(aes(KL_divs_6mer),binwidth = 0.0005,fill="#00BFC4")+
  geom_vline(xintercept = 0.018)+
  geom_vline(xintercept = 0.03)+
  scale_y_continuous(expand = c(0,0),limits = c(0,7))+
  scale_x_continuous(limits = c(0.011,0.049))

ggsave(p_selex_QC_stat,filename = "2024_WRKY/Me-selex/KL_divs_6mer_QC_1.pdf",width =9.3,height = 2.5 )
save(m_selex_files,file="2024_WRKY/Me-selex/KL_divs_6mer_QC_1.Rdata")


  selex_files= read_xlsx("2024_WRKY/WRKY_For_analysis.xlsx",sheet = "SELEX")
  reads_files=selex_files$clean_reads
  KL_divs=fjComm::mclapply(reads_files,function(file){seqs=readLines(file);SELEX_6mer_ic(seqs)},mc.cores = 20) %>% unlist()
  reads_num=fjComm::mclapply(reads_files,function(file){seqs=readLines(file);seqs=fjComm::rmdup(as.data.frame(seqs))[[1]];length(seqs)},mc.cores = 10) %>% unlist()
  selex_files %<>% mutate(KL_divs_6mer=KL_divs, total_reads=reads_num)
  p_selex_QC_stat=ggplot(selex_files)+geom_histogram(aes(KL_divs_6mer),binwidth = 0.0005,fill="#F8766D")+
    geom_vline(xintercept = 0.018)+
    geom_vline(xintercept = 0.03)+
    scale_y_continuous(expand = c(0,0),limits = c(0,7))+
    scale_x_continuous(limits = c(0.011,0.049))
  ggsave(p_selex_QC_stat,filename = "2024_WRKY/Me-selex/KL_divs_6mer_QC_2.pdf",width =9.3,height = 2.5 )
  save(selex_files,file = "2024_WRKY/Me-selex/KL_divs_6mer_QC_2.Rdata")




  KL_dist_cutoff=0
  load("~/2024_WRKY/Me-selex/KL_divs_6mer_QC_1.Rdata")
  load("~/2024_WRKY/Me-selex/KL_divs_6mer_QC_2.Rdata")
  selex_files <- selex_files[selex_files$KL_divs_6mer>KL_dist_cutoff,]
  m_selex_files <- m_selex_files[m_selex_files$KL_divs_6mer>KL_dist_cutoff,]
  
  name <- intersect(m_selex_files$TF,selex_files$TF)
  m_df <- m_selex_files[m_selex_files$TF%in%name &m_selex_files$KL_divs_6mer>KL_dist_cutoff,] %>% group_by(TF) %>% arrange(desc(KL_divs_6mer)) %>% dplyr::slice_head(n = 1) %>% ungroup()
  m_df$TF <- paste0(m_df$TF,"_m")
  s_df <- selex_files[selex_files$TF%in%name &selex_files$KL_divs_6mer>KL_dist_cutoff,]%>% group_by(TF) %>% arrange(desc(KL_divs_6mer)) %>% dplyr::slice_head(n = 1) %>% ungroup()
  s_df$TF <- paste0(s_df$TF,"_s")
  merge_df <- bind_rows(s_df,m_df)
  
  cl=makeCluster(20)
  registerDoParallel(cl)
  result <-   foreach(i=seq(nrow(merge_df)))%dopar%{
    library(readr)
    library(universalmotif)
    library(fjComm)
    gk <- fjComm::gkmer_enrich(merge_df$clean_reads[i] %>% read_lines(),len = 8)
  }
  stopImplicitCluster()
  stopCluster(cl)
  
  
  thrid <- lapply(result, function(x) x$counts) %>% do.call(cbind, .)
  colnames(thrid) <- merge_df$TF
  
  s_cor <- cor(thrid[,1:(nrow(merge_df)/2)] %>% as.matrix())#,method = "spearman")
  m_cor <- cor(thrid[,(nrow(merge_df)/2+1):(nrow(merge_df))] %>% as.matrix())#,method = "spearman")
  heatmap <- m_cor-s_cor
  heatmap %<>% set_rownames(rownames(heatmap) %>% str_replace("_m","")) %>% set_colnames(colnames(heatmap) %>% str_replace("_m",""))
  plot=ggheat(heatmap,clustering = "both")+
    scale_fill_gradientn(colours = colorspace::diverge_hcl(100,palette = "Tropic"),limits=c(-0.7,0.7),oob=scales::squish,name="Pearson's r\n(mSELEX-SELEX)")+
    ylab("")+xlab("")+gg_theme_Publication(7)+gg_axis_x_label_angle(45)
  gg_save_pdf(plot,8.7,6,filename = "plot/m_induced_diver_pearson")
  
  s_cor <- cor(thrid[,1:(nrow(merge_df)/2)] %>% as.matrix(),method = "spearman")
  m_cor <- cor(thrid[,(nrow(merge_df)/2+1):(nrow(merge_df))] %>% as.matrix(),method = "spearman")
  heatmap <- m_cor-s_cor
  heatmap %<>% set_rownames(rownames(heatmap) %>% str_replace("_m","")) %>% set_colnames(colnames(heatmap) %>% str_replace("_m",""))
  plot=ggheat(heatmap,clustering = "both")+
    scale_fill_gradientn(colours = colorspace::diverge_hcl(100,palette = "Tropic"),limits=c(-0.4,0.4),oob=scales::squish,name="Spearman's rho\n(mSELEX-SELEX)")+
    ylab("")+xlab("")+gg_theme_Publication(7)+gg_axis_x_label_angle(45)
  
  plot


