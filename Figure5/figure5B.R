fjComm::clear_()

env=environment()
KL_dist_cutoff=0.018

load("/wrk/jiangdingkun/2024_WRKY/Me-selex/KL_divs_6mer_QC_1.Rdata")
load("/wrk/jiangdingkun/2024_WRKY/Me-selex/KL_divs_6mer_QC_2.Rdata")
selex_files <- selex_files[selex_files$KL_divs_6mer>KL_dist_cutoff,]
m_selex_files <- m_selex_files[m_selex_files$KL_divs_6mer>KL_dist_cutoff,]

name <- intersect(m_selex_files$TF,selex_files$TF) 
m_df <- m_selex_files[m_selex_files$TF%in%name &m_selex_files$KL_divs_6mer>KL_dist_cutoff,] %>% group_by(TF) %>% arrange(desc(KL_divs_6mer)) %>% dplyr::slice_head(n = 1) %>% ungroup()
m_df$TF <- paste0(m_df$TF,"_m")
s_df <- selex_files[selex_files$TF%in%name &selex_files$KL_divs_6mer>KL_dist_cutoff,]%>% group_by(TF) %>% arrange(desc(KL_divs_6mer)) %>% dplyr::slice_head(n = 1) %>% ungroup()
s_df$TF <- paste0(s_df$TF,"_s")
merge_df <- bind_rows(s_df,m_df)


TF1="WRKY35"
wrky33 <- merge_df %>% dplyr::filter(TF %in% c(paste0(c(TF1),c("_s")),paste0(c(TF1),c("_m")))) %>% arrange(desc(TF))

  seqs_s=read_lines(wrky33$clean_reads[1]) %>% fjComm::length_adjust(101L)
pfm_s=fjComm::pfm_from_seed(seqs_or_file = seqs_s,seed1 = "AAGT",gapLen = 0,seed2 = "CAAC",two_strands = T,flankLen = 3,all_start_with_specified_gap = T,rmdup_fg = T)
  pfm_s %<>% apply(2, function(x)x/sum(x))
  seqs_m=read_lines(wrky33$clean_reads[2]) %>% fjComm::length_adjust(101L)
pfm_m=fjComm::pfm_from_seed(seqs_or_file = seqs_m,seed1 = "AAGT",gapLen = 0,seed2 = "CAAC",two_strands = T,flankLen = 3,all_start_with_specified_gap = T,rmdup_fg = T)
  pfm_m %<>% apply(2, function(x)x/sum(x))

  pseudo_freq=0.1
A_T_ratio_s=(pfm_s[1,]+pseudo_freq)/(pfm_s[4,]+pseudo_freq)
A_T_ratio_m=(pfm_m[1,]+pseudo_freq)/(pfm_m[4,]+pseudo_freq)
A_T_ratio=rbind(A_T_ratio_m,A_T_ratio_s)  %>% set_rownames(c("mSELEX","SELEX"))
plotdf=data.frame(A_T_ratio) %>% set_colnames(-3:10) %>% rownames_to_column() %>%  pivot_longer(!rowname) %>% set_colnames(c("group","pos","A_T_ratio")) %>% mutate(pos=as.integer(pos),group=factor(group,levels=c("SELEX","mSELEX")))
p_ratio=ggplot(plotdf)+geom_bar(aes(x=pos,y=A_T_ratio,fill=group),stat = "identity",position="dodge",width = 0.4)+
        theme(axis.line.x = element_blank(),axis.ticks.x = element_blank())+scale_x_continuous(breaks = -3:10) +gg_axis_y_noExp() + xlab("") +ylab("A/T ratio")




