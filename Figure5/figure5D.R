fjComm::clear_()
library(RColorBrewer)

TF=c("WRKY14","WRKY65","WRKY29")
tt=mclapply(TF, function(TF){
env=environment()
KL_dist_cutoff=0.018
load("~/KL_divs_6mer_QC_1.Rdata")
load("~/KL_divs_6mer_QC_1.Rdata")
selex_files <- selex_files[selex_files$KL_divs_6mer>KL_dist_cutoff,]
m_selex_files <- m_selex_files[m_selex_files$KL_divs_6mer>KL_dist_cutoff,]

name <- intersect(m_selex_files$TF,selex_files$TF) #%>% setdiff(too_weak) #除去信号太弱的
m_df <- m_selex_files[m_selex_files$TF%in%name &m_selex_files$KL_divs_6mer>KL_dist_cutoff,] %>% group_by(TF) %>% arrange(desc(KL_divs_6mer)) %>% dplyr::slice_head(n = 1) %>% ungroup()
m_df$TF <- paste0(m_df$TF,"_m")
s_df <- selex_files[selex_files$TF%in%name &selex_files$KL_divs_6mer>KL_dist_cutoff,]%>% group_by(TF) %>% arrange(desc(KL_divs_6mer)) %>% dplyr::slice_head(n = 1) %>% ungroup()
s_df$TF <- paste0(s_df$TF,"_s")
merge_df <- bind_rows(s_df,m_df)

Selex=merge_df %>% dplyr::filter(TF==paste0(env$TF,"_s")) %>% .$clean_reads
Me_selex=merge_df %>% dplyr::filter(TF==paste0(env$TF,"_m")) %>% .$clean_reads


Me_selex <-  read_lines(Me_selex)
Selex    <-  read_lines(Selex)

klen=8
p1=fjComm::kmerCntBit(Selex, k = klen, diffLen = T, collapse = T, asDf = T, all_possible_k = T, pseudo = 0)
p1_shuffle=fjComm::kmerCntBit(universalmotif::shuffle_sequences(Selex %>% DNAStringSet()) %>% as.character(), k = klen, diffLen = T, collapse = T, asDf = T, all_possible_k = T, pseudo = 0)
p1$counts=p1$counts/p1_shuffle$counts
p2=fjComm::kmerCntBit(Me_selex, k = klen, diffLen = T, collapse = T, asDf = T, all_possible_k = T, pseudo = 0)
p2_shuffle=fjComm::kmerCntBit(universalmotif::shuffle_sequences(Me_selex %>% DNAStringSet()) %>% as.character(), k = klen, diffLen = T, collapse = T, asDf = T, all_possible_k = T, pseudo = 0)
p2$counts=p2$counts/p2_shuffle$counts

merge=dplyr::left_join(p1, p2, by = "kmer")

merge <- merge %>%
  mutate(
    group = case_when(

      grepl("GTTGACT|AGTCAAC",merge$kmer) ~ "GTCAAC",
      grepl("GTTAACT|AGTTAAC",merge$kmer) ~ "GTTAAC",
      grepl("TTTGACT|AGTCAAA",merge$kmer) ~ "GTCAAA",
      grepl("TTTTACT|AGTAAAA",merge$kmer) ~ "GTAAAA",
      grepl("GTTTACT|AGTAAAC",merge$kmer) ~ "GTAAAC",
      !grepl("C|G",merge$kmer) ~ "NoC"

    )
  )
merge$group[which(is.na(merge$group))]="other"

df1=merge[which(merge$group=="GTCAAC"),]
df2=merge[which(merge$group=="GTAAAC"),]
df3=merge[which(merge$group=="GTCAAA"),]
df4=merge[which(merge$group=="GTAAAA"),]
df5=merge[which(merge$group=="NoC"),]
df6=merge[which(merge$group=="GTAAAC"),]
df7=merge[which(merge$group=="other"),]

p <- ggplot(merge)+
  geom_smooth(data = df1, mapping = aes(x=counts.x,y=counts.y), formula = y ~ x+0,method = "lm",se = FALSE,fullrange = T,color = "#377EB8")+
  geom_smooth(data = df3, mapping = aes(x=counts.x,y=counts.y), formula = y ~ x+0,method = "lm",se = FALSE,fullrange = T,color = "#4DAF4A")+
  geom_point(data = df3,mapping =   aes(x=counts.x,y=counts.y), color="#4DAF4A")+
  scale_color_manual(values  = brewer.pal(6,name = "Set1"))+
  gg_theme_Publication(base_size = 7)+gg_axis_x_noExp()+gg_axis_y_noExp()+
  xlab("8-mer enrichment (SELEX)")+
  ylab("8-mer enrichment (mSELEX)")+

  labs(title = TF)

},mc.cores = 3)


