
    {
    pwm <- read_jaspar("~/class11_pwm.txt") %>% convert_motifs("TFBSTools-PWMatrix") %>% do.call(PWMatrixList,.)
    pwm <- pwm[1:11]
    names(pwm@listData) <- sapply(1:11, function(x){pwm@listData[[x]]@name})
    
    selex_files <- read.xlsx("~/WRKY_For_analysis.xlsx",sheet = "SELEX")
    m_selex_files <- read.xlsx("~/WRKY_For_analysis.xlsx",sheet = "M_SELEX")
    m_selex_files$TF <- paste0(m_selex_files$TF,"_m")
    selex_files$TF <- paste0(selex_files$TF,"_s")
    selex_files <- bind_rows(selex_files,m_selex_files)
  }
  
  
  
  cl=makeCluster(30)
  registerDoParallel(cl)
  rich_in_101N <- foreach(i=seq(nrow(selex_files)))%dopar%{
    pseudo = 1000
    library(magrittr)
    library(stringr)
    library(motifmatchr)
    library(fjComm)
    library(Biostrings)
    seqs <- fjComm::getSeq_fqfachrFile(selex_files$for_motif_disc[i])
    seqlens <- seqs %>% nchar %>% table() %>% which.max() %>% names %>% as.integer()
    seqs <- fjComm::length_adjust(seqs, seqlens)
    seqs <- c(seqs, revComp(seqs)) %>% DNAStringSet()
    a <-  matchMotifs(pwm,seqs,out="matches",p.cutoff = 1e-05) 
    a <- a@assays@data@listData$motifMatches  %>% colSums()
    b <-  matchMotifs(pwm,seqs %>% shuffle_sequences(),out="matches",p.cutoff = 1e-05)
    b <- b@assays@data@listData$motifMatches  %>% colSums()
    counts <- (a+pseudo) / (b+pseudo)
    names(counts) <- names(pwm)
    counts
  }
  stopImplicitCluster()
  stopCluster(cl)
  names(rich_in_101N) <- paste0(selex_files$TF,"_",selex_files$repeats)

  df <- rich_in_101N %>% largeListToDf() %>% t() %>% `colnames<-`(c(names(pwm))[1:11])
  df <- df[c(1:99,df[100:139,] %>% rownames() %>% order()+99),]
  df <- apply(df, 1, function(x){log2(x+1)/max(log2(x+1))}) %>% t()
  {
    df1 <- df[1:99,]
    df2 <- df[100:139,]
    SELEX <- colMeans(df1)-2*colMeans(df1)
    MESELEX <- colMeans(df2)
    bar <- rbind(SELEX,MESELEX) %>% melt()
    bar <- ggplot(bar, aes(x = Var2, y = value, fill = Var1)) + 
      geom_col(position = "identity") +  # 不堆叠，直接按数值绘制
      coord_flip() + 
      scale_y_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1, 1)) +
      scale_fill_manual(values = c("SELEX" = "#F8766D", "MESELEX" = "#00BFC4")) +
      theme_bw()
    
    df1 <- melt(as.matrix(df1))
    p1 <- ggplot(df1,aes(x=Var2,y=Var1,fill=value))+geom_tile()+
      scale_fill_gradientn(colours = brewer.pal(n = 9,name = "PuBu"),name="lgo2x/log2max",limits = c(0, 1))+
      ylab("")+xlab("")+
      theme_bw()+theme(panel.grid=element_blank()) #
    p1
    
    df2 <- melt(as.matrix(df2))
    p2 <- ggplot(df2,aes(x=Var2,y=Var1,fill=value))+geom_tile()+
      scale_fill_gradientn(colours = brewer.pal(n = 9,name = "PuBu"),name="lgo2x/log2max",limits = c(0, 1))+
      ylab("")+xlab("")+
      theme_bw()+theme(panel.grid=element_blank()) #
    p2 
    
    plot <- plot_grid(p1,p2,bar,ncol = 3)

    
  }
  
  

  
  {
    
    selex_files <- read.xlsx("~/WRKY_For_analysis.xlsx",sheet = "SELEX")
    m_selex_files <- read.xlsx("~/WRKY_For_analysis.xlsx",sheet = "M_SELEX")
    m_selex_files$TF <- paste0(m_selex_files$TF,"_m")
    selex_files$TF <- paste0(selex_files$TF,"_s")
    selex_files <- bind_rows(selex_files,m_selex_files)

    subject <- c("RGTCAR","AAAGTC","TTTTCCAC","AAGTTTTC","ATCGGTAGCACGA","TACTGCGCTTAGT","TAAAGATTACTAATAGGAA") %>% DNAStringSet()
    subject@ranges@NAMES <- c("W-box","WK-box","WT-box","WRKY18","N-control","PRE","SURE")
    
    cl=makeCluster(20)
    registerDoParallel(cl)
    W_BOX_RICH <- foreach(i=seq(nrow(selex_files)))%dopar%{
      pseudo = 200 #200
      library(motifmatchr)
      library(fjComm)
      library(Biostrings)
      
      seqs <- fjComm::getSeq_fqfachrFile(selex_files$for_motif_disc[i])
      seqs <- c(seqs) %>% DNAStringSet()
      shuffled_seqs <- seqs %>% shuffle_sequences()
      a1 <-  vcountPDict(subject[1:4],subject =seqs,fixed = F,max.mismatch = 0)%>% rowSums()  # %>% apply( 1, function(x){ length(which(x>0)) })#%>% rowSums()  
      a2 <-  vcountPDict(subject[1:4],subject =shuffled_seqs,fixed =F,max.mismatch = 0)  %>% rowSums()  #%>% apply( 1, function(x){ length(which(x>0)) })# %>% rowSums()
      a <- (a1+pseudo)/(a2+pseudo)
      b1 <-  vcountPDict(subject[5:7],subject =seqs,fixed = F,max.mismatch = 3)%>% rowSums()  # %>% apply( 1, function(x){ length(which(x>0)) })#%>% rowSums()
      b2 <-  vcountPDict(subject[5:7],subject =shuffled_seqs,fixed = F,max.mismatch = 3) %>% rowSums()  #%>% apply( 1, function(x){ length(which(x>0)) })#%>% rowSums()
      b <- (b1+pseudo)/(b2+pseudo)
      
      counts <- c(a,b)
      names(counts) <- selex_files$TF[i]
      counts
    }
    stopImplicitCluster()
    stopCluster(cl)


    
    names(W_BOX_RICH) <- paste0(selex_files$TF,"_",selex_files$repeats)
    df <- W_BOX_RICH %>% largeListToDf() %>% t() #%>% `rownames<-`(selex_files$TF %>% rename_duplicates()) %>% `colnames<-`(subject@ranges@NAMES)
    df <- apply(df, 1, function(x){log2(x)/max(log2(x))}) %>% t()
    df <- apply(df, 2, function(x) ifelse(x < 0, 0, x))
    colnames(df) <-names(subject) 
    
