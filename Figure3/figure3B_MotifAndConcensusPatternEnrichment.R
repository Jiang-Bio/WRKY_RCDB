#################################################################################
#Figure3B Calculate the enrichment of 11 motifs in SELEX libraries
#################################################################################
    rename_duplicates <- function(strings) {
      counts <- table(strings)
      seen <- setNames(rep(1, length(counts)), names(counts))
      new_strings <- character(length(strings))
      
      for (i in seq_along(strings)) {
        current_string <- strings[i]
        if (counts[current_string] > 1) {
          seen_count <- seen[current_string]
          new_strings[i] <- paste(current_string, seen_count, sep = "_")
          seen[current_string] <- seen_count + 1
        } else {
          new_strings[i] <- paste(current_string, 1, sep = "_")
        }
      }
      
      return(new_strings)
    }

   #load matrix
    pwm <- read_jaspar("~/class11_pwm.txt") %>% convert_motifs("TFBSTools-PWMatrix") %>% do.call(PWMatrixList,.)
    pwm <- pwm[1:11]
    names(pwm@listData) <- sapply(1:11, function(x){pwm@listData[[x]]@name})
    
   #load file
   selex_files <- list.files("~/SELEX/",pattern = "gz$") %>% as.data.frame() %>% `colnames<-`(c("for_motif_disc"))
   selex_files$TF <- selex_files$for_motif_disc %>% basename() %>% str_extract("WRKY[0-9]+") %>% paste0("_s")
   m_selex_files <- list.files("~/SELEX/",pattern = "gz$") %>% as.data.frame() %>% `colnames<-`(c("for_motif_disc"))
   m_selex_files$TF <- m_selex_files$for_motif_disc %>% basename() %>% str_extract("WRKY[0-9]+") %>% paste0("_m")
   selex_files <- bind_rows(selex_files,m_selex_files)
   selex_files$repeats <- rename_duplicates(selex_files$TF)
  
  
  #Calculate the enrichment of motifs in each SELEX library
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
  #Normalize to visualization
  df <- apply(df, 1, function(x){log2(x+1)/max(log2(x+1))}) %>% t()

  {
    df1 <- df[1:99,]
    df2 <- df[100:139,]
    SELEX <- colMeans(df1)-2*colMeans(df1)
    MESELEX <- colMeans(df2)
    bar <- rbind(SELEX,MESELEX) %>% melt()
    bar <- ggplot(bar, aes(x = Var2, y = value, fill = Var1)) + 
      geom_col(position = "identity") +  
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


  #################################################################################
  #Calculate the enrichment of the WRKY consensus sequence in each SELEX library
  #################################################################################

  {
   #load file
   selex_files <- list.files("~/SELEX/",pattern = "gz$") %>% as.data.frame() %>% `colnames<-`(c("for_motif_disc"))
   selex_files$TF <- selex_files$for_motif_disc %>% basename() %>% str_extract("WRKY[0-9]+") %>% paste0("_s")
   m_selex_files <- list.files("~/SELEX/",pattern = "gz$") %>% as.data.frame() %>% `colnames<-`(c("for_motif_disc"))
   m_selex_files$TF <- m_selex_files$for_motif_disc %>% basename() %>% str_extract("WRKY[0-9]+") %>% paste0("_m")
   selex_files <- bind_rows(selex_files,m_selex_files)
   selex_files$repeats <- rename_duplicates(selex_files$TF)

    subject <- c("RGTCAR","AAAGTC","TTTTCCAC","AAGTTTTC","ATCGGTAGCACGA","TACTGCGCTTAGT","TAAAGATTACTAATAGGAA") %>% DNAStringSet()
    subject@ranges@NAMES <- c("W-box","WK-box","WT-box","WRKY18","N-control","PRE","SURE")
    
    cl=makeCluster(20)
    registerDoParallel(cl)
    W_BOX_RICH <- foreach(i=seq(nrow(selex_files)))%dopar%{
      pseudo = 200 
      library(motifmatchr)
      library(fjComm)
      library(Biostrings)
      
      seqs <- fjComm::getSeq_fqfachrFile(selex_files$for_motif_disc[i])
      seqs <- c(seqs) %>% DNAStringSet()
      shuffled_seqs <- seqs %>% shuffle_sequences()
      a1 <-  vcountPDict(subject[1:4],subject =seqs,fixed = F,max.mismatch = 0)%>% rowSums() 
      a2 <-  vcountPDict(subject[1:4],subject =shuffled_seqs,fixed =F,max.mismatch = 0)  %>% rowSums()  
      a <- (a1+pseudo)/(a2+pseudo)
      b1 <-  vcountPDict(subject[5:7],subject =seqs,fixed = F,max.mismatch = 3)%>% rowSums()
      b2 <-  vcountPDict(subject[5:7],subject =shuffled_seqs,fixed = F,max.mismatch = 3) %>% rowSums() 
      b <- (b1+pseudo)/(b2+pseudo)
      
      counts <- c(a,b)
      names(counts) <- selex_files$TF[i]
      counts
    }
    stopImplicitCluster()
    stopCluster(cl)
    
    names(W_BOX_RICH) <- paste0(selex_files$TF,"_",selex_files$repeats)
    df <- W_BOX_RICH %>% largeListToDf() %>% t()
    df <- apply(df, 1, function(x){log2(x)/max(log2(x))}) %>% t()
    df <- apply(df, 2, function(x) ifelse(x < 0, 0, x))
    colnames(df) <-names(subject) 
    
    df1 <- df[1:99,]
    df2 <- df[100:139,]
    df2 <- df2[order(rownames(df2) ),]
    SELEX <- colMeans(df1)-2*colMeans(df1)
    MESELEX <- colMeans(df2)
    bar <- rbind(SELEX,MESELEX) %>% melt()
    bar <- ggplot(bar, aes(x = Var2, y = value, fill = Var1)) + 
      geom_col(position = "identity") +  
      coord_flip() + 
      scale_y_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1, 1)) +
      scale_fill_manual(values = c("SELEX" = "#F8766D", "MESELEX" = "#00BFC4")) +
      theme_bw()
    bar
    
    df1 <- melt(as.matrix(df1))
    p1 <- ggplot(df1,aes(x=Var2,y=Var1,fill=value))+geom_tile()+
      scale_fill_gradientn(colours = brewer.pal(n = 9,name = "PuBu"),name="lgo2x/log2max",limits = c(0, 1))+
      ylab("")+xlab("")+
      theme_bw()+theme(panel.grid=element_blank()) 
    p1
    
    df2 <- melt(as.matrix(df2))
    p2 <- ggplot(df2,aes(x=Var2,y=Var1,fill=value))+geom_tile()+
      scale_fill_gradientn(colours = brewer.pal(n = 9,name = "PuBu"),name="lgo2x/log2max",limits = c(0, 1))+
      ylab("")+xlab("")+
      theme_bw()+theme(panel.grid=element_blank()) 
    p2 
    
    plot <- plot_grid(p1,p2,bar,ncol = 3)
