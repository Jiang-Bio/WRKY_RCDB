#####################################################################
##3G RNA####
#####################################################################
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
  
    #load gene raw count file
    df <- read.table("~/diff_tissue_RNASeq.count",header = T)
    count_table<- df[,c(7:16)]
    colnames(count_table) <- c("Cotyledon","Cotyledon","Flower","Flower","Silique","Silique","Stem","Stem","Root","Root" ) %>% rename_duplicates()
    #TPM calculation
    rnaseq_tpm <-  t(t(count_table)/colSums(count_table))*10^6
    rnaseq_tpm <- cbind(df[,1:6],rnaseq_tpm)
    rnaseq_tpm$Geneid <- gsub(".Araport11.447",replacement = "",rnaseq_tpm$Geneid)
    
    
    library(parallel)
    library(parallelly)
    library(foreach)
    library(doParallel)
    cl=makeCluster(20)
    registerDoParallel(cl)
    result <- list()
    result <- foreach(x=seq(1:10))%dopar%{
      library(SummarizedExperiment)
      library(universalmotif)
      library(dplyr)
      library(Biostrings)
      library(motifmatchr)
      library(fjComm)
      pseudo=1000
      top=0.10
      bottom=0.10
      upstream = 300
      downstream = 300
      up_gene <- rnaseq_tpm[order(rnaseq_tpm[,x+6],decreasing = T),1][1:round(nrow(rnaseq_tpm)*top)]
      down_gene <- rnaseq_tpm[order(rnaseq_tpm[,x+6],decreasing = F),1][1:round(nrow(rnaseq_tpm)*bottom)]
      library(TxDb.Athaliana.BioMart.plantsmart51)
      promoter <- promoters(TxDb.Athaliana.BioMart.plantsmart51,upstream = upstream,downstream = downstream) %>% GenomicRanges::trim() %>% {.[width(.)==downstream+upstream]}
      promoter$tx_name <- gsub("\\.[0-9]+","",promoter$tx_name)
      
      up_gene <- promoter[which(promoter$tx_name %in% up_gene)] 
      up_gene <- up_gene[!duplicated(up_gene$tx_name)]
      seqlevels(up_gene, pruning.mode="coarse") <- seqlevels(BSgenome.Athaliana.TAIR.TAIR9::Athaliana)
      up_seq <- getSeq(BSgenome.Athaliana.TAIR.TAIR9::Athaliana,up_gene)
      
      down_gene <- promoter[which(promoter$tx_name %in% down_gene)] 
      down_gene <- down_gene[!duplicated(down_gene$tx_name)]
      seqlevels(down_gene, pruning.mode="coarse") <- seqlevels(BSgenome.Athaliana.TAIR.TAIR9::Athaliana)
      down_seq <- getSeq(BSgenome.Athaliana.TAIR.TAIR9::Athaliana,down_gene)
      
      
      df <- motifmatchr::matchMotifs(pwms = pwm,subject = up_seq,out ="matches",p.cutoff=1e-05,bg="subject")
      match <- df@assays@data@listData$motifMatches %>% colSums()
      up_data <- match %>% as.data.frame() %>% `colnames<-`("value")
      up_data$"group" <-  colnames(rnaseq_tpm)[x+6]; up_data$"motif" <- names(pwm.l)#;up_data$"class" <- "up"

      
      df <- motifmatchr::matchMotifs(pwms = pwm,subject = down_seq,out ="matches",p.cutoff=1e-05,bg="subject")
      match <- df@assays@data@listData$motifMatches %>% colSums()
      down_data <- match %>% as.data.frame() %>% `colnames<-`("value")
      down_data$"group" <- colnames(rnaseq_tpm)[x+6]; down_data$"motif" <- names(pwm.l)#;down_data$"class" <- "down"
      up_data$value <-  up_data$value/down_data$value
      up_data
    }
    stopImplicitCluster()
    stopCluster(cl)
    
    dat <- bind_rows(result) 

    dat <- mutate(dat,motif = fct_relevel(motif, "class11", "class10", "class9", 
                                      "class8", "class7", "class6", 
                                      "class5", "class4","class3","class2","class1"))
    dat$group <- dat$group %>%  gsub(".[0-9]+_[0-9]+_","",.)
    
    dat <- mutate(dat,group  = fct_relevel(group , "Flower_1" , "Flower_2" , "Cotyledon_1"   ,
                                       "Cotyledon_2"  ,"Silique_1" ,"Silique_2",
                                       "Root_1"  ,"Root_2"   ,  "Stem_1"   , "Stem_2" ))
  
    
    col_draw <- colorspace::diverge_hcl(200,palette = "Blue-Red3",rev = F)#
    plot_RNA <- ggplot(data = dat,mapping = aes(x=group,y = as.factor(motif),fill=log2(value)))+
      geom_tile()+ 
      coord_fixed(ratio = 1)+
      scale_fill_gradientn(colors = col_draw,name = "log2 Enrichment",limits=c(-1,1))+
      theme(
        axis.text.x = element_text(angle = 30,size = 10,hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title =element_text(size=15)
      )    +labs(title = "10%  ")

    
  

#####################################################################
##3G ATI####
#####################################################################
  {

    list.files("jdk_data/1_ribosome/data2_GSE152938_kidney/",pattern = "gz$") %>% as.data.frame() %>% `colnames<-`(c(""))
    file <- read.table("~/diff_tissues_ATI.txt",header = TRUE)
    file$tissue <- rename_duplicates(file$tissue)
    cl=makeCluster(21)
    registerDoParallel(cl)
    result <- list()
    result <- foreach(x=seq(nrow(file)))%dopar%{
      library(SummarizedExperiment)
      library(universalmotif)
      library(dplyr)
      library(Biostrings)
      library(motifmatchr)
      library(fjComm)
      pseudo=1000
      seq <- readLines(file$clean_reads_30n[x]) %>% DNAStringSet()  
      
      df <- motifmatchr::matchMotifs(pwms = pwm,subject = seq,out ="matches",p.cutoff=1e-05,bg="subject")
      match <- df@assays@data@listData$motifMatches %>% colSums()
      
      df_s <- motifmatchr::matchMotifs(pwms = pwm,subject = seq %>% universalmotif::shuffle_sequences(),out ="matches",p.cutoff=1e-05,bg="subject")
      match_s <- df_s@assays@data@listData$motifMatches %>% colSums()
      
      heat_data <- (match+pseudo)/(match_s+pseudo) %>% as.data.frame() %>% `colnames<-`("value")
      heat_data$"group" <- file$tissue[x]; heat_data$"motif" <- names(pwm.l)
      heat_data
    }
    stopImplicitCluster()
    stopCluster(cl)
    
    
    save(result,file = "~/diff_tissue_ATI.Rdata")
  }

  result <- bind_rows(result)
  result$value <- log2(result$value)
  result <- mutate(result,motif = fct_relevel(motif, "class11", "class10", "class9", 
                                              "class8", "class7", "class6", 
                                              "class5", "class4","class3","class2","class1"))

  result <- mutate(result,group = fct_relevel(group,
                                              "flower_1"   , "flower_2"    ,"flower_3"   ,"flower_4"    ,"flower_5" ,
                                              "cotyledon_1", "cotyledon_2" ,"cotyledon_3", "cotyledon_4",
                                              "silique_1"  , "silique_2"   ,"silique_3"  ,"silique_4",
                                              "root_1"     , "root_2"      ,"root_3"     ,"root_4",
                                              "stem_1"     , "stem_2"      ,"stem_3"     ,"stem_4"))
  
  col_draw <- colorspace::diverge_hcl(200,palette = "Blue-Red3",rev = T)#
  plot_ATI <- ggplot(data = result,mapping = aes(x=group,y = as.factor(motif),fill=value))+
    geom_tile()+ 
    coord_fixed(ratio = 1)+
    scale_fill_gradientn(colors =  col_draw %>% rev(),name = "log2 Enrichment",limits=c(-5.1,5.1))+
    theme(
      axis.text.x = element_text(angle = 30,size = 10,hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title =element_text(size=15)
    )    +labs(title = "ATI pseudo=1000 ")
  


#####################################################################
##3G ATAC####
#####################################################################

  {
    #load file
    file <- list.files("~/ATAC/pre_processed/peak_calling/narrow",pattern = "bed$") %>% 
      paste0("~/ATAC/pre_processed/peak_calling/narrow",.) %>% as.data.frame()
    
    file$"tissues" <- file$. %>% basename() %>% gsub("_summits.bed","",.) %>% gsub("^.*_[0-9]+_","",.)
    colnames(file)[1] <- "file"
    pwm.l <- read_jaspar("~/class11_pwm.txt")
    pwm.l <- convert_motifs(motifs = pwm.l,class ="TFBSTools-PWMatrix" ) 
    for (i in 1:length(pwm.l)) {names(pwm.l)[i] <- pwm.l[[i]]@name    }
    pwm <- do.call(PWMatrixList,pwm.l)
    
    library(parallel)
    library(parallelly)
    library(foreach)
    library(doParallel)
    cl=makeCluster(15)
    registerDoParallel(cl)
    result <- list()
    result <- foreach(x=seq(nrow(file)))%dopar%{
      library(SummarizedExperiment)
      library(universalmotif)
      library(dplyr)
      library(Biostrings)
      library(motifmatchr)
      library(fjComm)
      library(rtracklayer)
      bed <- rtracklayer::import(file$file[x])
      pseudo=1000
      bed <- SummarizedExperiment::resize(bed , width(bed) + 100, fix = "center")
      bed <- bed[grep("^Chr", bed@seqnames),]
      
      seqlevels(bed, pruning.mode="coarse") <- seqlevels(BSgenome.Athaliana.TAIR.TAIR9::Athaliana)
      seq <- BSgenome::getSeq(BSgenome.Athaliana.TAIR.TAIR9::Athaliana,bed)
      df <- motifmatchr::matchMotifs(pwms = pwm,subject = seq,out ="matches",p.cutoff=1e-05,bg="subject")
      match <- df@assays@data@listData$motifMatches %>% colSums()

      match_s <- lapply(1:length(pwm), function(x){
        df_s <- motifmatchr::matchMotifs(pwms = pwm,subject = seq %>% universalmotif::shuffle_sequences(),out ="matches",p.cutoff=1e-05,bg="subject")
        match_s <- df_s@assays@data@listData$motifMatches %>% colSums()
        match_s
      }) %>% largeListToDf() %>% rowMeans()
      
      heat_data <- (match+pseudo)/(match_s+pseudo) %>% as.data.frame() %>% `colnames<-`("value")
      heat_data$"group" <- file$tissues[x]; heat_data$"motif" <- names(pwm.l)
      heat_data
    }
    stopImplicitCluster()
    stopCluster(cl)
  }
 #plot
  result <- bind_rows(result) 
  shuff_result <- bind_rows(shuff_result)
  shuff_result$group %<>%  paste0("Shuffle_",.)
  result <- rbind(shuff_result,result)
  result <- mutate(result,motif = fct_relevel(motif, 
                                              "class1", "class2", "class3", 
                                              "class4", "class5", "class6", 
                                              "class7", "class8","class9","class10","class11"))
  result <- mutate(result,group = fct_relevel(group,
                                              "flower1","flower2","flower3",
                                              "fruit1", "fruit2", "fruit3",
                                              "leaf1",  "leaf2",  "leaf3",
                                              "root1",  "root2",  "root3",
                                              "stem1",  "stem2",  "stem3",
                                              "Shuffle_flower1","Shuffle_flower2","Shuffle_flower3",
                                              "Shuffle_fruit1", "Shuffle_fruit2", "Shuffle_fruit3",
                                              "Shuffle_leaf1",  "Shuffle_leaf2",  "Shuffle_leaf3",
                                              "Shuffle_root1",  "Shuffle_root2",  "Shuffle_root3",
                                              "Shuffle_stem1",  "Shuffle_stem2",  "Shuffle_stem3"))
  
  plot_ATAC <- ggplot(data = result,mapping = aes(x=group,y = as.factor(motif),fill=log2(value)))+
    geom_tile()+ 
    coord_fixed(ratio = 1)+
    scale_fill_gradientn(colors =  colorspace::diverge_hcl(200,palette = "Blue-Red3",rev = F),name = "log2 Enrichment",limits=c(-1,1))+
    theme(
      axis.text.x = element_text(angle = 30,size = 10,hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title =element_text(size=15)
    )    +labs(title = "ATAC")
