##########################################################################
#Figure3A compares, clusters and plots all monomer motifs
##########################################################################
pwmlist <- read_jaspar("~/ALL_mono_P_S_T.pwm")
  #矩阵比对
  pwmlist_alig <- DNAmotifAlignment(threshold = 0,pfms = pwmlist %>% convert_motifs(class = "motifStack-pfm"),
                                    revcomp= rep(FALSE, length(pwmlist)))# %>% convert_motifs("TFBSTools-PWMatrix") %>% TF_logo.2(outfile = "2.txt.pdf")
  #裁切矩阵保留核心区域使其等长
  for (i in 1:length(pwmlist_alig)) {
    pwmlist_alig[[i]]@mat <- pwmlist_alig[[i]]@mat[,-c(1:2,16:17)]
  }
  c <- pwmlist_alig
  for (i in 1:length(c)) {
    c[[i]]@name <- c[[i]]@name %>%  gsub("P","_1",.) %>% gsub("S","_2",.) %>% gsub("T","_3",.) %>% gsub("F","_4",.)
    names(c)[i] <- c[[i]]@name 
    
  }
  #一维展开矩阵并计算距离进行聚类
  p <- data.frame()
  for (i in 1:length(c)) {
    p1 <- cbind(c[[i]]@mat[1,] %>% t(),c[[i]]@mat[2,] %>% t(),c[[i]]@mat[3,] %>% t(),c[[i]]@mat[4,] %>% t()) %>% t()
    if (i==1) {
      p <- rbind(p,p1)
    }else{  p <- cbind(p,p1)  }
  }
  colnames(p) <- names(c);rownames(c) <- NULL
  p <- cor(p,method = "spearman")

  a <- dist(x = p,method = "maximum")
  hc <- do.call("hclust", list(d = a))
  library(ggtree)
  phylog <- ade4::hclust2phylog(hc,add.tools = TRUE)

  p <- ggtree(phylog,layout = "circular",size=0.6)  +
    geom_text(aes(label=node),size=1.5)+
    geom_tiplab(hjust=-0.5,size=4.2,align = TRUE,offset = 0)
  
  # p <- flip(p, 319,241) #旋转两个节点
  p <- ggtree::rotate(p,220)#旋转单个节点
  # geom_text(aes(label=node),size=1)
  dt <- lapply(c, function(pwm){
    dt <- data.frame()
    for (i in seq(ncol(pwm@mat))) {
      id=pwm@name
      group=which.max(pwm@mat[,i]) %>% names()
      size=pwm@mat[which.max(pwm@mat[,i]),i]
      strain=i
      a <- c(id,group,size,strain)
      dt <- rbind(dt,a)
    }
    colnames(dt) <- c("id","group","size","strain")
    return(dt)
  }) %>% do.call(rbind,.)
  dt$size <- as.numeric(dt$size)
  dt$strain <- as.numeric(dt$strain)
  row.names(dt) <- NULL

  p1 <- names(c)[names(c) %>% grep("m",.)]
  p2 <-names(c)[-( names(c) %>% grep("m",.))]
  text1 <-paste0("geom_strip(\"", p1, "\", \"",p1,"\", offset = 0.06, barsize = 20, extend = 0.6, color = \"#00BFC4\", offset.text = 3)") %>% paste0(collapse = "+")#offset :颜色bar 的外移内移
  text2 <-paste0("geom_strip(\"", p2, "\", \"",p2,"\", offset = 0.06, barsize = 20, extend = 0.6, color = \"#F8766D\", offset.text = 3)") %>% paste0(collapse = "+")
  {
    #根据高度0.37进行砍树，将结果分为11类
    tree <- cutree(tree = hc,h = 0.37)
    class1 <- tree[tree==1] %>% names()
    class2 <- tree[tree==2] %>% names()
    class3 <- tree[tree==3] %>% names()
    class4 <- tree[tree==4] %>% names()
    class5 <- tree[tree==5] %>% names()
    class6 <- tree[tree==6] %>% names()
    class7 <- tree[tree==7] %>% names()
    class8 <- tree[tree==8] %>% names()
    class9 <- tree[tree==9] %>% names()
    class10 <- tree[tree==10] %>% names()
    class11 <- tree[tree==11] %>% names()

    class1  <- paste0("geom_strip(\"", class1 , "\", \"",class1 ,"\", offset = 0.15, barsize = 10, extend = 0.6, color = \"#FFFFB3\",angle = 10,hjust = \"center\",fontsize = 4, offset.text = 2)") %>% paste0(collapse = "+")
    class2  <- paste0("geom_strip(\"", class2 , "\", \"",class2 ,"\", offset = 0.15, barsize = 10, extend = 0.6, color = \"#E5C494\",angle = 10,hjust = \"center\",fontsize = 4, offset.text = 2)") %>% paste0(collapse = "+")
    class3  <- paste0("geom_strip(\"", class3 , "\", \"",class3 ,"\", offset = 0.15, barsize = 10, extend = 0.6, color = \"#B3DE69\",angle = 10,hjust = \"center\",fontsize = 4, offset.text = 2)") %>% paste0(collapse = "+")
    class4  <- paste0("geom_strip(\"", class4 , "\", \"",class4 ,"\", offset = 0.15, barsize = 10, extend = 0.6, color = \"#85CBBF\",angle = 10,hjust = \"center\",fontsize = 4, offset.text = 2)") %>% paste0(collapse = "+")
    class5  <- paste0("geom_strip(\"", class5 , "\", \"",class5 ,"\", offset = 0.15, barsize = 10, extend = 0.6, color = \"#F0766D\",angle = 10,hjust = \"center\",fontsize = 4, offset.text = 2)") %>% paste0(collapse = "+")
    class6  <- paste0("geom_strip(\"", class6 , "\", \"",class6 ,"\", offset = 0.15, barsize = 10, extend = 0.6, color = \"#80A9CB\",angle = 10,hjust = \"center\",fontsize = 4, offset.text = 2)") %>% paste0(collapse = "+")
    class7  <- paste0("geom_strip(\"", class7 , "\", \"",class7 ,"\", offset = 0.15, barsize = 10, extend = 0.6, color = \"#B6B2D2\",angle = 10,hjust = \"center\",fontsize = 4, offset.text = 2)") %>% paste0(collapse = "+")
    class8  <- paste0("geom_strip(\"", class8 , "\", \"",class8 ,"\", offset = 0.15, barsize = 10, extend = 0.6, color = \"#F4C5DD\",angle = 10,hjust = \"center\",fontsize = 4, offset.text = 2)") %>% paste0(collapse = "+")
    class9  <- paste0("geom_strip(\"", class9 , "\", \"",class9 ,"\", offset = 0.15, barsize = 10, extend = 0.6, color = \"#FF9D00\",angle = 10,hjust = \"center\",fontsize = 4, offset.text = 2)") %>% paste0(collapse = "+")
    class10  <- paste0("geom_strip(\"", class10, "\",\"",class10 ,"\", offset = 0.15, barsize = 10, extend = 0.6,color = \"#F35EAB\",angle = 10,hjust = \"center\",fontsize = 4, offset.text = 2)") %>% paste0(collapse = "+")
    class11  <- paste0("geom_strip(\"", class11, "\",\"",class11 ,"\", offset = 0.15, barsize = 10, extend = 0.6,color = \"#6A3D9A\",angle = 10,hjust = \"center\",fontsize = 4, offset.text = 2)") %>% paste0(collapse = "+")
   
  }
  
  plot_text=paste0("p + ",text1,"+",text2,"+",class1,"+",class2,"+",class3,"+",
                   class4,"+",class5,"+",class6,"+",class7,"+",class8,"+",class9,"+",class10,"+",class11)
  plot <-  eval(parse(text = plot_text))
  

  p1 <- plot +geom_fruit(
    data=dt,
    geom=geom_star,
    mapping=aes(x=strain,y=id, fill=group,size=size),
    starstroke=0,
    starshape = 13,
    offset =0.24,
    pwidth=0.30
  )+
    scale_size_continuous(range=c(0, 4),
                          limits=c(min(dt$size), max(dt$size)),breaks=c(1, 2, 3))+
    scale_fill_manual(values=c("#109648", "#255C99", "#F7B32B", "#D62839"),
                      guide=guide_legend(keywidth = 0.4, keyheight = 0.4, order=4,
                                         override.aes = list(starstroke=0.3)))+
    geom_tiplab(hjust=-0.5,size=4.2,align = TRUE,offset = 0)+
    theme(legend.title = element_text(size = 20),   
          legend.text = element_text(size = 20),    
          legend.key.size = unit(10, "lines")   )





#根据树的聚类结果生成11类代表性motif矩阵
    class1 <- tree[tree==1] %>% names()
    class2 <- tree[tree==2] %>% names()
    class3 <- tree[tree==3] %>% names()
    class4 <- tree[tree==4] %>% names()
    class5 <- tree[tree==5] %>% names()
    class6 <- tree[tree==6] %>% names()
    class7 <- tree[tree==7] %>% names()
    class8 <- tree[tree==8] %>% names()
    class9 <- tree[tree==9] %>% names()
    class10 <- tree[tree==10] %>% names() 
    class11 <- tree[tree==11] %>% names()
  
  class <- list(class1,class3,class4,class2,class7,class6,class11,class10,class9,class8,class5)
  #生成代表性motif，平均IC大于0.5的，赋0.8的权重，小于0.5的赋0.2的权重
  pwm.l <- lapply(1:length(class), function(i){
    c1 <- c[class[[i]]] %>% convert_motifs("universalmotif-universalmotif")
   for(x in 1:length(c1)){ if( 
     ((c1[[x]]@icscore )/ncol(c1[[x]]@motif) ) >0.5 ){
     c1[[x]]@motif <- c1[[x]]@motif *0.8
    }else{ c1[[x]]@motif <- c1[[x]]@motif *0.2}  }

    a3 <- paste0("c1[[",c(1:length(c1)),"]]@motif",collapse = "+") %>% parse(text = .) %>% eval()
    a3 <- sweep(a3, 2, colSums(a3), FUN = "/") %>% as.matrix() %>% convert_motifs("universalmotif-universalmotif")
    a3@name <- glue("class{i}")
    a3
  })

    write_jaspar(motifs = pwm.l,file = "~/class11_pwm.txt",overwrite = T)
    
