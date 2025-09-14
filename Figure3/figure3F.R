
    
    width=5000
    p.cutoff = 1e-4
    Txdb <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
    seqlevels(Txdb) <- c(paste0("Chr", 1:5),"Mt","Pt")
    seqlevels(Txdb,pruning.mode="coarse")<- seqlevels(Txdb)[str_detect(seqlevels(Txdb),"[cC][hH][rR][0-9]+")]
    genome=readDNAStringSet("~/aragenome/ara.fa")
    
    promoters_=promoters(Txdb,upstream = width,downstream = width) %>% GenomicRanges::trim() %>% {.[width(.)==2*width]}
    seqs=genome[promoters_]
    seqnum=length(seqs)
    
    
    class_motif <- read_jaspar("~/class11_pwm.txt" ) %>% convert_motifs("TFBSTools-PWMatrix")
    for (i in 1:length(class_motif)) {
      class_motif[[i]]@profileMatrix <- class_motif[[i]]@profileMatrix[,-c(11:13)]
    }
    motifs_pwm <- do.call(PWMatrixList,class_motif)
    
    
    
    
    motif.enrich.in.TSS <- function(motifs_pwm,width=5000){
      
      motifs_rev=lapply(motifs_pwm, function(motif){
        motif@profileMatrix %>% .[,ncol(.):1] %>% set_rownames(qw("A C G T")) %>%  universalmotif::create_motif() %>% {.@name=motif@name;.}
      })
      motifs_rev_pwm=map(motifs_rev, ~convert_motifs(.x,"TFBSTools-PWMatrix")) %>% {do.call(PWMatrixList,.)}
      motif_hits=motifmatchr::matchMotifs(motifs_pwm, seqs %>% as.character(),
                                          out = "positions",bg ="subject",p.cutoff = p.cutoff,
                                          ranges=GenomicRanges::GRanges(seqnames="Chr1",ranges=IRanges(start = rep(1,seqnum),width = rep(2*width,seqnum) ))  )

      results=lapply(seq_along(motifs_pwm),function(ind){
        (motif_hits[[ind]]@ranges@start) %>% table() %>% as.data.frame() %>% mutate(motif=motifs_pwm[[ind]]@name)   #name 是null
      })
      result_df=do.call(rbind,results) %>% set_colnames(qw("spacing freq motif"))
      result_df %<>% pivot_wider(names_from = spacing,values_from = freq)
      result_df %<>% column_to_rownames("motif")
      
   
      trans_dist=result_df %>% as.matrix() %>% {.[is.na(.)]=0;.} %>% sweep(1,rowSums(.),"/") %>%
        apply(1,function(x)zoo::rollmean(x,20)) %>% t

      col_draw=colorspace::sequential_hcl(200,palette = "PuBuGn",rev = T)
      ComplexHeatmap::Heatmap(trans_dist,col = col_draw, cluster_rows = F, show_row_dend = F, cluster_columns =F,use_raster = T,
                              row_names_gp = gpar(fontsize=12))

    motif.enrich.in.TSS(motifs_pwm,width = width)
    dev.off()
  }
