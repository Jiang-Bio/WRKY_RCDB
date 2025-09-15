library(fjComm)
bk_div_pseudo <- 50
max_mismatch  <- 0
pattern1      <- "GTMAA"
pattern2      <- "GTMAA"

  files <- c(
    "~/A10_10_WRKY4N_1__101n.gz",
    "~/A11_11_WRKY4N_2__101n.gz",
    "~/A4_4_WRKY4_1__101n.gz"   ,
    "~/A5_5_WRKY4_2__101n.gz"   ,
    "~/D1_37_WRKY4C_2__101n.gz"  ,
    "~/D2_38_WRKY4C_4__101n.gz"  ,
    "~/B1_13_WRKY4NC_1__101n.gz",
    "~/B2_14_WRKY4NC_2__101n.gz"
  )
#
cl <- makeCluster(10)
registerDoParallel(cl)
result <- foreach(x = seq(length(files))) %dopar% {
  library(fjComm)
  gap  <- 0:30
  seqs <- getSeq_fqfachrFile(files[x])
  seqs <- c(seqs, revComp(seqs))
  
  all_results <- fjComm::dimer_enrichment(
    seqs,
    direct        = c("DR", "IR", "ER"),
    half_pattern1 = pattern1,
    half_pattern2 = pattern2,
    max_mismatch  = max_mismatch,
    gap           = gap,
    bk_use_lm_fit = TRUE,
    bk_sub        = TRUE,
    bk_div        = TRUE,
    bk_div_pseudo = bk_div_pseudo
  )
}
stopImplicitCluster()
stopCluster(cl)

names(result) <- c(fileNames %>% rename_duplicates())

IR <- c()
DR <- c()
ER <- c()

for (i in 1:length(result)) {
  IR <- rbind(IR, result[[i]]$IR)
  DR <- rbind(DR, result[[i]]$DR)
  ER <- rbind(ER, result[[i]]$ER)
}

rownames(IR) <- names(result)
rownames(DR) <- names(result)
rownames(ER) <- names(result)

merge_m <- cbind(ER, IR, DR)
merge_m <- t(apply(merge_m, 1, function(row) {
  (row - min(row)) / (max(row) - min(row))
}))

ComplexHeatmap::Heatmap(
  matrix           = merge_m,
  cluster_rows     = TRUE,
  cluster_columns  = FALSE,
  border           = FALSE,
  column_split     = c(rep("ER", 31), rep("IR", 31), rep("DR", 31)),
  col              = colorspace::sequential_hcl(n = 50, "PuBu") %>% rev(),
  column_names_gp  = gpar(fontsize = 5),
  row_names_gp     = gpar(fontsize = 8),
  name             = "(enrich)_max",
  row_title        = "dimer",
  column_labels    = c(0:30, 0:30, 0:30)
)

