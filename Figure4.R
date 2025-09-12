# =========================================
# figure 4A 4D
# =========================================

# ------------------------
# 1. Define input files and extract sample names
# ------------------------
files <- c(
  "./A10_10_WRKY4N_1__101n.gz",
  "./A11_11_WRKY4N_2__101n.gz",
  "./A4_4_WRKY4_1__101n.gz",
  "./A5_5_WRKY4_2__101n.gz",
  "./D1_37_WRKY4C_2__101n.gz",
  "./D2_38_WRKY4C_4__101n.gz",
  "./B1_13_WRKY4NC_1__101n.gz",
  "·/B2_14_WRKY4NC_2__101n.gz"
)

# Extract WRKY sample names
fileNames <- str_extract(basename(files), "WRKY[0-9A-Za-z]+")
files <- files[!is.na(fileNames)]
fileNames <- fileNames[!is.na(fileNames)]

# ------------------------
# 2. Set analysis parameters
# ------------------------
bk_div_pseudo <- 50
max_mismatch <- 0
pattern1 <- "GTMAA"
pattern2 <- "GTMAA"

# ------------------------
# 3. Parallel processing of sequences
# ------------------------
library(parallel)
library(foreach)
library(doParallel)

cl <- makeCluster(10)
registerDoParallel(cl)

result <- foreach(x = seq_along(files)) %dopar% {
  library(fjComm)
  
  # Define gap range
  gap <- 0:30
  
  # Load sequences from FASTQ/fasta file
  seqs <- getSeq_fqfachrFile(files[x])
  
  # Include reverse complement sequences
  seqs <- c(seqs, revComp(seqs))
  
  # Calculate dimer enrichment
  all_results <- fjComm::dimer_enrichment(
    seqs,
    direct = c("DR", "IR", "ER"),
    half_pattern1 = pattern1,
    half_pattern2 = pattern2,
    max_mismatch = max_mismatch,
    gap = gap,
    bk_use_lm_fit = TRUE,
    bk_sub = TRUE,
    bk_div = TRUE,
    bk_div_pseudo = bk_div_pseudo
  )
  all_results
}

stopImplicitCluster()
stopCluster(cl)

# Assign names to results
names(result) <- fileNames %>% rename_duplicates()

# ------------------------
# 4. Combine ER, IR, DR matrices
# ------------------------
IR <- DR <- ER <- c()
for (i in seq_along(result)) {
  IR <- rbind(IR, result[[i]]$IR)
  DR <- rbind(DR, result[[i]]$DR)
  ER <- rbind(ER, result[[i]]$ER)
}

rownames(IR) <- rownames(DR) <- rownames(ER) <- names(result)

merge_m <- cbind(ER, IR, DR)

# Normalize each row to 0–1
merge_m <- t(apply(merge_m, 1, function(row) (row - min(row)) / (max(row) - min(row))))

# ------------------------
# 5. Plot heatmap
# ------------------------
pdf("~/WRKY4CN.pdf", width = 15, height = 12)
ComplexHeatmap::Heatmap(
  matrix = merge_m,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  border = FALSE,
  column_split = c(rep("ER", 31), rep("IR", 31), rep("DR", 31)),
  col = colorspace::sequential_hcl(n = 50, "PuBu") %>% rev(),
  column_names_gp = gpar(fontsize = 5),
  row_names_gp = gpar(fontsize = 8),
  name = "(enrich)_max",
  row_title = "dimer",
  column_labels = c(0:30, 0:30, 0:30)
)
dev.off()
