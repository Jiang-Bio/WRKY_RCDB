
# # 运行 meme
# meme <- c("/wrk/chenhao/meme/bin/meme ~/WRKY.pep.form.tair.1.24.fa -protein -oc ~/pep_structure/meme/ -nmotifs 10 -minw 6 -maxw 20")

data <- readLines("~/pep_structure/meme/meme.txt")
start_index <- grep("sites sorted by position p-value", data) + 4

motif <- list()
for (i in start_index) {
  end_index <- (grep("^-+$", data[i:length(data)]) + i - 2)[1]
  motif[[as.character(i)]] <- data[i:end_index]
}


df <- c()
for (i in seq_along(motif)) {
  a <- motif[[i]] %>%
    str_split(pattern = "\\s+") %>%
    largeListToDf()

  a <- a[1:2, ] %>% t()
  colnames(a) <- c("molecule", "start")

  a <- as.data.frame(a)
  a$start <- as.numeric(a$start)
  a$end   <- a$start + 20
  a$gene  <- as.character(i)

  df <- rbind(df, a)
}

df$gene     <- df$gene %>% paste0("motif", .)
df$molecule <- df$molecule %>% str_extract(pattern = "WRKY[0-9]+")
df          <- df[, c("molecule", "gene", "start", "end")]
df$strand   <- "forward"

data <- df


dummies <- make_alignment_dummies(
  data,
  aes(xmin = start, xmax = end, y = molecule, id = gene),
  on = "motif1"
)

pep_structure <- ggplot(data, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow(
    arrowhead_height = unit(3, "mm"),
    arrowhead_width  = unit(0.7, "mm")
  ) +
  geom_blank(data = dummies) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x  = element_blank(),  
    axis.line.x  = element_blank(),  
    axis.ticks.x = element_blank()
  )
pep_structure

