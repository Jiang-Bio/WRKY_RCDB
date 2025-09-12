# -----------------------------------------------------------
#fig3A####
# -----------------------------------------------------------
## -------------------------------
## Motif alignment, clustering, and circular tree visualization with stars
## -------------------------------

library(universalmotif)
library(motifStack)
library(ggtree)
library(ggnewscale)
library(ggplot2)
library(fjComm)
library(stringr)
library(dplyr)
library(ade4)

# Load PWM list
pwmlist <- read_jaspar("2024_WRKY/8_motif_net/ALL_mono_P_S_T.pwm")

# Align motifs
pwmlist_alig <- DNAmotifAlignment(
  threshold = 0,
  pfms = pwmlist %>% convert_motifs(class = "motifStack-pfm"),
  revcomp = rep(FALSE, length(pwmlist))
)

# Trim PWM matrices (remove flanking columns)
for (i in seq_along(pwmlist_alig)) {
  pwmlist_alig[[i]]@mat <- pwmlist_alig[[i]]@mat[, -c(1:2, 16:17)]
}

# Rename motifs for clarity
c <- pwmlist_alig
for (i in seq_along(c)) {
  c[[i]]@name <- c[[i]]@name %>%
    gsub("P", "_1", .) %>%
    gsub("S", "_2", .) %>%
    gsub("T", "_3", .) %>%
    gsub("F", "_4", .)
  names(c)[i] <- c[[i]]@name
}

# Construct motif correlation matrix for clustering
p <- data.frame()
for (i in seq_along(c)) {
  p1 <- t(rbind(c[[i]]@mat[1, ], c[[i]]@mat[2, ], c[[i]]@mat[3, ], c[[i]]@mat[4, ]))
  p <- if (i == 1) p1 else cbind(p, p1)
}
colnames(p) <- names(c)
p <- cor(p, method = "spearman")

# Hierarchical clustering
a <- dist(x = p, method = "maximum")
hc <- hclust(a)

# Convert to phylog object for circular tree plotting
phylog <- ade4::hclust2phylog(hc, add.tools = TRUE)
p <- ggtree(phylog, layout = "circular", size = 0.6) +
  geom_text(aes(label = node), size = 1.5) +
  geom_tiplab(hjust = -0.5, size = 4.2, align = TRUE, offset = 0) +
  ggtree::rotate(220) # Rotate tree

# Prepare star plot data for each PWM
dt <- lapply(c, function(pwm) {
  df <- data.frame()
  for (i in seq_len(ncol(pwm@mat))) {
    id <- pwm@name
    group <- which.max(pwm@mat[, i]) %>% names()
    size <- pwm@mat[which.max(pwm@mat[, i]), i]
    strain <- i
    df <- rbind(df, c(id, group, size, strain))
  }
  colnames(df) <- c("id", "group", "size", "strain")
  df$size <- as.numeric(df$size)
  df$strain <- as.numeric(df$strain)
  df
}) %>% do.call(rbind, .)
row.names(dt) <- NULL

# Separate motif names into two groups
p1 <- names(c)[str_detect(names(c), "m")]
p2 <- names(c)[!str_detect(names(c), "m")]

# Generate geom_strip strings
text1 <- paste0("geom_strip(\"", p1, "\", \"", p1, "\", offset = 0.06, barsize = 20, extend = 0.6, color = \"#00BFC4\", offset.text = 3)") %>% paste0(collapse = "+")
text2 <- paste0("geom_strip(\"", p2, "\", \"", p2, "\", offset = 0.06, barsize = 20, extend = 0.6, color = \"#F8766D\", offset.text = 3)") %>% paste0(collapse = "+")

# Generate class strips based on clustering
tree <- cutree(hc, h = 0.37)
colors <- c("#FFFFB3","#E5C494","#B3DE69","#85CBBF","#F0766D","#80A9CB","#B6B2D2","#F4C5DD","#FF9D00","#F35EAB","#6A3D9A")
class_strips <- lapply(seq_along(colors), function(i){
  class_nodes <- names(tree[tree == i])
  paste0("geom_strip(\"", class_nodes, "\", \"", class_nodes, "\", offset = 0.15, barsize = 10, extend = 0.6, color = \"", colors[i], "\", angle = 10, hjust = \"center\", fontsize = 4, offset.text = 2)") %>% paste0(collapse = "+")
}) %>% paste0(collapse = "+")

# Combine all strips
plot_text <- paste0("p + ", text1, "+", text2, "+", class_strips)
plot <- eval(parse(text = plot_text))

# Add stars to circular tree
p1 <- plot + geom_fruit(
  data = dt,
  geom = geom_star,
  mapping = aes(x = strain, y = id, fill = group, size = size),
  starstroke = 0,
  starshape = 13,
  offset = 0.24,
  pwidth = 0.30
) +
  scale_size_continuous(range = c(0, 4), limits = c(min(dt$size), max(dt$size)), breaks = c(1, 2, 3)) +
  scale_fill_manual(values = c("#109648", "#255C99", "#F7B32B", "#D62839"),
                    guide = guide_legend(keywidth = 0.4, keyheight = 0.4, order = 4,
                                         override.aes = list(starstroke = 0.3))) +
  geom_tiplab(hjust = -0.5, size = 4.2, align = TRUE, offset = 0) +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.key.size = unit(10, "lines")
  )

# Save plot
ggsave(filename = "2024_WRKY/8_motif_net/plot/fig3a.pdf", plot = p1, width = 20, height = 20)
