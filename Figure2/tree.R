# 构建系统发育树
 #bash ~/align_phylotree.sh \  -i ~/WRKY.pep.form.tair.1.24.fa \   -o . -n WRKY -t 30

library(treeio)
library(ggtree)

tree <- read.newick("~/WRKY_72.iqtree.treefile")

ggtree(tree, layout = "rectangular", branch.length = "none") +
  geom_tiplab() +
  geom_text(aes(label = node), hjust = 0.5, vjust = 0.5)   # label=node 便于旋转树

ggtree::fortify(tree)
data.tree <- as.tibble(tree)


tree[["tip.label"]] <- tree[["tip.label"]] %>%
  gsub("AT.*_AT", "", .) %>%
  gsub("AT", "", .)


ggtree(tree, layout = "rectangular", branch.length = "none") +
  geom_tiplab(hjust = -0.3) +
  geom_text(aes(label = node), hjust = 0.5, vjust = 0.5) +
  xlim(NA, 14)


subfamily <- read.table(
  file = "~/subfamily.txt",
  header = TRUE
)
subfamily <- split(subfamily, subfamily$Subfamily)

c1 <- c(); c2 <- c()
for (i in seq(length(subfamily))) {
  group <- subfamily[[i]]$TF
  for (ii in seq_len(length(group))) {
    c1 <- c(
      c1,
      data.tree[grep(glue("{group[ii]}$"), x = data.tree$label), 4]
    )
  }
  c2[[i]] <- as.character(c1)
  c1 <- c()
}

groupInfo <- list(
  GroupI   = c2[[1]],
  GroupIIa = c2[[2]],
  GroupIIb = c2[[3]],
  GroupIIc = c2[[4]],
  GroupIId = c2[[5]],
  GroupIIe = c2[[6]],
  GroupIII = c2[[7]]
)

tree1 <- groupOTU(tree, groupInfo)
phylotree <- ggtree(tree1, layout = "rectangular", branch.length = "none") +
  geom_tiplab(aes(color = group), hjust = -0.0) +
  geom_text(aes(label = node), hjust = 0.5, vjust = 0.5) +
  scale_color_manual(values = c(
    "black", "#E41A1C", "#4DAF3e", "blue",
    "#FF7F00", "#377EB8", "#A65628", "yellow"
  )) +
  theme(legend.position = "top") +
  xlim(NA, 14)


phylotree <- flip(phylotree, 22, 21)  # 旋转两个节点
phylotree <- ggtree::rotate(phylotree, 92)  # 旋转单个节点

phylotree
