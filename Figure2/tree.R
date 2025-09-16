# 构建系统发育树
 bash ~/2024_WRKY/all_species_wrky_phylotree/align_phylotree.sh \  -i ~/2024_WRKY/all_species_wrky_phylotree/72.WRKY.pep.form.tair.1.24.fa \   -o . -n WRKY_72 -t 30

library(treeio)
library(ggtree)

# 读取树
tree <- read.newick("~/2024_WRKY/7_Fig2/phylo_tree/WRKY_72.iqtree.treefile")

# 初步绘制树，显示节点编号
ggtree(tree, layout = "rectangular", branch.length = "none") +
  geom_tiplab() +
  geom_text(aes(label = node), hjust = 0.5, vjust = 0.5)   # label=node 便于旋转树

# 将树转为可读数据框
ggtree::fortify(tree)
data.tree <- as.tibble(tree)

# 整理 tip.label
tree[["tip.label"]] <- tree[["tip.label"]] %>%
  gsub("AT.*_AT", "", .) %>%
  gsub("AT", "", .)

# 绘制整理后的树
ggtree(tree, layout = "rectangular", branch.length = "none") +
  geom_tiplab(hjust = -0.3) +
  geom_text(aes(label = node), hjust = 0.5, vjust = 0.5) +
  xlim(NA, 14)

# 读取 subfamily 信息
subfamily <- read.table(
  file = "2024_WRKY/1_data.paths.for.analysis/subfamily.txt",
  header = TRUE
)
subfamily <- split(subfamily, subfamily$Subfamily)

# 生成 groupInfo
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

# 着色分组绘制
tree1 <- groupOTU(tree, groupInfo)
p <- ggtree(tree1, layout = "rectangular", branch.length = "none") +
  geom_tiplab(aes(color = group), hjust = -0.0) +
  geom_text(aes(label = node), hjust = 0.5, vjust = 0.5) +
  scale_color_manual(values = c(
    "black", "#E41A1C", "#4DAF3e", "blue",
    "#FF7F00", "#377EB8", "#A65628", "yellow"
  )) +
  theme(legend.position = "top") +
  xlim(NA, 14)

# 节点旋转
p <- flip(p, 22, 21)  # 旋转两个节点
p <- ggtree::rotate(p, 92)  # 旋转单个节点

# 保存结果
ggsave(
  filename = "2024_WRKY/7_Fig2/phylo_tree/plot/phylotree.pdf",
  plot = p,
  height = 15
)
