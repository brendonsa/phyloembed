library(TreeSim)
library(ape)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    stop(paste("Missing argument", flag))
  }
  args[[idx + 1]]
}

leaf <- as.integer(get_arg("--leaf"))
height <- as.numeric(get_arg("--height"))
speciation_rate <- as.numeric(get_arg("--speciation_rate"))
extinction_rate <- as.numeric(get_arg("--extinction_rate"))
seed <- as.integer(get_arg("--seed"))
out_nwk <- get_arg("--out_nwk")
out_nexus <- get_arg("--out_nexus")
out_times <- get_arg("--out_times")

dir.create(dirname(out_nwk), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_nexus), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_times), recursive = TRUE, showWarnings = FALSE)

set.seed(seed)
tree <- sim.bd.taxa.age(leaf, 1, speciation_rate, extinction_rate, 1, height, TRUE)[[1]]

depths <- node.depth.edgelength(tree)
tip_depths <- depths[seq_len(Ntip(tree))]
target_depth <- height
scale_factor <- target_depth / max(tip_depths)
tree$edge.length <- tree$edge.length * scale_factor

depths <- node.depth.edgelength(tree)
tip_depths <- depths[seq_len(Ntip(tree))]
eps <- 1e-8
for (i in seq_len(Ntip(tree))) {
  diff <- target_depth - tip_depths[i]
  if (abs(diff) > eps) {
    edge_idx <- which(tree$edge[, 2] == i)
    tree$edge.length[edge_idx] <- tree$edge.length[edge_idx] + diff
  }
}

options(digits = 16, scipen = 999)
nwk <- write.tree(tree, digits = 16)
nwk <- sub(";$", "", nwk)
nwk <- sub(":[0-9.eE+-]+$", "", nwk)
write.tree(tree, file = out_nwk, digits = 16)

tree$tip.label <- as.character(seq_along(tree$tip.label))
nwk_num <- write.tree(tree, digits = 16)
nwk_num <- sub(";$", "", nwk_num)
nwk_num <- trimws(nwk_num)
nwk_num <- sub(":[0-9.eE+-]+$", "", nwk_num)

cat(
  "#NEXUS\nBEGIN TREES;\nTREE sp_tree = ",
  nwk_num,
  ";\nEND;\n",
  file = out_nexus,
  sep = ""
)

times <- node.depth.edgelength(tree)
times_df <- data.frame(node = seq_along(times), time = times)
write.table(times_df, file = out_times, sep = "\t", quote = FALSE, row.names = FALSE)
