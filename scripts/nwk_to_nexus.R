suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript scripts/nwk_to_nexus.R <in.nwk> <out.nex>")

infile <- args[1]
outfile <- args[2]

tr <- read.tree(infile)

nwk <- write.tree(tr)

con <- file(outfile, open = "wt")
writeLines(c(
  "#NEXUS",
  "BEGIN TREES;",
  paste0("  TREE tree1 = ", nwk),
  "END;"
), con)
close(con)
