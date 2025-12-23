library("TreeSim")
library("ape")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) stop(paste("Missing", flag))
  args[idx + 1]
}

l <- as.integer(get_arg("--l"))
h <- as.numeric(get_arg("--h"))

br <- 10^(-7)
dr <- 0

st <- sim.bd.taxa.age(l, 1, br, dr, 1, h, TRUE)[[1]]

# SimPhy: root must NOT have branch length
st$root.edge <- NULL

nwk <- write.tree(st)
nwk <- gsub("[\r\n\t ]+", "", nwk)          # make it one token
if (!grepl(";$", nwk)) nwk <- paste0(nwk, ";")

cat(nwk)
