suppressPackageStartupMessages(library(TreeSim))
suppressPackageStartupMessages(library(ape))

parse_args <- function(x) {
  out <- list()
  i <- 1
  while (i <= length(x)) {
    key <- x[i]
    if (startsWith(key, "--")) {
      key <- sub("^--", "", key)
      out[[key]] <- x[i + 1]
      i <- i + 2
    } else stop(paste("Unexpected token:", x[i]))
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

n_tips <- as.integer(args[["n_tips"]])
tree_height <- as.numeric(args[["tree_height"]])
lambda <- as.numeric(args[["lambda"]])
mu <- as.numeric(args[["mu"]])
seed <- as.integer(args[["seed"]])
out <- args[["out"]]

set.seed(seed)

st <- sim.bd.taxa.age(
  n = n_tips,
  numbsim = 1,
  lambda = lambda,
  mu = mu,
  age = tree_height
)[[1]]


st$node.time <- node.depth.edgelength(st)


write.tree(st, file = out)
