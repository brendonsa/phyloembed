#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(phangorn)
  library(dplyr)
  library(TreeDist)
})

option_list <- list(
  make_option(c("--t1"), type = "character", help = "First tree"),
  make_option(c("--t2"), type = "character", help = "Second tree"),
  make_option(c("--out"), type = "character", help = "Output CSV")
)

opt <- parse_args(OptionParser(option_list = option_list))

BIG_PENALTY <- 1e9

metric_names <- c(
  "RF_raw", "RF_normalized", "Weighted_RF", "BranchScore_KF",
  "JRF", "JRF_normalized", "Nye", "Nye_normalized",
  "SMI", "SMI_normalized"
)

# ---- helpers ----
normalize_branch_lengths <- function(tr) {
  if (is.null(tr$edge.length)) return(tr)
  el <- tr$edge.length
  el[is.na(el)] <- 0
  el[el < 0] <- 0
  s <- sum(el)
  tr$edge.length <- if (s > 0) el / s else el
  tr
}

# per-metric safety wrapper: only this metric gets BIG_PENALTY on failure
safe_metric <- function(expr_fun) {
  tryCatch(
    {
      val <- expr_fun()
      if (length(val) == 0 || !is.finite(val[1])) BIG_PENALTY else val[1]
    },
    error = function(e) BIG_PENALTY
  )
}

# preprocess trees; if this fails, we can't compute *any* metric → all penalized
preprocess_trees <- function() {
  tryCatch(
    {
      t1 <- read.tree(opt$t1)
      t2 <- read.tree(opt$t2)

      common_tips <- intersect(t1$tip.label, t2$tip.label)
      if (length(common_tips) < 3) {
        stop("fewer than 3 common tips between trees – cannot compare.")
      }

      t1 <- keep.tip(t1, common_tips)
      t2 <- keep.tip(t2, common_tips)

      t1 <- unroot(multi2di(t1))
      t2 <- unroot(multi2di(t2))

      t1 <- normalize_branch_lengths(t1)
      t2 <- normalize_branch_lengths(t2)

      list(tree1 = t1, tree2 = t2)
    },
    error = function(e) {
      message(sprintf(
        "[!] Preprocessing failed for %s vs %s: %s",
        opt$t1, opt$t2, e$message
      ))
      NULL
    }
  )
}

prep <- preprocess_trees()

if (is.null(prep)) {
  # preprocessing error → big value for all metrics
  metrics <- setNames(as.list(rep(BIG_PENALTY, length(metric_names))), metric_names)
} else {
  tree1 <- prep$tree1
  tree2 <- prep$tree2

  rf_raw         <- safe_metric(function() RF.dist(tree1, tree2, normalize = FALSE)[1])
  rf_normalized  <- safe_metric(function() RF.dist(tree1, tree2, normalize = TRUE)[1])
  w_rf           <- safe_metric(function() wRF.dist(tree1, tree2)[1])
  kf             <- safe_metric(function() KF.dist(tree1, tree2))
  jrf            <- safe_metric(function() JaccardRobinsonFoulds(tree1, tree2, normalize = FALSE)[1])
  jrf_normalized <- safe_metric(function() JaccardRobinsonFoulds(tree1, tree2, normalize = TRUE)[1])
  ns             <- safe_metric(function() NyeSimilarity(tree1, tree2, normalize = FALSE)[1])
  ns_normalized  <- safe_metric(function() NyeSimilarity(tree1, tree2, normalize = TRUE)[1])
  smi            <- safe_metric(function() SharedPhylogeneticInfo(tree1, tree2)[1])
  smi_normalized <- safe_metric(function() SharedPhylogeneticInfo(tree1, tree2, normalize = TRUE)[1])

  metrics <- list(
    RF_raw         = rf_raw,
    RF_normalized  = rf_normalized,
    Weighted_RF    = w_rf,
    BranchScore_KF = kf,
    JRF            = jrf,
    JRF_normalized = jrf_normalized,
    Nye            = ns,
    Nye_normalized = ns_normalized,
    SMI            = smi,
    SMI_normalized = smi_normalized
  )
}

# ---- output ----
results <- tibble(
  metric = names(metrics),
  value  = unlist(metrics, use.names = FALSE)
)

write.csv(results, opt$out, row.names = FALSE)
cat(sprintf("[✓] Tree distances saved → %s\n", opt$out))
