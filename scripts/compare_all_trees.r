#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(phangorn)
  library(dplyr)
  library(TreeDist)
  library(tibble)
})

option_list <- list(
  make_option(c("--dataset-dir"), type = "character", help = "results/<dataset> directory"),
  make_option(c("--out"), type = "character", help = "Output CSV path"),
  make_option(
    c("--ref-name"),
    type = "character",
    default = "tree.nwk",
    help = "Reference tree filename inside dataset-dir (default: tree.nwk)"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$`dataset-dir`) || is.null(opt$out)) {
  stop("Missing required args: --dataset-dir and --out")
}

dataset_dir <- opt$`dataset-dir`
ref_name <- opt$`ref-name`
ref_path <- file.path(dataset_dir, ref_name)

if (!file.exists(ref_path)) {
  stop(sprintf("Reference tree not found: %s", ref_path))
}

BIG_PENALTY <- 1e9

metric_names <- c(
  "RF_raw", "RF_normalized", "Weighted_RF", "BranchScore_KF",
  "JRF", "JRF_normalized", "Nye", "Nye_normalized",
  "SMI", "SMI_normalized"
)

normalize_branch_lengths <- function(tr) {
  if (is.null(tr$edge.length)) return(tr)
  el <- tr$edge.length
  el[is.na(el)] <- 0
  el[el < 0] <- 0
  s <- sum(el)
  tr$edge.length <- if (s > 0) el / s else el
  tr
}

safe_metric <- function(expr_fun) {
  tryCatch(
    {
      val <- expr_fun()
      if (length(val) == 0 || !is.finite(val[1])) BIG_PENALTY else val[1]
    },
    error = function(e) BIG_PENALTY
  )
}

preprocess_trees <- function(t1_path, t2_path) {
  tryCatch(
    {
      t1 <- read.tree(t1_path)
      t2 <- read.tree(t2_path)

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
      message(sprintf("[!] Preprocessing failed for %s vs %s: %s", t1_path, t2_path, e$message))
      NULL
    }
  )
}

compare_pair <- function(ref_path, tree_path) {
  prep <- preprocess_trees(ref_path, tree_path)

  if (is.null(prep)) {
    metrics <- setNames(as.list(rep(BIG_PENALTY, length(metric_names))), metric_names)
  } else {
    tree1 <- prep$tree1
    tree2 <- prep$tree2

    metrics <- list(
      RF_raw         = safe_metric(function() RF.dist(tree1, tree2, normalize = FALSE)[1]),
      RF_normalized  = safe_metric(function() RF.dist(tree1, tree2, normalize = TRUE)[1]),
      Weighted_RF    = safe_metric(function() wRF.dist(tree1, tree2)[1]),
      BranchScore_KF = safe_metric(function() KF.dist(tree1, tree2)),
      JRF            = safe_metric(function() JaccardRobinsonFoulds(tree1, tree2, normalize = FALSE)[1]),
      JRF_normalized = safe_metric(function() JaccardRobinsonFoulds(tree1, tree2, normalize = TRUE)[1]),
      Nye            = safe_metric(function() NyeSimilarity(tree1, tree2, normalize = FALSE, similarity=FALSE)[1]),
      Nye_normalized = safe_metric(function() NyeSimilarity(tree1, tree2, normalize = TRUE, similarity=FALSE)[1]),
      SMI            = safe_metric(function() SharedPhylogeneticInfo(tree1, tree2)[1]),
      SMI_normalized = safe_metric(function() SharedPhylogeneticInfo(tree1, tree2, normalize = TRUE)[1])
    )
  }

  c(list(tree = basename(tree_path)), metrics)
}

nwk_files <- list.files(dataset_dir, pattern = "\\.nwk$", full.names = TRUE, recursive = FALSE)
others <- nwk_files[basename(nwk_files) != ref_name]

if (length(others) == 0) {
  out_df <- tibble(
    tree = character(0),
    RF_raw = numeric(0),
    RF_normalized = numeric(0),
    Weighted_RF = numeric(0),
    BranchScore_KF = numeric(0),
    JRF = numeric(0),
    JRF_normalized = numeric(0),
    Nye = numeric(0),
    Nye_normalized = numeric(0),
    SMI = numeric(0),
    SMI_normalized = numeric(0)
  )
} else {
  rows <- lapply(others, function(p) compare_pair(ref_path, p))
  out_df <- bind_rows(lapply(rows, as_tibble))

  cols <- c(
    "tree", "RF_raw", "RF_normalized", "Weighted_RF", "BranchScore_KF",
    "JRF", "JRF_normalized", "Nye", "Nye_normalized",
    "SMI", "SMI_normalized"
  )
  out_df <- out_df[, intersect(cols, names(out_df)), drop = FALSE]

  if ("RF_normalized" %in% names(out_df)) {
    out_df <- out_df %>% arrange(.data$RF_normalized)
  }
}

write.csv(out_df, opt$out, row.names = FALSE)
cat(sprintf("[✓] Tree distances saved → %s\n", opt$out))
