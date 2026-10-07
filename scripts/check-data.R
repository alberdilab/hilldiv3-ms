#!/usr/bin/env Rscript

# Check the staged files independently of the source archive.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript")
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
if (!requireNamespace("data.table", quietly = TRUE) ||
    !requireNamespace("ape", quietly = TRUE)) {
  stop("Install the R packages data.table and ape")
}

samples <- data.table::fread(
  file.path(root, "data", "sample_metadata.tsv"),
  data.table = FALSE, check.names = FALSE
)
reps <- data.table::fread(
  file.path(root, "data", "representative_metadata.tsv"),
  data.table = FALSE, check.names = FALSE
)
abundance <- data.table::fread(
  cmd = paste(
    "xz -dc",
    shQuote(file.path(root, "data", "mag_relative_abundance.tsv.xz"))
  ),
  check.names = FALSE, data.table = FALSE
)
tree <- ape::read.tree(file.path(root, "data", "representative_tree.nwk"))

stopifnot(
  nrow(samples) == 360L,
  nrow(reps) == 5518L,
  nrow(abundance) == 5518L,
  ncol(abundance) == 361L,
  identical(abundance$genome_id, reps$genome_id),
  identical(names(abundance)[-1L], samples$sample_id),
  !anyDuplicated(samples$sample_id),
  !anyDuplicated(reps$genome_id),
  !anyDuplicated(reps$secondary_cluster),
  length(table(samples$habitat_class)) == 8L,
  all(as.integer(table(samples$habitat_class)) == 45L),
  ape::Ntip(tree) == 5518L,
  setequal(tree$tip.label, reps$genome_id),
  !anyDuplicated(tree$tip.label),
  length(tree$edge.length) == nrow(tree$edge),
  all(is.finite(tree$edge.length) & tree$edge.length > 0)
)
values <- as.matrix(abundance[-1L])
storage.mode(values) <- "double"
totals <- colSums(values)
stopifnot(
  all(is.finite(values) & values >= 0),
  all(totals > 0.99 & totals < 1.01)
)
cat(
  "OK: ", nrow(values), " genomes, ", ncol(values),
  " samples, ", ape::Ntip(tree), " tree tips; profile sums ",
  sprintf("%.6f", min(totals)), "–", sprintf("%.6f", max(totals)),
  ".\n", sep = ""
)
