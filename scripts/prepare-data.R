#!/usr/bin/env Rscript

# Build the staged Microflora Danica tables from pinned public source files.
# Usage: Rscript scripts/prepare-data.R [MAG metadata] [abundance TSV.xz] [tree]

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L || length(args) > 3L) {
  stop("Usage: Rscript scripts/prepare-data.R [MAG metadata] [abundance TSV.xz] [tree]")
}
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
source_dir <- file.path(root, "sources")
paths <- c(
  file.path(source_dir, "mags_shallow_all.tsv"),
  file.path(source_dir, "MFD_SRnodrep_tax_relative_abundance.tsv.xz"),
  file.path(source_dir, "gtdb-shallow.tree")
)
if (length(args)) paths[seq_along(args)] <- args
if (!all(file.exists(paths))) stop("Missing source file; run scripts/fetch-source-data.sh")
if (!requireNamespace("data.table", quietly = TRUE) ||
    !requireNamespace("ape", quietly = TRUE)) {
  stop("Install the R packages data.table and ape")
}

samples <- data.table::fread(
  file.path(root, "data", "sample_metadata.tsv"), data.table = FALSE,
  check.names = FALSE
)
reps <- data.table::fread(
  file.path(root, "data", "representative_metadata.tsv"), data.table = FALSE,
  check.names = FALSE
)
if (nrow(samples) != 360L || nrow(reps) != 5518L ||
    anyDuplicated(samples$sample_id) ||
    anyDuplicated(samples$abundance_column) ||
    anyDuplicated(reps$genome_id) ||
    anyDuplicated(reps$secondary_cluster)) {
  stop("The locked sample or representative panel is invalid")
}

mags <- data.table::fread(
  paths[1], select = c("bin", "secondary_cluster"), data.table = FALSE
)
if (nrow(mags) != 19253L || anyDuplicated(mags$bin) ||
    !setequal(reps$genome_id, intersect(reps$genome_id, mags$bin))) {
  stop("The MAG metadata differs from the locked panel")
}

profiles <- data.table::fread(
  cmd = paste("xz -dc", shQuote(normalizePath(paths[2]))),
  select = c("clade_name", samples$abundance_column),
  check.names = FALSE
)
if (!identical(names(profiles)[-1L], samples$abundance_column)) {
  stop("Selected abundance columns differ from the sample metadata")
}
strain <- profiles[grepl("|t__", profiles$clade_name, fixed = TRUE)]
strain[, genome_id := sub("[.]fa$", "", sub("^.*[|]t__", "", clade_name))]
cluster_by_bin <- stats::setNames(mags$secondary_cluster, mags$bin)
strain[, secondary_cluster := unname(cluster_by_bin[genome_id])]
if (nrow(strain) != 19246L || anyNA(strain$secondary_cluster) ||
    anyDuplicated(strain$genome_id)) {
  stop("The published strain rows do not map uniquely to MAG clusters")
}

sample_columns <- samples$abundance_column
collapsed <- strain[
  , lapply(.SD, sum), by = secondary_cluster, .SDcols = sample_columns
]
index <- match(collapsed$secondary_cluster, reps$secondary_cluster)
if (anyNA(index) || anyDuplicated(index)) {
  stop("Collapsed profiles have unknown or duplicate species clusters")
}
abundance <- matrix(
  0, nrow = nrow(reps), ncol = nrow(samples),
  dimnames = list(reps$genome_id, samples$sample_id)
)
abundance[index, ] <- as.matrix(collapsed[, ..sample_columns]) / 100
if (any(!is.finite(abundance)) || any(abundance < 0) ||
    any(colSums(abundance) < 0.99 | colSums(abundance) > 1.01)) {
  stop("Relative abundances contain invalid values or incomplete profiles")
}

tree <- ape::read.tree(paths[3])
tree$tip.label <- gsub("(^'|'$)", "", tree$tip.label)
if (anyDuplicated(tree$tip.label) ||
    !all(reps$genome_id %in% tree$tip.label)) {
  stop("The source tree does not contain every representative exactly once")
}
tree <- ape::keep.tip(tree, reps$genome_id)
if (!setequal(tree$tip.label, reps$genome_id) ||
    is.null(tree$edge.length) ||
    any(!is.finite(tree$edge.length) | tree$edge.length <= 0)) {
  stop("The representative tree is invalid")
}

output <- file.path(root, "data", "mag_relative_abundance.tsv.xz")
con <- xzfile(output, open = "wt", compression = 9)
tryCatch(
  write.table(
    data.frame(genome_id = rownames(abundance), abundance, check.names = FALSE),
    file = con, sep = "\t", quote = FALSE, row.names = FALSE
  ),
  finally = close(con)
)
ape::write.tree(tree, file = file.path(root, "data", "representative_tree.nwk"))

cat(
  "Staged ", nrow(abundance), " representatives, ", ncol(abundance),
  " samples, and a ", ape::Ntip(tree), "-tip tree.\n", sep = ""
)
