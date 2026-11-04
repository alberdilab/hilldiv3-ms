#!/usr/bin/env Rscript

# Distil KEGG ortholog annotations into genome-inferred functional traits
# (GIFTs) with distillR, for the biosynthesis (B) and degradation (D) domains.
# Usage: Rscript scripts/distill-functional-traits.R kegg_hits_primary.tsv.gz
# The input has one row per annotated gene and columns mag, gene, ..., ko.

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L || length(args) != 1L) {
  stop("Usage: Rscript scripts/distill-functional-traits.R <KEGG hits TSV.gz>")
}
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
if (!requireNamespace("data.table", quietly = TRUE) ||
    !requireNamespace("distillR", quietly = TRUE)) {
  stop("Install the R packages data.table and distillR")
}

reps <- data.table::fread(
  file.path(root, "data", "representative_metadata.tsv"),
  select = "genome_id", data.table = FALSE
)
genomes <- reps$genome_id

hits <- data.table::fread(
  cmd = paste("gzip -dc", shQuote(args[1L])),
  select = c("mag", "ko", "ec"), data.table = FALSE
)
stopifnot(all(hits$mag %in% genomes))

# The GIFT database defines bundles with KO identifiers and EC numbers.
# distill() reads a KO as K##### and EC numbers only inside a KEGG-style
# "[EC:a b c]" tag, taking one match per cell. Give each KO a row for itself
# and a second row carrying its EC numbers in that tag.
pairs <- unique(hits[c("ko", "ec")])
has_ec <- pairs$ec != "-"
ids <- rbind(
  data.frame(ko = pairs$ko, id = pairs$ko),
  data.frame(
    ko = pairs$ko[has_ec],
    id = paste0("[EC:", gsub(";", " ", pairs$ec[has_ec], fixed = TRUE), "]")
  )
)
hits <- unique(merge(unique(hits[c("mag", "ko")]), ids, by = "ko")[c("mag", "id")])

# Keep the biosynthesis (B) and degradation (D) bundles; drop structure (S).
gift_db <- distillR::GIFT_db
gift_db <- gift_db[gift_db$Domain %in% c("Biosynthesis", "Degradation"), ]

bundles <- distillR::distill(
  hits, GIFT_db = gift_db, genomecol = 1, annotcol = 2,
  stats = TRUE, verbosity = FALSE
)
elements <- distillR::to.elements(bundles, GIFT_db = gift_db)
functions <- distillR::to.functions(elements, GIFT_db = gift_db)

# MAGs without any KEGG hit are absent from distill(); report them as zeros.
complete <- function(x) {
  x <- as.matrix(x)
  out <- matrix(0, nrow = length(genomes), ncol = ncol(x),
                dimnames = list(genomes, colnames(x)))
  out[rownames(x), ] <- x
  out
}
write_table <- function(x, name) {
  x <- data.frame(genome_id = rownames(x), x, check.names = FALSE,
                  row.names = NULL)
  path <- file.path(root, "data", name)
  con <- pipe(paste("xz -9 >", shQuote(path)), "w")
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
  invisible(path)
}
missing <- setdiff(genomes, rownames(bundles))
message(length(missing), " of ", length(genomes), " MAGs have no KEGG hits")

write_table(complete(bundles), "mag_gift_bundles.tsv.xz")
write_table(complete(elements), "mag_gift_elements.tsv.xz")
write_table(complete(functions), "mag_gift_functions.tsv.xz")
