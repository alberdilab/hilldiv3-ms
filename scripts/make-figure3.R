#!/usr/bin/env Rscript
# Fig. 3 of the hilldiv3 paper: Microflora Danica worked example (six panels).
#
# Reuses the data preparation and hilldiv3 calls of the webbook chapters
# (01-data_overview, 02-diversity, 03-pairwise_beta, 04-taxonomic_shares) and
# draws panels a-f as svglite SVG at their final width (87.5 mm), then composes
# and exports the figure with compose.py and Inkscape.
#
# Run from anywhere:  Rscript scripts/make-figure3.R
# Expensive results (hilldiv3 calls) are cached as RDS in data/cache/ (git-ignored);
# delete a cache file to recompute it.
# Environment overrides: FIG_DIR (default ../hilldiv3/paper/figures), INKSCAPE.

suppressMessages({
  library(ggplot2)
  library(svglite)
})
stopifnot(utils::packageVersion("hilldiv3") >= "3.1.0.9000")

# ---- paths (project-relative) ----------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(script)) normalizePath(file.path(dirname(script), "..")) else getwd()
setwd(root)
fig_dir <- Sys.getenv("FIG_DIR", normalizePath("../hilldiv3/paper/figures", mustWork = FALSE))
raw_dir <- file.path(fig_dir, "raw"); src_dir <- file.path(fig_dir, "src")
final_dir <- file.path(fig_dir, "final"); cache_dir <- "data/cache"
for (d in c(raw_dir, src_dir, final_dir, cache_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

cached <- function(name, expr) {
  f <- file.path(cache_dir, paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  message("computing ", name, " ...")
  val <- expr
  saveRDS(val, f)
  val
}
plain <- function(x) { class(x) <- "data.frame"; x }

# ---- data (as in 00-data_preparation.Rmd) -----------------------------------
sample_metadata <- data.table::fread("data/sample_metadata.tsv", data.table = FALSE)
genome_metadata <- data.table::fread("data/representative_metadata.tsv", data.table = FALSE)
genome_abundance <- data.table::fread(cmd = "xz -dc data/mag_relative_abundance.tsv.xz",
                                      data.table = FALSE, check.names = FALSE)
genome_tree <- ape::read.tree("data/representative_tree.nwk")
stopifnot(identical(genome_abundance$genome_id, genome_metadata$genome_id),
          identical(names(genome_abundance)[-1L], sample_metadata$sample_id),
          setequal(genome_tree$tip.label, genome_metadata$genome_id))

habitat_levels <- unique(sample_metadata$habitat_label)
overview_abundance <- as.matrix(genome_abundance[-1L])
storage.mode(overview_abundance) <- "double"
rownames(overview_abundance) <- genome_abundance$genome_id
sample_ids <- sample_metadata$sample_id
habitat_short <- c("Field soil", "Grassland soil", "Forest soil", "Urban greenspace",
                   "Bog/fen soil", "Freshwater sediment", "Saltwater sediment", "Wastewater")
stopifnot(length(habitat_short) == length(habitat_levels))
short_of <- function(h) factor(habitat_short[match(h, habitat_levels)], levels = habitat_short)

# ---- palettes ---------------------------------------------------------------
# Habitats: Okabe-Ito (blue, sky blue, pink, black, grey) + Tol (green, sand, wine),
# none close to the neutral/phylogenetic colours below.
hab_cols <- c("Field soil" = "#0072B2", "Grassland soil" = "#56B4E9",
              "Forest soil" = "#117733", "Urban greenspace" = "#DDCC77",
              "Bog/fen soil" = "#882255", "Freshwater sediment" = "#CC79A7",
              "Saltwater sediment" = "#000000", "Wastewater" = "#999999")
type_cols <- c(Neutral = "#66C2A5", Phylogenetic = "#FC8D62")
phylum_cols <- c(Pseudomonadota = "#332288", Actinomycetota = "#DDCC77",
                 Acidobacteriota = "#882255", Thermoproteota = "#88CCEE", Other = "#BBBBBB")

# ---- theme: 87.5 mm panels, 6 pt text ---------------------------------------
panel_w <- 87.5; panel_h <- 53   # mm
fam <- "Helvetica"
theme_fig <- function() {
  theme_classic(base_size = 6.5, base_family = fam) %+replace%
    theme(axis.line = element_line(linewidth = 0.25, colour = "grey20"),
          axis.ticks = element_line(linewidth = 0.25, colour = "grey20"),
          axis.text = element_text(size = 6, colour = "grey10"),
          axis.title = element_text(size = 6.5),
          strip.background = element_blank(),
          strip.text = element_text(size = 6.5, face = "bold"),
          legend.text = element_text(size = 6), legend.title = element_text(size = 6.5),
          legend.key.size = unit(2.6, "mm"),
          panel.spacing = unit(2.5, "mm"),
          plot.margin = margin(4.5, 2, 1, 1, "mm"))   # top margin leaves room for the panel letter
}
save_panel <- function(p, name, h = panel_h) {
  f <- file.path(raw_dir, paste0(name, ".svg"))
  ggsave(f, p, width = panel_w, height = h, units = "mm", device = svglite)
  message("wrote ", f)
}

# ---- A: map (01-data_overview.Rmd, overview_map; one map, coloured by habitat)
d <- sample_metadata
d$Habitat <- short_of(d$habitat_label)
stopifnot(all(is.finite(d$latitude)), all(is.finite(d$longitude)))
denmark <- map_data("world", region = "Denmark")
pA <- ggplot() +
  geom_polygon(data = denmark, aes(long, lat, group = group),
               fill = "#f0eeeb", colour = "#a9aaa8", linewidth = 0.2) +
  geom_point(data = d, aes(longitude, latitude, colour = Habitat), size = 0.8, alpha = 0.85, stroke = 0) +
  scale_colour_manual(values = hab_cols, name = "Habitat (panels a, c-e)") +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  # x range extended eastwards so that the shared habitat legend sits beside the map
  coord_quickmap(xlim = c(7.9, 19.2), ylim = c(54.5, 57.8), expand = FALSE) +
  scale_x_continuous(breaks = c(8, 10, 12, 14)) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)") +
  theme_fig() +
  theme(legend.position = c(0.99, 0.5), legend.justification = c(1, 0.5),
        legend.key.height = unit(2.8, "mm"), legend.spacing.y = unit(0, "mm"),
        legend.background = element_blank(), legend.key = element_blank(),
        axis.line = element_blank(), panel.border = element_rect(fill = NA, linewidth = 0.25))
save_panel(pA, "fig3a")

# ---- B: alpha diversity (02-diversity.Rmd, hilldiv) -------------------------
alpha_neutral <- cached("alpha_neutral",
  plain(hilldiv3::hilldiv(overview_abundance, q = c(0, 1, 2), type = "neutral")))
alpha_phylogenetic <- cached("alpha_phylogenetic",
  plain(hilldiv3::hilldiv(overview_abundance, q = c(0, 1, 2), tree = genome_tree,
                          type = "phylogenetic", reference = "pool")))
alpha_neutral$type <- "Neutral"; alpha_phylogenetic$type <- "Phylogenetic"
alpha <- rbind(alpha_neutral, alpha_phylogenetic)
alpha$habitat <- sample_metadata$habitat_label[match(alpha$sample, sample_ids)]
stopifnot(!anyNA(alpha$habitat), all(is.finite(alpha$value)), all(alpha$value > 0))
# Chosen for compactness: q = 1 only, one facet per diversity type (free x scale).
a1 <- alpha[alpha$q == 1, ]
a1$Habitat <- factor(habitat_short[match(a1$habitat, habitat_levels)], levels = rev(habitat_short))
a1$type <- factor(a1$type, levels = c("Neutral", "Phylogenetic"))
pB <- ggplot(a1, aes(value, Habitat, fill = type)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.25, alpha = 0.55,
               colour = "grey20", median.linewidth = 0.4) +
  geom_point(position = position_jitter(height = 0.15, width = 0, seed = 42),
             size = 0.35, alpha = 0.5, stroke = 0, colour = "grey25") +
  scale_fill_manual(values = type_cols) +
  facet_wrap(~ type, scales = "free_x") +
  labs(x = "Hill diversity, q = 1", y = NULL) +
  theme_fig() + theme(legend.position = "none")
save_panel(pB, "fig3b")

# ---- C: diversity profiles (02-diversity.Rmd, hillprof) ---------------------
habitat_profiles <- sapply(habitat_levels, function(h)
  rowMeans(overview_abundance[, sample_metadata$habitat_label == h, drop = FALSE]))
rownames(habitat_profiles) <- rownames(overview_abundance)
profile_neutral <- cached("profile_neutral",
  plain(hilldiv3::hillprof(habitat_profiles, q = seq(0, 2, by = 0.25), type = "neutral")))
profile_phylogenetic <- cached("profile_phylogenetic",
  plain(hilldiv3::hillprof(habitat_profiles, q = seq(0, 2, by = 0.25), tree = genome_tree,
                           type = "phylogenetic", reference = "pool")))
profile_neutral$type <- "Neutral"; profile_phylogenetic$type <- "Phylogenetic"
profiles <- rbind(profile_neutral, profile_phylogenetic)
stopifnot(nrow(profiles) == 2L * 9L * length(habitat_levels),
          all(is.finite(profiles$value)), all(profiles$value > 0))
profiles$Habitat <- short_of(profiles$sample)
profiles$type <- factor(profiles$type, levels = c("Neutral", "Phylogenetic"))
pC <- ggplot(profiles, aes(q, value, colour = Habitat)) +
  geom_line(linewidth = 0.45) +
  scale_colour_manual(values = hab_cols) +
  scale_y_log10() + annotation_logticks(sides = "l", linewidth = 0.2, short = unit(0.6, "mm"),
                                        mid = unit(1, "mm"), long = unit(1.4, "mm")) +
  facet_wrap(~ type, scales = "free_y") +
  labs(x = "Diversity order (q)", y = "Hill diversity (log scale)") +
  theme_fig() + theme(legend.position = "none")
save_panel(pC, "fig3c")

# ---- D: PCoA of hillpair, phylogenetic q = 1 (03-pairwise_beta.Rmd) --------
# Only q = 1 is computed here (values are identical to the q = c(0,1,2) call).
pair_phylo <- cached("pair_phylogenetic_q1_C",
  plain(hilldiv3::hillpair(overview_abundance, q = 1, metric = "C", tree = genome_tree,
                           type = "phylogenetic", out = "tibble")))
first <- match(pair_phylo$first, sample_ids); second <- match(pair_phylo$second, sample_ids)
stopifnot(nrow(pair_phylo) == choose(length(sample_ids), 2L), !anyNA(first), !anyNA(second),
          all(first < second), !anyDuplicated(paste(first, second)),
          all(is.finite(pair_phylo$value)), all(pair_phylo$value >= -1e-10 & pair_phylo$value <= 1 + 1e-10))
dist_m <- matrix(0, length(sample_ids), length(sample_ids), dimnames = list(sample_ids, sample_ids))
dist_m[cbind(first, second)] <- pair_phylo$value
dist_m[cbind(second, first)] <- pair_phylo$value
ord <- stats::cmdscale(stats::as.dist(dist_m), k = 2, eig = TRUE, add = TRUE)  # Cailliez correction
variance <- 100 * ord$eig[1:2] / sum(pmax(ord$eig, 0))
ordination <- data.frame(PCoA1 = ord$points[, 1], PCoA2 = ord$points[, 2],
                         Habitat = short_of(sample_metadata$habitat_label))
message(sprintf("PCoA axes: %.1f%%, %.1f%%", variance[1], variance[2]))
pD <- ggplot(ordination, aes(PCoA1, PCoA2, colour = Habitat)) +
  geom_point(size = 0.8, alpha = 0.8, stroke = 0) +
  scale_colour_manual(values = hab_cols) +
  labs(title = "Phylogenetic, q = 1 (C dissimilarity)",
       x = sprintf("PCoA 1 (%.1f%%)", variance[1]), y = sprintf("PCoA 2 (%.1f%%)", variance[2])) +
  theme_fig() + theme(legend.position = "none",
                      plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5))
save_panel(pD, "fig3d")

# ---- shared sets (04-taxonomic_shares.Rmd) ----------------------------------
phylum <- sub("^p__", "", genome_metadata$phylum)
focus_phyla <- c("Pseudomonadota", "Actinomycetota", "Acidobacteriota", "Thermoproteota")
contrast_habitats <- c("Agricultural field soil", "Natural grassland soil")
share_sets <- lapply(focus_phyla, function(n) genome_metadata$genome_id[phylum == n])
names(share_sets) <- focus_phyla
share_sets$Other <- genome_metadata$genome_id[!phylum %in% focus_phyla]
stopifnot(sum(lengths(share_sets)) == nrow(genome_metadata))
soil_samples <- sample_metadata$habitat_label %in% contrast_habitats

# ---- E: hillshare, Gaiellales vs Streptosporangiales vs Actinomycetota -------
actino_orders <- c("Gaiellales", "Streptosporangiales")
actino_order_sets <- lapply(actino_orders, function(n)
  genome_metadata$genome_id[phylum == "Actinomycetota" & genome_metadata$order == paste0("o__", n)])
names(actino_order_sets) <- actino_orders
actino_order_sets$Other <- setdiff(genome_metadata$genome_id, unlist(actino_order_sets, use.names = FALSE))
stopifnot(all(lengths(actino_order_sets[actino_orders]) > 0L),
          sum(lengths(actino_order_sets)) == nrow(genome_metadata))
# Phylum share: neutral q = 1 on the 90 soil samples (same values as the chapter's all-sample call).
phylum_shares <- cached("share_phylum_neutral_q1_soil",
  plain(hilldiv3::hillshare(overview_abundance[, soil_samples, drop = FALSE], q = 1,
                            type = "neutral", sets = share_sets)))
order_shares <- cached("share_orders_neutral_q1_soil",
  plain(hilldiv3::hillshare(overview_abundance[, soil_samples, drop = FALSE], q = 1,
                            type = "neutral", sets = actino_order_sets)))
for (x in list(phylum_shares, order_shares)) stopifnot(all(is.finite(x$share)))
phylum_shares$habitat <- sample_metadata$habitat_label[match(phylum_shares$sample, sample_ids)]
order_shares$habitat <- sample_metadata$habitat_label[match(order_shares$sample, sample_ids)]
actino <- rbind(phylum_shares[phylum_shares$set == "Actinomycetota", c("sample", "habitat", "set", "share")],
                order_shares[order_shares$set %in% actino_orders, c("sample", "habitat", "set", "share")])
actino$Taxon <- factor(actino$set, levels = c("Actinomycetota", actino_orders))
actino$Habitat <- factor(actino$habitat, levels = contrast_habitats, labels = c("Field", "Grassland"))
stopifnot(nrow(actino) == 3L * 2L * 45L, !anyNA(actino$Habitat))
actino_medians <- aggregate(share ~ set + habitat, actino, median)
print(actino_medians)
e_cols <- c(Field = unname(hab_cols["Field soil"]), Grassland = unname(hab_cols["Grassland soil"]))
pE <- ggplot(actino, aes(Habitat, share, fill = Habitat)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.55, linewidth = 0.25,
               median.linewidth = 0.4) +
  geom_point(position = position_jitter(width = 0.13, height = 0, seed = 42),
             size = 0.35, alpha = 0.5, stroke = 0) +
  scale_fill_manual(values = e_cols) +
  scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x),
                     expand = expansion(mult = c(0.03, 0.07))) +
  facet_wrap(~ Taxon, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = "Neutral share (q = 1),\nfield vs grassland soils") +
  theme_fig() + theme(legend.position = "none", strip.text = element_text(size = 5.5, face = "bold"))
save_panel(pE, "fig3e")

# ---- F: hillcontrib stacked bars --------------------------------------------
contrast_profiles <- habitat_profiles[, contrast_habitats, drop = FALSE]
contrib_neutral <- cached("contrib_neutral", plain(hilldiv3::hillcontrib(
  contrast_profiles, q = 1, metric = "C", type = "neutral", sets = share_sets)))
contrib_phylogenetic <- cached("contrib_phylogenetic", plain(hilldiv3::hillcontrib(
  contrast_profiles, q = 1, metric = "C", tree = genome_tree, type = "phylogenetic", sets = share_sets)))
stopifnot(abs(sum(contrib_neutral$share) - 1) < 1e-10, abs(sum(contrib_phylogenetic$share) - 1) < 1e-10)
contributions <- rbind(transform(contrib_neutral, Type = "Neutral"),
                       transform(contrib_phylogenetic, Type = "Phylogenetic"))
contributions$Type <- factor(contributions$Type, levels = c("Neutral", "Phylogenetic"))
contributions$set <- factor(contributions$set, levels = names(share_sets))
message(sprintf("C dissimilarity: neutral %.4f, phylogenetic %.4f",
                sum(contrib_neutral$contribution), sum(contrib_phylogenetic$contribution)))
print(contributions[c("Type", "set", "contribution", "share")])
tot <- aggregate(contribution ~ Type, contributions, sum)
pF <- ggplot(contributions, aes(Type, contribution, fill = set)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.2) +
  geom_text(data = tot, aes(Type, contribution, label = sprintf("%.3f", contribution), fill = NULL),
            hjust = -0.15, size = 6 / .pt, family = fam) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 0.34), breaks = seq(0, 0.3, by = 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(limits = rev(levels(contributions$Type))) +
  scale_fill_manual(values = phylum_cols) +
  labs(x = NULL, y = "Contribution to C dissimilarity (q = 1)\nfield vs grassland soil", fill = NULL) +
  guides(fill = guide_legend(nrow = 2)) +
  theme_fig() + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
                      legend.box.margin = margin(-1, 0, 0, 0), legend.key.size = unit(2.4, "mm"))
save_panel(pF, "fig3f")

# ---- compose and export -----------------------------------------------------
panels <- file.path(raw_dir, sprintf("fig3%s.svg", letters[1:6]))
composed <- file.path(src_dir, "fig3.svg")
status <- system2("python3", c("-I", Sys.getenv("COMPOSE_PY", "../hilldiv3/paper/figures/compose.py"), composed,
                               "--width", "180", "--cols", "2", panels))
stopifnot(status == 0L)
ink <- Sys.getenv("INKSCAPE", "/Applications/Inkscape.app/Contents/MacOS/inkscape")
if (file.exists(ink)) {
  system2(ink, c(composed, "--export-type=pdf", "--export-text-to-path",
                 paste0("--export-filename=", file.path(final_dir, "fig3.pdf"))), stderr = FALSE)
  system2(ink, c(composed, "--export-type=png", "--export-dpi=600", "--export-background=white",
                 paste0("--export-filename=", file.path(final_dir, "fig3.png"))), stderr = FALSE)
} else message("Inkscape not found; skipped PDF/PNG export")
message("done")
