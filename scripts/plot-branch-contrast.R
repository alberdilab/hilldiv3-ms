# Plot two four-tip soil clades on a common q = 1 branch-gap scale.
branch_contrast_example <- function(contrast_profiles, genome_tree, sets,
                                    contrib_neutral,
                                    contrib_phylogenetic) {
  focus_clades <- list(
    Actinomycetota = c("LIB-MJ151-D7_04.1", "LIB-MJ262-F7_03.4",
                      "LIB-MJ171-G4_04.8",
                      "LIB-MJ253-D2_01_v_LIB-MJ404-D6_03.4"),
    Pseudomonadota = c("LIB-MJ402-A8_04.4", "LIB-MJ137-F1_02.2",
                       "LIB-MJ201-G3_03.8", "LIB-MJ342-D11_02.11")
  )
  stopifnot(setequal(rownames(contrast_profiles), genome_tree$tip.label),
            all(colSums(contrast_profiles) > 0),
            all(vapply(names(focus_clades), function(name) {
              all(focus_clades[[name]] %in% sets[[name]])
            }, logical(1))))

  profiles <- sweep(contrast_profiles, 2L,
                    colSums(contrast_profiles), "/")
  n_tip <- ape::Ntip(genome_tree)
  edges <- genome_tree$edge
  abundance <- matrix(0, n_tip + genome_tree$Nnode, 2L)
  abundance[seq_len(n_tip), ] <-
    profiles[genome_tree$tip.label, , drop = FALSE]
  postorder <- ape::reorder.phylo(genome_tree, "postorder")
  for (i in seq_len(nrow(postorder$edge))) {
    parent <- postorder$edge[i, 1L]
    child <- postorder$edge[i, 2L]
    abundance[parent, ] <- abundance[parent, ] + abundance[child, ]
  }
  edge_abundance <- abundance[edges[, 2L], , drop = FALSE]
  tree_depth <- sum(genome_tree$edge.length * rowSums(edge_abundance))
  stopifnot(tree_depth > 0)
  values <- 2 * edge_abundance / tree_depth
  xlogx <- function(x) {
    out <- numeric(length(x))
    positive <- x > 0
    out[positive] <- x[positive] * log(x[positive])
    out
  }
  gap <- genome_tree$edge.length * (
    rowMeans(matrix(xlogx(values), ncol = 2L)) -
      xlogx(rowMeans(values))
  )
  gap <- pmax(0, gap)

  child_edges <- function(node) {
    rows <- which(edges[, 1L] == node)
    c(rows, unlist(lapply(edges[rows, 2L], child_edges),
                   use.names = FALSE))
  }
  tip_order <- function(node) {
    if (node <= n_tip) return(node)
    children <- edges[edges[, 1L] == node, 2L]
    unlist(lapply(children, tip_order), use.names = FALSE)
  }
  clade_layout <- function(name, ids) {
    node <- ape::getMRCA(genome_tree, ids)
    stem_edge <- which(edges[, 2L] == node)
    focus_edges <- c(stem_edge, child_edges(node))
    tip_nodes <- tip_order(node)
    stopifnot(length(stem_edge) == 1L, length(focus_edges) == 7L,
              setequal(genome_tree$tip.label[tip_nodes], ids))

    # Include the source-tree stem, which extract.clade() would omit.
    x <- rep(NA_real_, nrow(abundance))
    y <- rep(NA_real_, nrow(abundance))
    x[edges[stem_edge, 1L]] <- 0
    y[tip_nodes] <- rev(seq_along(tip_nodes))
    set_y <- function(k) {
      if (!is.na(y[k])) return(y[k])
      children <- edges[edges[, 1L] == k, 2L]
      child_y <- vapply(children, set_y, numeric(1))
      y[k] <<- mean(range(child_y))
      y[k]
    }
    set_y(node)
    y[edges[stem_edge, 1L]] <- y[node]
    for (i in focus_edges) {
      x[edges[i, 2L]] <- x[edges[i, 1L]] + genome_tree$edge.length[i]
    }

    horizontal <- data.frame(
      clade = name,
      x = x[edges[focus_edges, 1L]],
      xend = x[edges[focus_edges, 2L]],
      y = y[edges[focus_edges, 2L]],
      gap = 1e4 * gap[focus_edges]
    )
    internal_nodes <- unique(edges[focus_edges[-1L], 1L])
    vertical <- do.call(rbind, lapply(internal_nodes, function(k) {
      children <- edges[edges[, 1L] == k, 2L]
      data.frame(clade = name, x = x[k], ymin = min(y[children]),
                 ymax = max(y[children]))
    }))
    tips <- data.frame(
      clade = name, id = genome_tree$tip.label[tip_nodes],
      x = x[tip_nodes], y = y[tip_nodes],
      field = 100 * profiles[genome_tree$tip.label[tip_nodes], 1L],
      grassland = 100 * profiles[genome_tree$tip.label[tip_nodes], 2L]
    )
    totals <- 100 * colSums(profiles[ids, , drop = FALSE])
    list(horizontal = horizontal, vertical = vertical, tips = tips,
         totals = totals, stem_gap = gap[stem_edge],
         focus_edges = focus_edges)
  }
  parts <- Map(clade_layout, names(focus_clades), focus_clades)

  get_rows <- function(field) do.call(rbind, lapply(parts, `[[`, field))
  horizontal <- get_rows("horizontal")
  vertical <- get_rows("vertical")
  tips <- get_rows("tips")
  totals <- data.frame(
    clade = names(parts),
    field = vapply(parts, function(z) unname(z$totals[1L]), numeric(1)),
    grassland = vapply(parts, function(z) unname(z$totals[2L]), numeric(1))
  )
  share <- function(table, name) {
    100 * table$share[match(name, table$set)]
  }
  facet_labels <- vapply(names(parts), function(name) {
    sprintf("%s   |   phylum share: %.1f%% neutral → %.1f%% phylogenetic",
            name, share(contrib_neutral, name),
            share(contrib_phylogenetic, name))
  }, character(1))
  names(facet_labels) <- names(parts)
  for (name in c("horizontal", "vertical", "tips", "totals")) {
    obj <- get(name)
    obj$clade <- factor(facet_labels[as.character(obj$clade)],
                        levels = facet_labels)
    assign(name, obj)
  }
  headers <- data.frame(clade = factor(facet_labels, levels = facet_labels))
  label_x <- max(tips$x) + 0.035
  field_x <- label_x + 0.31
  grassland_x <- field_x + 0.18
  headers$field_x <- field_x
  headers$grassland_x <- grassland_x
  headers$y <- 4.75
  totals$label_x <- label_x
  totals$field_x <- field_x
  totals$grassland_x <- grassland_x
  totals$y <- 0.25

  plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = vertical,
      ggplot2::aes(x = x, xend = x, y = ymin, yend = ymax),
      colour = "#aaa9a5", linewidth = 0.7
    ) +
    ggplot2::geom_segment(
      data = horizontal,
      ggplot2::aes(x = x, xend = xend, y = y, yend = y, colour = gap),
      linewidth = 2.8, lineend = "round"
    ) +
    ggplot2::geom_text(
      data = tips,
      ggplot2::aes(x = label_x, y = y, label = sub("^LIB-", "", id)),
      hjust = 0, size = 3.3
    ) +
    ggplot2::geom_text(
      data = tips,
      ggplot2::aes(x = field_x, y = y, label = sprintf("%.3f", field)),
      hjust = 1, size = 3.3, colour = "#a52929"
    ) +
    ggplot2::geom_text(
      data = tips,
      ggplot2::aes(x = grassland_x, y = y,
                   label = sprintf("%.3f", grassland)),
      hjust = 1, size = 3.3, colour = "#1b6b75"
    ) +
    ggplot2::geom_text(
      data = headers,
      ggplot2::aes(x = field_x, y = y),
      label = "Field (%)", hjust = 1, fontface = "bold",
      colour = "#a52929"
    ) +
    ggplot2::geom_text(
      data = headers,
      ggplot2::aes(x = grassland_x, y = y),
      label = "Grassland (%)", hjust = 1, fontface = "bold",
      colour = "#1b6b75"
    ) +
    ggplot2::geom_text(
      data = totals, ggplot2::aes(x = label_x, y = y),
      label = "Clade total", hjust = 0, fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = totals,
      ggplot2::aes(x = field_x, y = y, label = sprintf("%.3f", field)),
      hjust = 1, fontface = "bold", colour = "#a52929"
    ) +
    ggplot2::geom_text(
      data = totals,
      ggplot2::aes(x = grassland_x, y = y,
                   label = sprintf("%.3f", grassland)),
      hjust = 1, fontface = "bold", colour = "#1b6b75"
    ) +
    ggplot2::scale_colour_gradient(
      low = "#d1d0cd", high = "#a52929",
      limits = c(0, max(horizontal$gap)),
      transform = scales::pseudo_log_trans(sigma = 0.01),
      breaks = c(0, 0.1, 1, 3), labels = c("0", "0.1", "1", "3"),
      name = expression(q == 1~"branch gap"~(x*10^{-4})),
      guide = ggplot2::guide_colourbar(
        barwidth = grid::unit(7, "cm"), ticks = TRUE)
    ) +
    ggplot2::facet_wrap(~clade, ncol = 1L) +
    ggplot2::coord_cartesian(xlim = c(0, grassland_x + 0.03),
                             ylim = c(0, 5.05), clip = "off") +
    ggplot2::labs(title = "Contrasting branch changes in two phyla",
                  subtitle = paste("Same field–grassland soil profiles;",
                                   "shared nonlinear colour scale")) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.spacing = grid::unit(1.1, "lines"),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background = ggplot2::element_rect(fill = "#f4f1ed", colour = NA),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(10, 15, 10, 10)
    )

  list(plot = plot, totals = totals, branch_gap = gap,
       stem_gap = vapply(parts, `[[`, numeric(1), "stem_gap"),
       focus_edges = lapply(parts, `[[`, "focus_edges"))
}
