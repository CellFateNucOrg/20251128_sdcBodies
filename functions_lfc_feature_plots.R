## Functions for plotting log2 fold changes by genomic features.
## Used by plot_lfc_by_features.R. Can also be pasted into functions_exploratory_analysis.R.
##
## Requires: GenomicRanges, rtracklayer, ggplot2, dplyr, ggpubr, rstatix
## and tableToGranges() from functions_exploratory_analysis.R

# helpers ------------------------------------------------------------------

# genes that were actually tested (non-zero LFC and an adjusted p value)
tested_genes <- function(df) {
  df[!is.na(df$padj) & df$log2FoldChange != 0, ]
}

# y axis limits depend on whether LFC shrinkage was used
ylimits_for <- function(runName, shrink = c(-0.5, 0.5), noShrink = c(-1.5, 1.5)) {
  if (grepl("noShrink", runName)) noShrink else shrink
}

# save a plot; an error gives a warning instead of stopping the loops. Returns TRUE if saved.
save_plot <- function(plot, file, width = 19, height = 13) {
  tryCatch({
    ggsave(file, plot, width = width, height = height, units = "cm")
    TRUE
  }, error = function(e) {
    warning("could not make ", basename(file), ": ", conditionMessage(e), call. = FALSE)
    FALSE
  })
}

# write a list of pages (each list(plot, width, height) in cm) to one multipage pdf;
# the page size is the largest of the pages
save_multipage <- function(pages, file) {
  w <- max(sapply(pages, `[[`, "width")) / 2.54
  h <- max(sapply(pages, `[[`, "height")) / 2.54
  pdf(file, width = w, height = h, onefile = TRUE)
  on.exit(dev.off())
  for (pg in pages) {
    tryCatch(print(pg$plot), error = function(e)
      warning("page skipped in ", basename(file), ": ", conditionMessage(e), call. = FALSE))
  }
  invisible(file)
}

# import a bed file and give it UCSC names and the genome's seqinfo (needed by trim())
load_sites <- function(file, si) {
  x <- import(file)
  seqlevelsStyle(x) <- "UCSC"
  seqlevels(x, pruning.mode = "coarse") <- seqlevels(si)
  seqinfo(x) <- si
  x
}

# genome tiles, flagged if they contain the centre of a target site
make_tiles <- function(seqlengths, binWidth, targets) {
  tiles <- tileGenome(seqlengths, tilewidth = binWidth, cut.last.tile.in.chrom = TRUE)
  tiles$target <- overlapsAny(tiles, resize(targets, width = 1, fix = "center"), ignore.strand = TRUE)
  tiles
}

# lacO sites carried by each contrast's strain, from a `lacO` column (columns id and lacO,
# e.g. in contrasts.csv). Values are site names from the name column of the bed file;
# several sites separated by ";"; empty = no site.
# Returns a list: contrast id -> site names, or NULL if there is no lacO column
# (then every site is treated as present in every contrast).
sites_per_contrast <- function(tab, targets) {
  if (!"lacO" %in% colnames(tab)) return(NULL)
  if (is.null(targets$name) || all(is.na(targets$name)))
    stop("there is a lacO column, but the target bed file has no site names (4th column)")
  s <- lapply(as.character(tab$lacO), function(x) {
    x <- trimws(unlist(strsplit(ifelse(is.na(x), "", x), ";")))
    x[x != ""]
  })
  names(s) <- tab$id
  unknown <- setdiff(unlist(s), targets$name)
  if (length(unknown) > 0)
    warning("lacO sites not in the bed file: ", paste(unknown, collapse = ", "), call. = FALSE)
  s
}

# id/lacO table for sites_per_contrast(), from a `lacO` column in the samplesheet: for each
# contrast, the lacO value(s) of the samples in its target group (column `variable`).
# Returns a table without a lacO column if the samplesheet has none.
lacO_from_samplesheet <- function(contrasts, samplesheet) {
  if (!"lacO" %in% colnames(samplesheet)) return(data.frame(id = contrasts$id))
  lac <- sapply(seq_len(nrow(contrasts)), function(i) {
    v <- contrasts$variable[i]
    if (!v %in% colnames(samplesheet)) stop("samplesheet has no column '", v, "' (contrast ", contrasts$id[i], ")")
    x <- as.character(samplesheet$lacO[as.character(samplesheet[[v]]) == as.character(contrasts$target[i])])
    if (length(x) == 0) warning("no samples for target '", contrasts$target[i], "' of ", contrasts$id[i], call. = FALSE)
    x <- unique(trimws(unlist(strsplit(ifelse(is.na(x), "", x), ";"))))
    x <- x[x != ""]
    if (length(unique(samplesheet$lacO[as.character(samplesheet[[v]]) == as.character(contrasts$target[i])])) > 1)
      warning("target samples of ", contrasts$id[i], " have different lacO values; using all: ",
              paste(x, collapse = ";"), call. = FALSE)
    paste(x, collapse = ";")
  })
  data.frame(id = contrasts$id, lacO = lac)
}

# per contrast: TRUE for genes in a bin that contains one of that strain's lacO sites
flag_target_bins <- function(df, targets, tiles, siteMap) {
  if (is.null(siteMap)) return(df)   # keep "any site" flag from annotate_features()
  hit <- findOverlaps(resize(targets, width = 1, fix = "center"), tiles, ignore.strand = TRUE, select = "first")
  siteBin <- setNames(paste(as.character(seqnames(targets)), start(tiles)[hit]), targets$name)
  allowed <- unlist(lapply(names(siteMap), function(g) paste(g, siteBin[siteMap[[g]]])))
  df$target <- paste(df$group, df$seqnames, df$binStart) %in% allowed
  df
}


# annotation ---------------------------------------------------------------

# GRanges of the results with chromosome region type, anchor status and genome bin
annotate_features <- function(results, chrRegions, anchors, tiles) {
  gr <- tableToGranges(results)
  seqlevelsStyle(gr) <- "UCSC"

  # arm / centre / tip domains (Rockman & Kruglyak 2009)
  gr$chr_region <- NA_character_
  ol <- findOverlaps(gr, chrRegions, ignore.strand = TRUE)
  gr$chr_region[queryHits(ol)] <- chrRegions$DomainType[subjectHits(ol)]

  # rex anchor sites (X only)
  gr$anchor <- "autosome"
  gr$anchor[as.character(seqnames(gr)) == "chrX"] <- "nonAnchor"
  gr$anchor[overlapsAny(gr, anchors, ignore.strand = TRUE)] <- "Anchor"
  gr$anchor <- factor(gr$anchor, levels = c("autosome", "nonAnchor", "Anchor"))

  # genome bin containing the gene centre
  hit <- findOverlaps(resize(gr, width = 1, fix = "center"), tiles, ignore.strand = TRUE, select = "first")
  gr$binStart <- start(tiles)[hit]
  gr$target <- tiles$target[hit]
  gr
}


# plots --------------------------------------------------------------------

# boxplot of LFC per level of `xvar` (e.g. "seqnames", "chr_region", "anchor"),
# with gene counts and FDR-adjusted Wilcoxon p values of each level vs all genes.
# One panel per contrast, side by side.
boxplot_lfc <- function(df, xvar, ylimits, title, xlab, subtitle = NULL) {
  df <- tested_genes(df)
  df <- df[!is.na(df[[xvar]]), ]
  df[[xvar]] <- droplevels(factor(df[[xvar]]))
  df$group <- droplevels(factor(df$group))

  obsCounts <- df %>%
    group_by(across(all_of(c("group", xvar)))) %>%
    summarise(count = n(), .groups = "drop")
  wilcoxt <- df %>%
    group_by(group) %>%
    wilcox_test(as.formula(paste("log2FoldChange ~", xvar)), ref.group = "all") %>%
    adjust_pvalue(method = "fdr") %>%
    p_format(new.col = TRUE, accuracy = 1e-32)

  ggplot(df, aes(x = .data[[xvar]], y = log2FoldChange)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    facet_wrap(~group, nrow = 1) +
    coord_cartesian(ylim = ylimits) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = title, subtitle = subtitle, x = xlab, y = "Log2 Fold Change") +
    geom_text(data = obsCounts, aes(label = count, y = ylimits[1]), color = "blue", angle = 90) +
    stat_pvalue_manual(data = wilcoxt, label = "p.adj.format", y.position = ylimits[2] * 0.20,
                       remove.bracket = TRUE, angle = 90, hide.ns = TRUE, vjust = 0.5)
}

# boxplot of LFC in genome bins along one chromosome; bins containing a target are red
# (per contrast if the target column was set with flag_target_bins())
plot_binned <- function(df, chrom, binWidth, ylimits, title, subtitle = NULL, labelEvery = 5) {
  sub <- tested_genes(df[as.character(df$seqnames) == chrom, ])
  if (nrow(sub) == 0) stop("no tested genes on ", chrom)
  sub$bin <- factor(sub$binStart, levels = sort(unique(sub$binStart)))
  sub$group <- droplevels(factor(sub$group))

  obsCounts <- sub %>% group_by(group, bin) %>% summarise(count = n(), .groups = "drop")
  wilcoxt <- sub %>%
    group_by(group) %>%
    wilcox_test(log2FoldChange ~ bin, ref.group = "all") %>%
    adjust_pvalue(method = "fdr") %>%
    p_format(new.col = TRUE, accuracy = 1e-32)

  brks <- levels(sub$bin)[seq(1, nlevels(sub$bin), by = labelEvery)]

  ggplot(sub, aes(x = bin, y = log2FoldChange)) +
    geom_boxplot(aes(fill = target), outlier.shape = NA, show.legend = FALSE) +
    facet_wrap(~group, ncol = 1) +
    coord_cartesian(ylim = ylimits) +
    scale_fill_manual(values = c(`FALSE` = "grey90", `TRUE` = "red")) +
    geom_hline(yintercept = 0, alpha = 0.3) +
    scale_x_discrete(name = paste0(chrom, " (Mb, bin start)"), breaks = brks,
                     labels = (as.numeric(brks) - 1) / 1e6) +   # tile starts are 1-based
    labs(title = title, subtitle = subtitle, y = "Log2 Fold Change") +
    geom_text(data = obsCounts, aes(label = count, y = ylimits[1] * 0.9), color = "blue", angle = 90, size = 4) +
    stat_pvalue_manual(data = wilcoxt, label = "p.adj.format", y.position = ylimits[2] * 0.5,
                       remove.bracket = TRUE, angle = 90, hide.ns = TRUE, size = 3)
}

# LFC of each gene (as a line from 0) in a window of regionSize around the targets on one
# chromosome. Targets are marked with a dashed line and a square. With siteMap (see
# sites_per_contrast) the square is red in contrasts whose strain carries the site and
# an empty frame otherwise; without it, squares are coloured by site name if the bed file has one.
plot_zoom <- function(gr, targets, chrom, regionSize, xlim = NULL, ylim = NULL, siteMap = NULL) {
  t <- targets[as.character(seqnames(targets)) == chrom]
  if (length(t) == 0) stop("no targets on ", chrom)
  win <- trim(resize(t, width = regionSize, fix = "center"))
  df <- data.frame(subsetByOverlaps(gr, win))
  tdf <- data.frame(t)

  p <- ggplot(df, aes(x = start / 1e6, ymin = 0, ymax = log2FoldChange, color = log2FoldChange > 0)) +
    geom_linerange() +
    scale_color_manual(values = c(`TRUE` = "royalblue", `FALSE` = "darkorange"),
                       labels = c(`TRUE` = "up", `FALSE` = "down"), name = "LFC") +
    geom_hline(yintercept = 0) +
    geom_vline(data = tdf, aes(xintercept = start / 1e6), linetype = "dashed", color = "grey") +
    facet_grid(group ~ .) +
    labs(x = paste0(chrom, " genomic position (Mb)"), y = "Log2 Fold Change")

  # site name at the top of the first panel, next to its dashed line
  if ("name" %in% colnames(tdf) && any(!is.na(tdf$name))) {
    groups <- levels(factor(gr$group))
    lab <- tdf[!is.na(tdf$name), ]
    lab$group <- factor(groups[1], levels = groups)
    p <- p + geom_text(data = lab, aes(x = start / 1e6, y = Inf, label = name), inherit.aes = FALSE,
                       hjust = -0.05, vjust = 1.3, size = 2.5, color = "grey30")
  }

  if (!is.null(siteMap)) {
    # one square per site and contrast: filled if that strain carries the site
    groups <- levels(factor(gr$group))
    sdf <- do.call(rbind, lapply(groups, function(g) {
      d <- tdf
      d$group <- g
      d$present <- d$name %in% siteMap[[g]]
      d
    }))
    sdf$group <- factor(sdf$group, levels = groups)
    p <- p + geom_point(data = sdf, aes(x = start / 1e6, y = 0, fill = present), inherit.aes = FALSE,
                        shape = 22, size = 2.5) +
      scale_fill_manual(values = c(`TRUE` = "firebrick3", `FALSE` = "transparent"),   # FALSE: frame only
                        labels = c(`TRUE` = "in strain", `FALSE` = "not in strain"),
                        name = "lacO", drop = FALSE)
  } else if ("name" %in% colnames(tdf) && any(!is.na(tdf$name))) {
    p <- p + geom_point(data = tdf, aes(x = start / 1e6, y = 0, fill = name), inherit.aes = FALSE,
                        shape = 22, size = 2.5) + labs(fill = "Target")
  } else {
    p <- p + geom_point(data = tdf, aes(x = start / 1e6, y = 0), inherit.aes = FALSE,
                        shape = 15, color = "firebrick3", size = 2)
  }
  p + coord_cartesian(xlim = xlim, ylim = ylim)
}

# one zoom plot per chromosome, all with the same scales: the same y range, and x windows
# of the same width (the widest window, e.g. a chromosome with two targets), centred on
# each chromosome's targets
plot_zoom_panels <- function(gr, targets, chroms, regionSize, siteMap = NULL) {
  wins <- lapply(chroms, function(ch)
    trim(resize(targets[as.character(seqnames(targets)) == ch], width = regionSize, fix = "center")))
  span <- max(sapply(wins, function(w) max(end(w)) - min(start(w))))
  inWin <- overlapsAny(gr, unlist(GRangesList(wins)))
  ylim <- range(0, gr$log2FoldChange[inWin], na.rm = TRUE)
  lapply(seq_along(chroms), function(i) {
    mid <- (min(start(wins[[i]])) + max(end(wins[[i]]))) / 2
    plot_zoom(gr, targets, chroms[i], regionSize,
              xlim = (mid + c(-span, span) / 2) / 1e6, ylim = ylim, siteMap = siteMap)
  })
}


# tables -------------------------------------------------------------------

# number of significant genes (all, up, down) per contrast and level of each feature
count_significant <- function(df, features = c(chromosome = "seqnames", chr_region = "chr_region",
                                                anchor = "anchor"), padjCutoff = 0.05) {
  sig <- df[!is.na(df$padj) & df$padj < padjCutoff, ]
  bind_rows(lapply(features, function(v) {
    sig[!is.na(sig[[v]]), ] %>%
      group_by(across(all_of(c("group", v)))) %>%
      summarise(n_significant = n(),
                n_up = sum(log2FoldChange > 0),
                n_down = sum(log2FoldChange < 0), .groups = "drop") %>%
      rename(level = all_of(v)) %>%
      mutate(level = as.character(level))
  }), .id = "feature") %>%
    mutate(padj_cutoff = padjCutoff)
}
