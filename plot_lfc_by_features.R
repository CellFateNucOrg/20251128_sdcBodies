## Plot log2 fold changes by chromosome, chromosome region type, rex anchor status,
## genome bins along the target chromosomes, and zoomed around the target sites.
## Loops over every paramset (runName) in params.yml; all contrasts of a run go in the same file.
## At the end, the plots of all runs are combined into one multipage pdf per plot type
## (<plot type folder>/<plot name>_allRuns.pdf, next to the runName folders).

library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)
library(ggplot2)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)
library(ggpubr)
library(rstatix)
library(yaml)

theme_set(
  theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = ggtext::element_markdown(size = 9),
          axis.title.x = ggtext::element_markdown(size = 9),
          title = ggtext::element_markdown(size = 9))
)

# settings -----------------------------------------------------------------
serverPath <- "/Volumes/meister.data"
#serverPath <- "Z:/MeisterLab"

workDir <- paste0(serverPath, "/sinem/sequencing/20261028_se_ill_rna_l3_wPM1556_autosomal_bmm1")
setwd(workDir)

source(file.path(workDir, "scripts/functions_exploratory_analysis.R"))   # tableToGranges()
source(file.path(workDir, "scripts/functions_lfc_feature_plots.R"))

paramsFile      <- file.path(workDir, "scripts/params.yml")
contrastsFile   <- file.path(workDir, "contrasts.csv")
chrRegionsFile  <- file.path(serverPath, "publicData/various/Rockman_Kruglyak_2009_armsVcenter/celegans_rockman_kruglyak_2009_table1_domains_ce11_bp.tsv")
anchorsFile     <- file.path(serverPath, "publicData/various/rexSites_Das2023/382_X.eigs_cis.vecs_37peaks_p0.65_correct.bed")
targetsFile     <- file.path(workDir, "scripts/lacO_sites.bed")   # targets for the binned and zoomed plots
# which lacO site(s) each strain carries: a lacO column in the samplesheet (site names from
# the 4th column of targetsFile, several separated by ";", empty = none). Each contrast uses
# the lacO of the samples in its target group. Without a lacO column, every site is shown
# as present in every contrast.
samplesheetFile <- file.path(workDir, "samplesheet.csv")

binWidth     <- 3e5
regionSize   <- 3e5
padjCutoff   <- 0.05
ylimitsShrink   <- c(-0.25, 0.25)   # y axis of all LFC plots for runs with LFC shrinkage
ylimitsNoShrink <- c(-0.75, 0.75)   # ... and for runs without shrinkage ("noShrink" in runName)
excludeChroms   <- c("chrM")      # genes on these chromosomes are left out of all plots and tables
targetChroms <- NULL   # NULL = every chromosome with a site in targetsFile, or e.g. c("chrI", "chrV")

# data shared by all runs ----------------------------------------------------
paramsFiles <- file.path(workDir, "scripts", c("params.yml", "params_raptor.yml"))

runNames <- unique(unlist(lapply(paramsFiles, function(f)
  sapply(read_yaml(f), `[[`, "paramset_name"))))

contrasts <- read.csv(contrastsFile, header = TRUE)

si <- seqinfo(Celegans)
chrRegions <- GRanges(read.delim(chrRegionsFile))
seqlevelsStyle(chrRegions) <- "UCSC"
anchors <- import(anchorsFile)
seqlevelsStyle(anchors) <- "UCSC"
targets <- load_sites(targetsFile, si)
tiles <- make_tiles(seqlengths(si), binWidth, targets)
siteMap <- sites_per_contrast(lacO_from_samplesheet(contrasts, read.csv(samplesheetFile, header = TRUE)), targets)
if (!is.null(siteMap)) message("lacO per contrast: ", paste(names(siteMap), sapply(siteMap, paste, collapse = ";"), sep = " = ", collapse = ", "))

if (is.null(targetChroms)) targetChroms <- as.character(unique(seqnames(targets)))
targetChroms <- setdiff(targetChroms, excludeChroms)
chromOrder <- setdiff(c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX", "chrM"), excludeChroms)

# make a plot, save it to <plotRoot>/<typeDir>/<runName>/<name>, and keep it for the
# combined all-runs pdf. An error in one plot gives a warning and the loop continues.
plotRoot <- file.path(workDir, "diffabund/custom/plots")
combined <- list()   # "<typeDir>/<name>" -> list of pages
output <- function(plot, typeDir, name, width, height) {
  p <- tryCatch(plot, error = function(e) {
    warning(runName, ": could not make ", name, ": ", conditionMessage(e), call. = FALSE)
    NULL
  })
  if (is.null(p)) return(invisible(FALSE))
  d <- file.path(plotRoot, typeDir, runName)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  if (save_plot(p, file.path(d, name), width, height)) {
    key <- file.path(typeDir, name)
    combined[[key]] <<- c(combined[[key]], list(list(plot = p, width = width, height = height)))
  }
  invisible(TRUE)
}

# loop over runs ---------------------------------------------------------------
for (runName in runNames) {
  message("== ", runName)
  rdsFile <- file.path(workDir, "diffabund/custom/rds", runName, "results_annotated.RDS")
  if (!file.exists(rdsFile)) {
    warning("missing, skipped: ", rdsFile, call. = FALSE)
    next
  }
  results <- readRDS(rdsFile)
  ylimits <- ylimits_for(runName, shrink = ylimitsShrink, noShrink = ylimitsNoShrink)

  gr <- annotate_features(results, chrRegions, anchors, tiles)
  gr <- gr[!as.character(seqnames(gr)) %in% excludeChroms]
  df <- as.data.frame(gr)
  df$seqnames <- factor(df$seqnames, levels = intersect(chromOrder, unique(as.character(df$seqnames))))

  dirTxt <- file.path(workDir, "diffabund/custom/txt", runName)
  dir.create(dirTxt, showWarnings = FALSE, recursive = TRUE)

  # number of significant genes per chromosome / region type / anchor status
  write.table(count_significant(df, padjCutoff = padjCutoff),
              file.path(dirTxt, "n_significant_byFeature.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  contrastIds <- intersect(contrasts$id, unique(as.character(df$group)))
  notFound <- setdiff(contrasts$id, contrastIds)
  if (length(notFound) > 0) warning(runName, ": contrasts not in results: ", paste(notFound, collapse = ", "), call. = FALSE)
  if (length(contrastIds) == 0) next

  # all contrasts in the same file, in contrasts.csv order
  df <- df[df$group %in% contrastIds, ]
  df$group <- factor(as.character(df$group), levels = contrastIds)
  gr <- gr[as.character(gr$group) %in% contrastIds]
  gr$group <- factor(as.character(gr$group), levels = contrastIds)
  df <- flag_target_bins(df, targets, tiles, siteMap)   # red bins only in strains with that site
  nC <- length(contrastIds)
  message("   contrasts: ", paste(contrastIds, collapse = ", "))

  # boxplots: contrasts side by side -> width grows with the number of contrasts
  output(boxplot_lfc(df, "seqnames", ylimits, "Differential expression by chromosome",
                     "Chromosome", runName),
         "byChromosome", "boxplot_byChr.pdf", width = 3 + 6 * nC, height = 13)

  output(boxplot_lfc(df, "chr_region", ylimits, "Differential expression by chromosome region type",
                     "Chromosome region type", runName),
         "byChrRegion", "boxplot_byChrRegionType.pdf", width = 3 + 5.5 * nC, height = 13)

  output(boxplot_lfc(df, "anchor", ylimits, "Differential expression, anchor vs non-anchor genes",
                     "Anchor status", runName),
         "byChrRegion", "boxplot_anchorVnotanchor.pdf", width = 3 + 5.5 * nC, height = 13)

  # binned and zoomed plots: contrasts stacked -> height grows with the number of contrasts
  for (chrom in targetChroms) {
    output(plot_binned(df, chrom, binWidth, ylimits, paste0(chrom, " ", binWidth / 1e3, "kb bins"), runName),
           "byChrRegion", paste0(chrom, "_", binWidth / 1e3, "kb.pdf"), width = 32, height = 5 + 7 * nC)
  }

  if (length(targetChroms) > 0) {
    output({
      zp <- plot_zoom_panels(gr, targets, targetChroms, regionSize, siteMap)   # same scales in all panels
      annotate_figure(ggarrange(plotlist = zp, ncol = length(zp), common.legend = TRUE, legend = "bottom"),
                      top = text_grob(runName, size = 9))
    }, "byChrRegion", paste0("zoomed_", regionSize / 1e6, "Mb_lineplots.pdf"),
    width = 4 + 7.5 * length(targetChroms), height = 5 + 6 * nC)
  }
}

# all runs combined: one multipage pdf per plot type, above the runName folders ---------
for (key in names(combined)) {
  f <- file.path(plotRoot, sub("\\.pdf$", "_allRuns.pdf", key))
  save_multipage(combined[[key]], f)
  message("combined ", length(combined[[key]]), " runs: ", f)
}
