
gg_logFC <- function(x, y=NULL, rg = range(xy, na.rm = T)*1.02,
                     l.breaks = log1p(c(0, 1, 5, 10, 50, 100, Inf)),
                     l.labels = c('1', '5', '10', '50', '100', '100+'),
                     nbins = 200,
                     xlab = "log2(FC) in x",
                     ylab = "log2(FC) in y",
                     main = "", get.r =T, add.vd = T,
                     DEgsel = NULL,
                     ...){
  # Make a 2d binned hexplot for showing logFC comparison
  require(ggplot2)
  require(viridisLite)
  if(is.null(y))
    xy <- x
  else
    xy <- cbind(x, y)
  xy <- as.data.frame(xy)
  colnames(xy) <- c("x", "y")
  g <- ggplot(data = xy, mapping = aes(x=x, y=y)) +
    stat_bin_hex(aes(fill = after_stat(cut(log1p(after_stat(count)),
                                           breaks = l.breaks,
                              labels = F, right = T, include.lowest = T))),
                 bins=nbins) +
    scale_fill_gradientn(colours = viridisLite::inferno(length(l.breaks)),
                      breaks=1:length(l.labels),
                      labels=l.labels,
                      name = 'count',
                      limits=c(1,length(l.labels))) +
    theme_classic() + xlim(rg) + ylim(rg) +
    xlab(xlab) + ylab(ylab) + ggtitle(main) + coord_fixed()
  if(get.r){
    if(any(is.na(xy))){
      message("Warning: removed NA values to compute correlation coefs")
      rmsel <- apply(xy, 1, function(r) any(is.na(r)))
      xy <- xy[!rmsel,]
      if(!is.null(DEgsel))
        DEgsel <- DEgsel[!rmsel]
    }
    cc <- cor(xy)[1,2]
    cc2 <- cor(xy, method='spearman')[1,2]
    cctxt <- paste0("r = ", round(cc, 3), '\nrho = ', round(cc2, 3))
    if(add.vd)
      cctxt <- paste0(cctxt, "\n", round(100*(cc^2), 1), "% VarDev")
    if(!is.null(DEgsel)){
      ccDE <- cor(xy[DEgsel,])[1,2]
      cctxt <- paste0(cctxt, "\nr(DE) = ", round(ccDE, 3))
      if(add.vd)
        cctxt <- paste0(cctxt, "\n", round(100*(ccDE^2), 1), "% VarDev(DE)")
    }
    g <- g + annotate(geom="text", hjust=1, vjust=0,
                      x = rg[2],
                      y = rg[1],
                      label = cctxt)
  }
  return(g)
}



run_DESeq2_age <- function(X, p){
  # Run DEseq2 wt vs. mutant (with age covariate)
  require(DESeq2)
  rownames(p) <- colnames(X)
  p$aes <- scale(p$ae) # scale age value
  dds <- DESeqDataSetFromMatrix(countData = X,
                                colData = p,
                                design = ~aes+strain)
  dds <- DESeq(dds, test = "Wald", fitType = "local")
  return(dds)
}


get_DEres <- function(dds, coefname=NULL,contrastname=NULL){
  # Get results table from deseq output, managing NAs
  if(!is.null(coefname)){
    res <- results(dds, name=coefname,alpha=0.05)
  } else {
    res <- results(dds, contrast=contrastname, alpha=0.05)
  }
  # manage NAs
  res$padj[is.na(res$padj)] <- 1
  res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
  return(res)
}



find_df <- function(ref, w, dfs=2:8,varExplained=0.99){
  # Compute ssq of spline fits to find optimal df
  w.idx <- seq(max(c(1L, which.min(abs(w[1]-ref$time))-1 )),
               min(c(length(ref$time), which.min(abs(w[2]-ref$time))+1)))
  # get time values of window
  ts <- ref$time[w.idx]
  # get PCs of reference window
  pc <- summary(stats::prcomp(t(ref$interpGE[,w.idx]), scale=F, center=T))
  # keep enough components for 99% var
  spc <- sum(pc$importance[3, ] < varExplained)+1
  pcfit <- pc$x[, 1L:spc]
  w <- pc$importance[2, 1L:spc]
  # compute ssq of fit residuals with different dfs
  # for each df, compute weighted ssq of spline fit on
  ssqs <- cbind(sapply(dfs, function(dfi){
    ssq <- sum(w * colSums(
      stats::residuals(stats::lm(pcfit~splines::ns(ts, df=dfi)))^2
    ))
  }))/(length(ts))
  return(ssqs)
}

log1ptpm_2rawcounts <- function(X, glengths, nreadbygl){
  # Transform log1p(tpm) to (artificial) raw counts
  # note : nreadbygl = colSums(rawcounts/genelengths)
  if(length(nreadbygl) != ncol(X))
    stop("nreadbygl != ncol(X)")
  if(length(glengths) != nrow(X))
    stop("glengths must be of length nrow(X)")
  X <- t( (t(exp(X) - 1)/1e6) * nreadbygl ) * glengths
  X[X<0] <- 0
  return(round(X))
}

ref_2counts <- function(ref, ae_values,
                        gltable = wormRef::Cel_genes[,c("wb_id",
                                                        "transcript_length")],
                        avg_librarysize = 4e6){
  # Get expression profiles of given age from a RAPToR reference as
  # (artificial) count data.
  # note : gltable must have WBids as col 1 and gene length as col 2.
  # ref expression profiles at given timepoints :
  rX <- RAPToR::get_refTP(ref, ae_values=ae_values, return.idx = F)
  # transform to counts
  gl <- gltable[match(rownames(rX), gltable[,1]), 2]
  rX <- log1ptpm_2rawcounts(rX, gl, nreadbygl = rep(avg_librarysize,
                                                    ncol(rX))/median(gl))
  return(rX)
}

run_DESeq2_ref <- function(X, p, formula, ref,
                           ae_values=NULL, window.extend=1, ns.df=3){
  # Run DEseq2 wt vs. mutant correcting for development with ref. data.
  # Do no not specify age in the formula, it is added directly by the function.
  # Age estimates should either be an 'ae' column of p or given as 'ae_values'.
  require(DESeq2)
  if(!any(colnames(p)=='ae') & is.null(ae_values)){
    stop("Age estimates should either be a column of p, or given to ae_values.")
  }
  if(!is.null(ae_values)){
    p$ae <- ae_values
  }
  ## Extract reference expression profiles in sample time window
  w.rg <- range(p$ae) + c(-window.extend, window.extend)
  w.idx <- seq(
    max(c(1, which.min(abs(w.rg[1]-ref$time))-1 )),
    min(c(length(ref$time), which.min(abs(w.rg[2]-ref$time))+1))
  )
  # ref window time values
  w.ts <- ref$time[w.idx]
  # ref window expression values
  w.GE <- ref_2counts(ref = ref, ae_values = w.ts)
  ## Join ref & sample data
  nX <- ncol(X)
  nR <- ncol(w.GE)
  # get overlapping genes & join expression data
  ovl <- RAPToR::format_to_ref(samp = X, refdata = w.GE, verbose = F)
  Xj <- as.matrix(cbind(ovl$refdata, ovl$samp))
  # get relevant fields from p
  f0 <- as.formula(formula)
  p2 <- p[, attr(terms(f0), "term.labels"), drop=F]
  # get 1st level of each variable
  lev0 <- lapply(p2, function(col) levels(col)[1])
  # join p data and add Time and ref/sample Batch
  pj <- cbind(
    Time = c(w.ts, p$ae),
    Batch = factor(rep(c('r', 's'), c(nR, nX))),
    rbind(do.call(cbind, lapply(lev0, rep, times=nR)), p2) # other terms
  )
  rownames(pj) <- colnames(Xj)
  # Estimate dispersions with samples only
  s0 <- pj$Batch=="s"
  dds0 <- DESeqDataSetFromMatrix(countData = Xj[,s0],
                                 colData = pj[s0,],
                                 design = f0)
  dds0 <- estimateSizeFactors(dds0)
  dds0 <- estimateDispersions(dds0, fitType = "local")
  dd <- dispersions(dds0) # store dispersions
  dd[is.na(dd)] <- 0 # remove NAs
  ## Build full DE model with reference
  # Add Time and ref/sample batch to model formula
  f1 <- update.formula(f0, substitute(
    ~ splines::ns(Time, df = ns.df) + Batch + .,
    list(ns.df=ns.df)
  ))
  dds <- DESeqDataSetFromMatrix(countData = Xj,
                                colData = pj,
                                design = f1)
  dds <- estimateSizeFactors(dds)
  # inject dispersions from sample-only model
  dispersions(dds) <- dd
  dds <- nbinomWaldTest(dds)
  return(dds)
}


