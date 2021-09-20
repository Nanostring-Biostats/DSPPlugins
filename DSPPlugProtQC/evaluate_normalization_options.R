# Evaluate Normalization Options #

# Produces plots to help select the normalization method 
# Supports: DSP-nCounter Protein, DSP-nCounter RNA
# Note: this script should be run only on the dataset after QC and BEFORE normalization/scaling
# Please do not use spaces, special characters, or numbers when adding factors
# in the DSPDA Annotation file

#        User Options        #
##############################

# Please enter additional factor(s) you would like to plot
plot_factor <- c("Enter factor Here" ,"etc")

##########################################################
#### End of User Inputs. PLEASE DO NOT CHANGE CODE BELOW HERE  ####
##########################################################
# MIT License
# Copyright 2020 Nanostring Technologies, Inc.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# Contact us: 
#   NanoString Technologies, Inc.
#   530 Fairview Avenue N
#   Seattle, WA 98109
# Tel: (888) 358-6266
##############################

##############################
#        Execution Code      # 
##############################

# main function called by DSP-DA:
main <- function(dataset, segmentAnnotations, targetAnnotations, outputFolder) {
  
  #### preliminaries ----------------------------------
  
  # which columns of segmentAnnotations to color plots by:
  colorby <- intersect(c("ScanName", "SegmentName", plot_factor),
                       colnames(segmentAnnotations))
  
  #  define the color scheme:
  cols <- assign_colors(annot = segmentAnnotations[, colorby, drop = FALSE])
  
  # identify control probes:
  igg.names <- targetAnnotations$TargetGUID[targetAnnotations$CodeClass == "Negative"]
  hk.names <- targetAnnotations$TargetGUID[targetAnnotations$CodeClass == "Control"]
  
  # compute normalization factors:
  normfactors <- compute_normalization_factors(
    segmentAnnotations = segmentAnnotations,
    targetAnnotations = targetAnnotations,
    dataset = dataset,
    igg.names = igg.names,
    hk.names = hk.names
  )
  
  
  #### plot IgG concordance: ------------------------
  if (length(igg.names) > 1) {
    pdf(file = file.path(outputFolder, "igg_concordance.pdf", fsep = .Platform$file.sep), width = 12)
    par(mar = c(4, 4, 4, 1))
    for (varname in names(cols)) {
      tempmat <- t(pmax(dataset[igg.names, ], 1))
      colnames(tempmat) <- targetAnnotations[match(igg.names, targetAnnotations$TargetGUID), "TargetName"]
      colnames(tempmat) <- paste0(colnames(tempmat), " counts")
      plot_concordance(
        mat = tempmat,
        col = cols[[varname]][as.character(segmentAnnotations[, varname])],
        collegend = cols[[varname]],
        legend.main = varname,
        main = "Negative controls"
      )
    }
    dev.off()
  }
  
  
  #### plot HK concordance: --------------------------
  if (length(hk.names) > 1) {
    pdf(file = file.path(outputFolder, "housekeeper_concordance.pdf", fsep = .Platform$file.sep), width = 12)
    par(mar = c(4, 4, 2, 1))
    for (varname in names(cols)) {
      tempmat <- t(pmax(dataset[hk.names, ], 1))
      colnames(tempmat) <- targetAnnotations[match(hk.names, targetAnnotations$TargetGUID), "TargetName"]
      colnames(tempmat) <- paste0(colnames(tempmat), " counts")
      plot_concordance(
        mat = tempmat,
        col = cols[[varname]][as.character(segmentAnnotations[, varname])],
        collegend = cols[[varname]],
        legend.main = varname,
        main = "Housekeepers"
      )
    }
    dev.off()
  }
  
  
  ##### plot concordance among norm factors: ------------------------
  if (ncol(normfactors) > 1) {
    pdf(file = file.path(outputFolder, "normalization_factor_concordance.pdf", fsep = .Platform$file.sep), width = 12)
    par(mar = c(4, 4, 2, 1))
    # pairs plots:
    for (varname in names(cols)) {
      tempmat <- pmax(normfactors, 1)
      colnames(tempmat)[colnames(tempmat) == "HK geomean"] <- "HK geomean\n(counts)"
      colnames(tempmat)[colnames(tempmat) == "Neg geomean"] <- "Neg geomean\n(counts)"
      colnames(tempmat)[colnames(tempmat) == "Nuclei"] <- "Nuclei\n(counts)"
      colnames(tempmat)[colnames(tempmat) == "Area"] <- "Area (microns squared)"
      
      plot_concordance(
        mat = tempmat,
        col = cols[[varname]][as.character(segmentAnnotations[, varname])],
        collegend = cols[[varname]],
        legend.main = varname,
        main = "Normalization factors"
      )
    }
    dev.off()
  }
  
  #### QC proteins' signal level ---------------------------
  
  pdf(
    file = file.path(outputFolder, "signal_qc.pdf", fsep = .Platform$file.sep),
    width = pmin(pmax(10, nrow(dataset) * 0.2), 24)
  )
  qc_protein_signal(raw = as.matrix(dataset), neg.names = igg.names, targetAnnotations = targetAnnotations)
  dev.off()
}



#' QC proteins
#'
#' boxplot of signal-to-noise data to look for proteins below background
#' @param raw Matrix of raw count values. Proteins in rows, segments in columns.
#' @param neg.names Names of the negative controls
#' @param neg.thresh Threshold below which to flag the neg control factor as unreliably low.
#' @param targetAnnotations targetAnnotations object from DSPDA export
#' @param qccols Vector of two color names
#' @return Draws a boxplot for evaluating whether proteins ever get above background
#' @export
qc_protein_signal <- function(raw, neg.names, targetAnnotations = NULL, neg.thresh = 20, qccols = c("#00008B80", "#FF000080")) {
  
  # estimate background:
  negfactor <- pmax(colMeans(raw[neg.names, , drop = FALSE]), 1)
  
  # calc snr
  snr <- sweep(raw, 2, negfactor, "/")
  
  igginds <- which(is.element(rownames(snr), neg.names))
  o <- c(igginds, setdiff(order(apply(snr, 1, median)), igginds))
  
  protnames <- rownames(snr)
  if (length(targetAnnotations) > 0) {
    protnames <- targetAnnotations[match(rownames(snr), targetAnnotations$TargetGUID), "TargetName"]
  }
  par(mar = c(11, 4, 2, 1))
  boxplot(t(log2(snr[o, ])),
          las = 2,
          outline = FALSE,
          ylim = range(log2(snr)),
          names = protnames[o],
          ylab = "Log2 signal-to-background ratio",
          cex.axis = .85 - 0.3 * (nrow(snr) > 60)
  )
  axis(2, at = 1, labels = 1, las = 2, cex = 0.5)
  points(jitter(rep(1:nrow(snr), ncol(snr))),
         log2(snr[o, ]),
         col = "#00008B80",
         # col = qccols[1 + (negfactor < neg.thresh)],
         pch = 16, cex = 0.5
  )
  abline(h = 0)
  abline(v = length(igginds) + 0.5, lty = 2)
  abline(h = 1, lty = 2)
  # legend("bottomright", pch = 16, col = qccols, legend = c(paste0(c("Neg mean >= ", "Neg mean < "), neg.thresh)))
}


#' Choose annotation columns that look interesting to plot
#'
#' Automatically identifies columns of the annotation that look worth coloring by in
#'  normaliation QC plots. (Selects )
#' @param annot The segment annotation data frame
#' @return A vector of column names
choose_annotation_columns <- function(annot) {
  
  # color by all columns with > 1 and < 20 unique values:
  nvalues <- c()
  for (varname in colnames(segmentAnnotations)) {
    nvalues[varname] <- length(unique(segmentAnnotations[, varname]))
  }
  colorby <- names(nvalues)[(nvalues > 1) & (nvalues <= 20)]
  
  return(colorby)
}

#' Assign colors to levels of factors
#'
#' Given a data frame of factors or character vectors, assigns a unique color to each unique
#'  value of each variable.
#' @param annot Annotation matrix or data frame
#' @return A named list of named color vectors. Each element in the list is a variable, and the
#'  named color vector it contains gives colors for each value the vector takes.
assign_colors <- function(annot) {
  # vector of colors to choose from:
  colvec <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
    "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
    "#A6CEE3", "#B2DF8A",
    "#FB9A99", "#FDBF6F",
    "#CAB2D6", "#FFFF99",
    "#33A02C", "#FF7F00", "#B15928",
    "#1F78B4", "#E31A1C", "#6A3D9A",
    sample(colors(), 100)
  )
  # subsidiary function to fade colors (works like the "alpha" function from the "scales" package)
  fadecols <- function(cols, fade = .5) {
    fcols <- c()
    for (i in 1:length(cols))
    {
      tmp <- as.vector(col2rgb(cols[i]) / 256)
      fcols[i] <- rgb(tmp[1], tmp[2], tmp[3], fade)
    }
    return(fcols)
  }
  colvec <- fadecols(colvec, 0.7)
  
  # assign colors:
  cols <- list()
  colorby <- colnames(annot)
  for (varname in colorby) {
    varlevels <- as.character(unique(annot[, varname]))
    cols[[varname]] <- colvec[1:length(varlevels)]
    names(cols[[varname]]) <- varlevels
    # remove the used colors from further consideration:
    colvec <- setdiff(colvec, cols[[varname]]) # (disabling this so the more bold colors are re-used)
  }
  
  return(cols)
}


#' Compute basic normalization factors
#'
#' Given standard DSP-DA export, computes all possible normalization factors: housekeepers, igg's, area, nuclei
#' @param segmentAnnotations segmentAnnotations export data frame
#' @param targetAnnotations targetAnnotations export data frame
#' @param dataset dataset export data frame
#' @param igg.names Character vector giving igg names
#' @param hk.names Character vector giving HK probe names
#' @return  A matrix of normalization factors, with segments in rows and factors in columns
compute_normalization_factors <- function(segmentAnnotations, targetAnnotations, dataset, igg.names, hk.names) {
  
  # igg and hk factors:
  if (length(igg.names) > 1) {
    igg.factor <- exp(colMeans(log(pmax(dataset[igg.names, , drop = FALSE], 1))))
  }
  if (length(hk.names) > 1) {
    hk.factor <- exp(colMeans(log(pmax(dataset[hk.names, , drop = FALSE], 1))))
  }
  
  # area and nuclei factors:
  if (any(colnames(segmentAnnotations) == "AOIArea")) {
    area.factor <- segmentAnnotations$AOIArea
  }
  if (any(colnames(segmentAnnotations) == "AOINucleiCount")) {
    nuclei.factor <- segmentAnnotations$AOINucleiCount
  }
  
  # matrix of all available factors:
  factornames <- c("igg.factor", "hk.factor", "area.factor", "nuclei.factor")
  names(factornames) <- c("Neg geomean", "HK geomean", "Area", "Nuclei")
  factornames <- factornames[is.element(factornames, ls())]
  
  factors <- c()
  for (fname in factornames) {
    factors <- cbind(factors, get(fname))
  }
  colnames(factors) <- names(factornames)
  
  return(factors)
}


#' Plot pairwise concordance
#'
#' Draws a pairs plot of concordance between multiple normalization factors/ probes
#' @param mat Matrix of values to be compared. Could be e.g. IgG probe counts,
#'  or different normalization factors. Each variable/factor is in a column.
#' @param col Vector of colors for points
#' @param collegend Named vector of colors, used to draw a legend.
#' @param legend.main Title for the color legend, used to name the variable of interest
#' @param main Plot title
#' @param pch point type argument passed to pairs()
#' @param cex point size argument passed to pairs()
#' @param ... Additional arguments passed to pairs()
#' @example
#' # simulate HKs:
#' x = pmax(rnorm(100, 50, 10), 0)
#' mat = sweep(matrix(rnorm(300, 0, 5), 100), 1, x, "+")
#' plot_concordance(mat = mat, col = rep(c("blue", "orange"), each = 50))
plot_concordance <- function(mat, col = rgb(0, 0, 0, 0.5),
                             collegend = NULL, legend.main = NULL,
                             main = "", pch = 16, cex = 1.5, ...) {
  
  # subsidiary function for printing the SD of a ratio, for use by pairs():
  print.sd.log.ratio <- function(x, y, digits = 2, prefix = "", cex.cor = 1, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    pairwise.stat <- sd(log2(pmax(x, 1)) - log2(pmax(y, 1)), na.rm = TRUE)
    txt <- round(pairwise.stat, 2)
    txt <- paste0(prefix, txt)
    legend("center", legend = txt, box.col = "white", cex = 1.5)
  }
  
  # draw pairs plot:
  par(mar = c(4, 4, 2, 1))
  
  pairs(mat,
        log = "xy",
        upper.panel = points,
        lower.panel = print.sd.log.ratio,
        oma = c(3, 3, 3, 35),
        col = col,
        pch = pch,
        cex = cex,
        labels = colnames(mat),
        main = main
  )
  # draw a legend:
  if (length(collegend) > 0) {
    legend("right",
           col = c(NA, collegend, NA, NA),
           legend = c(legend.main, names(collegend), "", "Numbers show SD(log2(ratios))"),
           pch = 16
    )
  }
}