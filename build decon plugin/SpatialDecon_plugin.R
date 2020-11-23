# Copyright (C) 2020, NanoString Technologies, Inc.
#    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#    You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us: # NanoString Technologies, Inc. # 530 Fairview Avenue N # Seattle, WA 98109 # Tel: (888) 358-6266


##########################################################
####      Arguments to be modified by user            ####
##########################################################

# Cell profile matrix filename:
cell_profile_filename <- "safeTME-for-tumor-immune.csv"

# annotation column giving pure tumor AOIs
# (this column must have the value "tumor" for indicate tumor AOIs)
pure_tumor_column_name <- "none"

# define variables to show in heatmaps:
variables_to_plot <- c("ScanName", "SegmentName")


# define coloring of annotations (optional):
# 1: set the next line to TRUE:
custom_annotation_coloring <- FALSE
# and then modify the example syntax below:
if (custom_annotation_coloring) {
  # example syntax:
  cols <- list(
    ScanName = c(
      "scan1" = "red",
      "scan2" = "blue",
      "scan3" = "darkblue",
      "scan4" = "forestgreen"
    ),
    SegmentName = c(
      "Tumor" = "chartreuse3",
      "TME" = "dodgerblue3"
    )
  )
  # for a list of R colors, see http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
}

# heatmap color scheme:
hmcols <- c("white", "darkblue")
# alternative heatmap colors:
# remove the "#" to enable either of these lines:
# viridis option B:
# hmcols = c("#000004FF", "#1B0C42FF", "#4B0C6BFF", "#781C6DFF", "#A52C60FF", "#CF4446FF", "#ED6925FF", "#FB9A06FF", "#F7D03CFF", "#FCFFA4FF")
# viridis option D:
# hmcols = c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF")


#### advanced arguments

heatmaptruncationlimit <- NULL
pdf_width <- 12
pdf_height <- 7
draw_svgs_instead_of_pdf <- FALSE
subset_of_cells_to_show <- NULL


##########################################################
#### end of arguments. DO NOT CHANGE CODE BELOW HERE  ####
##########################################################


library(logNormReg)
library(pheatmap)
library(viridis)
library(scales)
library(openxlsx)
library(dplyr)

main <- function(dataset, segmentAnnotations, targetAnnotations, outputFolder) {

  #### preliminaries ----------------------
  dataset <- as.matrix(dataset)

  # access cell profile matrix file:
  X <- as.matrix(read.csv(cell_profile_filename, header = TRUE, row.names = 1))


  # ARGUMENT (hidden): define cell types to be added together in the final result:
  # example syntax:
  # merges = list()
  # merges[["T"]] = c("CD8.T", "CD4.T")
  # merges[["myeloid"]] = c("macrophage", "monocyte", "DC")
  merges <- list()

  # parse merges:
  mergesFull <- NULL
  if (length(merges) > 0) {
    # initialize with 1:1 mapping:
    mergesFull <- list()
    for (name in colnames(X)) {
      mergesFull[[name]] <- name
    }
    # add merges:
    for (name in names(merges)) {
      # remove entries for cells specified by user and replace with their entries:
      mergesFull[merges[[name]]] <- NULL
      mergesFull[[name]] <- merges[[name]]
    }
  }


  # ARGUMENT (hidden): enter the name of the column giving nuclei counts
  nuclei_count_column_name <- "this_is_hidden_for_advanced_users" # "AOINucleiCount"

  # parse nuclei column
  cell_counts <- NULL
  if (is.element(nuclei_count_column_name, colnames(segmentAnnotations))) {
    cell_counts <- as.numeric(segmentAnnotations[, nuclei_count_column_name])
  }
  # if (!is.element(nuclei_count_column_name, colnames(segmentAnnotations))) {
  #  warning("The value entered for nuclei_count_column_name was not a column header in the segment annotations.
  #          Results will not be output on the scale of cell counts; just in abundance scores and proportions.")
  # }

  # parse pure tumor column
  is_pure_tumor <- NULL
  if (is.element(pure_tumor_column_name, colnames(segmentAnnotations))) {
    is_pure_tumor <- tolower(segmentAnnotations[, pure_tumor_column_name]) == "tumor"
    is_pure_tumor <- replace(is_pure_tumor, is.na(is_pure_tumor), FALSE)
  }
  if (!is.element(pure_tumor_column_name, colnames(segmentAnnotations)) & (pure_tumor_column_name != "none")) {
    warning("The value entered for pure_tumor_column_name was not a column header in the segment annotations.")
  }

  # format data for spatialdecon:
  norm <- dataset[targetAnnotations$TargetGUID, segmentAnnotations$segmentID]
  rownames(norm) <- targetAnnotations$TargetName
  if (all(is.element(c("ScanName", "ROIName", "SegmentName"), colnames(segmentAnnotations)))) {
    segmentAnnotations <- mutate(segmentAnnotations, 
                               segmentDisplayName = paste(ScanName, ROIName, SegmentName, sep=" | ")) 
    if (all(!duplicated(segmentAnnotations$segmentDisplayName))) {
      colnames(norm) <- segmentAnnotations$segmentDisplayName
      rownames(segmentAnnotations) <- segmentAnnotations$segmentDisplayName
    }
  }
  
  # calculate background:
  bg <- derive_GeoMx_background(
    norm = norm,
    probepool = targetAnnotations$ProbePool,
    negnames = targetAnnotations$TargetName[targetAnnotations$CodeClass == "Negative"]
  )



  #### run decon: ----------------------------------------
  # decon:
  res <- spatialdecon(
    norm = norm,
    bg = bg,
    X = X,
    is_pure_tumor = is_pure_tumor,
    cell_counts = cell_counts,
    cellmerges = mergesFull
  )


  #### write results files: ---------------------------------------------

  # initialize a workbook to write to excel:
  wb <- createWorkbook("decon")
  # write beta:
  addWorksheet(wb, "Abundance scores")
  writeData(wb,
    sheet = "Abundance scores",
    res$beta,
    colNames = TRUE, rowNames = TRUE
  )
  # write props:
  addWorksheet(wb, "Proportions")
  writeData(wb,
    sheet = "Proportions",
    res$prop_of_nontumor,
    colNames = TRUE, rowNames = TRUE
  )
  # write scaled beta:
  addWorksheet(wb, "Scaled abundance scores")
  writeData(wb,
    sheet = "Scaled abundance scores",
    sweep(res$beta, 1, apply(res$beta, 1, max), "/"),
    colNames = TRUE, rowNames = TRUE
  )

  # write cell counts if available:
  if (is.element("cell.counts", names(res))) {
    addWorksheet(wb, "Cell counts")
    writeData(wb,
      sheet = "Cell counts",
      res$cell.counts$cell.counts,
      colNames = TRUE, rowNames = TRUE
    )
  }
  # add segment annotations:
  addWorksheet(wb, "Segment annotations")
  writeData(wb,
    sheet = "Segment annotations",
    segmentAnnotations,
    colNames = TRUE, rowNames = FALSE
  )


  # save the excel workbook:
  saveWorkbook(wb,
    file = file.path(outputFolder, "spatialdecon_outputs.xlsx", fsep = .Platform$file.sep),
    overwrite = TRUE
  )



  # set to TRUE to save reverse decon results
  if (FALSE) {
    # reverse decon:
    rdecon <- reverseDecon(
      norm = norm,
      beta = res$beta,
      epsilon = NULL
    )
    write.csv(rdecon$resids, file = file.path(outputFolder, "reverse_decon_residuals.csv", fsep = .Platform$file.sep))
    # reverse decon summary stats of gene dependency on cell mixing:
    sumstats <- cbind(rdecon$cors, rdecon$resid.sd)
    colnames(sumstats) <- c("cor w cell mixing", "residual SD from cell mixing")
    write.csv(sumstats, file = file.path(outputFolder, "gene_dependence_on_cell_mixing.csv", fsep = .Platform$file.sep))
  }


  #### results figures: ---------------------------------------------

  # parse the argument for variables to plot:
  if (length(setdiff(variables_to_plot, colnames(segmentAnnotations))) > 0) {
    warning(paste0(
      "the variables_to_plot values",
      paste0(setdiff(variables_to_plot, colnames(segmentAnnotations)), collapse = ", "),
      " are not present in the segmentAnnotations"
    ))
  }
  variables_to_plot <- intersect(variables_to_plot, colnames(segmentAnnotations))
  heatmapannot <- NULL
  if (length(variables_to_plot) > 0) {
    heatmapannot <- segmentAnnotations[, variables_to_plot, drop = FALSE]
    rownames(heatmapannot) <- colnames(res$beta)
    # rownames(heatmapannot) <- segmentAnnotations$segmentDisplayName
  }

  # colors for variables_to_plot:
  if (!exists("cols")) {
    cols <- assign_colors(segmentAnnotations[, variables_to_plot])
  }

  # show just the original cells, not tumor abundance estimates derived from the is.pure.tumor argument:
  cells.to.plot <- intersect(rownames(res$beta), union(colnames(X), names(mergesFull)))

  ## show only a subset of cells if specified:
  if (length(subset_of_cells_to_show) >= 2) {
    cells.to.plot <- intersect(cells.to.plot, subset_of_cells_to_show)
  }


  # one pdf for all results:
  if (!draw_svgs_instead_of_pdf) {
    pdf(file = file.path(outputFolder, "spatialdecon_results.pdf", fsep = .Platform$file.sep), width = pdf_width, height = pdf_height)
  }

  #### heatmaps
  # abundances:
  if (length(heatmaptruncationlimit) == 1) {
    thresh <- heatmaptruncationlimit
  }
  if (length(heatmaptruncationlimit) == 0) {
    thresh <- signif(quantile(res$beta, 0.97), 2)
  }

  if (draw_svgs_instead_of_pdf) {
    svg(file = file.path(outputFolder, "spatialdecon_results_abundance_heatmap.svg", fsep = .Platform$file.sep), width = pdf_width, height = pdf_height)
  }
  p1 <- pheatmap(pmin(res$beta[cells.to.plot, ], thresh),
    col = colorRampPalette(hmcols)(100),
    fontsize_col = 4,
    angle_col= 90,
    annotation_col = heatmapannot,
    annotation_colors = cols,
    legend_breaks = c(round(seq(0, thresh, length.out = 5))[-5], thresh),
    legend_labels = c(round(seq(0, thresh, length.out = 5))[-5], paste0("Abundance scores,\ntruncated above at ", thresh))
    # main = paste0("Abundance scores, truncated above at ", thresh)
  )

  if (draw_svgs_instead_of_pdf) {
    dev.off()
  }
  # print(p1)

  # scaled abundances:
  epsilon <- min(res$beta[res$beta > 0])
  mat <- sweep(res$beta, 1, pmax(apply(res$beta, 1, max), epsilon), "/")

  if (draw_svgs_instead_of_pdf) {
    svg(file = file.path(outputFolder, "spatialdecon_results_scaled_abundance_heatmap.svg", fsep = .Platform$file.sep), width = pdf_width, height = pdf_height)
  }
  p3 <- pheatmap(mat[cells.to.plot, ],
    col = colorRampPalette(hmcols)(100),
    fontsize_col = 4,
    angle_col= 90,
    annotation_col = heatmapannot,
    annotation_colors = cols,
    legend_breaks = c(round(seq(0, 1, length.out = 5), 2)[-5], 1),
    legend_labels = c(round(seq(0, 1, length.out = 5), 2)[-5], "Scaled abundance\n(ratio to max)")
    # main = paste0("Abundance scores, truncated above at ", thresh)
  )
  if (draw_svgs_instead_of_pdf) {
    dev.off()
  }

  # proportions:
  props <- replace(res$prop_of_nontumor[cells.to.plot, ], is.na(res$prop_of_nontumor[cells.to.plot, ]), 0)

  if (draw_svgs_instead_of_pdf) {
    svg(file = file.path(outputFolder, "spatialdecon_results_proportions_heatmap.svg", fsep = .Platform$file.sep), width = pdf_width, height = pdf_height)
  }
  p2 <- pheatmap(props,
    col = colorRampPalette(hmcols)(100),
    fontsize_col = 4,
    angle_col= 90,
    annotation_col = heatmapannot,
    annotation_colors = cols,
    legend_breaks = round(seq(0, max(props) * 0.99, length.out = 5), 2),
    legend_labels = c(round(seq(0, max(props), length.out = 5), 2)[-5], "Proportion of all\nfitted populations")
  )
  # print(p2)

  if (draw_svgs_instead_of_pdf) {
    dev.off()
  }


  #### barplots setup:

  # use safeTME colors if the right cells are present:
  if (all(is.element(cells.to.plot, names(cellcols)))) {
    col <- cellcols[cells.to.plot]
  }
  if (!all(is.element(cells.to.plot, names(cellcols)))) {
    manycols <- c(
      "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
      "#B3DE69", "#FCCDE5", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
      "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#1B9E77", "#D95F02",
      "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
      sample(grDevices::colors(), 99)
    )
    col <- manycols[seq_len(length(cells.to.plot))]
    names(col) <- cells.to.plot
  }

  # choose an appropriate label size:
  namescex <- 1
  if (ncol(res$beta) > 20) {
    namescex <- 0.75
  }
  if (ncol(res$beta) > 40) {
    namescex <- 0.5
  }
  if (ncol(res$beta) > 80) {
    namescex <- 0.25
  }


  ### abundance barplot:

  if (draw_svgs_instead_of_pdf) {
    svg(file = file.path(outputFolder, "spatialdecon_results_abundance_barplot.svg", fsep = .Platform$file.sep), width = pdf_width, height = pdf_height)
  }
  layout(mat = matrix(c(1, 2, 3, 3), 2), widths = c(10, 3, 10, 3), heights = c(1, 8, 10))
  par(mar = c(0, 8.2, 0, 0.2))
  plot(p1$tree_col, labels = F, main = "", ylab = "", yaxt = "n")
  par(mar = c(15, 8, 0, 0))

  # data to plot:
  mat <- res$beta[cells.to.plot, p1$tree_col$order]
  # infer scale of negative y-axis for annotation colorbars
  ymin <- -max(colSums(mat)) * 0.15
  if (!is.finite(ymin)) {
    ymin <- 0
  }

  # draw barplot:
  bp <- barplot(mat,
    cex.lab = 1.5,
    col = col, border = NA,
    cex.names = namescex,
    las = 2, main = "", ylab = "Abundance scores",
    ylim = c(ymin, max(colSums(mat)))
  )
  # add color bars
  for (name in rev(variables_to_plot)) {
    yrange <- seq(ymin / 3, ymin, length.out = length(variables_to_plot) + 1)[match(name, variables_to_plot) + c(0, 1)]
    xwidth <- (bp[2] - bp[1]) / 2
    for (i in 1:ncol(mat)) {
      rect(bp[i] - xwidth, yrange[2], bp[i] + xwidth, yrange[1],
        border = NA, col = cols[[name]][segmentAnnotations[match(colnames(mat)[i], segmentAnnotations$segmentID), name]]
      )
    }
  }
  axis(2,
    at = seq(ymin / 3, ymin, length.out = length(variables_to_plot) + 2)[-c(1, length(variables_to_plot) + 2)],
    las = 2, labels = variables_to_plot, lty = 0, cex.axis = 0.75
  )
  # draw a legend:
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  frame()
  legendcols <- legendnames <- c()
  for (name in rev(names(cols))) {
    legendcols <- c(legendcols, NA, cols[[name]], NA)
    legendnames <- c(legendnames, name, names(cols[[name]]), NA)
  }
  legend("center",
    pch = 15,
    col = c(legendcols, rep(NA, 1), rev(col)),
    legend = c(legendnames, "Cell", rev(names(col)))
  )
  if (draw_svgs_instead_of_pdf) {
    dev.off()
  }




  ### proportion barplot:
  if (draw_svgs_instead_of_pdf) {
    svg(file = file.path(outputFolder, "spatialdecon_results_proportion_barplot.svg", fsep = .Platform$file.sep), width = pdf_width, height = pdf_height)
  }

  layout(mat = matrix(c(1, 2, 3, 3), 2), widths = c(10, 3, 10, 3), heights = c(1, 8, 10))
  par(mar = c(0, 8.2, 0, 0.2))
  plot(p1$tree_col, labels = F, main = "", ylab = "", yaxt = "n")
  par(mar = c(15, 8, 0, 0))

  # data to plot:
  mat <- res$prop_of_nontumor[cells.to.plot, p2$tree_col$order]
  mat <- replace(mat, is.na(mat), 0)
  # infer scale of negative y-axis for annotation colorbars
  ymin <- -0.15
  # draw barplot:
  bp <- barplot(mat,
    cex.lab = 1.5,
    col = col, border = NA,
    cex.names = namescex,
    las = 2, main = "", ylab = "Proportion of fitted cells",
    ylim = c(ymin, max(colSums(mat)))
  )
  # add color bars
  for (name in rev(variables_to_plot)) {
    yrange <- seq(ymin / 3, ymin, length.out = length(variables_to_plot) + 1)[match(name, variables_to_plot) + c(0, 1)]
    xwidth <- (bp[2] - bp[1]) / 2
    for (i in 1:ncol(mat)) {
      rect(bp[i] - xwidth, yrange[2], bp[i] + xwidth, yrange[1],
        border = NA, col = cols[[name]][segmentAnnotations[match(colnames(mat)[i], segmentAnnotations$segmentID), name]]
      )
    }
  }
  axis(2,
    at = seq(ymin / 3, ymin, length.out = length(variables_to_plot) + 2)[-c(1, length(variables_to_plot) + 2)],
    las = 2, labels = variables_to_plot, lty = 0, cex.axis = 0.75
  )
  # draw a legend:
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  frame()
  legendcols <- legendnames <- c()
  for (name in rev(names(cols))) {
    legendcols <- c(legendcols, NA, cols[[name]], NA)
    legendnames <- c(legendnames, name, names(cols[[name]]), NA)
  }
  legend("center",
    pch = 15,
    col = c(legendcols, rep(NA, 1), rev(col)),
    legend = c(legendnames, "Cell", rev(names(col)))
  )
  if (draw_svgs_instead_of_pdf) {
    dev.off()
  }





  # set to TRUE to draw a "spaceplot" for each cell type (showing abundace of the cell in space)
  if (FALSE) {
    # spaceplots: draw separately for each cell:
    # ARGUMENTS: define which variables specify xy coordinates:
    xpositionname <- "ROICoordinateX"
    ypositionname <- "ROICoordinateY"
    tissueidname <- "ScanName"


    if (all(is.element(c(xpositionname, ypositionname, tissueidname), colnames(segmentAnnotations)))) {
      if (all(!is.na(segmentAnnotations[, c(xpositionname, ypositionname, tissueidname)]))) {

        # stratify on SegmentName if it has multiple levels:
        # (assumption: SegmentName gives segment type, identifies the different kinds of AOIs that can share an ROI)
        aoitype <- rep("", nrow(segmentAnnotations))
        if (is.element("SegmentName", colnames(segmentAnnotations))) {
          if (length(unique(segmentAnnotations[, "SegmentName"])) > 1) {
            aoitype <- segmentAnnotations[, "SegmentName"]
            aoitype <- replace(aoitype, is.na(aoitype), "na")
          }
        }
        for (varname in variables_to_plot) {
          # abundance spaceplots
          for (atype in unique(aoitype)) {
            for (cell in cells.to.plot) {
              use <- aoitype == atype
              sp <- spaceplot(
                x = segmentAnnotations[use, xpositionname],
                y = segmentAnnotations[use, ypositionname],
                z = res$beta[cell, use],
                tissue = segmentAnnotations[use, tissueidname],
                tissue.order = NULL,
                tissuecols = NULL,
                tissuecols.alpha = 0.2,
                rescale = TRUE,
                cex = 3, col = cols[[varname]][segmentAnnotations[use, varname]],
                nrow = NULL, rowind = NULL, colind = NULL,
                expansion = 1.2,
                main = cell
              )
              for (name in names(sp$boundaries)) {
                text(
                  x = median(range(sp$boundaries[[name]]$x)),
                  y = max(sp$boundaries[[name]]$y),
                  name
                )
              }
            }
          }

          # proportion spaceplots
          for (atype in unique(aoitype)) {
            use <- aoitype == atype
            for (cell in cells.to.plot) {
              sp <- spaceplot(
                x = segmentAnnotations[use, xpositionname],
                y = segmentAnnotations[use, ypositionname],
                z = replace(res$prop_of_nontumor[cell, use], is.na(res$prop_of_nontumor[cell, use]), 0),
                tissue = segmentAnnotations[use, tissueidname],
                tissue.order = NULL,
                tissuecols = NULL,
                tissuecols.alpha = 0.2,
                rescale = TRUE,
                cex = 3, col = cols[[varname]][segmentAnnotations[use, varname]],
                nrow = NULL, rowind = NULL, colind = NULL,
                expansion = 1.2,
                main = cell
              )
              for (name in names(sp$boundaries)) {
                text(
                  x = median(range(sp$boundaries[[name]]$x)),
                  y = max(sp$boundaries[[name]]$y),
                  name
                )
              }
            }
          }
        }
      }
    } # end spaceplots
  }
  if (!draw_svgs_instead_of_pdf) {
    dev.off()
  }
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

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.




#' SpatialDecon mixed cell deconvolution algorithm
#'
#' Runs the generic SpatialDecon decon workflow, including:
#' \enumerate{
#' \item run deconvolution once
#' \item remove poorly-fit genes from first round of decon
#' \item re-run decon with cleaned-up gene set
#' \item compute p-values
#' }
#'
#' @param Y p-length expression vector or p * N expression matrix - the actual
#' (linear-scale) data
#' @param X p * K Training matrix.
#' @param bg Expected background counts. Provide a scalar to apply to all
#'  data points, or
#'  else a matrix/vector aligning with Y to provide more nuanced expected
#'   background.
#' @param weights The same as the weights argument used by lm
#' @param resid_thresh A scalar, sets a threshold on how extreme individual
#' data points' values
#'  can be (in log2 units) before getting flagged as outliers and set to NA.
#' @param lower_thresh A scalar. Before log2-scale residuals are calculated,
#' both observed and fitted
#'  values get thresholded up to this value. Prevents log2-scale residuals from
#'  becoming extreme in
#'  points near zero.
#' @param align_genes Logical. If TRUE, then Y, X, bg, and wts are row-aligned
#' by shared genes.
#' @param maxit Maximum number of iterations. Default 1000.
#' @return a list:
#' \itemize{
#' \item beta: matrix of cell abundance estimates, cells in rows and
#' observations in columns
#' \item sigmas: covariance matrices of each observation's beta estimates
#' \item p: matrix of p-values for H0: beta == 0
#' \item t: matrix of t-statistics for H0: beta == 0
#' \item se: matrix of standard errors of beta values
#' \item resids: a matrix of residuals from the model fit.
#' (log2(pmax(y, lower_thresh)) - log2(pmax(xb, lower_thresh))).
#' }
algorithm2 <- function(Y, X, bg = 0, weights = NULL,
                       resid_thresh = 3, lower_thresh = 0.5,
                       align_genes = TRUE, maxit = 1000) {

    # align genes:
    if (align_genes) {
        sharedgenes <- intersect(rownames(X), rownames(Y))
        Y <- Y[sharedgenes, ]
        X <- X[sharedgenes, ]
        if (is.matrix(bg)) {
            bg <- bg[sharedgenes, ]
        }
        if (is.matrix(weights)) {
            weights <- weights[sharedgenes, ]
        }
    }

    # format the data nicely:
    tidied <- tidy_X_and_Y(X, Y)
    X <- tidied$X
    Y <- tidied$Y
    if ((length(bg) > 0) & (is.vector(bg))) {
        bg <- matrix(bg, nrow = length(bg))
    }

    # select an epsilon (lowest non-zero value to use)
    epsilon <- min(Y[(Y > 0) & !is.na(Y)])


    # initial run to look for outliers:
    out0 <- deconLNR(
        Y = Y, X = X, bg = bg, weights = weights, epsilon = epsilon,
        maxit = maxit
    )
    # also get yhat and resids:
    out0$yhat <- X %*% out0$beta + bg
    out0$resids <- log2(pmax(Y, lower_thresh)) -
        log2(pmax(out0$yhat, lower_thresh))

    # ID bad genes:
    outliers <- flagOutliers(
        Y = Y,
        yhat = out0$yhat,
        wts = weights,
        resids = out0$resids,
        resid_thresh = resid_thresh
    )

    # remove outlier data points:
    Y.nooutliers <- replace(Y, outliers, NA)

    # re-run decon without outliers:
    out <- deconLNR(
        Y = Y.nooutliers,
        X = X,
        bg = bg,
        weights = weights,
        epsilon = epsilon
    )
    out$yhat <- X %*% out$beta + bg
    out$resids <- log2(pmax(Y.nooutliers, 0.5)) - log2(pmax(out$yhat, 0.5))

    # compute p-values
    tempbeta <- out$beta
    tempse <- tempp <- tempt <- tempbeta * NA
    for (i in seq_len(ncol(tempse))) {
        tempse[, i] <- suppressWarnings(sqrt(diag(out$sigmas[, , i])))
    }
    tempt <- (tempbeta / tempse)
    tempp <- 2 * (1 - stats::pt(tempt, df = nrow(X) - ncol(X) - 1))
    out$p <- tempp
    out$t <- tempt
    out$se <- tempse

    # structure of output: beta, hessians, yhat, resids
    return(out)
}



#' Function to format Y, X inputs for decon
#'
#' Takes user-supplied X and Y, checks for accuracy, aligns by dimnames, adds
#' dimnames if missing
#'
#' @param X X matrix
#' @param Y Data matrix
#' @return X and Y, both formatted as matrices, with full dimnames and aligned
#' to each other by dimname
tidy_X_and_Y <- function(X, Y) {

    # format as matrices:
    Ynew <- Y
    if (is.vector(Y)) {
        Ynew <- matrix(Y, nrow = length(Y), dimnames = list(names(Y), "y"))
    }
    Xnew <- X

    # check alignment:
    if (!identical(rownames(Y), rownames(X))) {
        warning("Rows (genes) of X and Y are mis-aligned.")
    }
    out <- list(X = Xnew, Y = Ynew)
}

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.


#' Collapse related cell types within a deconvolution result
#'
#' Given the input of an SpatialDecon result output and a list of which cell
#' types to combine,
#'  returns a reshaped deconvolution result object with the specified cell
#'  types merged.
#' @param fit The object (a list) returned by the SpatialDecon algorithm
#' @param matching A list object holding the mapping from beta's cell names to
#' official cell names.
#'  See str(safeTME.matches) for an example.
#' @return A reshaped deconvolution result object
#' @examples
#' data(mini_geomx_dataset)
#' # estimate background:
#' mini_geomx_dataset$bg <- derive_GeoMx_background(
#'   norm = mini_geomx_dataset$normalized,
#'   probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
#'   negnames = "NegProbe"
#' )
#' # run basic decon:
#' res0 <- spatialdecon(
#'   norm = mini_geomx_dataset$normalized,
#'   bg = mini_geomx_dataset$bg,
#'   X = safeTME
#' )
#' res1 <- collapseCellTypes(
#'     fit = res0,
#'     matching = safeTME.matches
#' )
#' @export
collapseCellTypes <- function(fit, matching) {

    # results object to hold the collapsed results:
    out <- fit

    # format matching list as a matrix to take a linear combination of beta:
    startingcellnames <- unlist(matching)
    A <- matrix(0, length(matching), nrow(fit$beta),
        dimnames = list(names(matching), rownames(fit$beta))
    )
    for (name in names(matching)) {
        cellnames <- matching[[name]]
        A[name, cellnames] <- 1
    }

    # apply A transformation to beta:
    for (name in c("beta", "prop_of_all", "prop_of_nontumor")) {
        if (is.element(name, names(fit))) {
            out[[name]] <- A[, startingcellnames] %*% fit[[name]][startingcellnames, ]
        }
    }

    # if Sigma provided, get vcov of beta2:
    if (is.element("sigmas", names(out))) {
        sigma <- fit$sigmas
        if (length(dim(sigma)) == 2) {
            out$sigmas <- A[, startingcellnames] %*%
                sigma[startingcellnames, startingcellnames, ] %*%
                t(A[, startingcellnames])
        }
        if (length(dim(sigma)) == 3) {
            out$sigmas <- array(NA,
                dim = c(nrow(A), nrow(A), dim(sigma)[3]),
                dimnames = list(rownames(A), rownames(A), dimnames(sigma)[[3]])
            )
            for (i in seq_len(dim(sigma)[3])) {
                out$sigmas[, , i] <- A %*% sigma[, , i] %*% t(A)
            }
        }
    }

    # re-calculate p, se, t:
    if (is.element("beta", names(out)) & is.element("sigmas", names(out))) {
        # compute p-values
        tempbeta <- out$beta
        tempse <- tempp <- tempt <- tempbeta * NA
        for (i in seq_len(ncol(tempse))) {
            tempse[, i] <- suppressWarnings(sqrt(diag(out$sigmas[, , i])))
        }
        out$se <- tempse
        out$t <- (tempbeta / tempse)
        if (is.element("X", names(out))) {
            out$p <- 2 * (1 - stats::pt(out$t, df = nrow(fit$X) - ncol(fit$X) - 1))
        }
    }

    return(out)
}

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.



#' Convert abundance measurements to cell counts
#'
#' Converts cell abundance scores to cell counts, under the assumption that the
#'  observation with
#'  the greatest sum of cell abundance scores is 100% modelled cells.
#' @param beta Matrix of cell abundance scores, with cells in rows and
#'  observations in columns. The
#'  assumption is that this matrix is from well-normalized data.
#' @param nuclei.counts Optional. A vector of total nuclei counts. If provided,
#' the function will
#'  output not only cells.per.100 but also total cells.
#' @param omit.tumor Logical. If FALSE, any rows of beta with "tumor" in their
#' name will be omitted.
#' @return A list with two elements, each a rescaled version of beta.
#' cells.per.100 gives estimated
#'  percents of total, and cell.counts is cells.per.100 * nuclei.counts.
convertCellScoresToCounts <- function(beta, nuclei.counts = NULL,
                                      omit.tumor = FALSE) {
    # strip tumor rows if called for:
    if (omit.tumor) {
        beta <- beta[!grepl("tumor", rownames(beta)), , drop = FALSE]
    }

    # calc max abundance scores:
    max.total.abundance <- max(colSums(beta))

    # calculate rescaled scores:
    out <- list()
    out$cells.per.100 <- beta / max.total.abundance * 100
    if (length(nuclei.counts) == ncol(beta)) {
        out$cell.counts <- sweep(out$cells.per.100, 2, nuclei.counts, "*") / 100
    }
    return(out)
}

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.



#' Convert results from an arbitrary training matrix to our standardized cell
#' types.
#'
#' Takes betas from a decon run, using any cell type names whatsoever, and maps
#'  them back to the
#' "official" cell types we use. Allows for multiple rows of beta to map to the
#'  same official cell type,
#' in which case those rows will be added up.
#'
#' @param beta K * n matrix of estimated beta values (cell type abundances)
#' @param matching A list object holding the mapping from beta's cell names to
#'  official cell names.
#'  See str(safeTME.matches) for an example.
#' @param stat The function used to combine related cell types. Defaults to sum.
#' @param na.rm Whether to ignore NAs while computing stat
#' @param sigma A list of covariance matrices of beta estimates, in the format
#'  output by spatialdecon.
#' @return A list with two elements:
#' \itemize{
#' \item beta: a matrix of cell abundances, with specified cell types added
#' together
#' \item sigma: an array of covariance matrices for each observation's beta
#' vector
#' }
convertCellTypes <- function(beta, matching, stat = sum,
                             na.rm = FALSE, sigma = NULL) {
    # format matching list as a matrix to take a linear combination of beta:
    A <- matrix(0, length(matching), nrow(beta),
        dimnames = list(names(matching), rownames(beta))
    )
    for (name in names(matching)) {
        cellnames <- matching[[name]]
        A[name, cellnames] <- 1
    }

    # apply A transformation to beta:
    beta2 <- A %*% beta

    # if Sigma provided, get vcov of beta2:
    if (length(sigma) > 0) {
        if (length(dim(sigma)) == 2) {
            sigma2 <- A %*% sigma %*% t(A)
        }
        if (length(dim(sigma)) == 3) {
            sigma2 <- array(NA,
                dim = c(nrow(A), nrow(A), dim(sigma)[3]),
                dimnames = list(rownames(A), rownames(A), dimnames(sigma)[[3]])
            )
            for (i in seq_len(dim(sigma)[3])) {
                sigma2[, , i] <- A %*% sigma[, , i] %*% t(A)
            }
        }
    }

    # if no Sigma, just return transformed beta:
    if (length(sigma) == 0) {
        return(beta2)
    }
    # if there is a sigma, return beta and the sigma:
    if (length(sigma) > 0) {
        out <- list(beta = beta2, sigma = sigma2)
        return(out)
    }
}

#' Default colors for the cell types in the safeTME matrix
#'
#' A named vector of colors, giving colors for the cell types of the safeTME
#'  matrix.
#'
#' @format A named vector
"cellcols"


#' Default colors for the cell types in the safeTME matrix
#'
#' A named vector of colors, giving colors for the cell types of the safeTME
#'  matrix.
#'
#' @format A named vector
"mini_geomx_dataset"


#' Mapping from granularly-defined cell populations to broaded cell populations
#'
#' Mapping from granularly-defined cell populations to broaded cell populations,
#'  for use by the convertCellTypes function.
#'
#' @format A list. Each element of the list contains the granular cell types
#'  that roll up
#'  to a single coarse cell type.
"safeTME.matches"


#' SafeTME matrix
#'
#' A matrix of expression profiles of 906 genes over 18 cell types.
#'
#' @format A matrix with 906 genes (rows) and 18 cell types (columns)
"safeTME"


#' Genes' biological variability in immune deconvolution from TCGA.
#'
#' Genes' biological SDs, as estimated from immune deconvolution from TCGA.
#' Used to weight genes in spatialdecon.
#'
#' @format A named vector giving SDs of 1179 genes.
"mean.resid.sd"

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.




#' Deconvolution using logNormReg package to run linear mean model and log error
#'  model
#'
#' Calls lognlm() to optimize the model.
#'
#' @param Y p-length expression vector or p * N expression matrix - the actual
#'  (linear-scale) data
#' @param X p * K Training matrix
#' @param bg scalar or matrix of expected background counts per data point.
#' @param weights The same as the weights argument used by lm
#' @param epsilon optional,  a very small non-zero number to use as a lower
#' threshold to make fits well-behaved
#' @param maxit Maximum number of iterations. Default 1000.
#' @return a list: beta (estimate), sigmas (covariance matrix of estimate,
#' derived by inverting the hessian from lognlm)
#' @import logNormReg
deconLNR <- function(Y, X, bg = 0, weights = NULL, epsilon = NULL,
                     maxit = 1000) {
    if (length(weights) == 0) {
        weights <- replace(Y, TRUE, 1)
    }
    if (is.vector(Y)) {
        Y <- matrix(Y, nrow = length(Y))
    }
    if (length(bg) == 1) {
        bg <- matrix(bg, nrow(Y), ncol(Y))
    }
    # choose "epsilon": a very small non-zero number to make fits well-behaved
    if (length(epsilon) == 0) {
        epsilon <- min(replace(Y, (Y == 0) & !is.na(Y), NA), na.rm = TRUE)
    }

    # matrix-like data for apply:
    mat <- rbind(Y, bg, weights)
    # fn to apply:
    fn <- function(zz) {
        # break into y, b, w:
        y <- zz[seq_len((length(zz)) / 3)]
        b <- zz[seq_len((length(zz)) / 3) + (length(zz) / 3)]
        wts <- zz[seq_len((length(zz)) / 3) + (length(zz) / 3) * 2]

        # remove NA data:
        use <- !is.na(y)
        y <- y[use]
        b <- b[use]
        Xtemp <- X[use, , drop = FALSE]
        wts <- wts[use]

        init <- rep(mean(y) / (mean(X) * ncol(X)), ncol(X))
        names(init) <- colnames(X)

        # run lognlm:
        fit <- lognlm(pmax(y, epsilon) ~ b + Xtemp - 1,
            lik = FALSE,
            weights = wts,
            start = c(1, init),
            method = "L-BFGS-B",
            lower = c(1, rep(0, ncol(Xtemp))),
            upper = c(1, rep(Inf, ncol(Xtemp))),
            opt = "optim",
            control = list(maxit = maxit)
        )
        fnout <- list(
            beta = fit$coefficients[-1],
            sigma = solve(fit$hessian)[-1, -1]
        )
        return(fnout)
    }
    # apply across all observations:
    fnlist <- apply(mat, 2, fn)
    # extract beta and sigmas:
    getbeta <- function(zz) {
        return(zz$beta)
    }
    getsigma <- function(zz) {
        return(zz$sigma)
    }

    beta <- vapply(X = fnlist, FUN = getbeta, FUN.VALUE = numeric(ncol(X)))
    rownames(beta) <- colnames(X)
    sigmas <- array(vapply(
        X = fnlist,
        FUN = getsigma,
        FUN.VALUE = numeric(ncol(X)^2)
    ),
    dim = c(ncol(X), ncol(X), ncol(Y)),
    dimnames = list(colnames(X), colnames(X), colnames(Y))
    )

    out <- list(beta = pmax(beta, 0), sigmas = sigmas)
    return(out)
}

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.



#' Derive background at the scale of the normalized data for GeoMx data
#'
#' Estimates per-datapoint background levels from a GeoMx experiment.
#' In studies with two or more probe pools, different probes will have different
#' background levels. This function provides a convenient way to account for
#' this phenomenon.
#'
#' @param norm Matrix of normalized data, genes in rows and segments in columns.
#'  Must include negprobes, and must have rownames.
#' @param probepool Vector of probe pool names for each gene, aligned to the
#' rows of "norm".
#' @param negnames Names of all negProbes in the dataset. Must be at least one
#'  neg.name within each probe pool.
#' @return A matrix of expected background values, in the same scale and
#'  dimensions as the "norm" argument.
#' @examples
#' data(mini_geomx_dataset)
#' # estimate background:
#' mini_geomx_dataset$bg <- derive_GeoMx_background(
#'   norm = mini_geomx_dataset$normalized,
#'   probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
#'   negnames = "NegProbe"
#' )
#' @export
derive_GeoMx_background <- function(norm, probepool, negnames) {

    # check data input:
    if (nrow(norm) != length(probepool)) {
        stop("nrow(norm) != length(probepool)")
    }

    # initialize:
    bg <- norm * 0


    # fill in expected background at scale of normalized data:
    for (pool in unique(probepool)) {

        # get the pool's negProbes:
        tempnegs <- intersect(negnames, rownames(norm)[probepool == pool])
        if (length(tempnegs) == 0) {
            stop(paste0(pool, " probe pool didn't have any negprobes specified"))
        }
        tempnegfactor <- colMeans(norm[tempnegs, , drop = FALSE])

        # fill in the corresponding elements of bg:
        bg[probepool == pool, ] <-
            sweep(bg[probepool == pool, ], 2, tempnegfactor, "+")
    }
    return(bg)
}

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.




#' Compute gene weights
#'
#' Compute gene weights from pre-defined biological noise estimates and from
#'  raw count/ error-model-defined technical noise estimates
#' @param norm Matrix of data used in decon. Used to define the gene list and
#' the
#' shape of the output "wts" matrix.
#' @param raw Matrix of raw data. If provided, used to define technical noise
#' @param error.model Which error model to use. Defaults to "dsp"
#' @param weight.by.TIL.resid.sd If TRUE, then genes are weighted in part based
#' on their
#'  biological variability as estimated by their residual SD from decon
#'   performed on TCGA.
#' @return A matrix of weights, in the same dimension as norm
deriveWeights <- function(norm, raw = NULL, error.model = "dsp",
                          weight.by.TIL.resid.sd = FALSE) {

    # get tech SDs if raw data provided:
    if (length(raw) == 0) {
        sds.tech <- matrix(0.1, nrow(raw), ncol(raw), dimnames = dimnames(raw))
    }
    if (length(raw) > 0) {
        sds.tech <- runErrorModel(
            counts = raw,
            platform = "dsp"
        )
    }

    # if the mean.resid.sd vector (which defines genes' biological SD) is in
    # the environment, get biological noise:
    if (!weight.by.TIL.resid.sd) {
        sds.bio <- matrix(0.1, nrow(raw), ncol(raw), dimnames = dimnames(raw))
    }
    if (weight.by.TIL.resid.sd) {
        sds.bio <- matrix(NA, nrow(raw), ncol(raw), dimnames = dimnames(raw))
        for (gene in intersect(
            names(mean.resid.sd),
            rownames(sds.bio)
        )) {
            sds.bio[gene, ] <- mean.resid.sd[gene]
        }
        sds.bio <- replace(sds.bio, is.na(sds.bio), mean(sds.bio, na.rm = TRUE))
    }

    # define total SD, and invert to get weights
    sds.tot <- sqrt(sds.tech^2 + sds.bio^2)
    wts <- 1 / sds.tech
    return(wts)
}

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.



#' Download a cell profile matrix
#'
#' Download a cell profile matrix from the online library
#'
#' @param matrixname A name
#' @return A cell profile matrix
#' @details Valid values for the matrixname argument include:
#' \itemize{
#' \item Airway_Epithelium
#' \item Atlas_Adult_Retina_10x
#' \item Census_Adult_Immune_10x
#' \item Census_Newborn_Blood_10x
#' \item Diff_Fetal_Neuron_SS2
#' \item FetalMaternal_Adult_Blood_10x
#' \item FetalMaternal_Adult_Blood_SS2
#' \item FetalMaternal_Adult_Decidua_10x
#' \item FetalMaternal_Adult_Decidua_SS2
#' \item FetalMaternal_Fetal_Placenta_10x
#' \item Human_brain
#' \item Human_Cell_Landscape
#' \item IBD_Adult_Colon_10x
#' \item Landscape_Adult_Liver_10x
#' \item Lung_plus_neutrophils
#' \item Mouse_Brain
#' \item Profiling_Adult_BoneMarrow_10x
#' \item Reprogram_Embryo_Dendritic_10x
#' \item Sensitivity_Adult_Esophagus_10x
#' \item Sensitivity_Adult_Lung_10x
#' \item Sensitivity_Adult_Spleen_10x
#' \item Somatic_Adult_Pancreas_SS2
#' \item SpatioTemporal_Adult_Kidney_10x
#' \item SpatioTemporal_Fetal_Kidney_10x
#' \item Tcell_Adult_Blood_10x
#' \item Tcell_Adult_BoneMarrow_10x
#' \item Tcell_Adult_Lung_10x
#' \item Tcell_Adult_LymphNode_10x
#' }
#' @examples
#' X <- download_profile_matrix(matrixname = "Human_brain")
#' head(X)
#' @export
download_profile_matrix <- function(matrixname) {

    # check formatting:
    if (length(matrixname) > 1) {
        stop("specify just one matrixname")
    }

    librarynames <- c(
        "Airway_Epithelium", "Atlas_Adult_Retina_10x", "Census_Adult_Immune_10x",
        "Census_Newborn_Blood_10x", "Diff_Fetal_Neuron_SS2",
        "FetalMaternal_Adult_Blood_10x", "FetalMaternal_Adult_Blood_SS2",
        "FetalMaternal_Adult_Decidua_10x", "FetalMaternal_Adult_Decidua_SS2",
        "FetalMaternal_Fetal_Placenta_10x", "Human_brain", "Human_Cell_Landscape",
        "IBD_Adult_Colon_10x", "Landscape_Adult_Liver_10x",
        "Lung_plus_neutrophils", "Mouse_Brain", "Profiling_Adult_BoneMarrow_10x",
        "Reprogram_Embryo_Dendritic_10x", "Sensitivity_Adult_Esophagus_10x",
        "Sensitivity_Adult_Lung_10x", "Sensitivity_Adult_Spleen_10x",
        "Somatic_Adult_Pancreas_SS2", "SpatioTemporal_Adult_Kidney_10x",
        "SpatioTemporal_Fetal_Kidney_10x", "Tcell_Adult_Blood_10x",
        "Tcell_Adult_BoneMarrow_10x", "Tcell_Adult_Lung_10x",
        "Tcell_Adult_LymphNode_10x"
    )
    if (!is.element(matrixname, librarynames)) {
        warning(paste0(matrixname, " is not an expected cell profile matrix name."))
    }

    X <- as.matrix(utils::read.csv(paste0(
        "https://raw.githubusercontent.com/patrickjdanaher/cell-profile-library/master/profile_matrices/",
        matrixname, ".csv"
    ), row.names = 1))

    X <- X[rowSums(X) > 0, ]
    return(X)
}

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.



#' Identify outlier genes in a decon result
#'
#' Analyses a decon result's residuals to identify poorly-fit genes and flag
#' them for removal.
#'  Rule: flag anything with
#'
#' @param Y p-length expression vector or p * N expression matrix - the actual
#' (linear-scale) data
#' @param yhat Expectation of Y given the decon fit
#' @param resids Log2-scale residuals of Y vs. yhat
#' (log2(pmax(Y, 0.5)) - log2(yhat))
#' @param wts Matrix of data point weights, aligned to dimensions of resids
#' @param resid_thresh A scalar, sets a threshold on how extreme individual
#' data points' values
#'  can be (in log2 units) before getting flagged as outliers and set to NA.
#' @return a vector of names of poorly-fit genes
flagOutliers <- function(Y, yhat, resids, wts, resid_thresh = 3) {

    # get weighted resids:
    if (length(wts) == 0) {
        wres <- resids
    }
    if (length(wts) > 0) {
        wres <- resids * wts
    }

    # flag bad genes:
    outlier_genes <- c() # <-- this line makes it so no outlier genes are filtered
    # flag bad data points: (not doing anything for now)
    outlier_data_points <- abs(resids) > resid_thresh
    return(outlier_data_points)
}

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.


#' Draw coxcomb plots as points in a graphics window
#'
#' Draws a scatterplot where each point is a circular barplot, intended to show
#' decon results
#'
#' @param x Vector of x coordinates
#' @param y Vector of y coordinates
#' @param b matrix or cell abundances, with columns aligned with the elements
#' of x and y
#' @param col vector of colors, aligned to the rows of b.
#' @param legendwindow Logical. If TRUE, the function draws a color legend in a
#'  new window
#' @param rescale.by.sqrt Logical, for whether to rescale b by its square root
#' to make value proportional to
#'  shape area, not shape length.
#' @param border Color of pie segment border, defauls to NA/none
#' @param add Logical. If TRUE, the function draws florets atop an existing
#' graphics device (TRUE) or call a new device (FALSE).
#' @param cex Floret size. Florets are scaled relative to the range of x and y;
#' this further scales up or down.
#' @param bty bty argument passed to plot()
#' @param xaxt xaxt argument passed to plot()
#' @param yaxt yaxt argument passed to plot()
#' @param xlab xlab, defaults to ""
#' @param ylab ylab, defaults to ""
#' @param ... additional arguments passed to plot()
#' @return Draws a coxcomb plot, returns no data.
#' @examples
#' data(mini_geomx_dataset)
#' # estimate background:
#' mini_geomx_dataset$bg <- derive_GeoMx_background(
#'   norm = mini_geomx_dataset$normalized,
#'   probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
#'   negnames = "NegProbe"
#' )
#' # run basic decon:
#' res0 <- spatialdecon(
#'   norm = mini_geomx_dataset$normalized,
#'   bg = mini_geomx_dataset$bg,
#'   X = safeTME
#' )
#' # draw florets:
#' florets(
#'   x = mini_geomx_dataset$annot$x,
#'   y = mini_geomx_dataset$annot$y,
#'   b = res0$beta, cex = 2
#' )
#' @export
florets <- function(x, y, b, col = NULL, legendwindow = FALSE,
                    rescale.by.sqrt = TRUE, border = NA, add = FALSE, cex = 1,
                    bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                    ...) {
    # rescale b by sqrt so magnitude is proportional to coxcomb area, not length
    if (rescale.by.sqrt) {
        b <- sqrt(b)
    }
    # make b a matrix:
    if (is.vector(b)) {
        b2 <- matrix(b, nrow = length(b))
        rownames(b2) <- names(b)
        b <- b2
        rm(b2)
    }
    # choose colors if not given:
    if ((length(col) == 0) &
        all(is.element(rownames(b), names(cellcols)))) {
        col <- cellcols[rownames(b)]
    }
    if ((length(col) == 0) &
        !all(is.element(rownames(b), names(cellcols)))) {
        manycols <- c(
            "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
            "#B3DE69", "#FCCDE5", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
            "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#1B9E77", "#D95F02",
            "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
            sample(grDevices::colors(), 99)
        )
        col <- manycols[seq_len(nrow(b))]
        names(col) <- rownames(b)
    }
    # convert colors to matrix of the same dimension as b:
    if (length(col) == 1) {
        col <- matrix(col, nrow = nrow(b), ncol = ncol(b))
    }
    if (is.vector(col)) {
        col <- matrix(col, nrow = nrow(b), ncol = ncol(b))
    }

    # get radians:
    angles <- seq(0, 2 * pi, length.out = nrow(b) + 1)

    # scale b based on the range of x and y:
    maxrange <- max(diff(range(x, na.rm = TRUE)), diff(range(y, na.rm = TRUE)))
    b <- b * maxrange / mean(b, na.rm = TRUE) * 0.007 * cex

    # draw plot:
    if (!add) {
        graphics::plot(x, y,
            col = 0, bty = bty, xaxt = xaxt, yaxt = yaxt,
            xlab = xlab, ylab = ylab, ...
        )
    }

    # draw florets:
    if (nrow(b) > 1) {
        for (i in seq_len(length(x))) {
            for (j in seq_len(nrow(b))) {
                tempangles <- seq(angles[j], angles[j + 1], length.out = 20)
                xt <- b[j, i] * cos(tempangles)
                yt <- b[j, i] * sin(tempangles)
                graphics::polygon(x[i] + c(0, xt), y[i] + c(0, yt),
                    col = col[j, i],
                    border = border, lwd = 0.5
                )
            }
        }
    }

    # if just one point, draw a full circle:
    if (nrow(b) == 1) {
        for (i in seq_len(length(x))) {
            for (j in seq_len(nrow(b))) {
                tempangles <- seq(angles[j], angles[j + 1], length.out = 20)
                xt <- b[j, i] * cos(tempangles)
                yt <- b[j, i] * sin(tempangles)
                graphics::polygon(x[i] + xt, y[i] + yt,
                    col = col[j],
                    border = border, lwd = 0.5
                )
            }
        }
    }

    # draw a legend:
    if (legendwindow) {
        graphics::plot(0, 0,
            col = 0, xlim = c(-1, 1), ylim = c(-1, 1), xaxt = "n",
            yaxt = "n", xlab = "", ylab = "", ...
        )
        for (j in seq_len(length(angles))) {
            graphics::lines(c(0, 0.75 * cos(angles[j])), c(0, 0.75 * sin(angles[j])),
                col = col[j], lwd = 2
            )
            graphics::text(0.85 * cos(angles[j]), 0.85 * sin(angles[j]),
                rownames(b)[j],
                cex = 1.4
            )
        }
    }
}

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.



#' Estimate a tumor-specific profile and merge it with the pre-specified cell
#'  profile matrix (X)
#'
#' Given the input of "tumor-only" AOI's, estimates an collection of
#'  tumor-specific
#' expression profiles and merges them with the immune cell expression
#' training matrix.
#' The process:
#' \enumerate{
#' \item log2/normalized data from tumor-only AOIs is clustered with hclust,
#' and cutree() is used to define clusters.
#' \item 2. Each cluster's geomean profile is merged into the immune cell
#' profile matrix.
#' }
#'
#' @param norm matrix of normalized data
#' @param bg matrix of expected background, on the scale of norm.
#' @param pure_tumor_ids Vector identifying columns of norm that are pure tumor.
#'  Can be indices, logicals or column names.
#' @param X The training matrix
#' @param K the number of clusters to fit
#' @return an updated X matrix with new columns, "tumor.1", "tumor.2", ...
#' @examples
#' data(mini_geomx_dataset)
#' mini_geomx_dataset$bg <- derive_GeoMx_background(
#'   norm = mini_geomx_dataset$normalized,
#'   probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
#'   negnames = "NegProbe"
#' )
#' safeTME.with.tumor <- mergeTumorIntoX(
#'   norm = mini_geomx_dataset$norm,
#'   bg = mini_geomx_dataset$bg,
#'   pure_tumor_ids = mini_geomx_dataset$annot$AOI.name == "Tumor",
#'   X = safeTME,
#'   K = 3
#' )
#' @export
mergeTumorIntoX <- function(norm, bg, pure_tumor_ids, X, K = 10) {

    # round up 0 values in norm:
    min.nonzero <- min(norm[norm > 0], na.rm = TRUE)
    norm <- pmax(norm, min.nonzero)

    # subset data to only the pure tumor IDs:
    norm <- norm[, pure_tumor_ids, drop = FALSE]
    bg <- bg[, pure_tumor_ids, drop = FALSE]

    # bg-subtract:
    norm <- pmax(norm - bg, min(norm) / 20)

    # fix K if too big:
    if (ncol(norm) < K) {
        K <- ncol(norm)
    }

    # case 1: want to use every column in norm as a separate profile:
    #  (includes case of just one column in norm)
    if (K == ncol(norm)) {
        tumorX <- norm
    }

    # case 2: if many tumor AOIs, get profiles for K clusters of data:
    if (K < ncol(norm)) {
        # cluster and cut:
        h <- stats::hclust(stats::dist(t(log2(norm))))
        cut <- stats::cutree(h, k = K)
        # get clusters' geomean profiles:
        tumorX <- c()
        for (cid in unique(cut)) {
            tumorX <- cbind(
                tumorX,
                exp(rowMeans(log(norm[, cut == cid, drop = FALSE])))
            )
        }
        colnames(tumorX) <- paste0("tumor.", seq_len(ncol(tumorX)))
    }

    # align tumorX with X:
    sharedgenes <- intersect(rownames(tumorX), rownames(X))
    tumorX <- tumorX[sharedgenes, ]
    X <- X[sharedgenes, ]

    # rescale tumor X:
    meanq90 <- max(mean(apply(X, 2, stats::quantile, 0.9)), 1e-3)
    tumorq90s <- apply(tumorX, 2, stats::quantile, 0.9)
    tumorX <- sweep(tumorX, 2, tumorq90s, "/") * meanq90

    # merge:
    out <- cbind(X, tumorX)
    return(out)
}

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.



#' SpatialDecon: A package for computating the notorious bar statistic.
#'
#' The SpatialDecon package estimates mixed cell type abundance in the regions
#' of spatially-resolved gene
#' expression studies, using the method of Danaher & Kim (2020), "Advances in
#'  mixed cell deconvolution enable
#' quantification of cell types in spatially-resolved gene expression data."
#' It is also appropriate to apply to bulk gene expression data.
#'
#' @section functions:
#' Functions to help set up deconvolution:
#' \itemize{
#'  \item derive_GeoMx_background Estimates the background levels from GeoMx
#'  experiments
#'  \item collapseCellTypes reformats deconvolution results to merge
#'  closely-related cell types
#'  \item download_profile_matrix Downloads a cell profile matrix.
#'  \item safeTME: a data object, a matrix of immune cell profiles for use in
#'   tumor-immune deconvolution.
#' }
#' Deconvolution functions:
#' \itemize{
#'  \item spatialdecon runs the core deconvolution function
#'  \item reverseDecon runs a transposed/reverse deconvolution problem, fitting
#'  the data as a function of cell abundance estimates.
#'   Used to measure genes' dependency on cell mixing and to calculate gene
#'    residuals from cell mixing.
#' }
#' Plotting functions:
#' \itemize{
#'  \item florets Plot cell abundance on a specified x-y space, with each point
#'   a cockscomb plot showing the cell abundances of that region/sample.
#'  \item TIL_barplot Plot abundances of tumor infiltrating lymphocytes (TILs)
#'   estimated from the safeTME cell profile matrix
#' }
#' @examples
#' data(mini_geomx_dataset)
#' # estimate background:
#' mini_geomx_dataset$bg <- derive_GeoMx_background(
#'   norm = mini_geomx_dataset$normalized,
#'   probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
#'   negnames = "NegProbe"
#' )
#' # run basic decon:
#' res0 <- spatialdecon(
#'   norm = mini_geomx_dataset$normalized,
#'   bg = mini_geomx_dataset$bg,
#'   X = safeTME
#' )
#' # run decon with bells and whistles:
#' res <- spatialdecon(
#'   norm = mini_geomx_dataset$normalized,
#'   bg = mini_geomx_dataset$bg,
#'   X = safeTME,
#'   cellmerges = safeTME.matches,
#'   cell_counts = mini_geomx_dataset$annot$nuclei,
#'   is_pure_tumor = mini_geomx_dataset$annot$AOI.name == "Tumor"
#' )
#' @docType package
#' @name SpatialDecon-package
NULL

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.


#' Reverse deconvolution
#'
#' Performs "reverse deconvolution", modelling each gene expression's ~
#' cell scores.
#' Returns a matrix of "fitted" expression values, a matrix of residuals,
#'  a matrix of
#' reverse decon coefficients for genes * cells.
#'
#' @param norm Matrix of normalized data, with genes in rows and observations
#' in columns
#' @param beta Matrix of cell abundance estimates, with cells in rows and
#' observations in columns.
#'  Columns are aligned to "norm".
#' @param epsilon All y and yhat values are thresholded up to this point when
#' performing decon.
#'  Essentially says, "ignore variability in counts below this threshold."
#' @return A list:
#' \itemize{
#' \item coefs, a matrix of coefficients for genes * cells, where element i,j is
#'  interpreted as
#' "every unit increase in cell score j is expected to increase expression of
#' gene i by _".
#' \item yhat, a matrix of fitted values, in the same dimension as norm
#' \item resids, a matrix of log2-scale residuals from the reverse decon fit,
#'  in the same
#'  dimension as norm
#' \item cors, a vector giving each gene's correlation between fitted and
#' observed expression
#' \item resid.sd, a vector of each gene's residual SD, a metric of how much
#' variability genes
#'  have independend of cell mixing.
#' }
#' @import logNormReg
#' @examples
#' data(mini_geomx_dataset)
#' # estimate background:
#' mini_geomx_dataset$bg <- derive_GeoMx_background(
#'   norm = mini_geomx_dataset$normalized,
#'   probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
#'   negnames = "NegProbe"
#' )
#' # run basic decon:
#' res0 <- spatialdecon(
#'   norm = mini_geomx_dataset$normalized,
#'   bg = mini_geomx_dataset$bg,
#'   X = safeTME
#' )
#' # run reverse decon:
#' rdecon <- reverseDecon(
#'   norm = mini_geomx_dataset$norm,
#'   beta = res0$beta
#' )
#' @export
reverseDecon <- function(norm, beta, epsilon = NULL) {

    # remove cell types with no SD:
    beta <- beta[apply(beta, 1, stats::sd) > 0, ]
    # remove cell types that get removed by lm()
    # (presumably removed due to linear dependence)
    lm1 <- stats::lm(norm[1, ] ~ t(beta))
    beta <- beta[!is.na(lm1$coef[setdiff(names(lm1$coef), "(Intercept)")]), ,
        drop = FALSE
    ]

    # run reverse decon for all genes:
    rd <- function(y) {
        fit <- suppressWarnings(
            logNormReg::lognlm(y ~ t(beta),
                lik = FALSE,
                method = "L-BFGS-B",
                lower = rep(0, ncol(beta) + 1),
                upper = rep(Inf, ncol(beta) + 1),
                opt = "optim",
                control = list(maxit = 1000)
            )
        )
        return(fit$coefficients)
    }
    coefs <- t(apply(norm, 1, rd))
    colnames(coefs)[-1] <- rownames(beta)

    # get yhat
    yhat <- norm * NA
    for (ind in seq_len(ncol(yhat))) {
        yhat[, ind] <- coefs %*% c(1, beta[, ind])
    }

    # auto-select a reasonable epsilon if not provided
    if (length(epsilon) == 0) {
        epsilon <- stats::quantile(norm[norm > 0], 0.01)
    }

    # get resids:
    resids <- log2(pmax(norm, epsilon)) - log2(pmax(yhat, epsilon))

    # get summary stats:
    cors <- suppressWarnings(diag(stats::cor(t(norm), t(yhat))))
    resid.sd <- apply(resids, 1, stats::sd)

    out <- list(
        coefs = coefs, yhat = yhat, resids = resids, cors = cors,
        resid.sd = resid.sd
    )
    return(out)
}

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.




#' Apply error model to estimate technical SD from raw counts
#'
#' Based on raw counts, uses past data to estimate each raw count's log-scale
#' SD from technical noise.
#' Specifies different error models for different platforms.
#'
#' @param counts vector or matrix of raw counts
#' @param platform String specifying which platform was used to create
#' "rawCounts". Default to "dsp".
#'  Other options include "ncounter", "rsem" and "quantile".
#' @return a matrix of log2-scale SDs
runErrorModel <- function(counts, platform = "general") {
    if (platform == "ncounter") {
        sds <- counts * 0 + 0.1
        sds <- replace(sds, counts < 200, 0.2)
        sds <- replace(sds, counts < 100, 0.3)
        sds <- replace(sds, counts < 75, 0.4)
        sds <- replace(sds, counts < 50, 0.5)
        sds <- replace(sds, counts < 40, 0.7)
        sds <- replace(sds, counts < 30, 1)
        sds <- replace(sds, counts < 20, 3)
    }


    if (platform == "rsem") {
        sds <- counts * 0 + 0.5930982
        sds <- replace(sds, log2(counts) < 9.5, 0.6458475)
        sds <- replace(sds, log2(counts) < 8.5, 0.7847597)
        sds <- replace(sds, log2(counts) < 7.5, 1.0576471)
        sds <- replace(sds, log2(counts) < 6.5, 1.2990917)
        sds <- replace(sds, log2(counts) < 5.5, 1.5061735)
        sds <- replace(sds, log2(counts) < 4.5, 1.6930872)
        sds <- replace(sds, log2(counts) < 3.5, 1.7894239)
    }

    if (platform == "dsp") {
        predictsd.dsp <- function(rawcounts) {
            m <- log2(pmax(rawcounts, 1e-3))
            meanvec <- c(-1e-6, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, Inf)
            sdvec <- c(
                1.5, 1.383, 1.191, 0.800, 0.48, 0.301, 0.301,
                0.301, 0.263, 0.235, 0.235
            )

            s <- replace(m, TRUE, sdvec[1])
            for (i in seq_len(length(meanvec) - 1)) {
                s <- replace(s, m >= meanvec[i], sdvec[i + 1])
            }
            return(s)
        }

        if (is.vector(counts)) {
            sds <- vapply(
                X = counts,
                FUN = predictsd.dsp,
                FUN.VALUE = numeric(length(counts))
            )
        }
        if (is.matrix(counts)) {
            sds <- predictsd.dsp(counts)
        }
    }


    if (platform == "quantile") {
        if (is.vector(counts)) {
            quantile <- rank(counts) / length(counts)
        }
        if (is.matrix(counts)) {
            quantile <- matrix(rank(counts), nrow(counts)) / length(counts)
        }

        sds <- quantile * 0 + 0.1
        sds <- replace(sds, quantile < 0.2, 0.2)
        sds <- replace(sds, quantile < 0.15, 0.3)
        sds <- replace(sds, quantile < 0.1, 0.4)
        sds <- replace(sds, quantile < 0.05, 0.5)
        sds <- replace(sds, quantile < 0.01, 1)
    }
    return(sds)
}

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.



#' Mixed cell deconvolution of spatiall-resolved gene expression data
#'
#' Runs the spatialdecon algorithm with added optional functionalities.
#' Workflow is:
#' \enumerate{
#' \item compute weights from raw data
#' \item Estimate a tumor profile and merge it into the cell profiles matrix
#' \item run deconvolution once
#' \item remove poorly-fit genes from first round of decon
#' \item re-run decon with cleaned-up gene set
#' \item combine closely-related cell types
#' \item compute p-values
#' \item rescale abundance estimates, to proportions of total, proportions of
#'  immune, cell counts
#' }
#'
#' @param norm p-length expression vector or p * N expression matrix - the
#' actual (linear-scale) data
#' @param bg Same dimension as norm: the background expected at each data point.
#' @param X Cell profile matrix. If NULL, the safeTME matrix is used.
#' @param raw Optional for using an error model to weight the data points.
#'  p-length expression vector or p * N expression matrix - the raw
#'  (linear-scale) data
#' @param wts Optional, a matrix of weights.
#' @param resid_thresh A scalar, sets a threshold on how extreme individual data
#'  points' values
#'  can be (in log2 units) before getting flagged as outliers and set to NA.
#' @param lower_thresh A scalar. Before log2-scale residuals are calculated,
#'  both observed and fitted
#'  values get thresholded up to this value. Prevents log2-scale residuals from
#'  becoming extreme in
#'  points near zero.
#' @param align_genes Logical. If TRUE, then Y, X, bg, and wts are row-aligned
#'  by shared genes.
#' @param is_pure_tumor A logical vector denoting whether each AOI consists of
#'  pure tumor. If specified,
#'  then the algorithm will derive a tumor expression profile and merge it with
#'  the immune profiles matrix.
#' @param cell_counts Number of cells estimated to be within each sample. If
#' provided alongside norm_factors,
#'  then the algorithm will additionally output cell abundance esimtates on the
#'  scale of cell counts.
#' @param cellmerges A list object holding the mapping from beta's cell names to
#'  combined cell names. If left
#'  NULL, then defaults to a mapping of granular immune cell definitions to
#'   broader categories.
#' @param n_tumor_clusters Number of tumor-specific columns to merge into the
#' cell profile matrix.
#'  Has an impact only when is_pure_tumor argument is used to indicate pure
#'   tumor AOIs.
#'  Takes this many clusters from the pure-tumor AOI data and gets the average
#'  expression profile in each cluster.  Default 10.
#' @param maxit Maximum number of iterations. Default 1000.
#' @return a list:
#' \itemize{
#' \item beta: matrix of cell abundance estimates, cells in rows and
#' observations in columns
#' \item sigmas: covariance matrices of each observation's beta estimates
#' \item p: matrix of p-values for H0: beta == 0
#' \item t: matrix of t-statistics for H0: beta == 0
#' \item se: matrix of standard errors of beta values
#' \item prop_of_all: rescaling of beta to sum to 1 in each observation
#' \item prop_of_nontumor: rescaling of beta to sum to 1 in each observation,
#' excluding tumor abundance estimates
#' \item cell.counts: beta rescaled to estimate cell numbers, based on
#' prop_of_all and nuclei count
#' \item beta.granular: cell abundances prior to combining closely-related
#' cell types
#' \item sigma.granular: sigmas prior to combining closely-related cell types
#' \item cell.counts.granular: cell.counts prior to combining closely-related
#' cell types
#' \item resids: a matrix of residuals from the model fit.
#'  (log2(pmax(y, lower_thresh)) - log2(pmax(xb, lower_thresh))).
#' \item X: the cell profile matrix used in the decon fit.
#' }
#' @examples
#' data(mini_geomx_dataset)
#' # estimate background:
#' mini_geomx_dataset$bg <- derive_GeoMx_background(
#'   norm = mini_geomx_dataset$normalized,
#'   probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
#'   negnames = "NegProbe"
#' )
#' # run basic decon:
#' res0 <- spatialdecon(
#'   norm = mini_geomx_dataset$normalized,
#'   bg = mini_geomx_dataset$bg,
#'   X = safeTME
#' )
#' # run decon with bells and whistles:
#' res <- spatialdecon(
#'   norm = mini_geomx_dataset$normalized,
#'   bg = mini_geomx_dataset$bg,
#'   X = safeTME,
#'   cellmerges = safeTME.matches,
#'   cell_counts = mini_geomx_dataset$annot$nuclei,
#'   is_pure_tumor = mini_geomx_dataset$annot$AOI.name == "Tumor"
#' )
#' @export
spatialdecon <- function(norm, bg, X = NULL,
                         raw = NULL, wts = NULL,
                         resid_thresh = 3, lower_thresh = 0.5,
                         align_genes = TRUE,
                         is_pure_tumor = NULL, n_tumor_clusters = 10,
                         cell_counts = NULL,
                         cellmerges = NULL,
                         maxit = 1000) {

    #### preliminaries ---------------------------------

    # check formatting:
    if (!is.matrix(norm)) {
        stop("norm should be a matrix")
    }
    if ((length(X) > 0) & (!is.matrix(X))) {
        stop("X should be a matrix")
    }
    if ((length(raw) > 0) & (!is.matrix(raw))) {
        stop("raw must be numeric")
    }
    if ((length(wts) > 0) & (!is.matrix(wts))) {
        stop("wts must be numeric")
    }
    if ((length(cell_counts) > 0) & (!is.numeric(cell_counts))) {
        stop("cell_counts must be numeric")
    }


    if (length(bg) == 1) {
        bg <- matrix(bg, nrow(norm), ncol(norm),
            dimnames = list(rownames(norm), colnames(norm))
        )
    }

    # prep training matrix:
    if (length(X) == 0) {
        X <- safeTME
    }
    sharedgenes <- intersect(rownames(norm), rownames(X))
    if (length(sharedgenes) == 0) {
        stop("no shared gene names between norm and X")
    }
    if (length(sharedgenes) < 100) {
        stop(paste0(
            "Only ", length(sharedgenes),
            " genes are shared between norm and X - this may not be enough
                to support accurate deconvolution."
        ))
    }

    # calculate weights based on expected SD of counts
    # wts = replace(norm, TRUE, 1)
    if (length(raw) > 0) {
        weight.by.TIL.resid.sd <-
            length(intersect(colnames(X), colnames(safeTME))) > 10
        wts <- deriveWeights(norm,
            raw = raw, error.model = "dsp",
            weight.by.TIL.resid.sd = weight.by.TIL.resid.sd
        )
    }

    #### if pure tumor AOIs are specificed, get tumor expression profile --------
    if (sum(is_pure_tumor) > 0) {

        # derive tumor profiles and merge into X:
        # (derive a separate profile for each tissue)
        X <- mergeTumorIntoX(
            norm = norm,
            bg = bg,
            pure_tumor_ids = is_pure_tumor,
            X = X[sharedgenes, ],
            K = n_tumor_clusters
        )

        sharedgenes <- intersect(rownames(norm), rownames(X))
    }


    #### Run decon  -----------------------------------
    res <- algorithm2(
        Y = norm[sharedgenes, ],
        bg = bg[sharedgenes, ],
        X = X[sharedgenes, ],
        weights = wts[sharedgenes, ],
        maxit = maxit
    )


    #### combine closely-related cell types ------------------------------------

    if (length(cellmerges) > 0) {
        tempconv <- convertCellTypes(
            beta = res$beta,
            matching = cellmerges,
            stat = sum,
            na.rm = FALSE,
            sigma = res$sigmas
        )
        # overwrite original beta with merged beta:
        res$beta.granular <- res$beta
        res$sigma.granular <- res$sigmas
        res$sigmas <- NULL
        res$beta <- tempconv$beta
        res$sigma <- tempconv$sigma
    }


    #### compute p-values -------------------------------------------
    tempbeta <- res$beta
    tempse <- tempp <- tempt <- tempbeta * NA
    for (i in seq_len(ncol(tempse))) {
        tempse[, i] <- suppressWarnings(sqrt(diag(res$sigma[, , i])))
    }
    tempt <- (tempbeta / tempse)
    tempp <- 2 * (1 - stats::pt(tempt, df = length(sharedgenes) - ncol(X) - 1))
    res$p <- tempp
    res$t <- tempt
    res$se <- tempse


    #### rescale abundance estimates --------------------------------
    # (to proportions of total, proportions of immune, cell counts)

    # proportions:
    res$prop_of_all <- sweep(res$beta, 2, colSums(res$beta), "/")
    nontumorcellnames <- rownames(res$beta)[!grepl("tumor", rownames(res$beta))]
    res$prop_of_nontumor <- sweep(
        res$beta[nontumorcellnames, ], 2,
        colSums(res$beta[nontumorcellnames, ]), "/"
    )

    # on scale of cell counts:
    if (length(cell_counts) > 0) {
        res$cell.counts <- convertCellScoresToCounts(
            beta = res$beta,
            nuclei.counts = cell_counts,
            omit.tumor = TRUE
        )
        if (exists("res$beta.granular") > 0) {
            res$cell.counts.granular <- convertCellScoresToCounts(
                beta = res$beta.granular,
                nuclei.counts = cell_counts,
                omit.tumor = TRUE
            )
        }
    }

    # add other pertinent info to res:
    res$X <- X[rownames(res$resids), ]
    return(res)
}

# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression
# data
#    This program is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation, either version 3 of the License, or (at your option)
#    any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#    more details.
#    You should have received a copy of the GNU General Public License along
#    with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.


#' Barplot of abundance estimates
#'
#' Draw barplot of the "betas" from a decon fit
#'
#' @param mat Matrix of cell proportions or abundances, in the same dimensions
#' output by spatialdecon
#'  (cells in rows, observations in columns). User is free to re-order
#'  columns/observations in
#'  whatever order is best for display.
#' @param draw_legend Logical. If TRUE, the function draws a legend in a new
#' plot frame.
#' @param main Title for barplot
#' @param col Vector of colors for cell types. Defaults to pre-set colors for
#' the safeTME cell types.
#' @param ... Arguments passed to barplot()
#' @return Draws a barplot.
#' @examples
#' data(mini_geomx_dataset)
#' # estimate background:
#' mini_geomx_dataset$bg <- derive_GeoMx_background(
#'   norm = mini_geomx_dataset$normalized,
#'   probepool = rep(1, nrow(mini_geomx_dataset$normalized)),
#'   negnames = "NegProbe"
#' )
#' # run basic decon:
#' res0 <- spatialdecon(
#'   norm = mini_geomx_dataset$normalized,
#'   bg = mini_geomx_dataset$bg,
#'   X = safeTME
#' )
#' # run barplot:
#' TIL_barplot(mat = res0$beta)
#' # run barplot and draw a color legend
#' TIL_barplot(mat = res0$beta, draw_legend = TRUE)
#' @export
TIL_barplot <- function(mat, draw_legend = FALSE, main = "", col = NULL, ...) {


    # infer colors:
    if (length(col) == 0) {
        # use safeTME colors if the right cells are present:
        if (all(is.element(rownames(mat), names(cellcols)))) {
            col <- cellcols[rownames(mat)]
        }
        else {
            manycols <- c(
                "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
                "#B3DE69", "#FCCDE5", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
                "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#1B9E77", "#D95F02",
                "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
                sample(grDevices::colors(), 99)
            )
            col <- manycols[seq_len(nrow(mat))]
            names(col) <- rownames(mat)
        }
    }


    #  usecells <- intersect(rownames(mat), names(cellcols))
    usecells <- rownames(mat)

    # draw barplot:
    graphics::barplot(mat[usecells, ],
        cex.lab = 1.5,
        col = col, border = NA,
        las = 2, main = main, ...
    )

    # draw a legend:
    if (draw_legend) {
        graphics::frame()
        graphics::legend("center",
            fill = rev(col),
            legend = rev(names(col))
        )
    }
}

#### functions for drawing points in space and polygons behind them

# contents:
# - spaceplot: a fn to draw expression values in space. Calls the below two functoins
# - getBoundary: fn to infer a polygonal shape around a set of points
# - makeTissueGrid: fn to arrange all the xy points from different tissues in a study in a non-overlapping way



# dev note: what's needed here:
# - 0 and tiny values get an empty pch = 1


#' Plot a gene or score in space
#'
#' Wrapper function for other spatial plotting functions.
#' Draws plots the x-y coords of tissues, with point size giving variable value.
#' @param x X coords
#' @param y Y coords
#' @param z A non-negative vector of expression levels, pathway/cell scores, etc.
#' @param tissue A vector of tissue IDs
#' @param tissue.order Optional, vector of tissue names giving the order in which to plot them.
#'  If NULL, then alphabetical order will be used
#' @param tissuecols Named vector of colors to use for each tissue's polygon
#' @param use Vector of logicals specifying which elements of the data to use
#' @param rescale Logical for whether to rescale to give z a mean of 1.
#' @param cex A scalar, controls point size. (Points are further scaled by the values of z.)
#' @param boundaries Optional, a list of x,y vectors defining polygons, in format created by getBoundary()
#' @param ... Arguments passed to plot()
#' @return Draws a plot in xy space. Returns a list:
#'  x: new x coords
#'  y: new y coords,
#'  outlines: list of polygonal tissue boundaries' x-y coords
#' @example
#' # sim data from 5 tissues:
#' set.seed(0)
#' x = rnorm(50)
#' y = rnorm(50)
#' z = runif(50)
#' tissue = rep(letters[1:5], each = 10)
#' spaceplot(x, y, z, tissue, tissuecols = NULL, use = TRUE, rescale = FALSE,
#'           cex = 1, col = "#00008B80",
#'           nrow = NULL, rowind = NULL, colind = NULL, expansion = 1.2)
spaceplot <- function(x, y, z, tissue, tissue.order = NULL, tissuecols = NULL, tissuecols.alpha = 0.2, use = TRUE, rescale = FALSE,
                      cex = 1, col = "#00008B80",
                      nrow = NULL, rowind = NULL, colind = NULL, expansion = 1.2, ...) {

  # subset everything:
  if (length(col) == length(x)) {
    col <- col[use]
  }
  x <- x[use]
  y <- y[use]
  z <- z[use]
  tissue <- tissue[use]

  # first, get new xy coords:
  newxy.custom <- makeTissueGrid(
    x = x, y = y,
    tissue = tissue, tissue.order = tissue.order,
    rowind = rowind, colind = colind,
    nrow = NULL, expansion = expansion
  )
  x <- newxy.custom$x
  y <- newxy.custom$y

  # also get polygon tissue boundaries:
  boundaries <- list()
  xrange <- yrange <- c(NA, NA)
  for (tiss in unique(tissue)) {
    tempuse <- tissue == tiss
    boundaries[[tiss]] <- getBoundary(x = x[tempuse], y = y[tempuse], marg = 0.4)
    yrange[1] <- min(yrange, boundaries[[tiss]]$y, na.rm = T)
    yrange[2] <- max(yrange, boundaries[[tiss]]$y, na.rm = T)
    xrange[1] <- min(xrange, boundaries[[tiss]]$x, na.rm = T)
    xrange[2] <- max(xrange, boundaries[[tiss]]$x, na.rm = T)
  }

  # transform z:
  if (rescale) {
    z <- z / max(z)
  }

  # plotting:
  plot(x, y,
    cex = z * cex,
    col = col,
    pch = 16,
    xlab = "", ylab = "",
    xaxt = "n", yaxt = "n",
    ylim = yrange, xlim = xrange, ...
  )

  # tissue boundaries:
  if (length(tissuecols) == 0) {
    tissuecols <- rep("grey50", length(unique(tissue)))
    names(tissuecols) <- unique(tissue)
  }
  for (tiss in names(boundaries)) {
    polygon(
      x = boundaries[[tiss]]$x,
      y = boundaries[[tiss]]$y,
      col = alpha(tissuecols[[tiss]], tissuecols.alpha), border = NA
    )
  }

  # return coordinates:
  out <- list(x = x, y = y, boundaries = boundaries)
  return(out)
}






#' Function to define a polygon boundary around a set of points in xy space
#'
#' Fits a convex hull polygon around the points, leaving a bit of margin beyond the points
#' @param x vector of x coords
#' @param y vector of y coords
#' @param marg Amount of extra margin to draw. Default 10%
#' @return a list: x: coords of the polygon. y: y coords on the polygon
#' @example
#'  x = rnorm(30)
#'  y = rnorm(30)
#'  plot(x, y)
#'  bound = getBoundary(x, y)
#'  polygon(bound, col = rgb(0,0,1,0.5))
getBoundary <- function(x, y, marg = 0.2) {
  # get center of shape:
  mx <- mean(x)
  my <- mean(y)
  # expand all points away from center
  x2 <- mx + (x - mx) * (1 + marg)
  y2 <- my + (y - my) * (1 + marg)

  # get convex hull around expanded points:
  ch <- chull(x2, y2)

  out <- list(
    x = x2[ch],
    y = y2[ch]
  )
  return(out)
}


#' Define non-overlapping x-y coords for plotting multiple tissues in the same frame
#'
#' Given xy coords for segments from several tissues, shifts each tissue so they don't overlap
#' @param x Vector of x coords
#' @param y Vector of y coords
#' @param tissue Vector of tissue IDs corresponding to x and y
#' @param tissue.order Optional, vector of tissue names giving the order in which to plot them.
#'  If NULL, then alphabetical order will be used
#' @param rowind Optional, integer vector of row assignments for tissues, aligned to tissue.order.
#'  Must be provided with colind to work.
#' @param colind Optional, integer vector of column assignments for tissues, aligned to tissue.order
#'  Must be provided with rowind to work.
#' @param expansion A constant scaling factor, for how much to expand the margins between tissues
#' @return A list: x = new x coords. y = new y coords.
#' @examples
#' # sim data from 5 tissues:
#' set.seed(0)
#' x <- rnorm(50)
#' y <- rnorm(50)
#' tissue <- rep(letters[1:5], each = 10)
#'
#' # arrange in default layout:
#' newxy <- makeTissueGrid(
#'   x = x, y = y, tissue = tissue,
#'   tissue.order = NULL, rowind = NULL, colind = NULL,
#'   nrow = NULL, expansion = 1.2
#' )
#' plot(newxy, pch = 16, col = 0)
#' text(newxy$x, newxy$y, tissue, col = as.numeric(as.factor(tissue)))
#'
#' # specify a layout:
#' newxy.custom <- makeTissueGrid(
#'   x = x, y = y, tissue = tissue,
#'   rowind = c(1, 1, 1, 2, 2), colind = c(1, 2, 3, 1, 3),
#'   nrow = NULL, expansion = 1.2
#' )
#' plot(newxy.custom, pch = 16, col = 0)
#' text(newxy.custom$x, newxy.custom$y, tissue, col = as.numeric(as.factor(tissue)))
makeTissueGrid <- function(x, y, tissue,
                           tissue.order = NULL, rowind = NULL, colind = NULL,
                           nrow = NULL, expansion = 1.2) {

  # define the tissue ids and their order:
  if (length(tissue.order) > 0) {
    tissues <- tissue.order
  }
  else {
    tissues <- sort(unique(tissue))
  }

  # get tissue spans and centers:
  xspans <- yspans <- xcenters <- ycenters <- c()
  for (tiss in tissue) {
    rx <- range(x[tissue == tiss])
    ry <- range(y[tissue == tiss])
    xspans[tiss] <- diff(rx)
    yspans[tiss] <- diff(ry)
    xcenters[tiss] <- median(rx)
    ycenters[tiss] <- median(ry)
  }

  # define how far apart tissues should be in x and y space:
  xmarg <- max(xspans) * expansion
  ymarg <- max(yspans) * expansion

  # if not specified, define a grid of tissue centers:
  if ((length(rowind) == 0) | (length(colind) == 0)) {
    # define number of rows and columns:
    if (length(nrow) == 0) {
      nrow <- round(sqrt(length(tissues)))
    }
    ncol <- ceiling(length(tissues) / nrow)

    # assign rows and columns to each tissue:
    rowind <- ceiling((1:length(tissues)) / ncol)
    colind <- ceiling(((1:length(tissues))) %% ncol)
    colind <- replace(colind, colind == 0, ncol)
    # plot(rowind ~ colind);text(colind, rowind, tissues)
  }
  rowind <- max(rowind) - rowind + 1
  names(rowind) <- names(colind) <- tissues

  # now get xy offsets for each tissue
  xnew <- ynew <- replace(x, TRUE, NA)
  for (tiss in tissues) {
    tempxoffset <- colind[tiss] * xmarg - xcenters[tiss]
    tempyoffset <- rowind[tiss] * ymarg - ycenters[tiss]
    xnew[tissue == tiss] <- x[tissue == tiss] + tempxoffset
    ynew[tissue == tiss] <- y[tissue == tiss] + tempyoffset
  }

  out <- list(x = xnew, y = ynew)
  return(out)
}


mean.resid.sd <- c(
  1.577481401, 1.063036699, 1.202461016, 1.367237565, 1.410840762, 0.67591047, 1.220860563, 1.718933941, 1.303709992,
  1.259962537, 0.97212259, 1.4643764, 1.077864546, 1.174264436, 1.249325739, 0.909423981, 0.920114445, 1.006490125,
  1.258907901, 1.007493076, 0.940172986, 1.029183764, 1.469675992, 1.36356, 1.11092534, 1.299970315, 1.449455023,
  1.417454962, 1.123040125, 0.796768287, 0.573463116, 0.971067741, 1.233350143, 0.51832697, 1.105882294, 0.902535376,
  1.1605874, 1.360926845, 0.958075117, 1.146307688, 1.277235724, 1.354120935, 1.044601231, 1.188103271, 1.437049058,
  1.172566966, 0.918462213, 1.204440293, 1.011769123, 1.472673472, 1.266894665, 1.170205867, 1.059595804, 1.228521064,
  1.292023877, 1.073624528, 1.39595114, 0.778459441, 1.134830119, 1.520495321, 0.992562959, 1.131868233, 1.867585237,
  1.517617727, 1.013708238, 1.249759445, 1.570486101, 1.436094471, 1.36992249, 0.940075956, 1.023223048, 1.396310081,
  1.312143259, 1.419502189, 1.721734486, 1.165552358, 1.236673916, 1.25121613, 1.305309903, 0.936961311, 1.193286365,
  1.605876775, 1.220204411, 1.288371649, 1.028220374, 0.909285523, 0.966030369, 1.238257068, 0.995385915, 1.821200536,
  1.801160326, 1.750147235, 1.666706108, 1.504455857, 1.62063494, 1.332204254, 1.418750092, 1.261096187, 0.705344539,
  1.311575021, 0.827287134, 1.211839577, 1.095681399, 1.34357463, 1.188069364, 1.529732287, 1.065342658, 1.074453771,
  1.207828007, 0.960043932, 1.178987777, 1.315995177, 1.115073259, 0.877854218, 0.930084886, 1.137969477, 1.39418916,
  1.300118927, 1.325520239, 1.147939484, 0.440804715, 1.02171245, 1.087993929, 0.94294609, 1.288487132, 1.558632504,
  1.307828215, 1.65515318, 1.259022228, 0.771041818, 1.449401969, 1.339161479, 1.71608494, 0.911048932, 1.022336526,
  1.007857951, 1.423608861, 1.246143898, 1.176044433, 1.255504774, 0.987505095, 0.948945156, 1.38641615, 1.372336595,
  1.128771429, 1.467254101, 1.106608205, 1.574170158, 1.447836026, 1.287672115, 1.51286675, 1.433942152, 0.932152842,
  1.476188723, 1.170114286, 1.259560876, 1.585972434, 1.161052339, 0.103941778, 1.451116848, 1.385200724, 1.584160788,
  1.154038999, 1.412598944, 1.459491243, 1.514566413, 1.394193217, 1.375287138, 1.149042768, 0.96594859, 1.205704511,
  1.273861002, 1.079220286, 1.539934467, 1.244810539, 1.433701928, 1.623277213, 1.13918148, 1.409251588, 1.227752489,
  1.165592745, 1.47525718, 1.558673118, 1.501158225, 1.586362384, 1.546050539, 1.287041858, 1.486709964, 1.290536177,
  1.424513068, 1.953222313, 1.829566154, 0.907721409, 1.16338451, 1.43575358, 1.749554381, 1.543601509, 1.310876646,
  1.430167941, 1.136064073, 1.436997821, 1.416257506, 1.005553273, 1.165344935, 0.90164889, 1.24758323, 1.623868165,
  1.227956847, 1.24865627, 1.271939041, 1.345270077, 0.956369079, 1.018222645, 1.663891084, 1.348318411, 1.437417157,
  1.112257144, 1.892520021, 1.359259026, 1.16101317, 1.51021356, 1.230470984, 1.210910186, 1.524757183, 0.931457738,
  1.223199376, 1.26924158, 1.126480343, 1.031287703, 1.090056102, 0.606324542, 0.815152144, 0.998383818, 1.342805208,
  0.655015479, 1.284263675, 0.9049824, 0.881454384, 1.09043499, 0.843058353, 1.078335332, 1.201966201, 1.045499158,
  0.774579699, 1.165981314, 1.569138771, 1.297158281, 1.74660299, 1.540376757, 0.946667343, 1.32335967, 1.247582385,
  1.878899702, 1.857748195, 1.124153579, 1.534142594, 1.471595983, 1.253995611, 1.620601252, 1.574379208, 1.367192971,
  1.692613554, 1.729455351, 1.228789587, 1.200688049, 1.135530197, 1.100890614, 1.792693005, 1.380770408, 1.332946406,
  1.572561258, 1.310933212, 1.500771444, 1.325872809, 1.146853753, 0.846105951, 1.144204732, 0.767875744, 1.137239249,
  0.579208085, 1.08707585, 0.798817484, 1.559755416, 1.131514863, 1.55475927, 1.12357963, 1.410184555, 0.707778154,
  1.477811854, 1.482372738, 1.415538296, 1.425226429, 1.383692385, 1.477613334, 1.547880045, 1.463920647, 0.752459302,
  1.352064842, 1.31124381, 1.201125323, 1.117438879, 1.337032615, 1.477481575, 1.50370472, 1.490162234, 1.526000759,
  1.177004769, 1.287306433, 1.419863948, 1.258770214, 1.540261413, 1.58548279, 1.014412096, 1.2766158, 1.438912667,
  1.57196279, 1.089151415, 1.147196886, 1.376921501, 1.510803409, 1.15027614, 1.833374163, 1.406579896, 1.098100092,
  1.607186259, 0.962778659, 1.032498982, 1.403550379, 1.022834304, 1.325091774, 1.249348461, 1.630445586, 0.123384014,
  1.618939371, 0.240668056, 1.404683402, 0.95988365, 1.445672394, 1.051814618, 1.226648125, 1.371625096, 0.689211763,
  0.900140605, 1.009772211, 0.153843592, 2.039415531, 1.332207716, 1.006521856, 1.295409395, 1.209472589, 1.42346539,
  0.846003445, 1.100764583, 1.493953404, 0.994403258, 1.24208501, 1.5224595, 0.022578704, 1.021366102, 1.059178705,
  0.818839138, 1.543761053, 1.115244943, 1.280373101, 1.968219747, 1.271223773, 1.421210049, 1.134746702, 1.443893508,
  0.720701016, 1.045430381, 1.029296811, 1.116675081, 1.412383093, 1.434077318, 1.391448963, 1.549486519, 1.459813865,
  1.144367152, 1.332856806, 0.271050034, 1.705735944, 1.049812647, 1.435167483, 1.457480651, 1.521251521, 1.359991208,
  1.440904544, 0.694766497, 0.976756481, 0.581001343, 1.447827858, 1.29399291, 1.245104293, 1.078213647, 1.072155778,
  1.319718642, 1.037593432, 1.145481921, 1.584909821, 0.965845913, 1.527301337, 1.25109383, 1.840007561, 1.388470998,
  1.277612783, 1.344677898, 1.507741408, 1.135947484, 1.230591925, 1.23395835, 1.144364851, 2.265100043, 1.23804809,
  1.200185611, 1.298425602, 1.075764125, 1.432838811, 1.113810745, 1.559414622, 1.086558905, 1.301430891, 1.005542817,
  1.068481914, 1.300816083, 1.793939839, 1.475971838, 1.212388065, 1.279015961, 1.56970427, 1.254101252, 1.020843557,
  1.466232175, 1.341406459, 1.258439356, 1.40667883, 1.180640778, 1.203687299, 1.163165856, 1.328086204, 0.630420704,
  1.301804228, 1.309262499, 1.211976157, 1.287208297, 1.127168436, 1.112587397, 1.036944079, 1.321083935, 1.052447574,
  1.69112741, 1.385896243, 1.101881312, 1.055103438, 1.166441445, 1.872989547, 1.354367303, 1.492762648, 1.531201505,
  1.138132923, 1.510602806, 1.625404871, 0.821564058, 1.179816169, 0.9591359, 1.407139192, 1.671672666, 0.916006939,
  0.675520268, 1.02241593, 1.132708363, 1.060304641, 1.240220069, 1.40236613, 1.509869919, 1.332627653, 0.759530807,
  0.808911307, 1.544755905, 1.576936996, 1.338739482, 1.474934915, 1.383508948, 1.479553148, 1.204968066, 1.313838142,
  1.440955782, 1.279711202, 1.199077415, 1.208623412, 1.003730272, 1.106365212, 1.238125095, 1.183015197, 1.217492108,
  1.192149398, 1.304421453, 1.780783184, 1.559611755, 1.52161717, 1.563930563, 1.463222316, 1.590179142, 1.607075521,
  1.515378809, 1.136817609, 1.263762911, 1.175958752, 1.128217248, 1.399204848, 1.173957204, 1.184271139, 1.328562863,
  1.477519675, 1.140488485, 0.79626397, 1.125185344, 1.066303306, 1.288303714, 0.652769896, 0.708756335, 1.007770671,
  0.844385465, 1.301452491, 1.174790441, 1.181517974, 1.250820462, 1.30952121, 1.177910564, 0.721665136, 1.435999561,
  1.496272233, 1.206510899, 1.419019876, 1.199710289, 1.313291149, 1.323831585, 0, 0.95654704, 0.985321196, 1.489696459,
  1.543696824, 1.309524419, 1.234886913, 0.666096565, 1.398400278, 0.657683252, 1.005077345, 1.437006841, 0.574438123,
  0.900771015, 0.978384545, 1.278829783, 1.34007597, 1.34900327, 1.308160763, 0.891544437, 1.409684603, 0.200052255,
  1.179369155, 1.012096952, 0.967175163, 0.678131312, 1.184524426, 1.606707232, 0.000743488, 0.994552928, 1.430118813,
  1.012117808, 1.025095397, 0.24746253, 1.273861188, 0.279371195, 0.934530921, 0.99539256, 1.278039958, 1.102525503,
  1.075692203, 0.000891661, 1.119159852, 1.374578929, 0.790976289, 1.496462308, 1.518098582, 1.095476966, 1.305841997,
  1.357200274, 1.400088025, 1.330094342, 1.371492539, 1.522335247, 1.599325549, 1.144175609, 1.47513834, 1.634296078,
  1.356598711, 1.381311453, 1.459403023, 1.021432628, 0.890233474, 1.150025995, 1.40250287, 0.71966334, 0.338766966,
  1.973299825, 1.253958199, 0.91108129, 1.232723797, 0.96541708, 1.148551398, 0.78881444, 1.333658887, 1.264839955,
  0.528948514, 1.433746784, 1.005490974, 1.189013376, 0.663268569, 0.888386238, 0.962839287, 0.809445751, 0.864213101,
  0.621004685, 0.975295559, 1.269876654, 0.493797164, 1.295510238, 1.463441881, 0.594323996, 0.821918641, 1.697295652,
  1.210231053, 1.001970192, 1.237578903, 0.91697838, 0.902605353, 1.06268569, 1.348828245, 1.088227853, 1.486155323,
  1.47915698, 1.444306977, 1.434168022, 1.833588689, 1.422601419, 1.333559497, 0.514011244, 0.378583663, 1.376516151,
  1.311029728, 1.211945468, 1.413281918, 1.382221184, 1.123824229, 1.445234158, 1.147692168, 1.290794163, 1.34627919,
  1.265222126, 1.401369053, 0.975066496, 1.594219809, 1.362889683, 1.024342818, 1.301959307, 1.295840025, 1.066760483,
  1.205317856, 1.371688142, 1.507421123, 1.109430777, 1.373377366, 1.39982915, 0.466293467, 1.201740346, 1.111950348,
  0.924308073, 1.096144503, 1.361079922, 1.385328274, 1.342209911, 1.438439964, 0.928412197, 0.817455764, 1.202329666,
  1.200727454, 1.301061312, 1.344473405, 1.135775382, 0.909849756, 1.394070097, 1.236398106, 1.371908025, 0.955862596,
  1.427157288, 1.076116357, 1.326035973, 0.958824185, 1.519226341, 0.736585199, 1.138869264, 1.531811377, 1.65906135,
  1.517853867, 1.367844954, 1.155384958, 1.127595872, 0.855370591, 1.749548476, 1.353419645, 1.440343026, 1.347469126,
  1.001473088, 1.050221098, 1.284617523, 1.305784338, 1.253091866, 1.144361213, 1.33922808, 1.617199209, 1.268984315,
  1.365271749, 1.364819484, 1.613913852, 0.310440453, 1.248213706, 0.989479806, 1.31870228, 1.252976163, 0.920110586,
  1.383331323, 0.820396128, 1.274823211, 1.29359616, 1.562315339, 1.539128283, 0.643410504, 1.431774584, 1.428214823,
  1.569972562, 1.291036516, 1.614370533, 1.197474555, 1.050779073, 1.29995724, 1.362339187, 1.406648065, 1.156968851,
  1.311236236, 0.878116712, 1.175788612, 1.369037448, 2.162292103, 1.156231228, 1.18227028, 0.541472959, 1.416708495,
  1.492238861, 1.502736378, 1.277544309, 1.823921699, 1.466019423, 1.520265098, 1.367054808, 1.485721962, 0.811806602,
  0.945161133, 1.573035872, 1.020568101, 0.944873068, 1.485238358, 1.209931815, 1.409095311, 0.74644762, 1.029340043,
  1.08552754, 1.434410307, 1.522232471, 0.813016111, 1.200377892, 1.723817157, 0.943120832, 1.434742195, 1.092067724,
  1.058371116, 1.369853895, 1.538637545, 1.265651011, 0.107480022, 1.155796089, 1.237376806, 1.023410346, 1.204403674,
  1.189911763, 1.162500529, 0.06243859, 1.217025019, 1.443056251, 1.303961854, 0.983134596, 1.184826382, 1.313532965,
  1.372349837, 1.132544997, 1.166087606, 1.317971723, 1.338316636, 0.465127659, 1.215627858, 1.179937199, 1.01655281,
  1.291846298, 1.306141869, 1.302892928, 1.360096324, 0.641790739, 0.164168353, 0.057212783, 0.016185264, 1.271329425,
  1.241497703, 1.107399734, 1.446697232, 1.137320591, 1.152030195, 1.223540916, 1.142117607, 1.140171561, 1.155785293,
  0.927607027, 0.982465535, 1.333815193, 1.188579046, 1.240768782, 1.293283345, 0.995490225, 1.104639371, 1.169332393,
  1.125122463, 0.953540919, 1.265015594, 1.045356102, 1.225689173, 0.833202731, 1.220719077, 0.349440197, 1.212574215,
  1.189037517, 1.138332515, 1.418385221, 1.58025569, 1.234757752, 1.612607345, 1.409348779, 1.262271852, 0.436246892,
  1.143619976, 1.267073976, 1.611569905, 1.281763354, 0.591007901, 1.528196347, 1.18105755, 1.274072586, 1.456414613,
  1.252924223, 1.209457759, 1.44522785, 1.441879705, 1.238706881, 1.353620665, 1.3329761, 1.159430884, 1.428472535,
  1.265393735, 0.982819503, 0.838056404, 1.12855582, 1.332958282, 1.510955599, 1.264692695, 1.353911195, 1.201963253,
  1.693418303, 1.306368474, 1.291699301, 1.663494712, 1.418095966, 1.234179266, 1.320174427, 1.42837007, 1.05949392,
  1.556656023, 1.257351062, 1.45131609, 1.102666302, 0.279579016, 0.957774391, 1.05721677, 1.145137703, 0.969945427,
  1.381495667, 0.983879626, 1.945898264, 1.311384408, 0.222371337, 0.862244001, 1.313331217, 1.221539766, 1.110910218,
  1.379871962, 1.268425629, 1.468890321, 1.209202687, 1.10144052, 1.501951276, 1.089811797, 1.234379399, 1.399153004,
  0.856466332, 1.452805747, 1.190501436, 1.522956836, 1.00693471, 1.135400296, 1.554398873, 1.21501437, 1.076157579,
  1.520314752, 1.468723643, 0.971232819, 1.341264131, 1.143292158, 1.164704337, 1.569261554, 1.174645, 0.714411054,
  1.107001704, 1.465338291, 1.150913164, 1.049530248, 1.565522541, 0.726678832, 0.240951537, 1.360739177, 1.083446025,
  1.350761968, 1.572869625, 1.735157967, 1.281740768, 1.294512687, 1.076462326, 1.288832236, 1.229090564, 1.177452137,
  1.322240746, 1.074367875, 1.040137072, 0.894988808, 1.490516853, 0.949317821, 1.454802129, 1.429282626, 1.321542121,
  1.597831515, 1.565658963, 1.12495439, 0.904515272, 1.333960636, 1.216617635, 1.171310876, 1.343586139, 1.269616355,
  1.235851264, 1.39049654, 1.514703193, 1.466958036, 1.21060628, 1.029852035, 1.187487213, 1.36160251, 1.553700303,
  1.153230271, 1.228684366, 1.267712423, 0.615244328, 0.687187971, 1.316030478, 1.475858223, 0.380069232, 1.211266117,
  1.129857044, 1.147789795, 1.127144661, 1.297137368, 1.062814088, 1.047621592, 1.402211534, 1.305584009, 1.42247839,
  1.103657995, 1.014401839, 1.446228695, 0.449043346, 1.645977793, 0.726855958, 1.408374855, 1.29863375, 1.184416424,
  1.217320353, 1.385110228, 1.042224382, 1.228805272, 0.615812854, 1.169721182, 1.262682449, 0.58844309, 0.690692816,
  1.031242804, 1.356045097, 1.307040924, 1.297633075, 1.596745976, 1.446620222, 1.748079456, 1.312807865, 1.175481263,
  1.277257121, 1.291211765, 0.497645497, 0.730442851, 1.123187256, 1.507044672, 1.114459831, 1.308762262, 1.118486626,
  1.27449641, 1.09128827, 1.300085438, 1.152211716, 1.133235618, 0.853660821, 1.226943958, 1.214199883, 1.232446724,
  1.380792331, 1.142672099, 1.326431906, 1.681357774, 1.029403954, 1.413336398, 1.868536013, 0.599821846, 1.171726955,
  1.834922691, 0.817892015, 1.016401732, 0.944856956, 0.502742097, 1.564155013, 1.236841144, 1.413174785, 1.398726385,
  0.704210416, 1.200996063, 1.393976109, 1.545005202, 0.865460761, 1.144687407, 0.78117357, 1.514531879, 1.187600593,
  1.127090256, 1.288539408, 0.901906135, 1.539465252, 1.101644671, 1.124562568, 1.36687231, 1.144914114, 1.168758156,
  1.171859033, 1.01962954, 1.481420608, 1.132076084, 1.17189542, 1.120789294, 1.441422372, 1.454898793, 1.433281696,
  0.96569623, 0.9825766, 1.101295801, 1.602552871, 0.969804365, 1.162879837, 1.081990206, 1.145283425, 1.566260064,
  1.184522729, 0.958635813, 0.839954129, 1.229393193, 1.118962291, 1.227607595, 0.983163496, 1.082826966, 1.566962848,
  1.288896439, 1.046988428, 2.025674177, 1.971279444, 1.395272462, 1.038846945, 1.208795519, 1.48638255, 1.232969155,
  0.041697637, 1.394818997, 0.835639479, 0.801059316, 1.077637062, 1.27518101, 1.219420104, 0.997673684, 1.220394896,
  1.222828209, 1.071373298, 1.418951857, 0.619981506, 0.437257568, 1.599148642, 0.794195865, 0.532997103, 1.379407422,
  0.953630602, 1.730973249, 1.203055994, 1.506598505, 0.953461349, 1.438971903, 1.024243766, 1.27172747, 1.301156998,
  1.175343707, 1.131869734, 1.395285429, 1.496392949, 1.355647452, 1.42479394, 1.365792406, 1.207052606, 1.334451777,
  1.171923167, 1.177774054, 1.004071888, 0.940374404, 0.980308254, 1.335352942, 1.48343804, 1.3204372, 1.275575808,
  0.926708838, 1.268629573, 1.200164859, 1.348717903, 0.801366481, 1.298998951, 0.966527148, 0.989377371, 0.691569899,
  0.680521257, 1.148354277, 0.957133629, 1.00042614, 1.2098885, 0.784638682, 1.107204552, 1.561208594
)
names(mean.resid.sd) <- c(
  "A2M", "ABCB1", "ABCB4", "ABCC3", "ACAP1", "ACER1", "ACHE", "ACP5", "ACPP", "ADAM12", "ADAM23", "ADAM28",
  "ADAM33", "ADAMDEC1", "ADAMTS1", "ADAMTS10", "ADAMTS3", "ADAMTS5", "ADM", "ADORA3", "ADRB2", "AGBL3", "AIF1",
  "AIM2", "AK5", "AKAP12", "ALOX15", "ALOX5", "ALPL", "AMPD1", "ANGPT4", "ANGPTL1", "ANGPTL4", "ANKRD31", "ANKRD34B",
  "ANKRD55", "ANO5", "ANXA3", "AOC2", "APLNR", "APOBEC3A", "APOBEC3G", "APOBEC3H", "APOL3", "AQP1", "AQP9", "ARG1",
  "ARHGAP22", "ARHGEF15", "ARHGEF4", "ARSJ", "ASGR1", "ASGR2", "ASPN", "ATP1A3", "ATP8B4", "AXL", "AZU1", "B3GALT2",
  "B3GNT7", "BAALC", "BACH2", "BANK1", "BARX2", "BATF3", "BCAT1", "BCL11B", "BCL2", "BCL2A1", "BEND6", "BHLHA15",
  "BHLHE41", "BIRC3", "BLK", "BLNK", "BMP2", "BMP3", "BMP4", "BMP5", "BMX", "BNC2", "BRSK2", "BST1", "BTLA", "BTNL8",
  "C10orf105", "C12orf74", "C1orf127", "C1orf54", "C1QA", "C1QB", "C1QC", "C1S", "C2", "C3", "C3AR1", "C5AR1", "C7",
  "C8B", "C8G", "C9orf47", "CA4", "CA6", "CA8", "CACHD1", "CACNA1D", "CACNA1E", "CACNA2D1", "CACNA2D2", "CACNA2D3",
  "CACNB2", "CADM1", "CALCRL", "CALN1", "CAMP", "CASP5", "CAV1", "CAV2", "CBX2", "CCDC102B", "CCL1", "CCL11", "CCL13",
  "CCL14", "CCL17", "CCL18", "CCL19", "CCL20", "CCL22", "CCL23", "CCL28", "CCL4", "CCL5", "CCL7", "CCL8", "CCNA1",
  "CCND2", "CCNJL", "CCR10", "CCR2", "CCR3", "CCR4", "CCR5", "CCR7", "CCRL2", "CD14", "CD160", "CD163", "CD177",
  "CD180", "CD19", "CD1A", "CD1B", "CD1C", "CD1D", "CD1E", "CD2", "CD200", "CD200R1L", "CD207", "CD209", "CD24",
  "CD244", "CD247", "CD248", "CD27", "CD28", "CD300A", "CD300C", "CD300E", "CD33", "CD34", "CD36", "CD37", "CD38",
  "CD3D", "CD3E", "CD3G", "CD4", "CD40", "CD40LG", "CD44", "CD5", "CD52", "CD6", "CD68", "CD69", "CD7", "CD70", "CD72",
  "CD79A", "CD79B", "CD80", "CD84", "CD86", "CD8A", "CD8B", "CD93", "CD96", "CDA", "CDC20", "CDH11", "CDH12", "CDH2",
  "CDH23", "CDH5", "CDHR1", "CDK14", "CDKN1C", "CDKN2B", "CEACAM1", "CEACAM3", "CEACAM8", "CEBPA", "CES1", "CFH", "CFP",
  "CHAD", "CHI3L1", "CHI3L2", "CHIT1", "CHL1", "CHST15", "CKB", "CLC", "CLCN4", "CLEC10A", "CLEC14A", "CLEC17A", "CLEC4A",
  "CLEC4C", "CLEC4D", "CLEC4E", "CLEC5A", "CLEC6A", "CLEC7A", "CLEC9A", "CLECL1", "CLIC2", "CMA1", "CMKLR1", "CMTM2",
  "CNR1", "CNR2", "CNTN4", "CNTNAP2", "COBL", "COCH", "COL12A1", "COL13A1", "COL15A1", "COL19A1", "COL1A1", "COL1A2",
  "COL24A1", "COL4A1", "COL4A3", "COL4A4", "COL5A1", "COL5A2", "COL5A3", "COL6A2", "COL6A3", "COL8A1", "COL8A2",
  "COLEC12", "CORO2B", "CPA3", "CPM", "CPNE5", "CPNE7", "CPXM1", "CR2", "CRABP2", "CREB5", "CRHBP", "CRIP3", "CRISP3",
  "CRISPLD1", "CRLF2", "CRTAM", "CRYBB1", "CRYM", "CSF1", "CSF1R", "CSF2", "CSF3R", "CSNK1A1L", "CST7", "CTGF", "CTHRC1",
  "CTLA4", "CTNNA2", "CTSG", "CTSK", "CTSW", "CUX2", "CX3CR1", "CXCL1", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL14",
  "CXCL3", "CXCL5", "CXCL6", "CXCL9", "CXCR1", "CXCR2", "CXCR3", "CXCR4", "CXCR5", "CXCR6", "CXorf57", "CYBB", "CYP11A1",
  "CYP1B1", "CYP24A1", "CYP27A1", "CYP2E1", "CYP2J2", "CYP2S1", "CYP46A1", "CYP4F3", "CYTL1", "CYYR1", "DAB2", "DACH1",
  "DAPK2", "DCHS1", "DCN", "DEFA4", "DEFB1", "DGKK", "DKK1", "DKK2", "DKK3", "DLL1", "DLL4", "DNAH6", "DNAH8", "DNAJC5B",
  "DNASE1L3", "DNTT", "DOK7", "DPEP1", "DPEP2", "DPP4", "DPT", "DPYSL3", "DSC1", "DSEL", "DSG3", "DTHD1", "DTNA", "DUSP2",
  "DYTN", "EBF1", "EBI3", "ECSCR", "EDAR", "EDN1", "EFNA5", "EGF", "EGFL7", "EGR1", "EGR2", "EHF", "ELANE", "ELOVL4",
  "EMCN", "EMID1", "ENG", "ENHO", "ENPP3", "EOMES", "EPHB3", "EPHB6", "EPHX4", "EPX", "EREG", "ERG", "ETS1", "ETV1",
  "ETV4", "F13A1", "FAM111B", "FAM124B", "FAM153A", "FAM153C", "FAM177B", "FAM198B", "FAP", "FASLG", "FAT4", "FBLN5",
  "FBLN7", "FBN2", "FBP1", "FCAR", "FCER1A", "FCER2", "FCGBP", "FCGR2B", "FCGR3B", "FCN1", "FCRL1", "FCRL2", "FCRL3",
  "FCRL5", "FCRL6", "FCRLA", "FES", "FEZ1", "FFAR2", "FGF2", "FIBCD1", "FIBIN", "FKBP10", "FLI1", "FLT1", "FLT3", "FLT4",
  "FLVCR2", "FN1", "FOLR2", "FOSB", "FOSL1", "FOXP3", "FPR1", "FPR2", "FPR3", "FREM1", "FSD1", "FSTL3", "FXYD6", "FZD2",
  "GADD45G", "GAL3ST4", "GALR1", "GBP1", "GFI1", "GFPT2", "GGT5", "GIMAP5", "GIMAP6", "GIMAP7", "GIPR", "GJA4", "GLB1L2",
  "GLDC", "GLIS3", "GNA14", "GNAL", "GNG4", "GNG7", "GNLY", "GPA33", "GPBAR1", "GPC4", "GPNMB", "GPR1", "GPR171", "GPR18",
  "GPR183", "GPR19", "GPR20", "GPR25", "GPR27", "GPR65", "GPR68", "GPRC5B", "GPRC5D", "GRAP2", "GRIP1", "GUCY2D", "GYPE",
  "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "HAL", "HAPLN3", "HAVCR1", "HCK", "HDC", "HESX1", "HEY1", "HGF", "HHEX", "HHIP",
  "HIC1", "HIST1H2AE", "HIST1H2BG", "HK3", "HKDC1", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA",
  "HMGA2", "HOXA1", "HOXA3", "HOXA6", "HOXA7", "HPCAL4", "HPGD", "HPGDS", "HPSE", "HRASLS2", "HRH1", "HRH4", "HRK",
  "HSD11B1", "HSPA6", "HTR1B", "HTR1F", "HTR2B", "HTR7", "ICAM1", "ICAM2", "ICAM4", "ICOS", "ID1", "IDO1", "IDO2", "IER3",
  "IFI27", "IFI44L", "IFI6", "IFIT1", "IFIT2", "IFIT3", "IFNA10", "IFNG", "IGF1", "IGFBP3", "IGFBP4", "IGLL1", "IGSF6",
  "IL12B", "IL12RB2", "IL13", "IL15", "IL16", "IL17A", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL1R2", "IL1RL1",
  "IL1RN", "IL21", "IL21R", "IL22RA2", "IL23R", "IL26", "IL2RA", "IL2RB", "IL3", "IL31RA", "IL32", "IL34", "IL3RA", "IL4",
  "IL4I1", "IL5", "IL5RA", "IL6", "IL6ST", "IL7", "IL7R", "IL9", "INHBA", "INPP4B", "IQCA1", "IRF4", "IRF8", "ISM1",
  "ITGA1", "ITGA2", "ITGA4", "ITGA5", "ITGA7", "ITGAL", "ITGB2", "ITGB3", "ITGB4", "ITGB7", "ITGB8", "ITK", "ITLN1",
  "JAKMIP1", "KCNA3", "KCNA5", "KCNG1", "KCNG2", "KCNH7", "KCNH8", "KCNIP3", "KCNJ15", "KCNK10", "KCNK13", "KCNMB1",
  "KCNN3", "KCTD12", "KDR", "KIAA0087", "KIAA1211", "KIF19", "KIF5C", "KIR2DL1", "KIR2DL3", "KIR2DL4", "KIR3DL1", "KIR3DL2",
  "KIR3DX1", "KIRREL3", "KIT", "KLF1", "KLF4", "KLHDC7A", "KLHL33", "KLHL4", "KLK1", "KLRB1", "KLRC1", "KLRC2", "KLRC3",
  "KLRC4", "KLRD1", "KLRF1", "KLRG1", "KLRK1", "KNDC1", "KRT1", "KRT14", "KRT23", "KRT5", "KRT7", "KRT72", "KRT73", "KRT81",
  "KRT86", "KY", "KYNU", "L1TD1", "LAG3", "LAIR1", "LAIR2", "LAMA2", "LAMC2", "LAMP3", "LAT", "LAYN", "LCK", "LCN10",
  "LDB2", "LEF", "LEKR1", "LEP", "LGALS12", "LGALS2", "LGALS3", "LGALS9B", "LGALS9C", "LGR6", "LHCGR", "LIF", "LIFR",
  "LILRA1", "LILRA2", "LILRA4", "LILRA5", "LILRB1", "LILRB2", "LINGO2", "LIX1", "LONRF3", "LOXL3", "LPAR1", "LRG1", "LRMP",
  "LRP1B", "LRP4", "LRRC32", "LRRC36", "LRRC4", "LRRC43", "LRRN3", "LST1", "LTA", "LTB", "LTC4S", "LTF", "LTK", "LUM", "LY6E",
  "LY86", "LY9", "LYNX1", "LYPD2", "LYZ", "MACROD2", "MAF", "MAFB", "MAGEA11", "MAK", "MAL", "MAMDC2", "MAN1A1", "MAN1C1",
  "MAP1A", "MAP4K1", "MAP9", "MAPT", "MARCO", "MAST1", "MBL2", "MDGA2", "MDS2", "ME1", "MECOM", "MEFV", "MELK", "MEOX1",
  "MEP1A", "MGAM", "MGST1", "MIAT", "MIXL1", "MME", "MMP1", "MMP12", "MMP17", "MMP2", "MMP25", "MMP28", "MMP3", "MMP7",
  "MMP9", "MMRN2", "MNDA", "MPL", "MPP2", "MRC1", "MS4A1", "MS4A14", "MS4A2", "MS4A3", "MS4A4A", "MS4A6A", "MS4A7", "MSC",
  "MSLN", "MSR1", "MUC1", "MX1", "MXRA8", "MYBPC2", "MYCT1", "MYO7A", "MYOZ3", "NAALADL1", "NAT8L", "NCAM1", "NCF2", "NCR1",
  "NCR3", "NECAB2", "NEFL", "NELL2", "NETO1", "NFATC1", "NFE2", "NGF", "NID1", "NIPSNAP3B", "NKAIN2", "NKD1", "NKG7", "NLRP3",
  "NMBR", "NMUR1", "NOD2", "NOG", "NOS3", "NOTCH4", "NOV", "NOX3", "NPAS1", "NPL", "NPM2", "NPR1", "NPTX2", "NR4A1", "NR4A3",
  "NRCAM", "NRG1", "NRP1", "NT5E", "NTN3", "NTN4", "NTNG1", "NTRK1", "NYNRIN", "OAS2", "OASL", "OLFML2B", "OR2A4", "OR6K3",
  "OR6N1", "OR6N2", "OSCAR", "OSM", "OTOA", "OTX1", "OXCT2", "P2RX1", "P2RX6", "P2RY10", "P2RY13", "P2RY14", "PADI4", "PALMD",
  "PAPSS2", "PAQR5", "PATL2", "PAX5", "PAX7", "PCDH12", "PCDH9", "PCDHA5", "PCSK5", "PCSK6", "PCYT1B", "PDCD1", "PDCD1LG2",
  "PDE3A", "PDE6C", "PDE9A", "PDGFB", "PDGFD", "PDGFRA", "PDGFRB", "PDK4", "PDZK1IP1", "PDZRN3", "PECAM1", "PGLYRP1", "PHEX",
  "PI16", "PI3", "PID1", "PKD2L2", "PKIB", "PLA1A", "PLA2G4A", "PLA2G7", "PLAG1", "PLCB1", "PLCH2", "PLD4", "PLEKHG1", "PLEKHG7",
  "PLEKHH2", "PLLP", "PLVAP", "PLXDC1", "PLXNA4", "PMCH", "PNOC", "PODN", "POSTN", "POU2AF1", "POU2F2", "PPARG", "PPBP", "PPFIA4",
  "PPP2R2B", "PRF1", "PRG2", "PRKAR2B", "PRLR", "PROC", "PROK2", "PROM1", "PRR5L", "PRSS23", "PRSS35", "PSG2", "PTCRA", "PTGDR",
  "PTGER2", "PTGER3", "PTGER4", "PTGIR", "PTPRCAP", "PTPRM", "PTX4", "PZP", "QPCT", "RAB27B", "RAMP3", "RASA3", "RASD1",
  "RASGRP1", "RASGRP2", "RASGRP3", "RBM11", "RBM24", "RCAN2", "REG4", "REN", "RENBP", "RFESD", "RGS1", "RGS13", "RGS17",
  "RHOBTB3", "RIPK3", "RNASE2", "RNASE6", "RNF157", "RNF165", "ROBO1", "ROBO4", "ROR2", "RORC", "RPH3A", "RPL10L", "RRAD",
  "RRM2", "RSAD2", "RTN1", "RUNX3", "RXFP1", "RXFP2", "RYR1", "S100A12", "S100A8", "S100A9", "S100B", "S1PR3", "S1PR5", "SAMD3",
  "SAMD9", "SAMSN1", "SARDH", "SCARA5", "SCARF1", "SCN9A", "SCT", "SDS", "SELE", "SELL", "SERPINB2", "SERPINE1", "SERPINF1",
  "SERPING1", "SETBP1", "SEZ6L", "SFRP2", "SFRP5", "SFTPD", "SGCD", "SGK1", "SH2D1A", "SH2D1B", "SHD", "SIGLEC1", "SIGLEC15",
  "SIGLEC5", "SIRPG", "SIT1", "SKAP1", "SLAMF1", "SLAMF6", "SLAMF8", "SLC12A1", "SLC12A3", "SLC15A3", "SLC16A10", "SLC17A3",
  "SLC18A2", "SLC1A3", "SLC1A7", "SLC24A4", "SLC2A6", "SLC2A9", "SLC35F3", "SLC37A2", "SLC38A11", "SLC46A2", "SLC4A10", "SLC7A10",
  "SLC7A11", "SLC7A3", "SLC7A5", "SLC8A3", "SLCO2B1", "SLCO4C1", "SLCO5A1", "SLFN13", "SLPI", "SNAI1", "SNCA", "SOAT2", "SOCS1",
  "SORBS2", "SORCS3", "SOST", "SOX5", "SP140", "SPIB", "SPINK5", "SPOCK2", "SPON2", "SPP1", "SPRY1", "SRPX", "SRPX2", "SSPN",
  "SSTR3", "SSX1", "ST3GAL6", "ST6GAL1", "ST6GAL2", "ST6GALNAC1", "ST8SIA1", "STAB1", "STAP1", "STEAP4", "STXBP6", "STYK1",
  "SULT1C4", "SYBU", "SYN3", "SYNM", "SYNPO", "SYT17", "TACSTD2", "TAS1R3", "TBX21", "TBXAS1", "TCL1A", "TCL1B", "TCL6", "TCN1",
  "TCTEX1D1", "TDRD6", "TEK", "TEX101", "TFCP2L1", "TFPI", "TGFBR2", "TGM3", "TGM5", "THBD", "THBS1", "THBS2", "THEM5", "THEMIS",
  "THSD7A", "THY1", "TIAM1", "TIE1", "TIGIT", "TIMD4", "TIMP3", "TLR1", "TLR10", "TLR2", "TLR3", "TLR4", "TLR7", "TLR8", "TM4SF1",
  "TMEM108", "TMEM150B", "TMEM156", "TMEM171", "TMEM176A", "TMEM176B", "TMIE", "TMIGD2", "TMTC1", "TNC", "TNFAIP6", "TNFRSF11A",
  "TNFRSF11B", "TNFRSF13B", "TNFRSF17", "TNFRSF4", "TNFRSF8", "TNFSF14", "TNFSF15", "TNFSF4", "TNFSF9", "TNIP3", "TNXB", "TPM2",
  "TPPP3", "TPRG1", "TPSAB1", "TPSB2", "TPSD1", "TRAT1", "TREM1", "TREM2", "TREML2", "TRIM64", "TRPM6", "TSHR", "TSHZ2", "TSPAN18",
  "TSPAN7", "TTC16", "TWIST2", "TXK", "TYR", "UBASH3A", "UBD", "UGT1A8", "UGT2B11", "UGT8", "UNC13C", "UNC45B", "UPK3A", "USHBP1",
  "UTY", "VCAM1", "VCAN", "VEGFC", "VILL", "VIT", "VMO1", "VNN1", "VNN2", "VNN3", "VPREB3", "VSIG4", "WDR63", "WIF1", "WNT16",
  "WNT2", "WNT5A", "WNT5B", "WNT7A", "XCL1", "XCL2", "XCR1", "XKRX", "ZAP70", "ZBED2", "ZBP1", "ZBTB16", "ZBTB32", "ZC3H12D",
  "ZEB2", "ZMAT4", "ZMYND15", "ZNF135", "ZNF365", "ZNF366", "ZNF385D", "ZNF391", "ZNF442", "ZNF556", "ZNF683", "ZNF80",
  "ZNF831", "ZNF860"
)


cellcols <- c("red", "firebrick", "#FF66FF", "#CC0000", "#FF0000", "#FF6633", "#FF9900", "#FFFF00", "darkblue", "#000099", "#0000FF", "#3399CC", "#00FFFF", "#00FFFF", "#006600", "#33CC00", "#66CC66", "#33CC00", "#00FF00", "#9966CC", "#FFFF00", "#999999", "#996633", "#333333")
names(cellcols) <- c("CD4.T.cells", "CD8.T.cells", "Treg", "T.CD4.naive", "T.CD4.memory", "T.CD8.naive", "T.CD8.memory", "NK", "B", "B.naive", "B.memory", "plasma", "pDC", "pDCs", "macrophages", "monocytes", "monocytes.C", "monocytes.NC.I", "mDCs", "neutrophils", "mast", "fibroblasts", "endothelial.cells", "tumor")

