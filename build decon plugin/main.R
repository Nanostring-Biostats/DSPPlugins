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





##########################################################
#### end of arguments. DO NOT CHANGE CODE BELOW HERE  ####
##########################################################


library(logNormReg)
library(pheatmap)
library(viridis)
library(scales)
library(openxlsx)

main <- function(dataset, segmentAnnotations, targetAnnotations, outputFolder) {

  #### preliminaries ----------------------
  dataset <- as.matrix(dataset)

  # access cell profile matrix file:
  X <- as.matrix(read.csv(cell_profile_filename, header = T, row.names = 1))


  # ARGUMENT (hidden): define cell types to be added together in the final result:
  # example syntax:
  # merges = list()
  # merges[["T"]] = c("CD8.T", "CD4.T")
  # merges[["myeloid"]] = c("macrophage", "monocyte", "DC")
  merges <- list()

  # parse merges:
  merges.full <- NULL
  if (length(merges) > 0) {
    # initialize with 1:1 mapping:
    merges.full <- list()
    for (name in colnames(X)) {
      merges.full[[name]] <- name
    }
    # add merges:
    for (name in names(merges0)) {
      # remove entries for cells specified by user and replace with their entries:
      merges.full[merges[[name]]] <- NULL
      merges.full[[name]] <- merges[[name]]
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
  # colnames(norm) <- segmentAnnotations$segmentDisplayName

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
    cellmerges = merges.full
  )

  # reverse decon:
  # rdecon <- reverseDecon(
  #  norm = norm,
  #  beta = res$beta,
  #  epsilon = NULL
  # )


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



  # reverse decon resids
  # write.csv(rdecon$resids, file = file.path(outputFolder, "reverse_decon_residuals.csv", fsep = .Platform$file.sep))

  ## reverse decon summary stats of gene dependency on cell mixing:
  # sumstats <- cbind(rdecon$cors, rdecon$resid.sd)
  # colnames(sumstats) <- c("cor w cell mixing", "residual SD from cell mixing")
  # write.csv(sumstats, file = file.path(outputFolder, "gene_dependence_on_cell_mixing.csv", fsep = .Platform$file.sep))

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
  cells.to.plot <- intersect(rownames(res$beta), union(colnames(X), names(merges.full)))

  # one pdf for all results:
  pdf(file = file.path(outputFolder, "spatialdecon_results.pdf", fsep = .Platform$file.sep), width = 12)

  #### heatmaps
  # abundances:
  thresh <- signif(quantile(res$beta, 0.97), 2)
  p1 <- pheatmap(pmin(res$beta[cells.to.plot, ], thresh),
    col = colorRampPalette(hmcols)(100),
    fontsize_col = 4,
    annotation_col = heatmapannot,
    annotation_colors = cols,
    legend_breaks = c(round(seq(0, thresh, length.out = 5))[-5], thresh),
    legend_labels = c(round(seq(0, thresh, length.out = 5))[-5], paste0("Abundance scores,\ntruncated above at ", thresh))
    # main = paste0("Abundance scores, truncated above at ", thresh)
  )
  # print(p1)

  # scaled abundances:
  mat <- sweep(res$beta, 1, apply(res$beta, 1, max), "/")
  mat <- replace(mat, is.na(mat), 0)
  p3 <- pheatmap(mat,
    col = colorRampPalette(hmcols)(100),
    fontsize_col = 4,
    annotation_col = heatmapannot,
    annotation_colors = cols,
    legend_breaks = c(round(seq(0, 1, length.out = 5), 2)[-5], 1),
    legend_labels = c(round(seq(0, 1, length.out = 5), 2)[-5], "Scaled abundance\n(ratio to max)")
    # main = paste0("Abundance scores, truncated above at ", thresh)
  )
  # print(p3)

  # proportions:
  props <- replace(res$prop_of_nontumor[cells.to.plot, ], is.na(res$prop_of_nontumor[cells.to.plot, ]), 0)
  p2 <- pheatmap(props,
    col = colorRampPalette(hmcols)(100),
    fontsize_col = 4,
    annotation_col = heatmapannot,
    annotation_colors = cols,
    legend_breaks = round(seq(0, max(props) * 0.99, length.out = 5), 2),
    legend_labels = c(round(seq(0, max(props), length.out = 5), 2)[-5], "Proportion of all\nfitted populations")
  )
  # print(p2)

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
  layout(mat = matrix(c(1, 2, 3, 3), 2), widths = c(10, 3, 10, 3), heights = c(1, 8, 10))
  par(mar = c(0, 8.2, 0, 0.2))
  plot(p1$tree_col, labels = F, main = "", ylab = "", yaxt = "n")
  par(mar = c(15, 8, 0, 0))

  # data to plot:
  mat <- res$beta[cells.to.plot, p1$tree_col$order]
  # infer scale of negative y-axis for annotation colorbars
  ymin <- -max(colSums(mat)) * 0.15

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




  ### proportion barplot:
  layout(mat = matrix(c(1, 2, 3, 3), 2), widths = c(10, 3, 10, 3), heights = c(1, 8, 10))
  par(mar = c(0, 8.2, 0, 0.2))
  plot(p1$tree_col, labels = F, main = "", ylab = "", yaxt = "n")
  par(mar = c(15, 8, 0, 0))

  # data to plot:
  mat <- res$prop_of_nontumor[cells.to.plot, p2$tree_col$order]
  # infer scale of negative y-axis for annotation colorbars
  ymin <- -max(colSums(mat)) * 0.15
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





  # spaceplots: draw separately for each cell:
  if (FALSE) {
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
  dev.off()
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
