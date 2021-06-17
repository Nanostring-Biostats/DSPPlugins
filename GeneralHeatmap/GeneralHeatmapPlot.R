# General Heatmap #

# Produces heatmap with clinical annotations
# Supports: DSP-NGS CTA
# Note: 
# Please do not use spaces, special characters, or numbers when adding factors
# in the DSPDA Annotation file

##############################
#        User Options        #
##############################

# users can modify following arguments - currently set to defaults


# define annotations to show in heatmap (optional):
annotations_to_show <- NULL #e.g. c("TissueType", "Response")

# define coloring of annotations (optional):
# 1: set the next line to TRUE:
custom_annotation_colors <- FALSE
# 2: then modify the example syntax below
if (custom_annotation_colors) {
  # example syntax
  annotation_colors <- list(
    TissueType = c(
      "Normal" = "green", 
      "Disease" = "red"), 
    Response = c(
      "NR" = "blue", 
      "R" = "orange"))
  } else {
    annotation_colors <- NULL
  }
  # for a list of R colors, see http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf


# define annotation to use for supervised clustering (optional):
sort_by <- NULL #e.g. "Response"
# define order levels should appear in (optional):
sort_order <- NULL #e.g. c("Unknown", "Responder", "Non-Responder")


#### Advanced User Inputs

# choose heatmap color scheme (optional):
heatmap_colors <- NULL #e.g. c("navy", "white", "firebrick3")

# set heatmap color and legend scale (absolute value):
scale_cutoff <- 3

# set whether sample and target names should be shown:
show_sample_names <- FALSE
show_target_names <- FALSE

# set clustering options for unsupervised clustering:
clustering_distance <- "euclidean" #options: "euclidean", "correlation" (Pearson's), "maximum", "manhattan", "canberra", "binary", or "minkowski"
scale <- "row" #options: "none", "row", "column"

fontsize <- 10
plot_title <- "log2 Change from Mean"

# set output file type for plot:
file_type <- "pdf" # other options:"svg", "png", "tiff"

# set aspect ratio of output file:
pdf_width <- 12 
pdf_height <- 7 



##########################################################
#### End of User Inputs. PLEASE DO NOT CHANGE CODE BELOW HERE  ####
##########################################################
# MIT License
# Copyright 2020 Nanostring Technologies, Inc.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# Contact us: 
#   NanoString Technologies, Inc.
#   530 Fairview Avenue N
#   Seattle, WA 98109
# Tel: (888) 358-6266
##############################

##############################
#        Execution Code      # 
##############################

# dependent libraries
library(pheatmap)
library(RColorBrewer)
library(dplyr)

# main function called by DSP-DA
main <- function(dataset, segmentAnnotations, targetAnnotations, outputFolder) {
  
  targetCountMatrix <- dataset
  # make unique segment identifiers instead of GUIDs
  segmentAnnotationsMod <- segmentAnnotations %>%
    mutate(segmentDisplayName = paste(ScanName, ROIName, SegmentName, sep=" | "))
  # update count matrix column names with new segment unique ID instead of GUID. 
  names(targetCountMatrix) <- segmentAnnotationsMod[match(names(targetCountMatrix), segmentAnnotationsMod[ , "segmentID"]), "segmentDisplayName"]  
  # update count matrix rownames with targetNames
  rownames(targetCountMatrix) <- targetAnnotations[match(rownames(targetCountMatrix), targetAnnotations[ , "TargetGUID"]), "TargetName"]
  
  # setup output file
  if (file_type == "pdf") {
    pdf(file = file.path(outputFolder, "GeneralHeatmap.pdf", fsep = .Platform$file.sep), width = pdf_width, height = pdf_height)
  } else if (file_type == "svg") {
    svg(file = file.path(outputFolder, "GeneralHeatmap.svg", fsep = .Platform$file.sep), width = pdf_width, height = pdf_height)
  } else if (file_type == "png") {
    png(file = file.path(outputFolder, "GeneralHeatmap.png", fsep = .Platform$file.sep), width = pdf_width, height = pdf_height, units = "in", res = 2000)
  } else if (file_type == "tiff") {
    tiff(file = file.path(outputFolder, "GeneralHeatmap.tiff", fsep = .Platform$file.sep), width = pdf_width, height = pdf_height, units = "in", res = 150)
  } else {
    stop("Error: file_type must be pdf, svg, png, or tiff")
  }
  
  # checks for user inputs
  if (!is.numeric(scale_cutoff)) {
    stop("Error: scale_cutoff must be numeric")
  }
  if (!is.numeric(fontsize)) {
    stop("Error: fontsize must be numeric")
  }
  if (!is.numeric(pdf_width)) {
    stop("Error: pdf_width must be numeric")
  }
  if (!is.numeric(pdf_height)) {
    stop("Error: pdf_height must be numeric")
  }
  
  # call plotting function
  draw_general_heatmap(data = targetCountMatrix,
                       heatmap_colors = heatmap_colors,
                       scale_cutoff = scale_cutoff,
                       annotations = segmentAnnotationsMod,
                       annotations_to_show = annotations_to_show,
                       annotation_colors = annotation_colors,
                       sort_by = sort_by,
                       sort_order = sort_order,
                       clustering_distance = clustering_distance,
                       scale = scale,
                       show_sample_names = show_sample_names,
                       show_target_names = show_target_names,
                       fontsize = fontsize,
                       plot_title = plot_title)
  dev.off()
}


#' Draw general heatmap
#' 
#' @param data data frame of data to plot.  targets in rows, samples in columns
#' @param heatmap_colors list of low, med, high colors used for heatmap colors. if NULL, uses pheatmap default
#' @param scale_cutoff value to set range of heatmap colors and legend
#' @param annotations data frame of sample annotations
#' @param annotations_to_show list of annotations to show as color bars on heatmap. must include sort_by variable if specified
#' @param annotation_colors named list of lists specifying colors for each level of annotations_to_show
#' @param sort_by annotation to use for supervised clustering
#' @param sort_order list of levels of annotation defined by sort_by in desired order of appearance. optionally defined if sort_by is defined.  If NULL, then default to alphabetical order
#' @param clustering_distance distance measure used in clustering rows and columns
#' @param scale passed to pheatmap.  character indicating if the values should be centered and scaled in either the row direction, column direction, or none
#' @param show_sample_names TRUE or FALSE
#' @param show_target_names TRUE or FALSE
#' @param fontsize passed to pheatmap. base fontsize for the plot
#' @param plot_title title of heatmap
#' 
#' @return Heatmap with annotations and legends
#' 
#' @examples 
#' @import pheatmap
#' @import RColorBrewer
#' draw_general_heatmap(data = targetCountMatrix,
#' heatmap_colors = c("navy", "white", "firebrick3"),
#' scale_cutoff = 3,
#' annotations = segmentAnnotations,
#' annotations_to_show = c("RIOName", "ScanName"),
#' annotation_colors = NULL,
#' sort_by = "ScanName",
#' sort_order = c("HST 18.1 F 2-9", "HST 18.1 G 2-9"),
#' clustering_distance = "euclidean",
#' scale = "row",
#' fontsize = 10)
#' 
#' @export draw_general_heatmap

draw_general_heatmap <- function(data = targetCountMatrix,
                                 heatmap_colors = NULL,
                                 scale_cutoff = 3L,
                                 annotations = segmentAnnotations,
                                 annotations_to_show = NULL,
                                 annotation_colors = NULL,
                                 sort_by = NULL,
                                 sort_order = NULL,
                                 clustering_distance = "euclidean", 
                                 scale = "row",
                                 show_sample_names = FALSE,
                                 show_target_names = FALSE,
                                 fontsize = 10L, 
                                 plot_title = "log2 Change from Mean", ...) {
  
  
  # subset annotations data
  if (!is.null(annotations_to_show)) {
    annotations <- subset(annotations, select = c("segmentDisplayName", annotations_to_show))
    annot_rownames <- annotations[,1L]
    annotations <- annotations[,-1L, drop = FALSE]
    rownames(annotations) <- annot_rownames
  } else {
    annotations <- NULL
  }
  
  # sort annotations and set up white spaces between levels
  if (!is.null(sort_by)) {
    if (!is.null(sort_order)) {
      annotations <- annotations[order(factor(annotations[,sort_by], levels = sort_order)), ,drop = FALSE]
    } else {
      annotations <- annotations[order(annotations[,sort_by]), ,drop = FALSE]
    }
    
    gaps_col <- match(unique(annotations[,sort_by]), annotations[,sort_by])
    gaps_col <- gaps_col[2L:length(gaps_col)] - 1L
  } else{
    gaps_col <- NULL
  }
  
  # shift zeros to one
  if (any(as.vector(data) == 0)) {
    for (i in seq_along(colnames(data))) {
      data[which(data[,i] == 0),i] <- 1
    }
  }
  
  # log2 transform data
  data <- data.frame(log2(data))
  colnames(data) <- rownames(annotations)
  
  # remove all zero rows
  data = data[rowSums(data[]) != 0,]
  
  # set heatmap color palette
  if (!is.null(heatmap_colors)) {
    heatmap_colors <- colorRampPalette(heatmap_colors)(100L)
  } else {
    heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 7L, name =
                                                       "RdYlBu")))(100L)
  }
  
  scale_cutoff <- abs(scale_cutoff)
  
  # plot pheatmap
  p <- pheatmap(mat = data,
                color = heatmap_colors, 
                cluster_cols = is.null(sort_by),
                clustering_distance_rows = clustering_distance,
                clustering_distance_cols = clustering_distance,
                legend_breaks = seq(-scale_cutoff, scale_cutoff, 1L),
                legend_labels = seq(-scale_cutoff, scale_cutoff, 1L),
                breaks = seq(-scale_cutoff, scale_cutoff, length.out = 100L),
                fontsize = fontsize,
                labels_row = rownames(data),
                cellheight = ifelse(nrow(data) < 63L, 11L, NA),
                cellwidth = ifelse(ncol(data) < 51L, 11L, NA),
                border_color = NA,
                show_colnames = show_sample_names,
                show_rownames = show_target_names,
                annotation_col = annotations,
                annotation_colors = annotation_colors,
                gaps_col = gaps_col,
                scale = scale, 
                main = plot_title, ...)
  
  return(p)
  
}
