# QC Heatmap #

# Produces heatmap with probe detection barplot
# Supports: DSP
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




### Advanced User Inputs

detection_thresh <- 2
proportion_detect_thresh <- .10


# choose heatmap color scheme and breaks (optional):
heatmap_color_breaks <- c(0, 2, 5, 10, 50)
heatmap_color_palette <- rev(viridis(5)) #c("white", "white", "cadetblue2", "cadetblue4", "darkblue")

column_detection_barplot <- FALSE
cluster_columns <- TRUE
row_detection_barplot <- TRUE

plot_title <- "Signal-To-Noise Ratio"
legend_title <- "SNR"

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
library(circlize)
library(ComplexHeatmap)
library(ggplot2)
library(purrr)
library(dplyr)
library(viridis)

# main function called by DSP-DA
main <- function(dataset, segmentAnnotations, targetAnnotations, targetCountMatrix, outputFolder) {
  

  # make unique sample identifiers instead of GUIDs
  targetCountMatrix <- dataset
  segmentAnnotationsMod <- segmentAnnotations %>%
    mutate(segmentDisplayName = paste(ScanName, ROIName, SegmentName, sep=" | "))
  names(targetCountMatrix) <- segmentAnnotationsMod[match(names(targetCountMatrix), segmentAnnotationsMod[ , "segmentID"]), "segmentDisplayName"]  
  rownames(targetCountMatrix) <- targetAnnotations[match(rownames(targetCountMatrix), targetAnnotations[ , "TargetGUID"]), "TargetName"]
  
  # calculate SNR matrix from background
  bg <- derive_GeoMx_background(norm = targetCountMatrix, probepool = targetAnnotations$ProbePool, negnames = targetAnnotations$TargetName[targetAnnotations$CodeClass == "Negative"])
  targetSNR <- targetCountMatrix/bg
  
  # setup output file
  if (file_type == "pdf") {
    pdf(file = file.path(outputFolder, "QCHeatmap.pdf", fsep = .Platform$file.sep), width = pdf_width, height = pdf_height)
  } else if (file_type == "svg") {
    svg(file = file.path(outputFolder, "QCHeatmap.svg", fsep = .Platform$file.sep), width = pdf_width, height = pdf_height)
  } else if (file_type == "png") {
    png(file = file.path(outputFolder, "QCHeatmap.png", fsep = .Platform$file.sep), width = pdf_width, height = pdf_height, units = "in", res = 2000)
  } else if (file_type == "tiff") {
    tiff(file = file.path(outputFolder, "QCHeatmap.tiff", fsep = .Platform$file.sep), width = pdf_width, height = pdf_height, units = "in", res = 150)
  } else {
    stop("Error: file_type must be pdf, svg, png, or tiff")
  }
  
  # checks for user inputs
  

  
  
  
  # call plotting function
  draw_detection_heatmap(SNR_data = t(targetSNR), 
                          detection_thresh = detection_thresh,
                          annotations = segmentAnnotationsMod,
                          annotations_to_show = annotations_to_show,
                          annotation_colors = annotation_colors,
                          heatmap_color_breaks = heatmap_color_breaks, 
                          heatmap_color_palette = heatmap_color_palette,
                          heatmap_height = heatmap_height,
                          column_detection_barplot = column_detection_barplot, 
                          proportion_detect_thresh = proportion_detect_thresh, 
                          cluster_columns = cluster_columns,
                          row_detection_barplot = row_detection_barplot, 
                          plot_title = plot_title,
                          legend_title = legend_title)
  dev.off()
}


#' Derive background at the scale of the normalized data for GeoMx data
#'
#' Estimates per-datapoint background levels from a GeoMx experiment.
#' In studies with two or more probe pools, different probes will have different
#' background levels. This function provides a convenient way to account for this phenomenon.
#'
#' @param norm Matrix of normalized data, genes in rows and segments in columns.
#'  Must include negprobes, and must have rownames.
#' @param probepool Vector of probe pool names for each gene, aligned to the rows of "norm".
#' @param negnames Names of all negProbes in the dataset. Must be at least one neg.name within each probe pool.
#' @return A matrix of expected background values, in the same scale and dimensions as the "norm" argument.
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
  
  
  if (all(is.na(probepool))) {
    probepool <- rep(1, length(probepool))
    warning("all probepool info missing.  assuming one probe pool")
  } else if (any(is.na(probepool))) {
    stop("probepool has missing values")
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
    bg[probepool == pool, ] <- sweep(bg[probepool == pool, ], 2, tempnegfactor, "+")
  }
  return(bg)
}



#' Generate a heatmap of Signal to Noise Ratios (SNR) stratified by clinical variables
#' 
#' @param SNR_data matrix of Signal to Noise Ration (SNR) with genes as columns
#' @param detection_thresh minimum threshold for detection
#' @param annotations matrix of annotations with clinical variables as columns
#' @param annotations_to_show list of annotation column names to be shown as color bars on left side of heatmap
#' @param annotation_colors named list of lists specifying colors for each level of annotations_to_show
#' @param heatmap_color_breaks optional list of break points for heatmap colors
#' @param heatmap_color_palette optional list of colors for heatmap.  Must be same length as heatmap_color_breaks
#' @param column_detection_barplot TRUE or FALSE to determine whether to show column barplot of % detected genes
#' @param proportion_detect_thresh minimum threshold for proportion of detection, entered as a decimal point
#' @param cluster_columns TRUE or FALSE. If FALSE and column_detection_barplots = TRUE, then columns ordered by % detected
#' @param row_detection_barplot TRUE or FALSE to determine whether to show row barplot of % detected samples
#' @param plot_title title for plot displayed on the top.  Default is "Signal-To-Noise Ratio".  Set as "" for no title
#' @param legend_title title for heatmap legend.  Default is "SNR".  Set as "" for no title
#' 
#' @return Stratified heatmap with annotations and legends
#' 
#' @examples 
#' @import plotrix
#' @import circlize
#' @import ComplexHeatmap
#' draw_detection_heatmap(SNR_data = data_mat, 
#' detection_thresh = 2,
#' annotations = annot_mat,
#' annotations_to_show = c("Treatment_Type", "Response", "Tissue_Type"),
#' annotation_colors = NULL,
#' heatmap_color_breaks = c(0, 10, 20, 50, 100), 
#' heatmap_color_palette = c("red4", "red4", "tomato", "white", "blue"),
#' column_detection_barplot = FALSE,
#' proportion_detect_thresh = .10,
#' cluster_columns = TRUE)
#' 
# #' @export


draw_detection_heatmap = function(SNR_data, 
                                   detection_thresh = 2,
                                   annotations = NULL,
                                   annotations_to_show = NULL,
                                   annotation_colors = NULL,
                                   heatmap_color_breaks = c(0, 2, 5, 10, 50), 
                                   heatmap_color_palette = rev(viridis(5)),
                                   column_detection_barplot = FALSE, 
                                   proportion_detect_thresh = .10, 
                                   cluster_columns = TRUE,
                                   row_detection_barplot = FALSE, 
                                   plot_title = "Signal-To-Noise Ratio",
                                   legend_title = "SNR", ...){
  
  mat <- SNR_data
  
  # make color function for heatmap
  
  breaks <- heatmap_color_breaks
  colors <- heatmap_color_palette
  
  col_fun <- colorRamp2(breaks = breaks, colors = colors)
  
  ## row annotations ##
  
  set.seed(555)
  
  if (is.null(annotations) || is.null(annotations_to_show)) {
    row_ha_clinical = NULL
  } else {
    anno_df <- as.data.frame(annotations[,annotations_to_show])
    colnames(anno_df) <- annotations_to_show
    
    # add custom annotation colors if defined
    if (!is.null(annotation_colors)){
      row_ha_clinical <- rowAnnotation(df = anno_df, 
                                       col = annotation_colors)
    } else {
      row_ha_clinical <- rowAnnotation(df = anno_df)
    }
  }
  
  
  # default to no additional legend
  lgd = NULL
  
  ## column barplot annotations ##
  if (column_detection_barplot == TRUE) {
    # make vector of % detected (> detection_thresh) in each column
    col_detect_vec <- colSums(mat > detection_thresh)/nrow(mat)
    
    # make vector of colors for columns with less than proportion_dect_thresh
    col_color_vec <- col_detect_vec < proportion_detect_thresh
    col_color_vec[col_color_vec == "TRUE"] <- "red"
    col_color_vec[col_color_vec == "FALSE"] <- "grey"
    
    # column barplot annotation
    col_ha <- HeatmapAnnotation(detected = anno_barplot(col_detect_vec, 
                                                        bar_width = 1, 
                                                        gp = gpar(fill = col_color_vec)))
    
    # add annotation legend
    lgd <- Legend(labels = c(paste0(">", proportion_detect_thresh*100, "%"), 
                             paste0("<", proportion_detect_thresh*100, "%")), 
                  title = "Proportion Detected", 
                  legend_gp = gpar(fill = c("grey", "red")))
    
    if (cluster_columns == FALSE) {
      column_order <- names(sort(col_detect_vec, decreasing = TRUE))
    } else {
      column_order <- NULL
    }
    
  } else {
    col_ha <- NULL
    column_order <- NULL
  }
  
  ## row barplot annotation ##
  if (row_detection_barplot == TRUE) {
    # make vector of % deteched in each row
    row_detect_vec <- rowSums(mat > detection_thresh)/ncol(mat)
    
    # make vector of colors for rows with less than proportion_detect_thresh
    row_color_vec = row_detect_vec < proportion_detect_thresh
    row_color_vec[row_color_vec == "TRUE"] <- "red"
    row_color_vec[row_color_vec == "FALSE"] <- "grey"
    
    # row barplot annotation
    samples_to_point <- names(row_color_vec[row_color_vec == "red"])
    
    if (length(samples_to_point > 0)){
      row_ha <- rowAnnotation(detected = anno_barplot(row_detect_vec, 
                                                      bar_width = 1, 
                                                      gp = gpar(fill = row_color_vec)),
                              low_detection = anno_mark(at = match(samples_to_point, rownames(mat)), labels = samples_to_point))
    } else {
      row_ha <- rowAnnotation(detected = anno_barplot(row_detect_vec, 
                                                      bar_width = 1, 
                                                      gp = gpar(fill = row_color_vec)))
    }
    
    # add annotation legend
    lgd <- Legend(labels = c(paste0(">", proportion_detect_thresh*100, "%"), 
                             paste0("<", proportion_detect_thresh*100, "%")), 
                  title = "Proportion Detected", 
                  legend_gp = gpar(fill = c("grey", "red")))
    
  } else {
    row_ha <- NULL
  }
  
  
  # create heatmap object
  ht <- Heatmap(mat, name = "SNR", 
                col = col_fun, 
                heatmap_legend_param = list(title = legend_title, at = breaks, legend_height = unit(4, "cm")),
                left_annotation = row_ha_clinical, 
                right_annotation = row_ha,
                top_annotation = col_ha,
                show_row_names = FALSE,
                show_column_names = FALSE,
                show_column_dend = FALSE,
                cluster_columns = cluster_columns,
                column_order = column_order,
                column_title = plot_title, ...)
  
  draw(ht, annotation_legend_list = lgd)
  
}
