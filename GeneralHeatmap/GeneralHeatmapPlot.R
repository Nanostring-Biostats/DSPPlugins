
### function to draw general heatmap plot

# dependent libraries
library(pheatmap)
library(RColorBrewer)

# main function called by DSP-DA
main <- function(dataset, segmentAnnotations, targetAnnotations, outputFolder) {
  
  #### Function Inputs ------------------------------------------------
  
  # users can modify following arguments - currently set to defaults
  
  heatmap_colors <- NULL #e.g. c("navy", "white", "firebrick3")
  scale_cutoff <- 3
  annotations_to_show <- NULL #e.g. c("TissueType", "Response")
  annotation_colors <- NULL #e.g. list(TissueType = c(Normal = "green", Disease = "red"), Response = c(NR = "blue", R = "orange"))
  sort_by <- NULL
  sort_order <- NULL
  clustering_distance <- "euclidean" #options: "euclidean", "correlation" (Pearson's), "maximum", "manhattan", "canberra", "binary", or "minkowski"
  scale <- "row" #options: "none", "row", "column"
  fontsize <- 10
  
  # set output file type for plot
  file_type <- "png" # other options:"png", "tiff", "jpeg", "bmp", svg", "pdf"
  
  # set aspect ratio of output file
  height <- 500 # recommend 10 for file_type "svg" or "pdf", and 500 for "png", "tiff", "jpeg", or "bmp"
  width <- 500 # recommend 10 for file_type "svg" or "pdf", and 500 for "png", "tiff", "jpeg", or "bmp"
  
  
  
  #### Do not modify below! ------------------------------------------------
  
  # set row names to target names
  targetCountMatrix <- dataset
  rownames(targetCountMatrix) <- targetAnnotations[match(rownames(targetCountMatrix), targetAnnotations[ , "TargetGUID"]), "TargetName"]
  
  # call plotting function
  match.fun(file_type)(file.path(outputFolder, paste0("GeneralHeatmap.", file_type), fsep = .Platform$file.sep), width = width, height = height)
  draw_general_heatmap(data = targetCountMatrix,
                       heatmap_colors = heatmap_colors,
                       scale_cutoff = scale_cutoff,
                       annotations = segmentAnnotations,
                       annotations_to_show = annotations_to_show,
                       annotation_colors = annotation_colors,
                       sort_by = sort_by,
                       sort_order = sort_order,
                       clustering_distance = clustering_distance,
                       scale = scale,
                       fontsize = fontsize)
  dev.off()
}


#' Draw general heatmap
#' 
#' @param data data frame of data to plot.  targets in rows, samples in columns
#' @param heatmap_colors list of low, med, high colors used for heatmap colors. if NULL, uses pheatmap default
#' @param scale_cutoff value to limit standard deviations of z-score in either direction
#' @param annotations data frame of sample annotations
#' @param annotations_to_show list of annotations to show as color bars on heatmap. must include sort_by variable if specified
#' @param annotation_colors named list of lists specifying colors for each level of annotations_to_show
#' @param sort_by annotation to use for supervised clustering
#' @param sort_order list of levels of annotation defined by sort_by in desired order of appearance. optionally defined if sort_by is defined.  If NULL, then default to alphabetical order
#' @param clustering_distance distance measure used in clustering rows and columns
#' @param scale passed to pheatmap.  character indicating if the values should be centered and scaled in either the row direction, column direction, or none
#' @param fontsize passed to pheatmap. base fontsize for the plot
#' 
#' @return Heatmap with annotations and legends
#' 
#' @examples 
#' @import pheatmap
#' @import RColorBrewer
#' draw_general_heatmap(data = targetCountMatrix,
#' heatmap_colors = c("navy", "white", "firebrick3"),
#' scale_cutoffs = 3,
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
                                 scale_cutoff = 3,
                                 annotations = segmentAnnotations,
                                 annotations_to_show = NULL,
                                 annotation_colors = NULL,
                                 sort_by = NULL,
                                 sort_order = NULL,
                                 clustering_distance = "euclidean", 
                                 scale = "row",
                                 fontsize = 10, ...) {
  
  
  # subset annotations data
  if (!is.null(annotations_to_show)) {
    annotations <- subset(annotations, select = c("segmentID", annotations_to_show))
    annot_rownames <- annotations[,1]
    annotations <- annotations[,-1, drop = FALSE]
    rownames(annotations) <- annot_rownames
  } else {
    annotations = NULL
  }
  
  # sort annotations and set up white spaces between levels
  if (!is.null(sort_by)) {
    if (!is.null(sort_order)) {
      annotations <- annotations[order(factor(annotations[,sort_by], levels = sort_order)), ,drop = FALSE]
    } else {
      annotations <- annotations[order(annotations[,sort_by]), ,drop = FALSE]
    }
    
    gaps_col <- match(unique(annotations[,sort_by]), annotations[,sort_by])
    gaps_col <- gaps_col[2:length(gaps_col)] - 1
  } else{
    gaps_col <- NULL
  }
  
  # log2 transform data
  data <- data.frame(log2(data))
  
  # cap log2 data at scale_cutoff
  scale_cutoff <- abs(scale_cutoff)
  data <- pmin(pmax(t(scale(t(data))), - scale_cutoff), scale_cutoff)
  colnames(data) <- rownames(annotations)
  
  # set heatmap color palette
  if (!is.null(heatmap_colors)) {
    heatmap_colors <- colorRampPalette(heatmap_colors)(100)
  } else {
    heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                                       "RdYlBu")))(100)
  }
  
  # plot pheatmap
  p <- pheatmap(mat = data,
                color = heatmap_colors, 
                cluster_cols = is.null(sort_by),
                clustering_distance_rows = clustering_distance,
                clustering_distance_cols = clustering_distance,
                legend_breaks = seq(-scale_cutoff, scale_cutoff, 1),
                legend_labels = seq(-scale_cutoff, scale_cutoff, 1),
                fontsize = 10,
                labels_row = rownames(data),
                cellheight = ifelse(nrow(data) < 63, 11, NA),
                cellwidth = ifelse(ncol(data) < 51, 11, NA),
                border_color = NA,
                show_colnames = FALSE,
                show_rownames = nrow(data) < 63,
                annotation_col = annotations,
                annotation_colors = annotation_colors,
                gaps_col = gaps_col,
                scale = scale, ...)
  
  return(p)
  
}
