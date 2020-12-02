
### function to draw general heatmap plot

# dependent libraries
library(pheatmap)
library(RColorBrewer)

# main function called by DSP-DA
main <- function(dataset, segmentAnnotations, targetAnnotations, outputFolder) {
  
  #### Function Inputs ------------------------------------------------
  
  # users can modify following arguments - currently set to defaults
  
  # supervised/unsupervised ROIs option
  sort_by = NULL# annotation to use for clustering (supervised), if NULL, then unsupervised clustering
  sort_order = NULL # list of levels of annotation defined by sort_by, NULL for unsupervised clustering.  should be defined if sort_by is defined.  If NULL, then default to alphabetical order
  
  # annotations input (up to 3)
  annotations_to_show = NULL #e.g. annotations_to_show = c("TissueType", "Response")
  
  # annotations colors input
  annotations_colors = NULL # named list of lists. e.g. annotation_colors = list(TissueType = c(Normal = "blue", Disease = "red"), Response = c(R = "green", NR = "orange"))
  
  # heatmap colors option (high, med, low)
  heatmap_colors = NULL # list of low, med, high colors. if NULL, uses pheatmap default. e.g. c("navy", "white", "firebrick3")
  
  scale = "row" # passed to pheatmap. options are "none", "row", "column"
  scale_cutoff = 3
  
  # clustering distance & linkage option
  clustering_distance = "euclidean" # other options: "correlation" (Pearson's), "maximum", "manhattan", "canberra", "binary", or "minkowski"
  
  
  #### Do not modify below! ------------------------------------------------
  
  # set row names to target names
  targetCountMatrix <- dataset
  rownames(targetCountMatrix) <- targetAnnotations[match(rownames(targetCountMatrix), targetAnnotations[ , "TargetGUID"]), "TargetName"]
  
  # set up height and width for output PDF
  width <- ncol(targetCountMatrix)*1.6
  if (ncol(targetCountMatrix) > 50) {
    width <- 1000
  } else if (ncol(targetCountMatrix < 15)) {
    width <- width*2.5
  }
  
  height <- max(nrow(targetCountMatrix) * 1.3, 10)
  if (ncol(targetCountMatrix) > 63) {
    height <- 800
  }
  
  
  pdf(file = file.path(outputFolder, "GeneralHeatmap.pdf", fsep = .Platform$file.sep), width = width, height = height)
  p=draw_general_heatmap(data = targetCountMatrix[1:50,],
                       heatmap_colors = NULL,
                       scale_cutoff = 3,
                       annotations = segmentAnnotations,
                       annotations_to_show = c("ROIName"),
                       annotation_colors = NULL,
                       sort_by = NULL,
                       sort_order = NULL,
                       clustering_distance = "euclidean",
                       scale = "row")
  dev.off()
}





  
  draw_general_heatmap <- function(data = targetCountMatrix,
                                   heatmap_colors = NULL,
                                   scale_cutoff = 3,
                                   annotations = segmentAnnotations,
                                   annotations_to_show = NULL,
                                   annotation_colors = NULL,
                                   sort_by = NULL,
                                   sort_order = NULL,
                                   clustering_distance = "euclidean", 
                                   scale = "row") {
    
    
    # subset annotations data
    if (!is.null(annotations_to_show)) {
      annotations <- subset(annotations, select = c("segmentID", annotations_to_show))
      annot_rownames <- annotations[,1]
      annotations <- annotations[,-1, drop = FALSE]
      rownames(annotations) <- annot_rownames
    }
    
    # log2 transform data
    data <- data.frame(log2(data))
    
    # cap log2 data at scale_cutoff
    scale_cutoff <- abs(scale_cutoff)
    data <- pmin(pmax(t(scale(t(data))), - scale_cutoff), scale_cutoff)
    colnames(data) <- rownames(annotations)
    
    # set heatmap color palette
    if (!is.null(heatmap_colors)) {
      heatmap_colors = colorRampPalette(heatmap_colors)(100)
    } else {
      heatmap_colors = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                         "RdYlBu")))(100)
    }
    
    # plot pheatmap
    ph <- pheatmap(mat = data,
                   color = heatmap_colors, 
                   cluster_cols = is.null(sort_by),
                   clustering_distance_rows = clustering_distance,
                   clustering_distance_cols = clustering_distance,
                   fontsize = 10,
                   labels_row = rownames(targetCountMatrix),
                   cellheight = ifelse(nrow(data) < 63, 11, NA),
                   cellwidth = ifelse(ncol(data) < 51, 11, NA),
                   border_color = NA,
                   show_colnames = FALSE,
                   show_rownames = nrow(data) < 63,
                   annotation_col = annotations,
                   annotation_colors = annotation_colors,
                   scale = scale
                   #main = name
                   )
    
    return(ph)
    
  }
  