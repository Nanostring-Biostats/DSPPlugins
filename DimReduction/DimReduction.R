# Type of Plot:
plot_type <- "tSNE"
# Options: tSNE, UMAP, PCA

# Plot Parameters
color_by = "CD68" # tag, factor, or target
shape_by = "SegmentName" # tag, or factor
plot_font = list(family = "sans", size = 15)
# shape & color can be set to NULL
# font families include sans, serif, mono and
# may include specifically named fonts, but 
# these may not always render properly.

# Plot color options
plot_colors = list("orange2", "gray", "darkblue")
color_levels = c("High", "Mid", "Low")
# color_levels must match values in the color_by 
# column in the annotations file when using a tag
# or a factor. For Targets use "High", "Low" and "Mid".
#   "Mid" is optional

#plot_color_theme = "RdBu"
#reverse_theme = TRUE # reverse palette color
# Set to plot_color_theme = NULL if you'd 
# prefer to set your colors manually below
#plot_color_theme = "RdBu" # color palette

##################
# Execution Code #
##################
# Libraries:
library(Rtsne)
library(umap)
library(stats)
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)

# main function
main <- function(dataset, segmentAnnotations, targetAnnotations, outputFolder) {
  # Step 1: format dataset
  # Make a table to check the passed dataframes from DSP
  targetCountMatrix <- dataset
  
  # Make unique segment identifiers instead of GUIDs
  segmentAnnotations <- segmentAnnotations %>%
    mutate(segmentDisplayName = paste(ScanName, ROIName, SegmentName, sep=" | "))
  
  # Update count matrix column names with new segment unique ID instead of GUID. 
  names(targetCountMatrix) <- segmentAnnotations[match(names(targetCountMatrix),
                                                       segmentAnnotations[ , "segmentID"]),
                                                 "segmentDisplayName"]  
  
  # Update count matrix rownames with targetNames
  rownames(targetCountMatrix) <- targetAnnotations[match(rownames(targetCountMatrix),
                                                         targetAnnotations[ , "TargetGUID"]),
                                                   "TargetName"]
  # gather color parameters
  names(plot_colors) <- color_levels
  params <- list("plot_type" = plot_type,
                 "color_by" = color_by,
                 "shape_by" = shape_by,
                 "plot_font" = plot_font,
                 "plot_colors" = plot_colors)
  
  # Step 2:
  # Calculate PCA, tSNE, or UMAP dimensions
  DR_ann <- calc_DR(plot_type = plot_type,
                    targetCountMatrix = targetCountMatrix,
                    segmentAnnotations = segmentAnnotations)
  segmentAnnotations <- DR_ann$annot
  if(plot_type == "PCA") {
    var_est <- summary(DR_ann$data)$importance[3, ]
  } else {
    var_est <- NULL
  }

  # Step 3: Graph data
  # if color or shape is a gene symbol add it to the annotations for plotting
  if(is.null(color_by)) {
    colType <- 'Null'
  } else if(color_by %in% rownames(targetCountMatrix)) {
    if(!all(color_levels %in% c("High", "Mid", "Low"))) {
      stop('Error: Please use color_levels "High", "Mid", "Low" with a Target coloring')
    }
    segmentAnnotations[, color_by] <- unlist(log2(targetCountMatrix[color_by, ]))
    colType <- 'Target'
  } else {
    segmentAnnotations
    if(!all(unique(segmentAnnotations[, color_by] %in% color_levels))) {
      stop('Error: Mismatch in levels used with colors')
    }
    colType <- 'Annot'
  }
  params$colType <- colType
  
  # Error catch - shape must be an annotation factor or tag
  if(!shape_by %in% colnames(segmentAnnotations)) {
    stop('Shape parameter not found in Segment Annotations\n')
  }
  
  # iterate through appropriate dimensions
  plt_list <- list()
  dims <- grep("Dim[0-9]", colnames(segmentAnnotations))
  dims <- colnames(segmentAnnotations)[dims]
  dims <- combn(dims, 2)
  for(i in 1:ncol(dims)) {
    plt_list[i] <- 
      plot_DR(targetCountMatrix = targetCountMatrix,
              segmentAnnotations = segmentAnnotations,
              dims = unlist(dims[,i]),
              params = params,
              var_est = var_est)
  }
  
  # Step 4: Save Files
  # Save File      
  ggsave(filename = paste0(plot_type, "_with_",
                           color_by, "_and_",
                           shape_by, ".png"),
         plot = plt,
         device = "png", # type of plot
         dpi = 300,      # print DPI for print quality
         units = "in",   # inches
         width = 8,      # default 8 inches
         height = 6,     # default 6 inches
         path = outputFolder)
  
  #Variance Estimate Plot
  if(plot_type == "PCA") {
    v_dat <- data.frame(Variance = 100*var_est,
                        PC = 1:length(var_est))
    vplt <- ggplot(v_dat[1:min(15, nrow(v_dat)), ], 
                   aes(x = PC, y = Variance, fill = Variance)) +
      geom_bar(stat = 'identity') +
      geom_hline(yintercept = c(25,50,75), lty = 'dashed', color = 'black') +
      theme_bw(base_size = 15) +
      scale_y_continuous(expand = expansion(mult = 0),
                         limits = c(0, 100)) +
      scale_x_continuous(expand = expansion(mult = 0), limits = c(0, 16), breaks = seq(1,15,1)) +
      scale_fill_gradient(low = 'darkgray', high = 'dodgerblue2', limits = c(0,100)) +
      labs(x = 'Principal Component #', y = 'Cumulative Variance Explained') +
      theme(legend.position = "none")
    ggsave(filename = "PCA_VarianceExplained.png",
           plot = vplt,
           device = "png", # type of plot
           dpi = 300,      # print DPI for print quality
           units = "in",   # inches
           width = 8,      # default 8 inches
           height = 6,     # default 6 inches
           path = outputFolder)
  }
  
  # Save CSV table with plot dimension data for re-graphing
  write.csv(segmentAnnotations,
            file = file.path(outputFolder, paste0("Annotations with", plot_type, " Dims.csv"),
                             fsep = .Platform$file.sep),
            row.names = FALSE)
}

##### Helper Functions #####

# calc_DR: calculate dimension reduction values
# inputs: plot_type - which plot to generate
#         targetCountMatrix - data to work from
#         segmentAnnotations - annotations to use
calc_DR <- function(plot_type = NULL,
                    targetCountMatrix = NULL,
                    segmentAnnotations = NULL) {
  if(plot_type == "UMAP") {
    set.seed(seed = 28502)
    dr_data <- umap(t(log2(targetCountMatrix)))
    segmentAnnotations$Dim1 <- dr_data$layout[, 1]
    segmentAnnotations$Dim2 <- dr_data$layout[, 2]
  } else if(plot_type == "tSNE") {
    # prevent perplexity errors by automatically setting it to near, but not at,
    # max possible
    tSNE_perplexity <- floor((ncol(targetCountMatrix)-1)*.2) 
    set.seed(seed = 28502)
    dr_data <- Rtsne(t(log2(targetCountMatrix)), perplexity = tSNE_perplexity)
    segmentAnnotations$Dim1 <- dr_data$Y[, 1]
    segmentAnnotations$Dim2 <- dr_data$Y[, 2]
  } else if(plot_type == "PCA") {
    dr_data <- prcomp(t(log2(targetCountMatrix)))
    segmentAnnotations$Dim1 <- dr_data$x[, 1]
    segmentAnnotations$Dim2 <- dr_data$x[, 2]
    segmentAnnotations$Dim3 <- dr_data$x[, 3]
  } else {
    stop('Error: Additional plot types not yet supported\n')
  }
  rtn <- list(annot = segmentAnnotations, data = dr_data)
  if(plot_type == "PCA") {
    rtn <- c(rtn, list(var = var_est))
  }
  return(rtn)
}

# plot_DR: plot dim reduction based on above variables
# inputs: targetCountMatrix - data to work from
#         segmentAnnotations - annotations to use
#         plot_type - which plot to generate
#         dims - names of dims to plot (x, y)
#         color_by - tag, factor, or target
#         shape_by - tag, or factor
#         plot_color_theme - color palette
#         reverse_theme = TRUE # reverse palette color
#         var_est - variance estimates for PCA

plot_DR <- function(targetCountMatrix = NULL,
                    segmentAnnotations = NULL,
                    dims = c("Dim1", "Dim2"),
                    params = params,
                    var_est = NULL) {
  # graph setup
  plt <- ggplot(segmentAnnotations,
                aes_string(x = dims[1], y = dims[2])) + 
    geom_point(aes_string(shape = params$shape_by,   # if either of these is null it will still graph
                          color = params$color_by),
               size = 3.5, alpha = 0.8) +
    theme_bw(base_size = params$plot_font$size) +
    theme(aspect.ratio = 1,
          text = element_text(family = params$plot_font$family)) 
  
  # Add variance estimates to PCA plot
  if(params$plot_type == "PCA") {
    plt <- plt + 
      labs(x = paste0("PC1 (Var = ", round(100*var_est[1], 1), "%)"),
           y = paste0("PC2 (Var = ", round(100*var_est[2], 1), "%)"))
  } else {
    plt <- plt +
      labs(x = paste0(params$plot_type, " Dimension 1"),
           y = paste0(params$plot_type, " Dimension 2"))
  }
  
  # Add color theme
  if(params$colType == "Annot") {
    plt <- plt +
      scale_color_manual(values = params$plot_colors)
  } else if(params$colType == "Target") {
    clr_name <- paste0(params$color_by, ',\nLog2 Counts')
    if(length(params$plot_colors) == 3) {
      plt <- plt +
        scale_color_gradient2(name = clr_name,
                              low = params$plot_colors$Low,
                              mid = params$plot_colors$Mid,
                              high = params$plot_colors$High,
                              midpoint = median(segmentAnnotations[, params$color_by]))
    } else {
      plt <- plt +
        scale_color_gradient(name = clr_name,
                             low = params$plot_colors$Low, 
                             high = params$plot_colors$High)
    }
  }
  return(plt)
}

# annot colors:
# if(is.null(params$plot_color_theme)) {
#   cols <- params$plot_colors
# } else {
#   n_cols <- length(unique(segmentAnnotations[, params$color_by]))
#   if(n_cols < brewer.pal.info[params$plot_color_theme, 1]) {
#     cols <- brewer.pal(n = max(3, n_cols),
#                        name = params$plot_color_theme)
#     cols <- cols[1:n_cols]
#   } else {
#     cols <- brewer.pal(n = brewer.pal.info[params$plot_color_theme, 1],
#                        name = params$plot_color_theme)
#     cols <- colorRampPalette(cols)(n_cols)
#   }
#   if(reverse_theme) {
#     cols <- rev(cols)
#   }
#   names(cols) <- names(params$plot_colors)
# }
#
# target colors:

# otherwise grab the palette, and determine where to replace the 'light'
# color with gray instead for visualization on a white background
# cols <- brewer.pal(n = brewer.pal.info[plot_color_theme, 1],
#                    name = plot_color_theme)
# if(reverse_theme) {
#   cols <- rev(cols)
# }
# if(brewer.pal.info[plot_color_theme, 2] == 'seq') {
#   # if sequential palette use gray -> color
#   # use 1 away from ends of palettes as these are brighter than the final colors
#   plt <- plt +
#     scale_color_gradient(name = clr_name,
#                          low = 'gray', high = cols[length(cols)-1])
#   
# } else {
#   # if divergent use color1 -> gray -> color2
#   plt <- plt +
#     scale_color_gradient2(name = clr_name,
#                           low = cols[2], mid = 'gray', high = cols[length(cols)-1],
#                           midpoint = median(segmentAnnotations[, color_by]))
# }