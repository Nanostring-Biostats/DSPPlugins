# Type of Plot:
plot_type <- "tSNE"
# Options: tSNE, UMAP, PCA

# Plot Parameters
plot_color = "CD68" # tag, factor, or target
plot_shape = "SegmentName" # tag, or factor
plot_font = list(family = "sans", size = 15)
# shape & color can be set to NULL
# font families include sans, serif, mono and
#  may include specifically named fonts, but 
#  these may not always render properly.

# Plot colors
plot_color_theme = "RdBu"
reverse_theme = TRUE # reverse palette color
# Set to plot_color_theme = NULL if you'd 
# prefer to set your colors manually below
# (color name or hexadecimal (e.g. "#ABABAB"))
plot_colors = c("orange2", "black", "purple2")
color_levels = c("High", "Mid", "Low")

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
  
  # Step 2:
  # Calculate PCA, tSNE, or UMAP dimensions
  DR_ann <- calc_DR(plot_type = plot_type,
                    targetCountMatrix = targetCountMatrix,
                    segmentAnnotations = segmentAnnotations)
  segmentAnnotations <- DR_ann$annot
  if(plot_type == "PCA") {
    var_est <- DR_ann$var
  } else {
    var_est <- NULL
  }
  rm("DR_ann")
  
  # Step 3: Graph data
  # if color or shape is a gene symbol add it to the annotations for plotting
  if(plot_color %in% rownames(targetCountMatrix)) {
    segmentAnnotations[, plot_color] <- unlist(log2(targetCountMatrix[plot_color, ]))
    colType <- 'Target'
  } else {
    colType <- 'Annot'
  }
  # Error catch - shape must be an annotation factor or tag
  if(!plot_shape %in% colnames(segmentAnnotations)) {
    stop('Shape parameter not found in Segment Annotations\n')
  }
  
  # graph setup
  plt <- ggplot(segmentAnnotations,
                aes(x = Dim1, y = Dim2)) + 
    geom_point(aes_string(shape = plot_shape,
                          color = plot_color),
               size = 3, alpha = 0.8) +
    theme_bw(base_size = plot_font$size) +
    theme(aspect.ratio = 1,
          text = element_text(family = plot_font$family)) 
  
  # Add variance estimates to PCA plot
  if(plot_type == "PCA") {
    plt <- plt + 
      labs(x = paste0("PC1 (Var = ", round(100*var_est[1], 1), "%)"),
           y = paste0("PC2 (Var = ", round(100*var_est[2], 1), "%)"))
  } else {
    plt <- plt +
      labs(x = paste0(plot_type, " Dimension 1"),
           y = paste0(plot_type, " Dimension 2"))
  }
  
  # Update colors if an annotation is used
  if(colType == "Annot") {
    n_cols <- length(unique(segmentAnnotations[, plot_color]))
    if(n_cols < brewer.pal.info[plot_color_theme, 1]) {
      cols <- brewer.pal(n = max(3, n_cols),
                         name = plot_color_theme)
    } else {
      cols <- brewer.pal(n = brewer.pal.info[plot_color_theme, 1],
                         name = plot_color_theme)
      cols <- colorRampPalette(cols)(n_cols)
    }
    if(reverse_theme) {
      cols <- rev(cols)
    }
    plt <- plt +
      scale_color_manual(values = cols)
  } else {
    # for target based coloring, if set to qualitative palette, use 
    #    darkblue -> gray -> orange instead
    clr_name <- paste0(plot_color, ',\nLog2 Counts')
    if(brewer.pal.info[plot_color_theme, 2] == 'qual') {
      plt <- plt +
        scale_color_gradient2(name = clr_name,
                              low = 'darkblue', mid = 'gray', high = 'orange2',
                              midpoint = median(segmentAnnotations[, plot_color]))
    } else {
      # otherwise grab the palette, and determine where to replace the 'light'
      # color with gray instead for visualization on a white background
      cols <- brewer.pal(n = brewer.pal.info[plot_color_theme, 1],
                         name = plot_color_theme)
      if(reverse_theme) {
        cols <- rev(cols)
      }
      if(brewer.pal.info[plot_color_theme, 2] == 'seq') {
        # if sequential palette use gray -> color
        # use 1 away from ends of palettes as these are brighter than the final colors
        plt <- plt +
          scale_color_gradient(name = clr_name,
                               low = 'gray', high = cols[length(cols)-1])
        
      } else {
        # if divergent use color1 -> gray -> color2
        plt <- plt +
          scale_color_gradient2(name = clr_name,
                                low = cols[2], mid = 'gray', high = cols[length(cols)-1],
                                midpoint = median(segmentAnnotations[, plot_color]))
        
      }
    }
  }
  
  # Step 4: Save Files
  # Save File      
  ggsave(filename = paste0(plot_type, "_with_",
                           plot_color, "_and_",
                           plot_shape, ".png"),
         plot = plt,
         device = "png", # type of plot
         dpi = 300,      # print DPI for print quality
         units = "in",   # inches
         width = 8,      # default 8 inches
         height = 6,     # default 6 inches
         path = outputFolder)
  
  #Variance Estimate Plot
  if(plot_type == "PCA") {
    v_dat <- data.frame(Variance = 100*cumsum(var_est),
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
    umap_data <- umap(t(log2(targetCountMatrix)))
    segmentAnnotations$Dim1 <- umap_data$layout[, 1]
    segmentAnnotations$Dim2 <- umap_data$layout[, 2]
  } else if(plot_type == "tSNE") {
    # prevent perplexity errors by automatically setting it to near, but not at,
    # max possible
    tSNE_perplexity <- floor((ncol(targetCountMatrix)-1)*.2) 
    set.seed(seed = 28502)
    tsne_data <- Rtsne(t(log2(targetCountMatrix)), perplexity = tSNE_perplexity)
    segmentAnnotations$Dim1 <- tsne_data$Y[, 1]
    segmentAnnotations$Dim2 <- tsne_data$Y[, 2]
  } else if(plot_type == "PCA") {
    pca_data <- prcomp(t(log2(targetCountMatrix)))
    segmentAnnotations$Dim1 <- pca_data$x[, 1]
    segmentAnnotations$Dim2 <- pca_data$x[, 2]
    segmentAnnotations$Dim3 <- pca_data$x[, 3]
    var_est <- summary(pca_data)$importance[2, ]
  } else {
    stop('Error: Additional plot types not yet supported\n')
  }
  rtn <- list(annot = segmentAnnotations)
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
#         plot_color - tag, factor, or target
#         plot_shape - tag, or factor
#         plot_color_theme - color palette
#         reverse_theme = TRUE # reverse palette color
#         var_est - variance estimates for PCA

plot_DR <- function(targetCountMatrix = NULL,
                    segmentAnnotations = NULL,
                    plot_type = NULL,
                    dims = c("Dim1", "Dim2"),
                    plot_color = NULL,
                    plot_shape = NULL,
                    plot_font = list(family = "Arial"),
                    plot_color_theme = NULL,
                    reverse_theme = FALSE,
                    var_est = NULL) {
  # pull from above
}