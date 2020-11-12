# Type of Plot:
plot_type <- "PCA"
# Options: tSNE, UMAP, PCA

# Plot Parameters
color_by = "SegmentName" # tag, factor, target or NULL
shape_by = "SegmentName" # tag, factor or NULL
plot_font = list(family = "sans", size = 15)
save_as = "pdf"
# font families include sans, serif, mono and
# may include specifically named fonts, but 
# these may not always render properly.

# Plot color options
plot_colors = list("green3", "cyan3")
color_levels = c("PanCK-pos", "PanCK-neg")
# color_levels must match values in the color_by 
# column in the annotations file when using a tag
# or a factor. For Targets use "High", "Low" and "Mid".
#   "Mid" is optional
# the first entry in plot_colors may be set to a palette 
# if coloring using annotations. the palette order shall
# be alphabetical unless all values of color_by are provided
# in color_levels, in which case that order shall be used

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
library(openxlsx)

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
    var_est <- summary(DR_ann$data)$importance[2, ]
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
    lvls <- unique(segmentAnnotations[, color_by])
    # if(!all(lvls %in% color_levels)) {
    #   stop('Error: Mismatch in levels used with colors')
    # }
    colType <- 'Annot'
    # allow users to pass no or incomplete color levels for annotations
    if(is.null(color_levels)) {
      color_levels <- lvls
    } else if(!all(color_levels %in% lvls)) {
      new_lvls <- lvls[!lvls %in% color_levels]
      color_levels <- c(color_levels[color_levels %in% lvls],
                        new_lvls)
      if(length(color_levels) > length(plot_colors) & 
         (!plot_colors[[1]] %in% rownames(brewer.pal.info))) {
        plot_colors <- c(plot_colors,
                         extend_palette(n = length(new_lvls)))
      }
      # need to think through the logic on this.
      # one option: set color_levels -> NULL to let it 
      #  figure it out on it's own
    }
  }
  
  # Error catch - shape must be an annotation factor or tag
  # if(!is.null(shape_by) & !shape_by %in% colnames(segmentAnnotations)) {
  #   stop('Shape parameter not found in Segment Annotations\n')
  # }
  
  # # check color levels, grab new ones if they aren't
  # # in the appropriate order
  # lvls <- unique(segmentAnnotations[, color_by])
  # if(!all(color_levels %in% lvls)) {
  #   color_levels <- c(color_levels[color_levels %in% lvls],
  #                     lvls[!lvls %in% color_levels])
  # }
  # need to think through the logic on this.
  # one option: set color_levels -> NULL to let it 
  #  figure it out on it's own
  
  # gather color parameters
  names(plot_colors) <- color_levels[1:length(plot_colors)]
  params <- list("plot_type" = plot_type,
                 "color_by" = color_by,
                 "shape_by" = shape_by,
                 "plot_font" = plot_font,
                 "plot_colors" = plot_colors,
                 "color_levels" = color_levels,
                 "colType" = colType)
  
  # iterate through appropriate dimensions
  plt_list <- list()
  dims <- grep("Dim[0-9]", colnames(segmentAnnotations))
  dims <- colnames(segmentAnnotations)[dims]
  dims <- combn(dims, 2)
  for(i in 1:ncol(dims)) {
    plt_list[[paste0(dims[1,i],dims[2,i])]] <- 
      plot_DR(targetCountMatrix = targetCountMatrix,
              segmentAnnotations = segmentAnnotations,
              dims = unlist(dims[,i]),
              params = params,
              var_est = var_est)
  }
  
  # If PCA - create variance explained plot:
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
  }
  
  # Step 4: Save Files
  # Save File. PDF is first case
  if(tolower(save_as) == 'pdf') {
    pdf(file = file.path(outputFolder, 
                         paste0(plot_type, "_with_",
                                color_by, "_and_",
                                shape_by, ".pdf"),
                         fsep = .Platform$file.sep),
        width = 8,
        height = 6)
    print(plt_list)
    # Also save Variance Estimate Plot - PCA only
    if(plot_type == "PCA") {
      print(vplt)
    }
    dev.off()
  # Save individual figures for all graphs if not using a PDF
  } else {
    for(plt in names(plt_list)) {
      ggsave(filename = file.path(outputFolder, 
                                  paste0(plot_type, "_", plt, "_with_",
                                         color_by, "_and_",
                                         shape_by, ".", tolower(save_as)),
                                  fsep = .Platform$file.sep),
             plot = plt_list[[plt]],
             device = tolower(save_as),
             units = "in",
             dpi = 300,
             width = 8,
             height = 6)
    }
    if(plot_type == "PCA") {
      ggsave(filename = file.path(outputFolder, 
                                  paste0("PCAVariancePlot.", tolower(save_as)),
                                  fsep = .Platform$file.sep),
             plot = vplt,
             device = tolower(save_as), 
             units = "in",
             dpi = 300,
             width = 8,
             height = 6)
    }
  }
  
  # Save XLSX table with plot dimension data for re-graphing
  wb <- createWorkbook("DimensionReduction")
  addWorksheet(wb, "Segment Annotations")
  writeData(wb,
            sheet = "Segment Annotations",
            segmentAnnotations,
            colNames = TRUE, rowNames = FALSE)
  if(plot_type == "PCA") {
    addWorksheet(wb, "PC Loadings")
    writeData(wb,
              sheet = "PC Loadings",
              DR_ann$data$rotation,
              colNames = TRUE, rowNames = TRUE)
    addWorksheet(wb, "Variance Estimates")
    writeData(wb,
              sheet = "Variance Estimates",
              as.data.frame(var_est),
              colNames = FALSE, rowNames = TRUE)
  }
  saveWorkbook(wb,
               file = file.path(outputFolder,
                                paste0(plot_type, " Data.xlsx"),
                                fsep = .Platform$file.sep),
               overwrite = TRUE)
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
    dr_data <- umap(t(log2(targetCountMatrix)), random.state = 28502)
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
#         dims - names of dims to plot (x, y)
#         params - list of parameters for the plot
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
    d1 <- as.numeric(gsub("Dim", "", dims[1]))
    d2 <- as.numeric(gsub("Dim", "", dims[2]))
    plt <- plt + 
      labs(x = paste0("PC", d1,
                      " (Var = ", round(100*var_est[d1], 1), "%)"),
           y = paste0("PC", d2,
                      " (Var = ", round(100*var_est[d2], 1), "%)"))
  } else {
    plt <- plt +
      labs(x = paste0(params$plot_type, " Dimension 1"),
           y = paste0(params$plot_type, " Dimension 2"))
  }
  
  # Add colors if provided
  if(params$colType == "Annot") {
    # use a palette if provided, & check that all levels are provided
    if(params$plot_colors[[1]] %in% rownames(brewer.pal.info)) {
      lvls <- unique(segmentAnnotations[, params$color_by])
      n_cols <- length(lvls)
      # if # levels < palette max, use default colors
      if(n_cols < brewer.pal.info[params$plot_colors[[1]], 1]) {
        cols <- brewer.pal(n = max(3, n_cols),
                           name = params$plot_colors[[1]])
        cols <- cols[1:n_cols]
      # otherwise extrapolate colors to larger space
      } else {
        cols <- brewer.pal(n = brewer.pal.info[params$plot_colors[[1]], 1],
                           name = params$plot_colors[[1]])
        cols <- colorRampPalette(cols)(n_cols)
      }
      # if all levels in the color_levels variable
      if(all(lvls %in% color_levels)) {
        color_levels <- color_levels[color_levels %in% lvls]
        names(cols) <- color_levels
      }
    plt <- plt + 
      scale_color_manual(values = cols)
    # if a list of colors is provided use those
    } else {
      plt <- plt +
        scale_color_manual(values = params$plot_colors)
    }
  # coloring by Target value
  } else if(params$colType == "Target") {
    clr_name <- paste0(params$color_by, ',\nLog2 Counts')
    # 3 point gradient
    if(all(c("High","Mid","Low") %in% names(params$plot_colors))) {
      plt <- plt +
        scale_color_gradient2(name = clr_name,
                              low = params$plot_colors$Low,
                              mid = params$plot_colors$Mid,
                              high = params$plot_colors$High,
                              midpoint = median(segmentAnnotations[, params$color_by]))
    } else if(all(c("High","Low") %in% names(params$plot_colors))) {
      # 2 point gradient
      plt <- plt +
        scale_color_gradient(name = clr_name,
                             low = params$plot_colors$Low, 
                             high = params$plot_colors$High)
    } else {
      stop("Error: Unexpected labels for plot levels, please use 'High', 'Mid', 'Low'")
    }
  }
  return(plt)
}

# extend_palette: create a brewer.pal color list
#   not limited by size of palette
extend_palette <- function(palette = 'Set1',
                           n = 9) {
  pal_n <- brewer.pal.info[palette, 1]
  if(n <= pal_n) {
    pal <- brewer.pal(n, palette)
  } else {
    pal <- brewer.pal(pal_n, palette)
    pal <- colorRampPalette(pal)(n)
  }
  return(pal)
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