#        User Options        #
##############################
plot_type = "PCA"
# Options: tSNE, UMAP, PCA

# Plot Parameters
color_by = "SegmentName"
# tag, factor, target or NULL
shape_by = "SegmentName"
# tag, factor or NULL
size_by = "CD68"
# target only or NULL

plot_font = list(family = "sans",
                 size = 15)
# font families include sans, serif, mono and
# may include specifically named fonts, but 
# these may not always render properly.

save_as = "pdf"
# options: pdf, jpeg, tiff, png, bmp, or svg

plot_colors = list("green3", "cyan3")
color_levels = c("PanCK-pos", "PanCK-neg")
# color_levels must match values in the
# color_by column in the annotations file 
# when using a tag or a factor. color_levels 
# may be set to NULL, and the plugin will 
# automatically assign levels to colors
#
# When coloring by a target use "High", "Low" 
# and "Mid". "Mid" is optional.
#
# the first entry in plot_colors may be set 
# to a r color palette if coloring using 
# annotations. the palette order shall be 
# alphabetical unless values of color_by 
# are provided in color_levels, in which 
# case that order shall be used

##############################
# Copyright (C) 2020, NanoString Technologies, Inc.
#    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#    You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us: 
#   NanoString Technologies, Inc.
#   530 Fairview Avenue N
#   Seattle, WA 98109
# Tel: (888) 358-6266
##############################

##############################
#        Execution Code      # 
##############################
# Libraries:
library(Rtsne)        # for t-SNE
library(umap)         # for UMAP
library(stats)        # for princomp
library(ggplot2)      # for ggplot
library(reshape2)     # for melt/dcast
library(dplyr)        # for melt/dcast
library(RColorBrewer) # for brewer.pal
library(openxlsx)     # for *Workbook
library(testthat)     # for fail

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
  # if color is a gene symbol add it to the annotations for plotting
  if(is.null(color_by)) {
    colType <- 'Null'
  } else if(!color_by %in% rownames(targetCountMatrix) &
            !color_by %in% colnames(segmentAnnotations)) {
    fail(message = 'Color not found. Please confirm that your color feature is either a column in Segment Properties or a Target name from your target count matrix')
  } else if(color_by %in% rownames(targetCountMatrix)) {
    if(!all(color_levels %in% c("High", "Mid", "Low"))) {
      fail(message = 'Incorrect color level definition. Please use color_levels = c("High", "Mid", "Low") or c("High", "Low") when using a Target for coloring')
    }
    segmentAnnotations[, color_by] <- unlist(log2(targetCountMatrix[color_by, ]))
    colType <- 'Target'
  } else {
    lvls <- unique(segmentAnnotations[, color_by])
    colType <- 'Annot'
    # allow users to pass no or incomplete color levels for annotations
    if(is.null(color_levels)) {
      color_levels <- lvls
      
    } else if(!all(color_levels %in% lvls)) {
      new_lvls <- lvls[!lvls %in% color_levels]
      color_levels <- c(color_levels[color_levels %in% lvls],
                        new_lvls)
    }
    if(length(color_levels) > length(plot_colors) & 
       (!plot_colors[[1]] %in% rownames(brewer.pal.info))) {
      plot_colors <- c(plot_colors,
                       extend_palette(n = length(new_lvls)))
    }
  }
  
  # Size by calculation:
  if(!is.null(size_by)) {
    if(size_by %in% rownames(targetCountMatrix)) {
      if(isTRUE(size_by == color_by)) {
        trg <- unlist(targetCountMatrix[size_by, ])
        size_by <- paste0(size_by,'_linear')
        segmentAnnotations[[size_by]] <- trg
      } else {
        segmentAnnotations[[size_by]] <- unlist(targetCountMatrix[size_by, ])
      }
    } else {
      fail(message = 'Size not found. Please check to ensure that the target name specified for size is in the target count matrix')
    }
  }

  # gather color parameters
  names(plot_colors) <- color_levels[1:length(plot_colors)]
  params <- list("plot_type" = plot_type,
                 "color_by" = color_by,
                 "shape_by" = shape_by,
                 "size_by" = size_by,
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
  save_params <- c(plot_type, "with", color_by, shape_by, size_by)
  save_params <- save_params[!is.null(save_params)]
  if(length(save_params) == 2) {
    save_params <- plot_type
  } else {
    save_params <- paste(save_params, collapse = '_')
  }
  if(tolower(save_as) == 'pdf') {
    pdf(file = file.path(outputFolder, 
                         paste0(save_params, ".pdf"),
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
                                  paste0(save_params, ".",
                                         tolower(save_as)),
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
                                  paste0("PCAVariancePlot.",
                                         tolower(save_as)),
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
  wb <- createWorkbook(title = paste("Dimension Reduction,", plot_type), 
                       creator = "NanoString DSP Plugin")
  addWorksheet(wb = wb, "Segment Annotations")
  writeData(wb = wb,
            sheet = "Segment Annotations",
            x = segmentAnnotations,
            colNames = TRUE, rowNames = FALSE)
  if(plot_type == "PCA") {
    addWorksheet(wb, "Principal Components (All)")
    writeData(wb = wb,
              sheet = "Principal Components (All)",
              x = DR_ann$data$x,
              colNames = TRUE, rowNames = TRUE)
    setColWidths(wb = wb, sheet = "Principal Components (All)", cols = 1, widths = "auto")
    addWorksheet(wb, "PC Loadings")
    writeData(wb = wb,
              sheet = "PC Loadings",
              x = DR_ann$data$rotation[order(abs(DR_ann$data$rotation[, 1]), decreasing = TRUE), ],
              colNames = TRUE, rowNames = TRUE)
    setColWidths(wb = wb, sheet = "PC Loadings", cols = 1, widths = "auto")
    addWorksheet(wb, "Variance Estimates")
    writeData(wb = wb,
              sheet = "Variance Estimates", 
              x = t(summary(DR_ann$data)$importance),
              colNames = TRUE, rowNames = TRUE)
    setColWidths(wb = wb, sheet = "Variance Estimates", cols = 1:4, widths = "auto")
  }
  saveWorkbook(wb = wb,
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
    fail('Plot type not found. Please use "PCA", "tSNE", or "UMAP"')
  }
  rtn <- list(annot = segmentAnnotations, data = dr_data)
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
    theme_bw(base_size = params$plot_font$size) +
    theme(aspect.ratio = 1,
          text = element_text(family = params$plot_font$family)) 
  
  # Add size
 if(!is.null(params$size_by)) {
   sz_ttl <- paste0(gsub("_linear", "", params$size_by),
                    ",\nCounts")
   plt <- plt +
      geom_point(aes_string(shape = params$shape_by,   # if any of these is null it will still graph
                            color = params$color_by,
                            size = params$size_by),
                 alpha = 0.8) +
      scale_size_continuous(name = sz_ttl, range = c(2,7))
  } else {
    plt <- plt +
      geom_point(aes_string(shape = params$shape_by,   # if any of these is null it will still graph
                            color = params$color_by),
                 size = 3.5, alpha = 0.8)
  }
  
  # Keep legend symbol size consistent
  if(params$colType == "Target") {
    plt <- plt +
      guides(shape = guide_legend(override.aes = list(size = 4)))
  } else {
    plt <- plt +
      guides(shape = guide_legend(override.aes = list(size = 4)),
             color = guide_legend(override.aes = list(size = 4)))
  }
  
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
      rng <- range(segmentAnnotations[, params$color_by])
      mdpt <- rng[2] - (rng[2]-rng[1])/2
      plt <- plt +
        scale_color_gradient2(name = clr_name,
                              low = params$plot_colors$Low,
                              mid = params$plot_colors$Mid,
                              high = params$plot_colors$High,
                              midpoint = mdpt,
                              guide = "colorbar")
    } else if(all(c("High","Low") %in% names(params$plot_colors))) {
      # 2 point gradient
      plt <- plt +
        scale_color_gradient(name = clr_name,
                             low = params$plot_colors$Low, 
                             high = params$plot_colors$High,
                             guide = "colorbar")
    } else {
      fail('Color Endpoints not defined, please ensure "High" and "Low" levels are included in color_levels.')
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