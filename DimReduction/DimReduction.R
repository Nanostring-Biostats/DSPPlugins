# Dimension Reduction #
# Version 1.2 #

# Produces PCA, tSNE, or UMAP visualizations
# Supports: DSP-nCounter Protein, DSP-nCounter RNA, DSP-NGS CTA, DSP-NGS WTA (mouse & human)
# Note: this script should be run only on a dataset AFTER normalization
# Please do not use spaces, special characters, or numbers when adding factors
# in the DSPDA Annotation file

#        User Options         #
###############################

plot_type = "PCA"
# Options: tSNE, UMAP, PCA

# Plot Parameters

color_by = "ScanName"
# select factor OR target of interest OR NULL
shape_by = "SegmentName"
# select factor of interest OR NULL
size_by = NULL
# to size by a target, replace NULL with target of interest (ex: "CD68")

plot_font = list(family = "sans",
                 size = 15)
# font families include sans, serif, mono and
# may include specifically named fonts, but 
# these may not always render.

save_as = "pdf"
# options: pdf, jpeg, tiff, png, bmp, or svg

# customize colors for datapoints
plot_colors = NULL
color_levels = NULL
# color_levels *must* match the names in the
# color_by factor column in the DSPDA annotation file 
# when using the factor of interest. color_levels 
# may be set to NULL, and the plugin will 
# automatically assign colors to the levels.
#
# When coloring by a target use "High", "Low" 
# and "Mid". "Mid" is optional.
# 
# Valid named color choices for plot_colors can be  
# found in the vignette or standard hexidecimal
# colors (e.g. '#3311FF')
# 
# for example, if we: color_by = "CD68", then
# plot_colors = list("red", "white", "blue")
# color_levels = c("High", "Mid", "Low")
#
# The first entry in plot_colors may also be set 
# to an *r color palette* when coloring using 
# annotations. For valid palette choices, refer
# to the vignette.
# The palette order shall be 
# alphabetical unless values of color_by 
# are provided in color_levels, in which 
# case that order shall be used.

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
  # check that requested output type is valid option
  if(!save_as %in% c('pdf', 'jpeg', 'tiff', 'png', 'bmp', 'svg')) {
    fail(message = 'Invalid filetype chosen. Please choose from: pdf, jpeg, tiff, png, bmp, or svg')
  }
  
  # if color is a gene symbol add it to the annotations for plotting, catch NA or NULL
  if(length(color_by) > 1) {
    fail(message = 'Color criteria must be only one value! Please pass in only one string value for color_by.')
  } else if(sum(is.null(color_by)) > 0 | sum(is.na(color_by)) > 0) {
    colType <- 'Null'
  } else if(!color_by %in% rownames(targetCountMatrix) &
            !color_by %in% colnames(segmentAnnotations)) {
    fail(message = 'Color_by value not found. Please confirm that your color feature is either a column in Segment Properties or a Target name from your target count matrix')
  } else if (color_by %in% rownames(targetCountMatrix) &
             color_by %in% colnames(segmentAnnotations)) {
    fail(message = 'Color choice found as a Target *and* in a column in Segment Properties. Please alter the column name or drop the gene to resolve the conflict.')
  } else if(color_by %in% rownames(targetCountMatrix)) {
    colType <- 'Target'
    segmentAnnotations[, color_by] <- unlist(log2(targetCountMatrix[color_by, ]))
    if(sum(is.null(color_levels)) > 0 | sum(is.na(color_levels)) > 0) {
      all_lvls <- color_levels <- c("High", "Low") #enforced default
    } else if(!all(color_levels %in% c("High", "Mid", "Low"))) {
      fail(message = 'Incorrect color level definition. Please use color_levels = c("High", "Mid", "Low") or c("High", "Low") when using a Target for coloring')
    } else if(!all(c("High", "Low") %in% color_levels)) {
      fail('Color Endpoints not defined, please ensure "High" and "Low" levels are included in color_levels.')
    } else {
      all_lvls <- color_levels
    }
  } else if (color_by %in% colnames(segmentAnnotations)) {
    colType <- 'Annot'
    all_lvls <- as.character(unique(segmentAnnotations[, color_by]))
    # allow users to pass no or incomplete color levels for annotations
    if(sum(is.null(color_levels)) > 0 | sum(is.na(color_levels)) > 0) {
      color_levels <- all_lvls # override if NA or NULL detected
    } else if(!all(color_levels %in% all_lvls)) {
      fail(paste0('Invalid level(s) chosen in color_levels not found in Segment Properties: ', 
                  toString(color_levels[which(!color_levels %in% all_lvls)])))
    } else if(!all(all_lvls %in% color_levels)) {
      color_levels <- c(color_levels, all_lvls[!all_lvls %in% color_levels])
    }
  }
  
  # color checking, lengthen palette to assign to all levels 
  # catch case of any NULL or NA value, as well as no color
  if (colType != 'Null') {
    if(sum(is.null(plot_colors)) == 0 & sum(is.na(plot_colors)) == 0) {
      # test for valid colors, overridden if first value is a valid palette
      if(sum(!(are_valid_colors(plot_colors) | 
                (plot_colors[[1]] %in% rownames(brewer.pal.info)))) > 0) {
        fail(message = paste0('Invalid color choice(s): ',
                              toString(plot_colors[which((are_valid_colors(plot_colors) | 
                                 (plot_colors[[1]] %in% rownames(brewer.pal.info))) == 0)]),
                              '. Please use an RBrewer palette, hexidecimal colors, or valid color. Find links in plug-in vignette.'))
      } else if(plot_colors[[1]] %in% rownames(brewer.pal.info)) { # use provided palette
        plot_colors <- c(extend_palette(palette = plot_colors[[1]], n = length(all_lvls)))
      } else if(length(color_levels) > length(plot_colors)) { # fill in more colors if needed
        plot_colors <- c(plot_colors, extend_palette(n = length(all_lvls[!all_lvls %in% color_levels])))
      }
      # cull if necessary
      if(length(plot_colors) > length(color_levels)) plot_colors <- plot_colors[1:length(color_levels)]
    } else {
      plot_colors <- extend_palette(n = length(all_lvls)) # override if NA or NULL detected
    }
  }

  # Size by calculation:
  if(!is.null(size_by)) {
    if (!is.na(size_by)) {
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
  } else if(tolower(save_as) == 'svg') {
    for(plt in names(plt_list)) {
      svg(file = file.path(outputFolder, 
                           paste0(save_params, "_", plt, ".svg"),
                           fsep = .Platform$file.sep),
          width = 8,
          height = 6)
      print(plt_list[[plt]])
      dev.off()
      if(plot_type == "PCA") {
        svg(file = file.path(outputFolder, 
                             "PCAVariancePlot.svg",
                             fsep = .Platform$file.sep),
            width = 8,
            height = 6)
        print(vplt)
        dev.off()
      }
    }
    # Save individual figures for all graphs if not using a PDF
  } else {
    for(plt in names(plt_list)) {
      ggsave(filename = file.path(outputFolder, 
                                  paste0(save_params, '_', plt, ".",
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
    dr_data <- prcomp(t(log2(targetCountMatrix)), scale. = TRUE)
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
      plt <- plt +
        scale_color_manual(values = params$plot_colors)
    # coloring by Target value
  } else if(params$colType == "Target") {
    clr_name <- paste0(params$color_by, ',\nLog2 Counts')
    # 3 point gradient
    if(all(c("High","Mid","Low") %in% names(params$plot_colors))) {
      rng <- range(segmentAnnotations[, params$color_by])
      mdpt <- rng[2] - (rng[2]-rng[1])/2
      plt <- plt +
        scale_color_gradient2(name = clr_name,
                              low = params$plot_colors[['Low']],
                              mid = params$plot_colors[['Mid']],
                              high = params$plot_colors[['High']],
                              midpoint = mdpt,
                              guide = "colorbar")
    } else if(all(c("High","Low") %in% names(params$plot_colors))) {
      # 2 point gradient
      plt <- plt +
        scale_color_gradient(name = clr_name,
                             low = params$plot_colors[['Low']], 
                             high = params$plot_colors[['High']],
                             guide = "colorbar")
    }
  }
  return(plt)
}

# extend_palette: create a brewer.pal color list
#   not limited by size of palette
# allow for extension with < 3 requested avoiding warning
extend_palette <- function(palette = 'Set1',
                           n = 9) {
  pal_n <- brewer.pal.info[palette, 1]
  if(n <= pal_n) {
    pal <- brewer.pal(max(n,3), palette)
  } else {
    pal <- brewer.pal(pal_n, palette)
    pal <- colorRampPalette(pal)(n)
  }
  return(pal[1:n])
}

# test supplied colors for validity
are_valid_colors <- function(clrs) {
  sapply(clrs, function(x) {
    tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE)
  })
}
