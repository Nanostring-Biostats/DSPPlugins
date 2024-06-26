# Combined labeled Volcano Plot
# Version 1.0 #

# Produces a Labeled Combined Volcano Plot from 2 DSPDA volcano plot inputs
# Supports: DSP-nCounter Protein, DSP-nCounter RNA, DSP-NGS CTA, DSP-NGS WTA (mouse & human)
# Note: this script can be run only after a DSPDA 'Statistical Tests' is performed 
# Please do not use spaces, special characters, or numbers when naming factors
# that are used in a statistical test

#        User Options        #
##############################
# Plot Parameters

output_format <- "png"
# options: png, jpg, tiff, svg, pdf, & bmp

######################## LABELING ########################  
plot_title <- "Enter Title Here"
# Volcano Plot Title

negative_label <- "Left Label"
# Negative (left) label for Fold Change from volcano plot

positive_label <- "Right Label"
# Positive (right) label for Fold Change from volcano plot

n_genes <- 10
# Number of top genes to label, gene_list overrides this variable if set

n_genes_per_dataset <- TRUE
# Set to TRUE to label top n_genes for each dataset, gene_list overrides this variable if set

gene_list <-NULL #c("PLOD3", "TMEM132E", "PITRM1", "DPP10", "OPA1", "CD4", "Wilms Tumor Protein", "IL-6R", "Nanog", "c-Maf", "Androgen Receptor")
# targets to label (genes or proteins) no matter where they are in the volcano plot. This overrides n_genes

# Optional Labels

show_legend <- TRUE
# Show color legend on figure


####################### THRESHOLDS #######################
pval_thresh <- 0.05
# P-value threshold, must set fdr_thresh to NULL to use

fdr_thresh <- 0.01
# FDR threshold, default threshold over pval_thresh

# Optional thresholds
fc_thresh <- 0.5
# Fold Change threshold, if set coloring options will change

label_fc <- FALSE
# Should genes below fc_thresh be labeled
######################### FONTS ##########################
font_size <- 8
# Font Size

label_size <- 2
# Label Font Size

font_family <- "sans"
# options include: serif, sans, mono

####################### PLOT SIZE ########################
plot_width <- 6
# Plot Width in inches

plot_height <- 4
# Plot Height in inches

####################### COLORING #########################
default_color <- "grey65"
# Default point color for points not in target groups or
# above thresholds

# Color options for plot. Must have more colors than
#   number of target groups 
color_options <-  c("#3A6CA1", "#FFD861", "#CF4244", 
                    "#47BAB4", "#474747", "#EB739E", 
                    "#318026", "#A66293", "#F28E2B", 
                    "#8F6954")


# Optional Color parameters
fc_color <- "grey30"
# Color of points below fc_thresh but above pval and/or 
#   fdr thresholds. change to the same as default_color
#   if you don't want these targets called out

# Optional Target labeling
target_groups <- "Protein_SNR5_filtered"
  #NULL #c("Hemostasis", "DNA Repair")
# Specific target groups to color in the plot
#   If variable is set, genes are colored by given target 
#     group; else genes colored if they are above thresholds

# Option for dataset coloring for Combined Plots
color_by_dataset <- TRUE
# If TRUE, color by dataset_name

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
library(ggplot2)
library(ggrepel)
library(plotly)
library(testthat)
library(scales)
library(stats)
library(stringr)
library(readxl)

# main function with GeoMxSet wrapper
main <- function(obj1, obj2, obj3, obj4){
  if(class(obj1) == "NanoStringGeoMxSet"){
    dataset <- exprs(obj1)
    segmentAnnotations <- pData(obj1)
    targetAnnotations <- fData(obj1)
    outputFolder <- obj3
  }else{
    dataset <- obj1
    segmentAnnotations <- obj2
    targetAnnotations <- obj3
    outputFolder <- obj4
  }
  volcanoPlot(dataset = dataset,
              segmentAnnotations = segmentAnnotations,
              targetAnnotations = targetAnnotations, 
              outputFolder = outputFolder)
}

# volcano plot function 
volcanoPlot <- function(dataset, segmentAnnotations, targetAnnotations, outputFolder){
  # List all xlsx and txt files in the current directory
  xlsx_files <- list.files(pattern = "\\.xlsx$")
  txt_files <- list.files(pattern = "\\.txt$")
  
  # Ensure all files have the same extension
  if(length(xlsx_files) > 0 && length(txt_files) > 0){
    stop("Files must have the same extension, either all .xlsx or all .txt")
  }
  
  # Determine the file type to process
  if(length(xlsx_files) >= 2){
    files_to_read <- xlsx_files
    read_function <- function(file) custom_read_excel(file)
  } else if(length(txt_files) >= 2){
    files_to_read <- txt_files
    read_function <- function(file) read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  } else {
    stop("Not enough files to process. Ensure there are at least two .xlsx or .txt files.")
  }
  
  # Read the volcano plot files
  data_list <- lapply(files_to_read, function(file) {
    data <- read_function(file)
    data$source <- file
    return(data)
  })
  
  # Combine the data (assuming all have the same structure)
  combined_results <- do.call(rbind, data_list)
  
  # Ensure these are numeric
  # Geographical number formatting check
  combined_results$Log2 <- as.numeric(gsub(x = combined_results$Log2, 
                                           pattern = ",", replacement = "."))
  combined_results$Pvalue  <- as.numeric(gsub(x = combined_results$Pvalue, 
                                              pattern = ",", replacement = "."))
  combined_results$X.log10.pvalue <- as.numeric(gsub(x = combined_results$X.log10.pvalue, 
                                                     pattern = ",", replacement = "."))
  
  combined_results$FDR <- as.numeric(gsub(x = combined_results$FDR, 
                                          pattern = ",", replacement = "."))
  
  #ensure consistent capitalization with Target.Name column name
  colnames(combined_results)[which(tolower(colnames(combined_results)) == "target.name")] <- "Target.Name"
  
  # test for valid input variables
  testVariableFormats(de=combined_results)
  
  # create volcano plot
  volcanoPlot_results <- plotVolcano(de=combined_results)
  
  if(output_format == "svg"){
    svg(filename=paste0(outputFolder, "/volcano_plot_", plot_title, ".", output_format), 
        width=plot_width, height=plot_height)
    print(volcanoPlot_results$plot)
    dev.off()
#    }else if(output_format == "html"){
#      htmlwidgets::saveWidget(volcanoPlot_results$plot, file=paste0(outputFolder, "/volcano_plot_", plot_title, ".html"))
    }else{
    ggsave(filename=paste0("volcano_plot_", plot_title, ".", output_format),
           plot=volcanoPlot_results$plot,
           device=output_format,
           path=outputFolder,
           height=plot_height,
           width=plot_width)
  }
  write.csv(x=volcanoPlot_results$gene_labels, 
              file=file.path(outputFolder,
                             paste0("labeled_genes_", plot_title, ".csv"), 
                             fsep=.Platform$file.sep),
               quote=FALSE, row.names=FALSE)
} 

#' plotVolcano
#'
#' create volcano plot figure using user input variables and DE results from DSPDA
#' @param de data frame of DE results from DSPDA
#' @return ggFigure ggplot of volcano plot
#' @export
plotVolcano <- function(de){
  
  # determine highest x point to make volcano plot equal on both sides
  maxFC <- max(abs(de$Log2))
  
  maxPval <- min(de$Pvalue)
  
  # find closest FDR value to fdr_thresh and use that pvalue to add y axis cutoff line
  if(!is.null(fdr_thresh)){
    fdr_pval <- mean(de$Pvalue[which(abs(de$FDR - fdr_thresh) == 
                                       min(abs(de$FDR - fdr_thresh)))])
  }else{
    fdr_pval <- NULL
  }
  
  # create basic volcano plot with correct formatting
  ggFigure <- ggplot(de, aes(x=Log2, y=Pvalue, 
                             text = paste("Target:", Target.Name)))+
    geom_point(color=default_color)+
    labs(y="Pvalue",
         x=paste(negative_label, "<-", "log2(FC)", "->", positive_label, sep=" "),
         title=plot_title)+
    theme_bw(base_size = font_size) +
    theme(text = element_text(family = font_family))+
    scale_x_continuous(limits=c(-maxFC, maxFC))
  
  # this makes for easier testing if not running in DSPDA, will flip yaxis of graph
  # scale_y_continuous(trans=change_axis_revlog_trans(base=10),
  #                    labels=function(x) format(x, trim=TRUE, digits=4,
  #                                              scientific=ifelse(maxPval < 0.0001, TRUE, FALSE), 
  #                                              drop0trailing=TRUE)) 
  
  if(show_legend == FALSE){
    ggFigure <- ggFigure + theme(legend.position="none")
  }
  
  if(color_by_dataset == TRUE){
    dataset_groups <- unique(de$dataset_name)
    if(length(dataset_groups) > length(color_options)){
      stop("Not enough colors provided in color_options")
    }
    
    color_options <- color_options[1:length(dataset_groups)]
    names(color_options) <- dataset_groups
    
    color_label <- "Dataset"
    gene_coloring <- de
    gene_coloring$Target_coloring <- gene_coloring$dataset_name
  }else if(!is.null(target_groups)){   # subset de to only include genes either in specified target groups or above pval/fdr threshold
    probe_groups <- strsplit(de$Target.group.membership.s, split=", ")
    de_subset_list <- lapply(X=target_groups, 
                             FUN=function(x){
                               #return de rows with specified target group named in column 
                               return(cbind(de[grep(x=de$Target.group.membership.s, pattern=x),], x))
                             })
    names(de_subset_list) <- target_groups
    gene_coloring <- as.data.frame(do.call(rbind, de_subset_list))
    names(gene_coloring)[names(gene_coloring) == "x"] <- "Target_coloring"
    gene_coloring$Target_coloring <- str_wrap(gene_coloring$Target_coloring, width=45)
    if(length(target_groups) > length(color_options)){
      stop("Not enough colors provided in color_options")
    }
    
    color_options <- color_options[1:length(target_groups)]
    names(color_options) <- target_groups
    color_label <- "Target Group\nMembership"
  }else{
    color_options <- color_options[1:2]
    if(!is.null(fdr_thresh) & !is.null(pval_thresh)){
      if(fdr_pval < pval_thresh){
        gene_coloring <- de[which(de$Pvalue < pval_thresh),]
        label_thresh_low <- paste("pval <", pval_thresh)
        label_thresh_high <- paste("FDR <", fdr_thresh)
        high_thresh <- fdr_pval
      }else{
        gene_coloring <- de[which(de$Pvalue < fdr_pval),]
        label_thresh_low <- paste("FDR <", fdr_thresh)
        label_thresh_high <- paste("pval <", pval_thresh)
        high_thresh <- pval_thresh
      }
      
      # label points as positive or negative FC
      gene_coloring$Target_coloring <- ifelse(test=gene_coloring$Log2 < 0, 
                                              yes=negative_label, 
                                              no=positive_label)
      
      # label points as above pval or FDR threshold
      gene_coloring$Target_coloring <- ifelse(test=gene_coloring$Pvalue >= high_thresh, 
                                              yes=paste(label_thresh_low, gene_coloring$Target_coloring), 
                                              no=paste(label_thresh_high, gene_coloring$Target_coloring))
      
      # make color options have muted colors for higher threshold
      color_options <- c(color_options, muted(color_options, l=80))
      names(color_options) <- c(paste(label_thresh_high, negative_label), 
                                paste(label_thresh_high, positive_label),
                                paste(label_thresh_low, negative_label), 
                                paste(label_thresh_low, positive_label))
    }else{
      if(is.null(fdr_thresh)){
        gene_coloring <- de[which(de$Pvalue < pval_thresh),]
        label_thresh <- paste("pval <", pval_thresh)
      }else{
        gene_coloring <- de[which(de$FDR < fdr_thresh),]
        label_thresh <- paste("FDR <", fdr_thresh)
      }
      gene_coloring$Target_coloring <- ifelse(test=gene_coloring$Log2 < 0, 
                                              yes=paste(label_thresh, negative_label), 
                                              no=paste(label_thresh, positive_label))
      
      names(color_options) <- c(paste(label_thresh, negative_label), 
                                paste(label_thresh, positive_label))
    }
    
    color_label <- "Significance:"
    
    # color by fc_thresh if fc_thresh is not NULL
    if(!is.null(fc_thresh)){
      gene_coloring$Target_coloring[abs(gene_coloring$Log2) < fc_thresh] <- paste("FC <", fc_thresh)
      color_options <- c(color_options, fc_color)
      names(color_options)[length(color_options)] <- paste("FC <", fc_thresh)
    }
  }
  if(color_by_dataset == FALSE){ #for other points without colors specified  
  color_options <- c(color_options, default_color)
  names(color_options)[length(color_options)] <- "Not Specified"
  }
  # add coloring to ggplot
  ggFigure <- ggFigure + geom_point(data=gene_coloring, aes(x=Log2, y=Pvalue, color=Target_coloring))+
    labs(color=color_label)+
    scale_color_manual(values=color_options)
  
  # add threshold line values to y axis
  yaxis <- data.frame(brk=as.numeric(pretty_breaks(n=4)(0:max(de$X.log10.pvalue))))
  yaxis$brk <- 10^-(yaxis$brk)
  
  # keep scientific notation if small enough pvalues when changing to character
  yaxis$label <- format(yaxis$brk, trim=TRUE, digits=4,
                        scientific=ifelse(maxPval < 0.0001, TRUE, FALSE), 
                        drop0trailing=TRUE)
  
  # add threshold lines if thresholds are not NULL
  if(!is.null(fc_thresh)){
    ggFigure <- ggFigure + geom_vline(xintercept=fc_thresh, linetype="dotted")+
      geom_vline(xintercept=-fc_thresh, linetype="dotted")+
      annotate("text", x=fc_thresh+0.35, y=1,family=font_family, size=font_size/3,
               label=paste0("FC=", round(2^fc_thresh, digits=2)))
  }
  if(!is.null(pval_thresh)){
    ggFigure <- ggFigure + geom_hline(yintercept=pval_thresh, linetype="dotted")
    yaxis <- rbind(yaxis, c(pval_thresh, paste0('pval=',pval_thresh)))
  }
  if(!is.null(fdr_thresh)){
    ggFigure <- ggFigure + geom_hline(yintercept=fdr_pval, linetype="dotted")
    yaxis <- rbind(yaxis, c(fdr_pval, paste0('FDR=',fdr_thresh)))
  }
  
  # order yaxis in increasing value 
  yaxis$brk <- as.numeric(yaxis$brk)
  yaxis <- yaxis[order(yaxis$brk, decreasing=F),]
  
  # subset de to only contain genes to label on plot, either by user specified genes or top n_genes by pval
  if(!is.null(gene_list)){
    gene_labels <- subset(de, subset=Target.Name %in% gene_list)
  }else{
    # remove gene labels for genes below lowest pvalue cutoff
    gene_labels <- de[which(de$Pvalue < min(fdr_pval, pval_thresh, na.rm = T)),]
    
    # only label genes above fc_thresh if set by user else only look at pvalue
    if(!is.null(fc_thresh) & label_fc == FALSE){
      gene_labels <- gene_labels[which(abs(gene_labels$Log2) > fc_thresh),]
    }
    
    # only keep top # of genes by pvalue
    if (n_genes_per_dataset == TRUE){
      gene_labels <- do.call(rbind, lapply(split(gene_labels, gene_labels$dataset_name), function(x) head(x[order(x$Pvalue), ], n_genes)))
        } else{
      gene_labels <- gene_labels[head(order(gene_labels$Pvalue, decreasing=FALSE), n=n_genes),]  
    }
  }  
  interactive=ifelse(output_format=="html", TRUE, FALSE)
  if (interactive==FALSE){
  # add gene labels to ggplot
  ggFigure <- ggFigure + geom_text_repel(data=gene_labels, aes(x=Log2, y=Pvalue, label=Target.Name), 
                                         family=font_family, force=5, fontface="bold", min.segment.length=0.1,
                                         size=label_size)
  }
  # add y axis labels 
  ggFigure <- ggFigure + scale_y_continuous(trans=change_axis_revlog_trans(base=10), breaks=as.numeric(yaxis$brk),
                                            labels=yaxis$label)
  
  colnames(gene_labels) <- c("Target tag", "Target group memership/s", "Target Name", "Log2", 
                             "Pvalue", "Adjusted pvalue", "-log10 pvalue", "-log10 adjusted pvalue",
                             "FDR")
  
  volcanoPlot <- list(ggFigure, gene_labels)
  names(volcanoPlot) <- c("plot", "gene_labels")
  
  if(interactive){
    ggplotlyFigure <- ggplotly(ggFigure)
    volcanoPlot <- list(ggplotlyFigure, gene_labels)
    names(volcanoPlot) <- c("plot", "gene_labels")
  }else {
    volcanoPlot <- list(ggFigure, gene_labels)
    names(volcanoPlot) <- c("plot", "gene_labels")
  }
  return(volcanoPlot)
}

#' testAreColors
#'
#' checks if all colors in a vector are valid color names
#' @param colors color names 
#' @return TRUE/FALSE statement on valid color status
#' @export
testAreColors <- function(colors) {
  return(sapply(colors, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error=function(e) FALSE)
  }))
}

#' testIdenticalClass
#'
#' raises error if given object is not the assumed variable class
#' @param object object to determine variable class
#' @param object_name name of object for error message
#' @param class_name expected variable class
#' @return None, errors out if class is not expected
#' @export
testIdenticalClass <- function(object, object_name, class_name){
  expect_identical(object=class(object), expected=class_name, 
                   label=paste(object_name, "variable must be",
                               class_name,"\n You've supplied a", class(object)))
}

#' testVariableFormats
#'
#' check all user input variables for class and other validity markers
#' @param de data frame of DE results from DSPDA
#' @return None, errors out if input is not valid
#' @export
testVariableFormats <- function(de=de_results){
  
  # check user input de file for number of columns and column name matching
  expected_col_names <- c("Target.tag", "Target.group.membership.s", "Target.Name", "Log2", 
                          "Pvalue", "Adjusted.pvalue", "X.log10.pvalue", "X.log10.adjusted.pvalue", 
                          "source", "dataset_name", "FDR"   )
  
  expect_identical(object=ncol(de), expected=length(expected_col_names),
                   label="Number of columns in given volcano plot tab delimited file do not match expected. 
                   Make sure file is tab delimited")
  
  expect_identical(object=colnames(de), expected=expected_col_names, 
                   label="Column names in given volcano plot tab delimited file do not match expected.")
  
  # check that all target_groups have at least one gene if not NULL 
  if(!is.null(target_groups)){
    invisible(lapply(X=target_groups, FUN=function(x){
      genes_in_group <- grep(x=de$Target.group.membership.s, pattern=x)
      if(length(genes_in_group) == 0){
        fail(message=paste(x, "is not a valid probe group. Please check spelling or remove from target_groups before continuing"))
      }
    }))
  }
  
  ############################### USER DEFINED VARIABLE CLASS CHECKS ###############################
  numeric_variables <- c("n_genes", "pval_thresh", "fdr_thresh", "fc_thresh", 
                         "font_size", "label_size", "plot_width", "plot_height")
  
  character_variables <- c("plot_title", "gene_list", "target_groups", 
                           "negative_label", "positive_label")
  
  logical_variables <- c("show_legend", "label_fc")
  
  for(v in 1:length(numeric_variables)){
    if(!is.null(eval(parse(text=numeric_variables[v])))){
      testIdenticalClass(object=eval(parse(text=numeric_variables[v])), 
                         object_name=numeric_variables[v], 
                         class_name="numeric")
    }
    
  }
  
  for(v in 1:length(character_variables)){
    if(!is.null(eval(parse(text=character_variables[v])))){
      testIdenticalClass(object=eval(parse(text=character_variables[v])), 
                         object_name=character_variables[v], 
                         class_name="character")
    }
  }
  
  for(v in 1:length(logical_variables)){
    if(!is.null(eval(parse(text=logical_variables[v])))){
      testIdenticalClass(object=eval(parse(text=logical_variables[v])), 
                         object_name=logical_variables[v], 
                         class_name="logical")
    }
  }
  
  # check that either n_genes or gene_list is not NULL
  if(is.null(gene_list) & is.null(n_genes)){
    fail(message="Either n_genes or gene_list must not be NULL
         both n_genes and gene_list are currently NULL")
  }
  
  # check that either pval_thresh or fdr_thresh is not NULL
  if(is.null(pval_thresh) & is.null(fdr_thresh)){
    fail(message="Either fdr_thresh or pval_thresh must not be NULL
         both fdr_thresh and pval_thresh are currently NULL")
  }
  
  # check that gene_list only contains genes in de results
  if(!is.null(gene_list)){
    expect_true(object=all(gene_list %in% de$Target.Name), 
                label="At least one gene in gene_list does not match genes in volcano plot file results")
  }
  
  ################################# USER DEFINED VARIABLE  CHECKS ################################## 
  
  allowed_fonts <- c("serif", "sans", "mono")
  
  # check that font_family is an allowable font
  if(!font_family %in% allowed_fonts){
    fail(message=paste(font_family, "is not a valid font. Allowed fonts are",
                       paste(allowed_fonts, collapse=", ")))
  }
  
  # check that default_color is an allowable color
  if(!testAreColors(colors=default_color)){
    fail(message=paste(paste(default_color, collapse=", "), "is not a valid color"))
  }
  
  # check that fc_color is an allowable color
  if(!testAreColors(colors=fc_color)){
    fail(message=paste(paste(fc_color, collapse=", "), "is not a valid color"))
  }
  
  # check that color_options are all allowable colors
  if(!all(testAreColors(colors=color_options))){
    error_colors <- color_options[which(!testAreColors(colors=color_options))]
    fail(message=paste(paste(error_colors, collapse=", "), "is/are not valid color(s)"))
  }
  
  # check that enough colors are given for labeled target groups
  expect_gte(object=length(color_options), expected=max(length(target_groups), 2), 
             label=paste("Not enough colors were given./n", max(length(target_groups), 2), 
                         "expected, only", length(color_options), "given"))
  
  expected_output_format <- c("png", "jpg", "tiff", "svg", "pdf", "bmp", "html")
  
  if(!output_format %in% expected_output_format){
    fail(message=paste("Output format not in expected list of formats.\n", output_format, 
                       "given\n expected", paste(expected_output_format, collapse=", ")))
  }
}

#' change_axis_revlog_trans
#'
#' reverse log transform axis; used to return pvalue rather than -log10(pvalue) on yaxis
#' @param base base in which logs are computed
#' @return revlog_trans reverse log transformation
#' @export
change_axis_revlog_trans <- function(base=exp(1)){
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  revlog_trans <- trans_new(name=paste0("revlog-", base), 
                            transform=trans, 
                            inverse=inv, 
                            breaks=log_breaks(base=base),
                            domain=c(1e-100, Inf))
  
  return(revlog_trans)
}



custom_read_excel <- function(file) {
  # Read the entire Excel file once
  full_data <- read_excel(file, col_names = FALSE)
  
  # Locate the row and column of the #Dataset Name cell
  dataset_name_row <- 3
  dataset_name_col <- which(full_data[dataset_name_row, ] == "#Dataset Name")
  
  # Check if #Dataset Name was found
  if(length(dataset_name_col) == 0) {
    stop(paste("#Dataset Name cell not found in row 3 of file:", file))
  }
  
  # Get the value of the cell to the right of #Dataset Name
  dataset_name <- full_data[dataset_name_row, dataset_name_col + 1, drop = TRUE]
  
  # Skip the first 6 rows to get the rest of the data
  data <- full_data[-(1:6), ]
  
  # Rename the columns based on the first non-skipped row
  colnames(data) <- full_data[7, ]
  
  # Remove the row with column names
  data <- data[-1, ]
  
  # Add the source column
  data$source <- file
  data$dataset_name <- dataset_name
  
  # Formatting
  colnames(data) <- gsub(pattern="\\W", replacement=".", colnames(data))
  if(any(startsWith(colnames(data), prefix="."))){
    starts_with_num <- which(startsWith(colnames(data), prefix="."))
    colnames(data)[starts_with_num] <- paste0("X", colnames(data)[starts_with_num])
  }
  
  # calculate FDR and add it as a column
  data$FDR <- as.numeric(p.adjust(data$Pvalue, method="fdr"))
  
  return(data)
}

# main(obj1 = NULL, obj2 = NULL, obj3 = NULL, obj4 = NULL)  # Replace with actual objects if necessary


