### functions to create labeled volcano plots:
# note: this code assumes the input of volcano plot

# TESTING DATASET #

dataset <- readxl::read_xlsx(path="Q3 Norm.xlsx", sheet="TargetCountMatrix")
segmentAnnotations <- as.data.frame(readxl::read_excel(path="Q3 Norm.xlsx", 
                                                       sheet="SegmentProperties"))
targetAnnotations <- as.data.frame(readxl::read_excel(path="Q3 Norm.xlsx", 
                                                      sheet="TargetProperties"))
targetAnnotations$TargetGroup <- as.character(targetAnnotations$TargetGroup)
outputFolder <- "test_output/"

##########################################################
####      Arguments to be modified by user            ####
##########################################################

# volcano plot results (tab delimited):
de_results_filename <- "VOLCANO PLOT.txt"

# output format for volcano plot, 
#   options include: png, jpg, tiff, svg
output_format <- "png"

######################## LABELING ######################## 
# Volcano Plot Title
plot_title <- "MSI vs MSS"

# Negative (left) label for Fold Change from volcano plot
negative_label <- "MSS"

# Positive (right) label for Fold Change from volcano plot
positive_label <- "MSI"

# Number of genes to label
#   gene_list overrides this variable if set
n_genes <- 25

# Gene list to label. These genes will get labeled no matter 
#   where they are in the volcano plot. This variable is 
#   the default for labeling genes over n_genes
gene_list <- NULL #c("IL2RG", "GLUL", "SPIB", "C2")

####################### THRESHOLDS #######################
# P-value threshold, default threshold over fdr_thresh
pval_thresh <- NULL

# FDR threshold, must set pval_thresh to NULL to use
fdr_thresh <- 0.05

# Fold Change threshold 
fc_thresh <- 0.75

# Should threshold lines be added to plot
thresh_lines <- TRUE

######################### FONTS ##########################
# Font Size
font_size <- 12

# Font Family
font_family <- "mono"

####################### PLOT SIZE ########################
# Plot Width in inches
plot_width <- 7

# Plot Height in inches
plot_height <- 7

####################### COLORING #########################
# Specific target groups to color in the plot
#   If variable is set, genes are colored by given target 
#     group; else genes colored if they are above thresholds
target_groups <- NULL #c("Hemostasis", "DNA Repair")

# Color options for plot. Must have more colors than
#   number of target groups 
color_options <- c("#3A6CA1", "#FFD861", "#CF4244", "#47BAB4", 
                   "#474747", "#EB739E", "#318026", "#A66293", 
                   "#F28E2B", "#8F6954")

##########################################################
#### end of arguments. DO NOT CHANGE CODE BELOW HERE  ####
##########################################################

library(ggplot2)
library(ggrepel)
library(testthat)
library(stats)
library(stringr)

# main function called by DSP-DA:
main <- function(dataset, segmentAnnotations, 
                 targetAnnotations, outputFolder){
  
  # access volcano plot file:
  de_results <- read.table(de_results_filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # test for valid input variables
  test_variable_formats(de = de_results)
  
  # calculate FDR and add it as a column
  de_results$FDR <- p.adjust(de_results$Pvalue, method = "fdr")
  
  # create volcano plot
  volcano_plot(de = de_results)
} 

#' volcano_plot
#'
#' create volcano plot figure using user input variables and DE results from DSPDA
#' @param de data frame of DE results from DSPDA
#' @return Draws a volcano plot
#' @export
volcano_plot <- function(de){
  
  # determine highest x point to make volcano plot equal on both sides
  maxFC <- max(abs(de$Log2))
  
  # create basic volcano plot with correct formatting
  gp <- ggplot(de, aes(x = Log2, y = X.log10.pvalue))+
    geom_point(color = "grey60")+
    labs(y = "Significance, -log10(pval)", 
         x = paste(negative_label, "<-", "log2(FC)", "->", positive_label, sep = " "),
         title = plot_title)+
    theme(aspect.ratio = 1, text=element_text(size=font_size,  family=font_family))+
    scale_x_continuous(limits = c(-maxFC, maxFC))
    #scale_y_continuous(sec.axis = sec_axis(trans = ~ 10^-., name = "pval")) 
    # attempt to put pval on a second y axis, 
  
  # subset de to only include genes either in specified target groups or above pval/fdr threshold
  if(!is.null(target_groups)){
    probe_groups <- strsplit(de$Target.group.membership.s, split = ", ")
    de_subset_list <- lapply(X = target_groups, 
                            FUN = function(x){
                              #return de rows with specified target group named in column 
                              return(cbind(de[grep(x = de$Target.group.membership.s, pattern = x),], x))
                              })
    names(de_subset_list) <- target_groups
    
    gene_coloring <- as.data.frame(do.call(rbind, de_subset_list))
    names(gene_coloring)[names(gene_coloring) == "x"] <- "Target_coloring"
    gene_coloring$Target_coloring <- str_wrap(gene_coloring$Target_coloring, width = 45)
    
    color_label <- "Target Group\nMembership"
  }else{
    if(is.null(pval_thresh)){
      gene_coloring <- de[which(de$FDR < fdr_thresh),]
      label_thresh <- paste("FDR <", fdr_thresh)
    }else{
      gene_coloring <- de[which(de$Pvalue < pval_thresh),]
      label_thresh <- paste("pval <", pval_thresh)
    }
    
    gene_coloring$Target_coloring <- ifelse(test = gene_coloring$Log2 < 0, 
                                            yes = paste(label_thresh, negative_label), 
                                            no = paste(label_thresh, positive_label))
    
    color_label <- "Significance"
  }
  
  # add coloring to ggplot
  gp <- gp + geom_point(data = gene_coloring, aes(x = Log2, y = X.log10.pvalue, color = Target_coloring))+
    labs(color = color_label)+
    scale_color_manual(values = color_options)
  
  # add threshold lines if specified
  if(thresh_lines == TRUE){
    gp <- gp + geom_vline(xintercept = fc_thresh, linetype = "dotted")+
      geom_vline(xintercept = -fc_thresh, linetype = "dotted")
    
    if(is.null(pval_thresh)){
      # fine closest FDR value to fdr_thresh and use that pvalue to add y axis cutoff line
      gp <- gp + geom_hline(yintercept = -log10(mean(de$Pvalue[which(abs(de$FDR - fdr_thresh) == 
                                                                       min(abs(de$FDR - fdr_thresh)))])), 
                            linetype = "dotted")
    }else{
      gp <- gp + geom_hline(yintercept = -log10(pval_thresh), linetype = "dotted")
    }
  }
  
  # subset de to only contain genes to label on plot, either by user specified genes or top n_genes by pval
  if(!is.null(gene_list)){
    gene_labels <- subset(de, subset = Target.Name %in% gene_list)
  }else{
    gene_labels <- de[head(order(abs(de$X.log10.pvalue), decreasing = TRUE), n = n_genes),]
  }  
  
  # add gene labels to ggplot
  gp <- gp + geom_text_repel(data = gene_labels, aes(label=Target.Name), force = 3.5)
  
  # save ggplot 
  ggsave(filename = paste0(outputFolder, plot_title, "_volcano_plot.", output_format), 
         plot = gp, device = output_format, width = plot_width, height = plot_height)
  
}

#' areColors
#'
#' checks if all colors in a vector are valid color names
#' @param x color names 
#' @return TRUE/FALSE statement on valid color status
#' @export
areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  })
}

#' identical_class
#'
#' raises error if given object is not the assumed variable class
#' @param object object to determine variable class
#' @param object_name name of object for error message
#' @param class_name expected variable class
#' @return None, errors out if class is not expected
#' @export
identical_class <- function(object, object_name, 
                            class_name){
  expect_identical(object = class(object), expected = class_name, 
                   label = paste(object_name, "variable must be", 
                                 class_name,"\n You've supplied a", class(object)))
}

#' test_variable_formats
#'
#' check all user input variables for class and other validity markers
#' @param de data frame of DE results from DSPDA
#' @return None, errors out if input is not valid
#' @export
test_variable_formats <- function(de = de_results){
  
  # check user input de file for number of columns and column name matching
  expected_col_names <- c("Target.tag", "Target.group.membership.s", "Target.Name", "Log2", 
                          "Pvalue", "Adjusted.pvalue", "X.log10.pvalue", "X.log10.adjusted.pvalue")
  
  expect_identical(object = ncol(de), expected = length(expected_col_names),
                   label = "Number of columns in given volcano plot tab delimited file do not match expected. Make sure file is tab delimited")
  
  expect_identical(object = colnames(de), expected = expected_col_names, 
                   label = "Column names in given volcano plot tab delimited file do not match expected.")
  
  expected_output_format <- c("png", "jpg", "tiff", "svg")
  
  if(!output_format %in% expected_output_format){
    fail(message = paste("Output format not in expected list of formats.\n", output_format, 
                         "given\n expected", paste(expected_output_format, collapse = ", ")))
  }
  
  ############################################ LABELING ############################################ 
  # check that n_genes is a numeric or NULL
  if(!is.null(plot_title)){
    identical_class(object = plot_title, object_name = "plot_title", class_name = "character")
  }
  
  # check that either n_genes or gene_list is not NULL
  if(is.null(gene_list) & is.null(n_genes)){
    fail(message = "Either n_genes or gene_list must not be NULL \n both n_genes and gene_list are currently NULL")
  }
  
  # check that n_genes is a numeric or NULL
  if(!is.null(n_genes)){
    identical_class(object = n_genes, object_name = "n_genes", class_name = "numeric")
  }
  
  # check that gene_list is a character list or NULL
  if(!is.null(gene_list)){
    identical_class(object = gene_list, object_name = "gene_list", class_name = "character")
    expect_true(object = all(gene_list %in% de$Target.Name), 
                label = "At least one gene in gene_list does not match genes in volcano plot file results")
  }
  
  ########################################### THRESHOLDS ###########################################
  # check that either pval_thresh or fdr_thresh is not NULL
  if(is.null(pval_thresh) & is.null(fdr_thresh)){
    fail(message = "Either fdr_thresh or pval_thresh must not be NULL \n both fdr_thresh and pval_thresh are currently NULL")
  }
  
  # check that pval_thresh is a numeric or NULL
  if(!is.null(pval_thresh)){
    identical_class(object = pval_thresh, object_name = "pval_thresh", class_name = "numeric")
  }
  
  # check that fdr_thresh is a numeric or NULL
  if(!is.null(fdr_thresh)){
    identical_class(object = fdr_thresh, object_name = "fdr_thresh", class_name = "numeric")
  }
  
  # check that fc_thresh is a numeric or NULL
  if(!is.null(fc_thresh)){
    identical_class(object = fc_thresh, object_name = "fc_thresh", class_name = "numeric")
  }
  
  # check that thresh_lines is a logical or NULL
  if(!is.null(thresh_lines)){
    identical_class(object = thresh_lines, object_name = "thresh_lines", class_name = "logical")
  }
  
  ############################################# FONTS ##############################################
  # check that font_size is a numeric or NULL
  if(!is.null(font_size)){
    identical_class(object = font_size, object_name = "font_size", class_name = "numeric")
  }
  
  # check that font_family is an allowable font
  if(!font_family %in% names(windowsFonts())){
    fail(message = paste("Given font_family is not a valid font. Allowed fonts are", 
                         paste(names(windowsFonts()), collapse = ", ")))
  }
  
  ########################################### PLOT SIZE ############################################
  # check that plot_width is a numeric or NULL
  if(!is.null(plot_width)){
    identical_class(object = plot_width, object_name = "plot_width", class_name = "numeric")
  }
  
  # check that plot_width is a numeric or NULL
  if(!is.null(plot_height)){
    identical_class(object = plot_height, object_name = "plot_height", class_name = "numeric")
  }
  
  ########################################### COLORING #############################################
  # check that target_groups is a character or NULL
  if(!is.null(target_groups)){
    identical_class(object = target_groups, object_name = "target_groups", class_name = "character")
  }
  
  # check that color_options are all allowable colors
  if(!all(areColors(x = color_options))){
    error_colors <- color_options[which(!areColors(x = color_options))]
    fail(message = paste(paste(error_colors, collapse = ", "), "is/are not valid color(s)"))
  }
  
  # check that enough colors are given for labeled target groups
  expect_gte(object = length(color_options), expected = max(length(target_groups), 2), 
             label = paste("Not enough colors were given./n", max(length(target_groups), 2), 
                           "expected, only", length(color_options), "given"))
}

main(dataset = dataset, 
     segmentAnnotations = segmentAnnotations, 
     targetAnnotations = targetAnnotations, 
     outputFolder = outputFolder)
