### functions to create labeled volcano plots:
# note: this code assumes the input of VOLCANO PLOT 
#       without header and tab delimited


##########################################################
####      Arguments to be modified by user            ####
##########################################################

# volcano plot results (tab delimited):
de_results_filename <- "VOLCANO PLOT.txt"

# output format for volcano plot, 
#   options include: png, jpg, tiff, svg, pdf, & bmp
output_format <- "jpg"

######################## LABELING ######################## 
# Volcano Plot Title
plot_title <- "AOI Surface Area Comparisons"

# Negative (left) label for Fold Change from volcano plot
negative_label <- "50 um Circle"

# Positive (right) label for Fold Change from volcano plot
positive_label <- "200 um Circle"

# Show color legend on figure
show_legend <- TRUE

# Number of genes to label
#   gene_list overrides this variable if set
n_genes <- 25

# Gene list to label. These genes will get labeled no matter 
#   where they are in the volcano plot. This variable is 
#   the default for labeling genes over n_genes
gene_list <- NULL #c("IL2RG", "GLUL", "SPIB", "C2")

####################### THRESHOLDS #######################
# P-value threshold, must set fdr_thresh to NULL to use
pval_thresh <- 0.05

# FDR threshold, default threshold over pval_thresh
fdr_thresh <- 0.01

# Fold Change threshold, if set coloring options will change
fc_thresh <- 0.5

# Should genes below fc_thresh be labeled
label_fc <- FALSE

######################### FONTS ##########################
# Font Size
font_size <- 8

# Label Font Size
label_size <- 2

# Font Family
#   options include: serif, sans, mono
font_family <- "sans"

####################### PLOT SIZE ########################
# Plot Width in inches
plot_width <- 6

# Plot Height in inches
plot_height <- 4

####################### COLORING #########################
# Default point color for points not in target groups or
#   above thresholds
default_color <- "grey65"

# Color of points below fc_thresh but above pval and/or 
#   fdr thresholds. change to the same as default_color
#   if you don't want these targets called out
fc_color <- "grey30"

# Specific target groups to color in the plot
#   If variable is set, genes are colored by given target 
#     group; else genes colored if they are above thresholds
target_groups <- NULL #c("Hemostasis", "DNA Repair")

# Color options for plot. Must have more colors than
#   number of target groups 
color_options <-  c("#3A6CA1", "#FFD861", "#CF4244", 
                    "#47BAB4", "#474747", "#EB739E", 
                    "#318026", "#A66293", "#F28E2B", 
                    "#8F6954")

##########################################################
#### end of arguments. DO NOT CHANGE CODE BELOW HERE  ####
##########################################################

library(ggplot2)
library(ggrepel)
library(testthat)
library(scales)
library(stats)
library(stringr)

# main function called by DSP-DA:
main <- function(dataset, segmentAnnotations, targetAnnotations, outputFolder){
  
  if(!file.exists(de_results_filename)){
    fail(message="Given volcano plot results file does not exist. 
         Please check file name and if the volcano plot file was uploaded to DSPDA")
  }
  
  # access volcano plot file:
  de_results <- read.table(de_results_filename, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
  #check if header is removed from volcano plot file
  if(!"target.name" %in% tolower(colnames(de_results))){
    
    #if header isn't removed in de_results_file then remove header and set column names
    header <- grep(pattern="target name", x=tolower(de_results))
    
    if(length(header) > 0){
      colnames(de_results) <- de_results[header, ]
      de_results <- de_results[-c(1:header), ]
      
      colnames(de_results) <- gsub(pattern="\\W", replacement=".", colnames(de_results))
      if(any(startsWith(colnames(de_results), prefix="."))){
        starts_with_num <- which(startsWith(colnames(de_results), prefix="."))
        colnames(de_results)[starts_with_num] <- paste0("X", colnames(de_results)[starts_with_num])
      }
    }else{
      fail(message="Too many rows were removed from VOLCANO PLOT.xlsx before running script. 
           The row with Target Name in the first column must be kept.")
    }
    
    de_results$Log2 <- as.numeric(de_results$Log2)
    de_results$Pvalue <- as.numeric(de_results$Pvalue)
  }
  
  #ensure consistent capitalization with Target.Name column name
  colnames(de_results)[which(tolower(colnames(de_results)) == "target.name")] <- "Target.Name"
  
  # test for valid input variables
  testVariableFormats(de=de_results)
  
  # calculate FDR and add it as a column
  de_results$FDR <- p.adjust(de_results$Pvalue, method="fdr")
  
  # create volcano plot
  volcanoPlot_results <- plotVolcano(de=de_results)
  
  if(output_format == "svg"){
    svg(filename=paste0(outputFolder, "/volcano_plot_", plot_title, ".", output_format), 
        width=plot_width, height=plot_height)
    print(volcanoPlot_results$plot)
    dev.off()
  }else{
    ggsave(filename=paste0("volcano_plot_", plot_title, ".", output_format),
           plot=volcanoPlot_results$plot,
           device=output_format,
           path=outputFolder,
           height=plot_height,
           width=plot_width)
  }
  
  write.table(x=volcanoPlot_results$gene_labels, 
              file=file.path(outputFolder,
                               paste0("labeled_genes_", plot_title, ".txt"), 
                               fsep=.Platform$file.sep),
              sep="\t", quote=FALSE, row.names=FALSE)
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
  ggFigure <- ggplot(de, aes(x=Log2, y=Pvalue))+
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
  
  # subset de to only include genes either in specified target groups or above pval/fdr threshold
  if(!is.null(target_groups)){
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
  
  color_options <- c(color_options, default_color)
  names(color_options)[length(color_options)] <- "Not Specified"
  
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
    gene_labels <- gene_labels[head(order(gene_labels$Pvalue, decreasing=FALSE), n=n_genes),]
  }  
  
  # add gene labels to ggplot
  ggFigure <- ggFigure + geom_text_repel(data=gene_labels, aes(x=Log2, y=Pvalue, label=Target.Name), 
                             family=font_family, force=5, fontface="bold", min.segment.length=0.1,
                             size=label_size)
  
  # add y axis labels 
  ggFigure <- ggFigure + scale_y_continuous(trans=change_axis_revlog_trans(base=10), breaks=as.numeric(yaxis$brk),
                          labels=yaxis$label)

  colnames(gene_labels) <- c("Target tag", "Target group memership/s", "Target Name", "Log2", 
                             "Pvalue", "Adjusted pvalue", "-log10 pvalue", "-log10 adjusted pvalue",
                             "FDR")
  
  volcanoPlot <- list(ggFigure, gene_labels)
  names(volcanoPlot) <- c("plot", "gene_labels")
  
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
                          "Pvalue", "Adjusted.pvalue", "X.log10.pvalue", "X.log10.adjusted.pvalue")
  
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
  
  expected_output_format <- c("png", "jpg", "tiff", "svg", "pdf", "bmp")
  
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

