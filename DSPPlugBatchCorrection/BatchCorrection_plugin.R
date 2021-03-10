# DSP-NGS Batch Correction #

# Runs PCA on samples
# Perform Batch Correction (BC)
# Runs PCA on BC samples
# Returns a BC data.frame


### ###############
### User Inputs ###
### ###############

## Statistical parameters
# The batching factor (required)
# batching_factor <- "SlideName"
# factors of biological interest
# to check for collinearity (or NULL)
# factors_of_interest <- c("SegmentName", "Tumor") # NULL
# The tolerance (1 or greater)
# vif_threshold <- 5

## Plot parameters
# color plots by one of these three:
# - "batching_factor"
# - "factors_of_interest"
# - NULL
color_by = "batching_factor"
# shape plots by one of these three:
# - "batching_factor"
# - "factors_of_interest"
# - NULL
shape_by = "factors_of_interest"
# size points by "PC3" or NULL
size_by = NULL
# font family and axes and label
# size
plot_font = list(
  family = "sans",
  size = 15
  )

### ###################
### End User Inputs ###
### ###################

### ####### ###
### License ###
### ####### ###
# MIT License
# Copyright 2020 Nanostring Technologies, Inc.
# Permission is hereby granted, free of charge, 
# to any person obtaining a copy of this software
# and associated documentation files (the “Software”),
# to deal in the Software without restriction, 
# including without limitation the rights to use, 
# copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit 
# persons to whom the Software is furnished to do so,
# subject to the following conditions:
# The above copyright notice and this permission notice
# shall be included in all copies or substantial portions
# of the Software.
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY
# OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO 
# EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
# AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
# USE OR OTHER DEALINGS IN THE SOFTWARE.
# Contact us: 
#   NanoString Technologies, Inc.
#   530 Fairview Avenue N
#   Seattle, WA 98109
# Tel: (888) 358-6266
### ########### ###
### End License ###
### ########### ###




### ############# ###
### Main Function ###
### ############# ###

# Main function used by DSP-DA
main <- function(dataset, segmentAnnotations, 
                 targetAnnotations, outputFolder) {

  # Load packages
  load_packages()
  
  # Check user inputs
  check_user_input(batching_factor, factors_of_interest, color_by,
                   shape_by, size_by, plot_font)

  # Run batch correction
  bc <- run_batch_correction(
    df=dataset,
    annots=segmentAnnotations,
    batching_factor=batching_factor
  )
  
  # Run the main QC function
  run_qc(
    df=dataset,
    bdf=bc,
    annots=segmentAnnotations,
    batching_factor=batching_factor,
    factors_of_interest=factors_of_interest
  )
  
  # Send bc data to disk
  write.table(bc, 
    file = file.path(outputFolder, 
         paste0("dataset.tsv"),
         fsep = .Platform$file.sep), 
    sep="\t", col.names=TRUE,
    row.names=TRUE, quote=FALSE)
              
  return(bc)
  
}

### ################ ###
### Helper functions ###
### ################ ###

# Loads packages
load_packages <- function(){

  # List of packages for session
  .packages = c(
    "lme4", # mixed models
    "parallel", # multiple processors
    "rlang", # tunneling data-variables
    "openxlsx" # for spreadsheets
  )

  # Load packages into session 
  null <- lapply(.packages, require, 
                 character.only=TRUE)
}

# Checks user input
# @param batching_factor is character provided by user
# @param factors_of_interest is NULL or a vector specifying
#        the column names of factors to compare with 
#        batch.
check_user_input <- function(batching_factor, factors_of_interest, color_by,
                             shape_by, size_by, plot_font){
  
  pass <- TRUE # passing flag that cascades down the function
  msg <- c() # for an error message to give user
  
  # Make sure the batching_factor 
  # is in segmentAnnotations <error>.
  if(!(bactching_factor %in% colnames(segmentAnnotations))){
    pass <- FALSE # hard fail.
    msg <- c(
      msg, 
      paste0("The batching factor provided, ",
      batching_factor, 
      " was not found. Please check spelling.\n"
      ))
  }
  # Check if batching_factor, if present, is numeric
  # and warn user if so.
  if(pass & inherits(
      x=segmentAnnotations[,bactching_factor], what=
      "numeric")
    ){
      # still could be meaningful but warn user.
      msg <- c(
        msg, 
        paste0("Warning! The batching factor, ", 
               bactching_factor, 
               ", is numeric. Will attempt ", 
               "to coerce to a factor.\n")
      )
  }

  # if factors of interest are specified, 
  # go through each one and verify there's
  # a column name in segmentAnnottaions.
  if(!is.null(factors_of_interest)){
    for(i in factors_of_interest){
      if(!(i %in% colnames(segmentAnnotations))){
        pass <- FALSE # hard fail.
        msg <- c(
          msg, 
          paste0("The factor \'",
                 i, 
                 "\' was not found. Please check spelling or set factors_of_interest to NULL.\n"
          ))
      }
    }
  }

  # Check that color_by is a valid
  if(!is.null(color_by)){
    if(!color_by %in% c("batching_factor", "factors_of_interest")){
      pass <- FALSE
      msg <- c(
        msg, 
        paste0("The color_by parameter given, \'",
               color_by, 
               "\', is not NULL or \'batching_factor'", 
               " or \'factors_of_interest\'. Please specify one of these three.\n"
        ))
    }
  }
  
  # Check that shape_by is a valid
  if(!is.null(shape_by)){
    if(!shape_by %in% c("batching_factor", "factors_of_interest")){
      pass <- FALSE
      msg <- c(
        msg, 
        paste0("The shape_by parameter given, \'",
               shape_by, 
               "\', is not NULL or \'batching_factor'", 
               " or \'factors_of_interest\'. Please specify one of these three.\n"
        ))
    }
  }
  
  # Check that size_by is a valid
  if(!is.null(size_by)){
    if(!size_by %in% c("PC3")){
      pass <- FALSE
      msg <- c(
        msg, 
        paste0("The size_by parameter given, \'",
               size_by, 
               "\', is not NULL or \'PC3\'.", 
               " Please specify NULL or \'PC3\'.\n"
        ))
    }
  }

  ## Check that the plot_font object is valid
  # It must be a list of length two with names
  # 'family' and 'size'. 'family must be a valid 
  # font family and size must be numeric and >0.
  if(!inherits(plot_font, "list")){
    pass <- FALSE
    msg <- c(
      msg, 
      paste0("The plot_font parameter must be a list.\n"
      ))
  } else if(length(plot_font)!=2){
    # error: it must be of length two: family, size
    pass <- FALSE
    msg <- c(
      msg, 
      paste0("plot_font must be a list of length two: family, size.\n"
      ))
  } else if(!all(names(plot_font) %in% c("family", "size"))){
    # error: names need to be family and size only.
    pass <- FALSE
    msg <- c(
      msg, 
      paste0("plot_font must be a list with names \'family\' and \'size\'.\n"
      ))
  }
  
  # size points by "PC3" or NULL
  size_by = NULL
  # font family and axes and label
  # size
  plot_font = list(
    family = "sans",
    size = 15
  )
  
  # If there were any messages, send to disk.
  if(length(msg)>0){
    write(msg, 
      file = file.path(outputFolder, 
            paste0("log.txt"),
            fsep = .Platform$file.sep))
  }
 
  # if user input did not pass, 
  # fail and send user the log.
  if(!pass){
    stop(paste0(
      "Errors found. ", 
      "Please see log.txt for message(s)."))
  } else {
    # return nothing if passed.
    return(NULL)
  }
}


# Main batch correction function
# @param df is dataset
# @param annots is the annotation data.frame
# @param batching_factor is the batching factor
#        specified by the user.
# @param n_process is the number of processors
#        to use. Defaults to 2 unless specified.
# @return a data.frame of the same dimensions
#        as dataset with batch corrected 
#        results.
run_batch_correction <- function(
  df, annots, batching_factor, 
  n_process=2){
  
  # advanced: 
  # set n_process to n_processors if
  # in the global environment.
  if("n_processors" %in% ls()){
    n_process <- n_processors
  }
  
  # Create clusters and export
  cl <- makeCluster(n_process)
  clusterExport(cl=cl, 
    varlist=c("df", "annots", 
      'batching_factor'),
    envir=environment())
  
  # Go through each feature by index.
  resid_list <- parLapply(
    cl, 
    1:nrow(df), function(i){
    require(lme4) # loaded in cluster
    require(tibble)
    
    # Pivot data for feature i
    # and combine with annot's batching_factor
    model_df <- t(log2(df[i,]))
    colnames(model_df)[1] <- "feature_i" # readability
    model_df <- add_column(
      as_tibble(model_df),
      "segmentID"=row.names(model_df), 
      .before=1)
    model_df <- base::merge(
      model_df, 
      annots[,c("segmentID", 
        batching_factor)],
      by="segmentID")
    
    # Run model
    the_formula = paste0(
      'feature_i ~ -1 + (1|', 
      batching_factor,
      ')')
    mod <- lme4::lmer(
      as.formula(the_formula), 
      data=model_df)
    
    # Get residuals and return
    out <- matrix(
      as.numeric(resid(mod)), nrow=1)
    colnames(out) <- model_df$segmentID
    row.names(out) <- row.names(df)[i]
    return(out)
  })
  
  # Stop cluster and rbind
  stopCluster(cl)
  resid_df <- do.call(rbind, resid_list)
  
  # Make sure rows and columns 
  # are aligned compared to 
  # df, make into data.frame, and return.
  resid_df <- resid_df[rownames(df),]
  resid_df <- resid_df[,colnames(df)]
  resid_df <- as.data.frame(resid_df)
  return(resid_df)
}

#' @title run_qc
#' @description: main qc wrapper function. Calls lower-level qc functions.
#' @param: df = dataset
#' @param: bdf = batch corrected data from run_batch_correction()
#' @param: annots = segmentAnnotations
#' @param: batching_factor = the factor used for batch effects
#' @param: factors_of_interest either NULL of a vector of 
#'         possibly collinear independent variables for which
#'         QC comparison with batching_factor are to be made.
#' @return: an excel workbook that comprises the QC metrics.

run_qc <- function(df, bdf, annots, batching_factor, factors_of_interest){
  require(openxlsx)
  # Instantiate spreadsheet
  wb <- createWorkbook(title = paste("Batch Correction QC"), 
                       creator = "NanoString DSP Plugin")
  
  ## Non-batch-corrected (nbc) data
  # PCA data and write results
  pca_nbc <- compute_pca(exp_data=df, log2_transform=TRUE)
  wb <- write_pca_results(the_wb=wb, pca_data=pca_nbc, the_name="Before BC")
  # Plot data
  wb <- plot_pca_data(the_wb=wb, pca_data=pca_nbc, annots=annots, the_name="Before BC", 
                      batching_factor=batching_factor, 
                      factors_of_interest=factors_of_interest)
  
  
  # PCA on batch-corrected (bc) data and write results
  pca_bc <- compute_pca(exp_data=df, log2_transform=FALSE)
  wb <- write_pca_results(the_wb=wb, pca_data=pca_bc, the_name="After BC")
  
  
  # Save workbook to disk
  saveWorkbook(wb = wb,
    file = file.path(outputFolder,
      paste0("batch_correction_QC.xlsx"),
      fsep = .Platform$file.sep), overwrite = TRUE)
  
  # return the workbook
  return(wb)
  
}

#' @title compute_pca
#' @description lower-level function to run PCA.
#' @param exp_data a data.frame of expression data with column names equal to samples.
#' @param log2_transform logical whether to log2 transform the data (default) or not. 
#' @return and object of class prcomp giving the PCA results.
compute_pca <- function(exp_data, log2_transform=TRUE){
  if(log2_transform){
    dr_data <- stats::prcomp(t(log2(exp_data)), center=TRUE, scale. = TRUE)
  } else {
    dr_data <- stats::prcomp(t(exp_data), center=TRUE, scale. = TRUE)
  }
  return(dr_data)
}

#' @title write_pca_results
#' @description writes PCA results given PCA data and name
#' @param the_wb an object of class Workbook 
#' @param pca_data the PCA data
#' @param the_name a char string given the name for prepending to sheet name
#' @details This function adds PCs (one for each sample), PC loadings
#' @return an object of class Workbook with the new sheets
write_pca_results <- function(the_wb, pca_data, the_name){
  
  # Add All PCs for each sample
  pcs_name <- paste0(the_name, " - PCs (All)")
  addWorksheet(the_wb, pcs_name)
  writeData(wb = the_wb,
            sheet = pcs_name,
            x = pca_data$x,
            colNames = TRUE, rowNames = TRUE)
  setColWidths(wb = the_wb, sheet = pcs_name, cols = 1, widths = "auto")

    # Add PC loadings
  loadings_name <- paste0(the_name, " - PCs Loadings")
  addWorksheet(the_wb, loadings_name)
  writeData(wb = the_wb,
            sheet = loadings_name,
            x = pca_data$rotation[order(abs(pca_data$rotation[, 1]), decreasing = TRUE), ],
            colNames = TRUE, rowNames = TRUE)
  setColWidths(wb = the_wb, sheet = loadings_name, cols = 1, widths = "auto")
  
  # Add variance estimates
  var_name <- paste0(the_name, " - Var. Est.")
  addWorksheet(the_wb, var_name)
  writeData(wb = the_wb,
            sheet = var_name, 
            x = t(summary(pca_data)$importance),
            colNames = TRUE, rowNames = TRUE)
  setColWidths(wb = the_wb, sheet = var_name, cols = 1:4, widths = "auto")
  
  return(the_wb)
}


#' @title plot_pca_data
#' @description plots PCA results and overlays batching_factor and factors_of_interest, if applicable
#' @param the_wb an object of class Workbook 
#' @param pca_data the PCA data
#' @param annots the segmentAnnotations
#' @param the_name a char string given the name for prepending to sheet name
#' @details Uses interaction in ggplot2 to color the points
#' @return an object of class Workbook
plot_pca_data <- function(the_wb, pca_data, annots, the_name, 
                    batching_factor, factors_of_interest){
  
  # to do. build out function to plot PC values and compute VIF values
  # goal: write the data for making the figure to a sheet, insert the plot, 
  #       and add a sheet summarizing the VIFs.
  
  # Pivot and Merge first 3 values to annotations' batch
  #   factor and factors of interest (if applicable)
  plot_df <- pca_data$x[,1:3]
  plot_df <- add_column(
    as_tibble(plot_df),
    "segmentID"=row.names(plot_df), 
    .before=1)
  plot_df <- base::merge(
    annots[,c("segmentID", batching_factor, eval(factors_of_interest))],
    plot_df,
    by="segmentID"
  )
  return(dim(plot_df))
  
  
  
  
}

print(p1)
insertPlot(wb=the_wb, pcs_name, width = 5, height = 3.5, fileType = "png", units = "in")



### #####
### Notes
### #####

# dataset has features (rows) by (samples)
# columns. These features and columns
# are referenced by alphanumeric strings
# and need may need to be translated.

# Column names in the `dataset` object 
# correspond to rows in the `SegmentID`
# column of `segmentAnnotations` object.

# Features. The targetAnnotations object
# links TargetGUID with TargetName.

