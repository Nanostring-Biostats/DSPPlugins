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
# cor_threshold <- 0.01

## Plot parameters
# color plots by one of these three:
# - "batching_factor"
# - "factors_of_interest"
# - NULL
# color_by = "batching_factor"
# shape plots by one of these three:
# - "batching_factor"
# - "factors_of_interest"
# - NULL
# shape_by = "factors_of_interest"
# size points by "PC3" or NULL
# size_by = NULL
# font family and axes and label
# size
# plot_font = list(
#   family = "sans",
#   size = 15
#   )

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
    "openxlsx", # for spreadsheets
    "ggplot2", # for plotting
    "cowplot" # side-by-side plots
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
  if(!(batching_factor %in% colnames(segmentAnnotations))){
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
      x=segmentAnnotations[,batching_factor], what=
      "numeric")
    ){
      # still could be meaningful but warn user.
      msg <- c(
        msg, 
        paste0("Warning! The batching factor, ", 
               batching_factor, 
               ", is numeric.\n")
      )
  }
  # Make sure there are at least two unique observation in batching_factor.
  if(pass & length(unique(segmentAnnotations[,batching_factor]))<2){
    # still could be meaningful but warn user.
    pass <- FALSE
    msg <- c(
      msg, 
      paste0("Error! The batching factor, ", 
             batching_factor, 
             ", has only ", 
             length(unique(segmentAnnotations[,batching_factor])), 
             " unique observations. There needs to be two or more levels.\n")
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
  
  # Make sure that shape_by has 6 or fewer values (default for ggplot)
  if(!is.null(shape_by)){
    if(length(factors_of_interest)==1){
      the_values <- unique(segmentAnnotations[,factors_of_interest])
    } else {
      the_values <- unique(apply(segmentAnnotations[,factors_of_interest], 1, paste, collapse="."))
    }
    if(length(the_values)>6){
      pass <- FALSE
      msg <- c(
        msg, 
        paste0("The shape_by parameter specified has ",
               length(the_values), 
               " combinations of values which is more than the maximum allowed.\n")
        )
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
  fonts_allowed <- c("sans") # acceptable fonts
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
  } else if(!plot_font$family %in% fonts_allowed){
    # error: the font specified is not in the tested fonts
    pass <- FALSE
    msg <- c(
      msg, 
      paste0("plot font family must be in the following list: ",
      paste(fonts_allowed, collapse=", "), 
      ".\n"
      ))
  } else if(!inherits(plot_font$size, "numeric")){
    # error: font size needs to be numeric
    pass <- FALSE
    msg <- c(
      msg, 
      paste0("plot font size must be a number >0.\n"
      ))   
  } else if(!plot_font$size > 0){
    # error: font size needs to be numeric
    pass <- FALSE
    msg <- c(
      msg, 
      paste0("plot font size must be larger than 0.\n"
      ))   
  }
  
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
#' @return: NULL. Results sent to disk.
run_qc <- function(df, bdf, annots, batching_factor, factors_of_interest){
  
  # Instantiate spreadsheet
  wb <- createWorkbook(title = paste("Batch Correction QC"), 
                       creator = "NanoString DSP Plugin")
  
  # compare independent factors
  ind_list <- compare_independent_factors(
    the_wb = wb,
    annots=annots, 
    batching_factor=batching_factor,
    factors_of_interest=factors_of_interest)
  wb <- ind_list[[1]]
  factors_of_interest <- ind_list[[2]] # removed 100% confounded ind. variables.
  
  ## Non-batch-corrected (nbc) data
  # PCA data and write results
  pca_nbc <- compute_pca(exp_data=df, log2_transform=TRUE)
  wb <- write_pca_results(the_wb=wb, pca_data=pca_nbc, the_name="Before BC")
  # Plot data
  nbc_plot_list <- plot_pca_data(the_wb=wb, pca_data=pca_nbc, annots=annots, the_name="Before BC", 
                      batching_factor=batching_factor, 
                      factors_of_interest=factors_of_interest)
  wb <- nbc_plot_list[[1]]
  nbc_pca_plot <- nbc_plot_list[[2]]
  
  ## Batch corrected (bc) data
  # PCA on batch-corrected (bc) data and write results
  pca_bc <- compute_pca(exp_data=df, log2_transform=FALSE)
  wb <- write_pca_results(the_wb=wb, pca_data=pca_bc, the_name="After BC")
  # Plot data
  bc_plot_list <- plot_pca_data(the_wb=wb, pca_data=pca_bc, annots=annots, the_name="After BC", 
                      batching_factor=batching_factor, 
                      factors_of_interest=factors_of_interest)
  wb <- bc_plot_list[[1]]
  bc_pca_plot <- bc_plot_list[[2]]
  
  # Compare nbc and bc datasets
  wb <- compare_batch_correction(
    the_wb=wb,
    annots=annots,
    pca_nbc=pca_nbc,
    pca_bc=pca_bc,
    batching_factor=batching_factor,
    factors_of_interest=factors_of_interest,
    nbc_pca_plot=nbc_pca_plot,
    bc_pca_plot=bc_pca_plot
  )
  
  # Save workbook to disk
  saveWorkbook(wb = wb,
    file = file.path(outputFolder,
      paste0("batch_correction_QC.xlsx"),
      fsep = .Platform$file.sep), overwrite = TRUE)
  
  # Sends all QC to disk
  return(NULL)
  
}


#' @title compare_independent_factors
#' @description Looks at independent factors and generates report
#' @param the_wb the Workbook
#' @param annots the annotations
#' @param batching_factor the batching_factor
#' @parma factors_of_interest the factors_of_interest (NULL or 1 or more)
compare_independent_factors <- function(the_wb, annots,batching_factor,factors_of_interest){
  
  sheet_name <- paste0("Compare Ind. Vars.")
  addWorksheet(the_wb, sheet_name)
  
  # Links the factors with their associates level names in the model
  # and gets the correlations from the alias function
  model_list <- tranlate_model_matrix(
    df=annots[,which(colnames(annots) %in% c(batching_factor, factors_of_interest))], 
    batching_factor=batching_factor, 
    factors_of_interest=factors_of_interest)
  names_link <- model_list[[1]]
  cors <- model_list[[2]]
  
  row_offset = 1 # Cascades down the logic to place results in Workbook.
  
  # Check if any factors are completely confounded with batching_factor. 
  # If any are, tell the user and drop those factors of interest.
  if("Complete" %in% attributes(cors)$names){
    complete_cores <- cors$Complete
    batching_levels <- as.character(names_link[which(names_link$the_factors %in% batching_factor),'the_levels_names'])
    
    factors_to_remove <- c()
    
    for(i in 1:nrow(complete_cores)){
      for(j in 1:ncol(complete_cores)){
        if(complete_cores[i,j]>0){
          # This combination is confounded.
          # Check if it affects a batching level.
          row_i_factor <- names_link$the_factors[which(as.character(names_link$the_levels_names)==rownames(complete_cores)[i])]
          col_j_factor <- names_link$the_factors[which(as.character(names_link$the_levels_names)==colnames(complete_cores)[j])]
          if(row_i_factor %in% batching_factor){
            # remove the associated column's factor
            factors_to_remove <- c(factors_to_remove, col_j_factor)
          }            
          if(col_j_factor %in% batching_factor){
            # remove the associated rows's factor
            factors_to_remove <- c(factors_to_remove, row_i_factor)
          } 
          
        }
      }
    }
    # tell user and remove factors:
    if(length(factors_to_remove)>0){
      msg <- paste0("There were ", 
                    length(factors_to_remove), " factor(s) (", 
                    paste(factors_to_remove, collapse=", "),
                    ") that were perfectly correlated with the batching factor \'",
                    batching_factor, 
                    "\'. Removing from list.")
      writeData(wb = the_wb, sheet = sheet_name, x = msg, startRow=row_offset,
                colNames = FALSE, rowNames = FALSE)
      row_offset <- row_offset+3
      factors_of_interest <- setdiff(factors_of_interest, factors_to_remove)
    }
  }
  
  # Flag if any of the correlations are greater than cor_threshold
  if("Partial" %in% attributes(cors)$names){
    if(any(abs(cors$Partial[-1,-1][upper.tri(cors$Partial[-1,-1])]) > cor_threshold)){
      msg <- paste0("One or more correlations are greater in magnitude than the cutoff of ", cor_threshold)
      writeData(wb = the_wb, sheet = sheet_name, x = msg, startRow=row_offset,
                colNames = FALSE, rowNames = FALSE)
      row_offset <- row_offset+3
    }
  }
  
  # Return the partial correlations
  if("Partial" %in% attributes(cors)$names){
    writeData(wb = the_wb, sheet = sheet_name, x = cors$Partial, startRow=row_offset,
              colNames = TRUE, rowNames = TRUE)
    row_offset <- row_offset+3
  }
  
  # Call out if no factors of interest remain
  if(length(factors_of_interest)==0){
    msg <- paste0("No factors of interest remain for comparison.")
    writeData(wb = the_wb, sheet = sheet_name, x = msg, startRow=row_offset,
              colNames = TRUE, rowNames = TRUE)
  }
  
  return(list(the_wb, factors_of_interest))
  
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
#' @details Note: this function calls several objects in global scope:
#'          
#' @return list of Workbook and plot
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
  
  # Copy of "simplified" plot_df to return to Workbook
  plot_df_to_return <- plot_df
  
  # Add static column names, if applicable, to be referenced by 
  # ggplot directly.
  if(!is.null(color_by)){
    if(color_by=="batching_factor"){
      plot_df$the_color <- plot_df[,batching_factor]
    } else if(color_by=="factors_of_interest"){
      # This can be one or more factors.
      if(length(factors_of_interest)==1){
        plot_df$the_color <- plot_df[,factors_of_interest]
      } else {
        plot_df$the_color <- apply(plot_df[,factors_of_interest], 1, paste, collapse=".")
      }
    } else {
      stop("Unexpected color_by found. Error 1.")
    }
  }
  # Same as above but for shape
  if(!is.null(shape_by)){
    if(shape_by=="batching_factor"){
      plot_df$the_shape <- plot_df[,batching_factor]
    } else if(shape_by=="factors_of_interest"){
      # This can be one or more factors.
      if(length(factors_of_interest)==1){
        plot_df$the_shape <- plot_df[,factors_of_interest]
      } else {
        plot_df$the_shape <- apply(plot_df[,factors_of_interest], 1, paste, collapse=".")
      }
    } else {
      stop("Unexpected color_by found. Error 2.")
    }
  } 
  # Similar to above but for size_by
  if(!is.null(size_by)){
    plot_df$the_size <- plot_df[,size_by]
  }  
  
  # Basic plot
  plt <- ggplot(data=plot_df, 
    aes(x=PC1, y=PC2)) 
  # Add conditional geom_point layer
  plt <- get_condiontal_geom_point(p=plt,
                            col_logic=ifelse(is.null(color_by), FALSE, TRUE), 
                            shp_logic=ifelse(is.null(shape_by), FALSE, TRUE),
                            siz_logic=ifelse(is.null(size_by), FALSE, TRUE))
  # Add variation explained for axes and dress plot
  var_explained <- t(summary(pca_data)$importance)[1:2,'Proportion of Variance']
  var_explained <- paste0(round(as.numeric(var_explained)*100,1), "%")
  plt <- plt + 
    xlab(paste0("PC1 (", var_explained[1], ")")) + 
    ylab(paste0("PC2 (", var_explained[2], ")")) + 
    theme_bw() + 
    theme(axis.title.x = element_text(angle=0, family=plot_font$family, size=plot_font$size),
          axis.title.y = element_text(angle=90, family=plot_font$family, size=plot_font$size)) + 
    ggtitle(label=the_name)
  
  # Add data and figure to Workbook
  sheet_name <- paste0(the_name, " - PC Plot")
  addWorksheet(the_wb, sheet_name)
  writeData(wb = the_wb,
            sheet = sheet_name, 
            x = plot_df_to_return, # i.e., simplified
            colNames = TRUE, rowNames = FALSE)
  setColWidths(wb = the_wb, sheet = sheet_name, cols = 1:ncol(plot_df), widths = "auto")
  
  # Add plot
  print(plt) # must print
  insertPlot(wb=the_wb, sheet = sheet_name, width = 5, height = 3.5, fileType = "png", units = "in")
  
  # return the Workbook object
  return(list(the_wb, plt))
  
}


#' @title get_condiontal_geom_point
#' @description Get the geom_point given conditionals
#' @param p the basic ggplot object to build upon
#' @param col_logic logical whether color is an aes.
#' @param shp_logic logical whether shape is an aes.
#' @param siz_logic logical whether size is an aes.
#' @details This is a lower level function called within plot_pca_data. 
#'          There are 2^3 different plotting possibilities. This function
#'          returns one of the 8 geom_point combinations.
#' @seealso \code{plot_pca_data}
#' @return a ggplot object with plotting aes.
get_condiontal_geom_point <- function(p, col_logic, shp_logic, siz_logic){
  if(col_logic){
    if(shp_logic){
      if(siz_logic){
        out <- p + geom_point(aes(color=the_color, shape=the_shape, size=the_size), alpha=0.7) +
          labs(color = paste0(get(color_by), collapse="."), shape=paste0(get(shape_by), collapse="."), size=size_by)
      } else {
        out <- p + geom_point(aes(color=the_color, shape=the_shape), alpha=0.7) + 
          labs(color = paste0(get(color_by), collapse="."), shape=paste0(get(shape_by), collapse="."))
      }
    } else {
      if(siz_logic){
        out <- p + geom_point(aes(color=the_color, size=the_size), alpha=0.7) +
          labs(color = paste0(get(color_by), collapse="."), size=size_by)
      } else {
        out <- p + geom_point(aes(color=the_color), alpha=0.7) +
          labs(color = paste0(get(color_by), collapse="."))
      }
    }
  } else {
    if(shp_logic){
      if(siz_logic){
        out <- p + geom_point(aes(shape=the_shape, size=the_size), alpha=0.7) +
          labs(shape=paste0(get(shape_by), collapse="."), size=size_by)        
      } else {
        out <- p + geom_point(aes(shape=the_shape), alpha=0.7) +
          labs(shape=paste0(get(shape_by), collapse="."))
      }
    } else {
      if(siz_logic){
        out <- p + geom_point(aes(size=the_size), alpha=0.7) +
          labs(size=size_by)
      } else {
        out <- p + geom_point(alpha=0.7)
      }
    }      
  }
  
  return(out)
}


#' @title compare_batch_correction
#' @description compares VIFs between the two models
#' @param the_wb the Workbook to append to
#' @param annots the annotations data.frame
#' @param pca_nbc the prcomp object without batch correction data
#' @param pca_bc the prcomp object with batch correction data
#' @param batching_factor the batching_factor
#' @param factors_of_interest the factors_of_interest
#' @param nbc_pca_plot a ggplot2 object of the non-batch corrected data
#' @param bc_pca_plot a ggplot2 object of the batch corrected data
#' @details Compares the variance inflation factors for each model. 
#'          This uses the vif_threshold in the global environment to flag.
#' @return the updated Workbook with the comparisons sheet (if applicable)
compare_batch_correction <- function(the_wb, annots, pca_nbc, 
        pca_bc, batching_factor, factors_of_interest,
        nbc_pca_plot, bc_pca_plot){
  
  # Add a new sheet
  sheet_name <- "Batch Correction Evaluation"
  addWorksheet(the_wb, sheet_name)
    
  ## Transform the two datasets to compare
  # Non-batch-corrected
  nbc_df <- pca_nbc$x[,1:3]
  nbc_df <- add_column(
    as_tibble(nbc_df),
    "segmentID"=row.names(nbc_df), 
    .before=1)
  nbc_df <- base::merge(
    annots[,c("segmentID", batching_factor, eval(factors_of_interest))],
    nbc_df,
    by="segmentID"
  )
  # Batch-corrected
  bc_df <- pca_bc$x[,1:3]
  bc_df <- add_column(
    as_tibble(bc_df),
    "segmentID"=row.names(bc_df), 
    .before=1)
  bc_df <- base::merge(
    annots[,c("segmentID", batching_factor, eval(factors_of_interest))],
    bc_df,
    by="segmentID"
  )
  
  # Go through each PC and compare batch
  outer <- do.call(rbind, lapply(paste0("PC", 1:3), function(pci){
    inner <- rbind(
        get_coeffs(df=nbc_df, response=pci, batching_factor=batching_factor, the_name="Non-Batch Corrected"),
        get_coeffs(df=bc_df, response=pci, batching_factor=batching_factor, the_name="Batch Corrected")
      )
    return(inner)
  }))
  
  writeData(wb = the_wb,
            sheet = sheet_name,
            x = outer,
            colNames = TRUE, rowNames = FALSE)
  setColWidths(wb = the_wb, sheet = sheet_name, cols = 1:8, widths = "auto")
  
  
  print(cowplot::plot_grid(nbc_pca_plot, bc_pca_plot)) # must print
  insertPlot(wb=the_wb, sheet = sheet_name, width = 12, height = 3.5, 
             fileType = "png", units = "in", startRow = 1, startCol = 10)
  
  return(the_wb)
  
}

#' @title get_coeffs
#' @description lower-level function for computing and returning linear model coefficients
#' @param df a data.frame with PC scores and batching_factor
#' @param batching_factor the name of the batching factor
#' @param the_name a name that will be used to give the model
#' @return a tibble giving the coefficients
get_coeffs <- function(df, response, batching_factor, the_name){
  the_formula <- as.formula(paste0(response, " ~ ", batching_factor))
  the_model <- lm(the_formula, data=df)
  the_summary <- summary(the_model)
  coeffs <- the_summary$coefficients
  coeffs <- add_column(as_tibble(coeffs),"Coeff"=row.names(coeffs), .before=1)
  if("(Intercept)" %in% coeffs$Coeff){
    coeffs <- coeffs[-c(which(coeffs$Coeff=="(Intercept)")),]
  }
  coeffs$adj.r.squared <- rep(the_summary$adj.r.squared, nrow(coeffs))
  coeffs <- add_column(coeffs, "PC"=rep(response, nrow(coeffs)), .before=1)
  coeffs <- add_column(coeffs, "Type"=rep(the_name, nrow(coeffs)), .before=1)
  
  return(coeffs)
}
  

#' @title tranlate_model_matrix
#' @description Takes a model of batch and 1 or more factors of interest and links the level names to the factor levels.
#' @param df a data.frame with batching_factors and one or more factors_of_interest
#' @param batching_factor the name of the batching factor
#' @param factors_of_interest the names of the factors of interest
#' @details This is a lower-level function to make it easier to parse out the alias results downstream.
#' @return a list giving the data.frame linking the factors of the model with all their level names and the output of the alias function
tranlate_model_matrix <- function(
  df, batching_factor, factors_of_interest){
  
  df$Y <- 1 # we're interested only in the independent variables
  
  # The model formula
  the_formula <- as.formula(
    paste0('Y ~ ', batching_factor, " + ", paste(factors_of_interest, collapse=" + "))
  )
  
  # get model structure and parse out 
  # level_names and factor assignments
  structure <- model.matrix(the_formula, data=df)
  the_levels_names <- attributes(structure)$dimnames[[2]]
  the_factor_assignments <- attributes(structure)$assign + 1 # 1 indexed with intercept = 1
  
  # fit model
  the_model <- lm(the_formula, data = df)
  
  # Get the model terms
  the_factors <- c("(Intercept)", attributes(terms(x = the_model))$term.labels)
  
  out <- data.frame("the_levels_names"=the_levels_names, 
                    "the_factor_assignments"=the_factor_assignments)
  out$the_factors <- the_factors[match(out$the_factor_assignments, 1:length(the_factors))]
  linker <- out[,3:1]
  
  # get the correlations:
  cors <- alias(the_model, partial=TRUE)
  
  return(list(linker, cors))
}


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

