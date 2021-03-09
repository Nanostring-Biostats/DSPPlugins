# DSP-NGS Batch Correction #

# Runs PCA on samples
# Perform Batch Correction (BC)
# Runs PCA on BC samples
# Returns a BC data.frame


### ###############
### User Inputs ###
### ###############

# batching_factor <- "SlideName"



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
  
  # Check user input
  check_user_input(batching_factor)
  
  # Run batch correction
  bc <- run_batch_correction(
    df=dataset,
    annots=segmentAnnotations,
    batching_factor=batching_factor
  )
  
  # Send dataset to disk
  x <- dataset
  write.table(x, 
    file = file.path(outputFolder, 
         paste0("dataset.tsv"),
         fsep = .Platform$file.sep), 
    sep="\t", col.names=TRUE,
    row.names=TRUE, quote=FALSE)
              
  return(x)
  
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
    "rlang" # tunneling data-variables
  )

  # Load packages into session 
  null <- lapply(.packages, require, 
                 character.only=TRUE)
}

# Checks user input
# @param batching_factor is character provided by user
# @min_observed is numeric and is the minimum 
#      number of observation per batch required for
#      lme4 to work.
check_user_input <- function(bactching_factor, min_observed=3){
  
  pass <- TRUE # passing flag
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
  cl <- makeCluster(n_processors)
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
  # df and return.
  resid_df <- resid_df[rownames(df),]
  resid_df <- resid_df[,colnames(df)]
  return(resid_df)
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

