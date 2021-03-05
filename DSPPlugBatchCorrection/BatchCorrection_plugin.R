# DSP-NGS Batch Correction #

# Runs PCA on samples
# Perform Batch Correction (BC)
# Runs PCA on BC samples
# Returns a BC data.frame


### ###############
### User Inputs ###
### ###############




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
main <- function(dataset, segmentAnnotations, targetAnnotations, outputFolder) {

  # Load packages
  load_packages()
  
  # Send dataset to disk
  x <- dataset
  write.table(x, 
              file = file.path(outputFolder, 
                    paste0("dataset.tsv"),
                    fsep = .Platform$file.sep), 
              sep="\t", col.names=TRUE, row.names=TRUE,
              quote=FALSE)
  return(x)
  

  
  
}

### ################ ###
### Helper functions ###
### ################ ###

# loads packages
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

### #####
### Notes
### #####

# dataset has features (rows) by (samples)
# columns. These features and columns
# are referenced by alphanumeric strings
# and need may need to be tranlated.

# Column names in the `dataset` object 
# correspond to rows in the `SegmentID`
# column of `segmentAnnotations` object.

# Features. The targetAnnotations object
# links TargetGUID with TargetName.

