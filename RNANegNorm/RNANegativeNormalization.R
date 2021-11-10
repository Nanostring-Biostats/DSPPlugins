# DSP-NGS RNA Negative Normalization #
# Version 1.1 #

# Performs background normalization on NGS RNA count data collapsed by target
# for multi-panel analyses (e.g. with spike-in(s))
# OUTPUT: a dataset within the DSPDA that contains background normalized counts
# Supports: DSP-NGS CTA, DSP-NGS WTA (mouse & human)
# The script should be run on the Biological Probe QC dataset

##########################################################
#### No User Inputs. PLEASE DO NOT CHANGE CODE BELOW HERE  ####
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
  
  normData <- negNorm(dataset = dataset,
                      segmentAnnotations = segmentAnnotations,
                      targetAnnotations = targetAnnotations, 
                      outputFolder = outputFolder)
  
  if(class(obj1) == "NanoStringGeoMxSet"){
    assayDataElement(obj1, "exprs") <- normData
    return(obj1)
  }else{
    return(normData)
  }
}

# negative normalization function
negNorm <- function(dataset, segmentAnnotations, targetAnnotations, outputFolder) {
  # Check for required columns
  cols_check <- c("ProbePool", "CodeClass", "TargetGUID")
  if(!all(cols_check %in% colnames(targetAnnotations))) {
    stop("Error: Required target annotation columns missing.")
  }
  if(any(is.na(targetAnnotations[["ProbePool"]]))) {
    stop("Error: Missing ProbePool designations. Make sure this is NGS RNA data.")
  }
  pools <- unique(targetAnnotations[["ProbePool"]])
  pool_neg_norm <- lapply(pools, 
    function(pool) {
      # Get pool and corresponding target counts
      pool_neg <- 
        targetAnnotations[targetAnnotations[["CodeClass"]] == "Negative" & 
                            targetAnnotations[["ProbePool"]] == pool, "TargetGUID"]
      if(length(pool_neg) < 1) {
        stop(paste0("Error: No negative could be located for probe pool ", 
                      pool, "."))
      }
      if(length(pool_neg) > 1) {
        stop(paste0("Error: More than one negative was located for probe pool ", 
                      pool, "."))
      }
      pool_targets <- 
        targetAnnotations[targetAnnotations[["ProbePool"]] == pool, "TargetGUID"]

      # Calculate normalization factor and normalized counts
      pool_neg_factors <- 
        unlist(dataset[pool_neg, ] / 
                 exp(mean(log(as.numeric(dataset[pool_neg, ])))))
      pool_counts <-
        as.matrix(dataset[pool_targets, ]) %*% diag(1 / pool_neg_factors)
   })

  #Collapse data back into one data frame
  neg_norm_df <- data.frame(do.call(rbind, pool_neg_norm))
  colnames(neg_norm_df) <- colnames(dataset)
  neg_norm_df <- neg_norm_df[rownames(dataset), ]
  return(neg_norm_df)
}
