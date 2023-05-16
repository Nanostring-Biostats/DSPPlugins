# TCR Analysis#
# Version 1.0 #

# Produces Gini Coefficients and Shannon entropy scores
# Supports: DSP-NGS CTA, DSP-NGS WTA (mouse & human) with TCR spike-in
# Note: this script should be run only on a dataset that has not been normalized
# Please do not use spaces, special characters, or numbers when adding factors
# in the DSPDA Annotation file

#        User Options        #
##############################

# option to filter out segments that fail QC
filter_segments = FALSE

# option to filter targets that occur in less than x% of AOIs (TCR genes will not be filtered)
filter_targets = FALSE

if (filter_targets) {
  percent_thresh = .1 #e.g. 10%
}

# used to calculate detection of TCR probes over background
LOQ_thresholds = c(2, 2.5, 3)

# subtract background from TCR probes?
TCR_probes_bg_subtraction = FALSE

# if background subtraction is turned on, should background be set at geomean 
# of negative controls or LOQ set at user threshold
if (TCR_probes_bg_subtraction) {
  background_method = "geomean" # or "LOQ"
  
  if (background_method == "LOQ") {
    bg_LOQ_thresh = 2
  }
}



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

library(Biobase)
library(GeomxTools)
library(openxlsx)

# main function with GeoMxSet wrapper

main <- function(obj1, obj2, obj3, obj4){
  if(class(obj1) == "NanoStringGeoMxSet") {
    dataset <- obj1
    outputFolder <- obj3
  } else {
    # change dataframes to GeoMxSet object
    assayData <- assayDataNew(exprs = as.matrix(obj1))
    
    rownames(obj2) <- obj2[, 1]
    phenoData <- Biobase::annotatedDataFrameFrom(as.matrix(obj2), byrow = FALSE)
    phenoData@data <- obj2
    
    rownames(obj3) <- obj3[, 1]
    featureData <- Biobase::annotatedDataFrameFrom(as.matrix(obj3), byrow = TRUE)
    featureData@data <- obj3
    
    protocol_mat <- matrix(NA, nrow = 94, ncol = 14, dimnames = list(c(rownames(obj2)), c("FileVersion", "SoftwareVersion",
                                                                                        "SampleID", "Plate_ID", "SeqSetId", "tamperedIni",
                                                                                        "Raw", "Trimmed", "Stitched", "Aligned", "FumiQ30",
                                                                                        "rtsQ30", "Date", "Well")))
    protocolData = Biobase::annotatedDataFrameFrom(protocol_mat, byrow = FALSE)
    protocolData@data <- as.data.frame(protocol_mat)
    protocolData@data[["SampleID"]] <- phenoData@data[["segmentID"]]
    
    dataset <- NanoStringGeoMxSet(
      assayData = assayData,
      phenoData = phenoData,
      featureData = featureData,
      experimentData = Biobase::MIAME(),
      annotation = character(),
      protocolData = protocolData,
      dimLabels = c("TargetGUID", "segmentID"),
      design = NULL,
      featureType = "Probe",
      analyte = "RNA",
      check = FALSE)
    
    other <- as.data.frame(t(rep(NA, 6)))
    colnames(other) <- c("PKCFileName" , "PCKFileVersion", "shiftedByOne", "shiftedByOneLogic", "DSPDAExport", "DSPDAGeoMxSetVersion")
    
    pkcs <- colnames(pData(dataset))[which(grepl("LOT", colnames(pData(dataset))))]
    other$PKCFileName <- paste(pkcs, collapse = ', ')
    other$shiftedByOne <- !any(exprs(dataset) == 0)
    dataset@experimentData@other <- other
    
    outputFolder <- obj4
  }
  
  tcrData <- tcr(dataset = dataset, outputFolder = outputFolder)
  
  if (class(obj1) == "NanoStringGeoMxSet") {
    pData(obj1) <- pData(tcrData)
    return(obj1)
  } else {
    exprsData=as.data.frame(exprs(tcrData))
    rownames(exprsData) = fData(tcrData)[, "TargetGUID"]
    colnames(exprsData) = pData(tcrData)[, "segmentID"]
    return(exprsData)
  }
}

# TCR function

tcr <- function(dataset, outputFolder) {
  
  #### preliminaries ----------------------
  # update target and segments names for assayData
  pData(dataset)[, "segmentDisplayName"] <- paste(pData(dataset)[, "ScanName"], pData(dataset)[, "ROIName"], pData(dataset)[, "SegmentName"], sep=" | ")
  
  tmp <- assayData(dataset)[["exprs"]]
  colnames(tmp) <- pData(dataset)[, "segmentDisplayName"]
  rownames(tmp) <- fData(dataset)[, "TargetName"]
  assayDataElement(dataset, "exprs", validate = FALSE) <- tmp
  
  rownames(pData(dataset)) <- pData(dataset)[, "segmentDisplayName"]
  rownames(fData(dataset)) <- fData(dataset)[, "TargetName"]
  rownames(pData(protocolData(dataset))) <- pData(dataset)[, "segmentDisplayName"]
  
  ### checks/warnings
  # error if not WTA/CTA
  if (!grepl("Atlas", dataset@experimentData@other[["PKCFileName"]], fixed = TRUE) & !grepl("Human NGS Whole Transcriptome Atlas", dataset@experimentData@other[["PKCFileName"]], fixed = TRUE)) {
    stop("Base panel must be WTA or CTA to run this script")
  }
  
  # error if not TCR spike-in
  if (!grepl("TCR", dataset@experimentData@other[["PKCFileName"]], fixed = TRUE)) {
    stop("Panel must include TCR spike-in to run this script")
  }
  
  
  #### QC using GeoMxTools
  
  # shift counts by one if not already done
  if (dataset@experimentData@other[["shiftedByOne"]] == FALSE) {
    dataset <- shiftCountsOne(dataset, useDALogic=TRUE)
  }
  
  # need to do before running segment QC
  pData(dataset)[, "Module"] <- dataset@experimentData@other[["PKCFileName"]]
  
  # run segment QC
  dataset <- setSegmentQCFlags(dataset) 
  
 
  #### Filter segments and targets
  
  # filter segments
  
  if (filter_segments) {
    lowNegatives <- which(pData(protocolData(dataset))[, "QCFlags"]$LowNegatives)
    dataset = dataset[, -lowNegatives]
  }
  
  # filter targets
  
  # get list of TCR relevant probes
  tcr_probes <- fData(dataset)$TargetName[grepl('TR[A/B/D/G][C/J/V]', fData(dataset)$TargetName)]
  
  if (filter_targets) {
    
    expression_threshold <- dim(dataset)[2]*percent_thresh
    low_targets <- which(rowSums(exprs(dataset) > 1) < expression_threshold)
    targets_to_remove <- low_targets[!(low_targets %in% which(dimnames(dataset)[[1]] %in% tcr_probes))]
    
    if (length(targets_to_remove) > 0) {
      dataset = dataset[-targets_to_remove, ]
    }
  }
  
  
  #### Normalize using GeoMxTools
  
  dataset <- normalize(dataset, norm_method="quant", 
                               desiredQuantile = .9, toElt = "q_norm")
  
  
  #### calculate detection of TCR probes over background (LOQ#)
  # output logical matrix with every probe being noted as TRUE or FALSE 
  # for being above (>) LOQ at the required thresholds specified by LOQ_thresholds

  # calculate the negative geomean for each module
  negativeGeoMeans <- 
    esBy(negativeControlSubset(dataset), 
         GROUP = "Module", 
         FUN = function(x) { 
           assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
         }) 
  pData(dataset)[, "NegGeoMean"] <- negativeGeoMeans
  
  # calculate the negative geo standard deviation for each module
  negativeGeoSD <- 
    esBy(negativeControlSubset(dataset), 
         GROUP = "Module", 
         FUN = function(x) { 
           assayDataApply(x, MARGIN = 2, FUN = ngeoSD, elt = "exprs") 
         }) 
  pData(dataset)[, "NegGeoSD"] <- negativeGeoSD
  
  # initialize array to store results
  LOQ_results <- array(NA, dim = c(dim(dataset), length(LOQ_thresholds)), dimnames = list(dimnames(dataset)[[1]], dimnames(dataset)[[2]], paste0("LOQ", LOQ_thresholds)))
  
  # determine which TCR targets are above LOQ thresh and save to results array
  for (i in 1:length(LOQ_thresholds)) {
    
    n = LOQ_thresholds[i]
    LOQ <- data.frame(row.names = colnames(dataset))
  
    vars <- c("NegGeoMean", "NegGeoSD")
    LOQ <- pData(dataset)[, vars[1]] * pData(dataset)[, vars[2]] ^ n
    
    above_LOQ <- t(esApply(dataset, MARGIN = 1, 
                           FUN = function(x) {
                             x > LOQ
                           }))
    LOQ_results[, , i] <- as.matrix(above_LOQ)
   
    pData(dataset)[ , paste0("LOQ", n)] <- LOQ
    pData(dataset)[ , paste0("Targets_Above_LOQ", n)] <- apply(above_LOQ, 2, sum)
    
  }
  
  
  
  #### subtract background from TCR probes if TCR_probes_bg_subtraction = TRUE
  if (TCR_probes_bg_subtraction == TRUE) {
    
    if (background_method == "geomean") {
      
      bg <- pData(dataset)[, "NegGeoMean"]
      
    } else if (background_method  == "LOQ") {
      
      if (bg_LOQ_thresh %in% LOQ_thresholds) {
        
        bg <- pData(dataset)[, paste0("LOQ", bg_LOQ_thresh)]
        
      } else {
        
          vars <- c("NegGeoMean", "NegGeoSD")
          bg <- pData(dataset)[, vars[1]] * pData(dataset)[, vars[2]] ^ bg_LOQ_thresh
          
      }
      
    } else {
      stop("Invalid background_method specified.  Must be either 'geomean' or 'LOQ'")
    }
    
    # subtract background
    assayDataElement(dataset, "bg_exprs") <- sweep(exprs(dataset), 2, bg, "-")
    
    # floor at zero so there are no negative values
    assayDataElement(dataset, "bg_exprs")[assayDataElement(dataset, "bg_exprs") < 0] <- 0
    
  }
  
  
  #### calculate the following summary statistics from TCR probes and append to GeoMxSetObject as new segment annotations
  
  # 1. Diversity statistics - Gini Coefficient
  pData(dataset)$TCR_Gini <- apply(exprs(dataset)[tcr_probes, ], 2, calculate_gini_coefficient)
  pData(dataset)$TCR_Gini_qnorm <- apply(assayDataElement(dataset, elt = "q_norm")[tcr_probes, ], 2, calculate_gini_coefficient)
  
  if (TCR_probes_bg_subtraction) {
    pData(dataset)$TCR_Gini_bgexprs <- apply(assayDataElement(dataset, elt = "bg_exprs")[tcr_probes, ], 2, calculate_gini_coefficient)
  }
  
  
  # 2. Diversity statistics - Shannon Entropy
  
  H <- calculate_shannon_diversity(t(exprs(dataset)[tcr_probes, ]))
  pData(dataset)$ShannonH <- H
  
  
  
  ### Output results - excel file with multiple sheets
  
  # initialize a workbook to write to excel
  wb <- createWorkbook("TCR_results")
  
  # write LOQ results
  for (i in 1:dim(LOQ_results)[3]) {
    
    addWorksheet(wb, dimnames(LOQ_results)[[3]][i])
    writeData(wb,
              sheet = dimnames(LOQ_results)[[3]][i],
              LOQ_results[tcr_probes, , i],
              colNames = TRUE, rowNames = TRUE
    )
  }
  
  # write Gini coefficient results
  addWorksheet(wb, "TCR_Gini")
  writeData(wb,
            sheet = "TCR_Gini",
            pData(dataset)["TCR_Gini"],
            colNames = TRUE, rowNames = TRUE)
  
  addWorksheet(wb, "TCR_Gini_qnorm")
  writeData(wb,
            sheet = "TCR_Gini_qnorm",
            pData(dataset)["TCR_Gini_qnorm"],
            colNames = TRUE, rowNames = TRUE)
  
  if (TCR_probes_bg_subtraction) {
    addWorksheet(wb, "TCR_Gini_bgexprs")
    writeData(wb,
              sheet = "TCR_Gini_bgexprs",
              pData(dataset)["TCR_Gini_bgexprs"],
              colNames = TRUE, rowNames = TRUE)
  }
  
  # write Shannon entropy results
  addWorksheet(wb, "ShannonH")
  writeData(wb,
            sheet = "ShannonH",
            pData(dataset)["ShannonH"],
            colNames = TRUE, rowNames = TRUE)
  
  # save excel workbook
  saveWorkbook(wb,
               file = file.path(outputFolder, "TCR_outputs.xlsx", fsep = .Platform$file.sep),
               overwrite = TRUE)
  
  return(dataset)
  
} 
  
  
#' calculate_gini_coefficient
#'
#' Function to calculate gini coefficient
#' @param x a numeric vector containing non-negative values
#' @return gini_coeff
#' @export
calculate_gini_coefficient <- function(x) {
  
  weights <- rep(1/length(x), length(x))
  
  x <- x[order(x)]
  p <- cumsum(weights)
  n <- length(x)
  nu <- cumsum(weights*x)
  n <- length(nu)
  nu <- nu/nu[n]
  
  gini_coeff <- sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
  
  return(gini_coeff)
}

#' calculate_shannon_diversity
#' 
#' Function to calculate shannon diversity
#' @param x a numeric matrix containing non-negative values
#' @return H vector of shannon diversity scores
#' @export
calculate_shannon_diversity <- function(x) {
  
  total_counts <- apply(x, 1, sum)
  x <- sweep(x, 1, total_counts, "/")
  x <- -x * log(x, exp(1))
  H <- apply(x, 1, sum, na.rm = TRUE)
  
  return(H)
}
