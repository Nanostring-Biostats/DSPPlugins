# title: RNA Negative Normalization
# description: Performs negative normalization on RNA count data collapsed by target
#              This script supports normalization for multi-panel analyses.
# input: QC'd RNA counts collapsed by target
# return: A dataset object that contains negative normalized counts

# Main function used by DSP-DA
main <- function(dataset, segmentAnnotations, targetAnnotations, outputFolder) {
  # Check for required columns
  cols_check <- c("ProbePool", "CodeClass", "TargetGUID")
  if(!all(cols_check %in% colnames(targetAnnotations))) {
    stop("Error: Required target annotation columns missing.")
  }
  if(any(is.na(targetAnnotations[["ProbePool"]]))) {
    stop("Error: Missing ProbePool designations. Make sure this is RNA data.")
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
