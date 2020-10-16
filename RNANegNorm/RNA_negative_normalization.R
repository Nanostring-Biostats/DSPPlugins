### functions to perform negative normalization:
# note: this code assumes the input is QC'd RNA counts collapsed by target



# main function called by DSP-DA:
main <- function(dataset, segmentAnnotations, targetAnnotations, outputFolder) {
  pools <- unique(targetAnnotations[["ProbePool"]])
  pool_neg_norm <- lapply(pools, 
    function(pool) {
      # Get pool and corresponding target counts
      pool_neg <- 
        targetAnnotations[targetAnnotations[["CodeClass"]] == "Negative" & 
                            targetAnnotations[["ProbePool"]] == pool, "TargetGUID"]
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
