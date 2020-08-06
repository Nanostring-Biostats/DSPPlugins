# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression data
# Copyright (C) 2020, NanoString Technologies, Inc.
#    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#    You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.
# Contact us:
# NanoString Technologies, Inc.
# 530 Fairview Avenue N
# Seattle, WA 98109
# Tel: (888) 358-6266
# pdanaher@nanostring.com


##########################################################
####      Arguments to be modified by user            #### 
##########################################################


# Specify the filename of the cell profile matrix:
# Instructions: specifying a new cell profile matrix: 
#   Use a .csv file with cell-type-specific expression profiles, with genes in 
#   rows and cell types in columns. Row names should be in HUGO gene names.
# Instructions: downloading a pre-specified cell profile matrix:
#   You can supply your own cell profile matrix, or you can download a pre-specified one from
#   https://github.com/Nanostring-Biostats/Extensions/cell-profile-library
# Instructions: uploading to DSPDA:
# ______________
cell_profile_filename = "safeTME-for-tumor-immune.csv"

# if you have segments selected to be pure tumor, free of immune and stroma, 
# then specify a column name giving the identities of those segments. The 
# code expects this column to have"tumor" entered for those segments, and other
# values elsewhere. 
pure_tumor_column_name = "none"

# enter the name of the column giving nuclei counts
nuclei_count_column_name = "AOINucleiCount"

# define cell types to be added together in the final result:
# example syntax:
# merges = list()
# merges[["T"]] = c("CD8.T", "CD4.T")
# merges[["myeloid"]] = c("macrophage", "monocyte", "DC")
merges = list()


##########################################################
#### end of arguments. DO NOT CHANGE CODE BELOW HERE  ####
##########################################################

library(logNormReg)
library(pheatmap)
library(viridis)

main <- function(dataset, segmentAnnotations, targetAnnotations, outputFolder) {
  
  #### preliminaries ----------------------
  dataset = as.matrix(dataset)
  
  # access cell profile matrix file:
  X = as.matrix(read.csv(cell_profile_filename, header = T, row.names = 1))
  
  # parse merges:
  merges.full = NULL
  if (length(merges) > 0) {
    # initialize with 1:1 mapping:
    merges.full = list()
    for (name in colnames(X)) {
      merges.full[[name]] = name
    }
    # add merges:
    for (name in names(merges0)) {
      # remove entries for cells specified by user and replace with their entries:
      merges.full[merges[[name]]] = NULL
      merges.full[[name]] = merges[[name]]
    }
  }
  
  # parse nuclei column
  cell_counts = NULL
  if (is.element(nuclei_count_column_name, colnames(segmentAnnotations))) {
    cell_counts = as.numeric(segmentAnnotations[, nuclei_count_column_name])
  }
  if (!is.element(nuclei_count_column_name, colnames(targetAnnotations))) {
    warning("The value entered for nuclei_count_column_name was not a column header in the target annotations. 
            Results will not be output on the scale of cell counts; just in abundance scores and proportions.")
  }
    
  # parse pure tumor column
  is_pure_tumor = NULL
  if (is.element(pure_tumor_column_name, colnames(targetAnnotations))) {
    is_pure_tumor = targetAnnotations[, pure_tumor_column_name] == "tumor"
    is_pure_tumor = replace(is_pure_tumor, is.na(is_pure_tumor), FALSE)
  }
  if (!is.element(pure_tumor_column_name, colnames(targetAnnotations)) & (pure_tumor_column_name != "none")) {
    warning("The value entered for pure_tumor_column_name was not a column header in the target annotations.")
  }
  
  # format data for spatialdecon:
  norm = dataset[targetAnnotations$TargetGUID, segmentAnnotations$segmentID]
  rownames(norm) = targetAnnotations$TargetName
  colnames(norm) = segmentAnnotations$segmentDisplayName

  # calculate background:
  bg = derive_GeoMx_background(norm = norm,
                               probepool = targetAnnotations$ProbePool,
                               negnames = targetAnnotations$TargetName[targetAnnotations$CodeClass == "Negative"])
  
  
  
  #### run decon: ----------------------------------------
  # decon:
  res = spatialdecon(norm = norm,
                     bg = bg, 
                     X = X, 
                     is_pure_tumor = is_pure_tumor, 
                     cell_counts = cell_counts, 
                     cellmerges = merges.full)
  
  # reverse decon:
  rdecon = reverseDecon(norm = norm, 
                        beta = res$beta, 
                        epsilon = NULL)
  
  
  #### write results files: ---------------------------------------------
  write.csv(res$beta, file = file.path(outputFolder, "cell_abundance_scores.csv", fsep = .Platform$file.sep))
  write.csv(res$p, file = file.path(outputFolder, "cell_pvalues.csv", fsep = .Platform$file.sep))
  if (is.element("cell.counts", names(res))) {
    write.csv(res$cell.counts, file = file.path(outputFolder, "cell_count_estimates.csv", fsep = .Platform$file.sep))
  }  
  if (is.element("beta.granular", names(res))) {
    write.csv(res$beta.granular, file = file.path(outputFolder, "cell_abundance_scores_granular.csv", fsep = .Platform$file.sep))
  }
  if (is.element("cell.counts.granular", names(res))) {
    write.csv(res$cell.counts.granular, file = file.path(outputFolder, "cell_count_estimates_granular.csv", fsep = .Platform$file.sep))
  }   
  
  # reverse decon resids
  write.csv(rdecon$resids, file = file.path(outputFolder, "reverse_decon_residuals.csv", fsep = .Platform$file.sep))

  # reverse decon summary stats of gene dependency on cell mixing:
  sumstats = cbind(rdecon$cors, rdecon$resid.sd)
  colnames(sumstats) = c("cor w cell mixing", "residual SD from cell mixing")
  write.csv(sumstats, file = file.path(outputFolder, "reverse_decon_residuals.csv", fsep = .Platform$file.sep))
  
  #### results figures: ---------------------------------------------
  
  # show just the original cells, not tumor abundance estimates derived from the is.pure.tumor argument:
  cells.to.plot = intersect(rownames(res$beta), union(colnames(X), names(merges.full)))
  
  # heatmap
  pdf(file = file.path(outputFolder, "cell_abundance_heatmap.pdf", fsep = .Platform$file.sep), width = 12)
  pheatmap(res$beta[cells.to.plot, ],
           col = viridis_pal(option = "B")(100),
           fontsize_col = 4)
  dev.off()
 
  # barplot
  pdf(file = file.path(outputFolder, "cell_abundance_barplot.pdf", fsep = .Platform$file.sep), width = 12)
  par(mar = c(15,5,2,1))
  TIL_barplot(mat = res$beta[cells.to.plot, ],
              draw_legend = TRUE)
  dev.off()
  
}



#### below here: all the functions of the SpatialDecon package ####