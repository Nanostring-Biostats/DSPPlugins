# build spatialdecon plugin
# workflow:
# 1 copies over the plugin main() function and arguments/header code
# 2 copies text from all files in InSituSort/R/
#   (and removes duplicate license text)
# 3 copies safeTME matrix into in-line code
# 4 copies the cell matches info into the R code
# The end result is a huge .R file, atop which the main() function is appended.


### TO DO:
#- test it out (temporarily overwriting the load X with a load safeTME command)
#- 
library("styler")


rm(list = ls())

# directory of spatialdecon library:
sdpath <- "C:\\Users\\pdanaher\\Documents\\Extensions\\SpatialDecon"

# directory of cell profiles library:
libpath <- "C:\\Users\\pdanaher\\Documents\\Extensions\\cell-profile-library"

# empty text:
zz <- list()


#### copy main() function (includes arguments):
temp <- readChar("main.R", file.info("main.R")$size)
zz[[length(zz) + 1]] <- temp


#### copy all R functions from spatialdecon

# point to the spatialdecon library:
sdecondir <- "C:\\Users\\pdanaher\\Documents\\Extensions\\SpatialDecon\\R"
# and copy all its R functions:
for (name in dir(sdecondir)) {
  fname <- paste0(sdecondir, "/", name)
  temp <- readChar(fname, file.info(fname)$size)

  # remove boilerplate
  temp <- gsub("\\# SpatialDecon: mixed cell deconvolution for spatial and/or bulk gene expression data\r\n", "", temp)
  temp <- gsub("\\# Copyright \\(C\\) 2020, NanoString Technologies, Inc\\.\r\n", "", temp)
  temp <- gsub("\\#    This program is free software\\: you can redistribute it and\\/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or \\(at your option\\) any later version\\.\r\n", "", temp)
  temp <- gsub("\\#    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE\\.  See the GNU General Public License for more details\\.\r\n", "", temp)
  temp <- gsub("\\#    You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.\r\n\\# Contact us:\r\n# NanoString Technologies, Inc\\.\r\n", "", temp)
  temp <- gsub("\\# 530 Fairview Avenue N\r\n", "", temp)
  temp <- gsub("\\# Seattle, WA 98109\r\n", "", temp)
  temp <- gsub("\\# Tel: \\(888\\) 358-6266\r\n", "", temp)
  temp <- gsub("\\# pdanaher@nanostring.com", "", temp)

  # remove package references:
  temp <- gsub("SpatialDecon::", "", temp)

  zz[[length(zz) + 1]] <- temp
}

#### copy data:
temp <- readChar("spatialdecon data.R", file.info("spatialdecon data.R")$size)
zz[[length(zz) + 1]] <- temp



cat("", file = "SpatialDecon_plugin.R")
for (i in 1:length(zz)) {
  cat(zz[[i]], file = "SpatialDecon_plugin.R", append = T, sep = " ", fill = TRUE)
}

# file.copy(from = "SpatialDecon_plugin.R", "final_script/SpatialDecon_plugin.R")

# setwd("final_script/")
# format:
#use_tidy_style()
