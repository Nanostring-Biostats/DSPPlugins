rm(list = ls())
source("SpatialDecon_plugin.txt")

# load test data:
load("test data/save Rdata result - 2020-05-04 18_14_53.79 export.RData")

# run it:
main(dataset = dataset, segmentAnnotations = segmentAnnotations, targetAnnotations = targetAnnotations, outputFolder = "testresults") 
