rm(list = ls())
source("SpatialDecon_plugin.txt")

# load test data:
load("test data/save Rdata result - 2020-05-04 18_14_53.79 export.RData")
load("test data/2.0.0.89 Dev1 M-275 NGS CTAx4.RData")

# run it:
main(dataset = dataset, segmentAnnotations = segmentAnnotations, targetAnnotations = targetAnnotations, outputFolder = "testresults") 
