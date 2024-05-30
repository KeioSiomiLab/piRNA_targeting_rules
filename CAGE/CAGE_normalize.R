.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2/")

library(CAGEr)

library(magrittr)
library(stringr)

ggsave <- function(...) {
  ggplot2::ggsave(...)
  invisible()
}
Covdir <- "/media/hermione/CAGE/figures/"
ctssdir <- "/media/hermione/CAGE/ctss/"
setwd(Covdir)
file_Cov <- dir(ctssdir, pattern = ".ctss$", full.names = TRUE)
file_Cov2 <- file_Cov  %>%
  str_replace_all(".ctss$", "") %>%
  str_replace(paste0(ctssdir, "/"), "")
ce <- CAGEexp(genomeName = "BSgenome.Melanogaster.dm6.flybase",
              inputFiles = file_Cov,
              inputFilesType = "ctss",
              sampleLabels = file_Cov2)

getCTSS(ce)
CTSStagCountSE(ce)
CTSScoordinatesGR(ce)
CTSStagCountDF(ce)

mergeSamples(ce, mergeIndex = c(1,1,1,2,2,2),
             mergedSampleLabels = c("siEGFP","siPiwi"))
ppi <- 300
png(paste0(Covdir,"cumulative.png"), width = 6*ppi, height = 6*ppi, res = ppi)
plotReverseCumulatives(ce, fitInRange = c(5,500000), onePlot = TRUE,xlim = c(1, 1e+06))
dev.off()

normalizeTagCount(ce,method = "powerLaw", fitInRange = c(5,500000),alpha = 1.05, T = 400000)
ce[["tagCountMatrix"]]

clusterCTSS(object = ce, threshold = 0.1, thresholdIsTpm = TRUE, nrPassThreshold = 0.1,
            method = "distclu",maxDist = 20, removeSingletons = TRUE, keepSingletonsAbove = 5)

tagClusters(ce, samples = "siEGFP")

cumulativeCTSSdistribution(ce, clusters = "tagClusters", useMulticore = T)
quantilePositions(ce, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)


png(paste0(Covdir,"interquantile.png"), width = 6*ppi, height = 6*ppi, res = ppi)
plotInterquantileWidth(ce, clusters = "tagClusters", tpmThreshold = 3, qLow = 0.1, qUp = 0.9)
dev.off()

aggregateTagClusters(ce, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)
ce$outOfClusters / ce$librarySizes * 100

consensusClustersGR(ce)



exportCTSStoBedGraph(ce, values = "normalized",format = "bedGraph", oneFile = TRUE)
exportCTSStoBedGraph(ce, values = "normalized",format = "BigWig", oneFile = FALSE)
exportToBed(object = ce, what = "tagClusters", qLow = 0.1, qUp = 0.9, oneFile = TRUE)
exportToBed(object = ce, what = "tagClusters", qLow = 0.1, qUp = 0.9, oneFile = FALSE)
exportToBed(object = ce, what = "tagClusters", oneFile = FALSE)


