library(TrenaMultiScore)
library(TrenaProjectErythropoiesis)
library(RUnit)
library(factoextra)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tmse")) {
   message(sprintf("--- creating instance of TrenaMultiScore"))
   tpe <- TrenaProjectErythropoiesis()
   tmse <- TrenaMultiScore(tpe, "TAL1");
   }
#------------------------------------------------------------------------------------------------------------------------
build.model <- function(targetGene, fimoThresholdAsNegativeExponent=5)
{
   printf("=========== building model for %s", targetGene)

   results.subDirectory <- sprintf("fimo%d", fimoThresholdAsNegativeExponent)
   filename <- sprintf("%s.RData", targetGene)
   if(!file.exists(results.subDirectory))
      dir.create(results.subDirectory)

   tms.tg <- TrenaMultiScore(tpe, targetGene);
   getGeneHancerRegion(tms.tg)
   findOpenChromatin(tms.tg)
   fimoThreshold <- 10^(-fimoThresholdAsNegativeExponent)
   findFimoTFBS(tms.tg, fimo.threshold=fimoThreshold)
   scoreMotifHitsForConservation(tms.tg)
   scoreMotifHitsForGeneHancer(tms.tg)
   addDistanceToTSS(tms.tg)

   mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28")
   addGeneExpressionCorrelations(tms.tg, mtx)
   addGenicAnnotations(tms.tg)
   addChIP(tms.tg)

   tbl <- getMultiScoreTable(tms.tg)
   tbl$cor[which(is.na(tbl$cor))] <- 0
   tbl$motifScore <- round(-log10(tbl$p.value), 2)
   tbl$targetGene <- targetGene

   dim(tbl)
   save(tbl, file=file.path(results.subDirectory, filename))
   invisible(tbl)

} # build.model
#------------------------------------------------------------------------------------------------------------------------
buildAll <- function()
{
  if(!exists("haney.erythropoiesis.tfs"))
     source("~/github/regulatoryGenomePaper/demos/common.R")

  tfs.oi <- c("GATA1", "GATA2", "FLI1", "SPI1")
  tbls.all <- lapply(tfs.oi, function(targetGene) build.model(targetGene, 6))
  names(tbls.all) <- tfs.oi

} # buildAll
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   buildAll()
