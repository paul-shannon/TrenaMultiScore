library(TrenaMultiScore)
library(TrenaProjectAD)
library(RUnit)
library(factoextra)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tp")) {
   message(sprintf("--- creating instance of TrenaProjectAD"))
   tp <- TrenaProjectAD()
   }
#------------------------------------------------------------------------------------------------------------------------
build.model <- function(targetGene, fimoThresholdAsNegativeExponent=5)
{
   printf("=========== building model for %s", targetGene)

   results.subDirectory <- sprintf("fimo%d", fimoThresholdAsNegativeExponent)
   filename <- sprintf("%s.RData", targetGene)
   if(!file.exists(results.subDirectory))
      dir.create(results.subDirectory)

   tms.tg <- TrenaMultiScore(tp, targetGene);
   getGeneHancerRegion(tms.tg)
   findOpenChromatin(tms.tg)

   findFimoTFBS(tms.tg, fimo.threshold=10^(-fimoThresholdAsNegativeExponent))
   scoreMotifHitsForConservation(tms.tg)
   scoreMotifHitsForGeneHancer(tms.tg)
   addDistanceToTSS(tms.tg)

   #mtx <- getExpressionMatrix(tp, "temporalCortex.15167x264")
   mtx <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/mayo.tcx.16969x262.covariateCorrection.log+scale.RData"))
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
  goi <- c("TREM2", "MEF2C", "APOE", "PILRA", "BIN1", "VGF", "TYROBP", "ATRX", "MAPT", "APP")

  tbls.all <- lapply(goi, function(targetGene) build.model(targetGene, fimoThresholdAsNegativeExponent=3))
  names(tbls.all) <- goi

} # buildAll
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   buildAll()
