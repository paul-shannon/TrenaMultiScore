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
build.model <- function(targetGene)
{
   printf("=========== building model for %s", targetGene)

   tms.tg <- TrenaMultiScore(tpe, targetGene);
   getGeneHancerRegion(tms.tg)
   findOpenChromatin(tms.tg)
   findFimoTFBS(tms.tg, fimo.threshold=1e-2)
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
   save(tbl, file=sprintf("fimo2/%s.RData", targetGene))
   invisible(tbl)


} # build.model
#------------------------------------------------------------------------------------------------------------------------
buildAll <- function()
{
  if(!exists("haney.erythropoiesis.tfs"))
     source("~/github/regulatoryGenomePaper/demos/common.R")

  tfs.oi <- c(haney.erythropoiesis.tfs(), "TBX15", "SIX1", "KLF1", "ZBTB7A", "TMCC2", "NR1H2")
  tbls.all <- lapply(tfs.oi, function(targetGene) build.model(targetGene))
  names(tbls.all) <- tfs.oi

} # buildAll
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   buildAll()
