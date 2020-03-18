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
goi <- function()
{
   tbl <- read.table("G2vsG3GenesUP.txt", sep="\t", as.is=TRUE, header=TRUE, nrow=100)
   dim(tbl)
   genes <- tbl$GeneName
   additional.genes <- read.table(file="additionalGenes.txt", stringsAsFactors=FALSE)$V1
   length(genes)
   length(additional.genes)
   all.goi <- sort(unique(c(genes, additional.genes)))
   length(all.goi)
   library(org.Hs.eg.db)
   tbl.ref <- select(org.Hs.eg.db, key=all.goi, columns=c("SYMBOL", "ENTREZID"), keytype="SYMBOL")
   successes <- tbl.ref[which(!is.na(tbl.ref$ENTREZID)),]$SYMBOL
   failures <- tbl.ref[which(is.na(tbl.ref$ENTREZID)),]$SYMBOL
   tbl.ref2 <- select(org.Hs.eg.db, key=failures, columns=c("SYMBOL", "ALIAS", "ENTREZID"), keytype="ALIAS")
   successes.round2 <- tbl.ref2[which(!is.na(tbl.ref2$SYMBOL)), "SYMBOL"]
   goi.final <- sort(unique(c(successes, successes.round2)))
   length(goi.final)  # 112

   return(goi.final)

} # goi
#------------------------------------------------------------------------------------------------------------------------
buildAll <- function()
{
  if(!exists("haney.erythropoiesis.tfs"))
     source("~/github/regulatoryGenomePaper/demos/common.R")

  tfs.oi <- c("GATA1", "GATA2", "FLI1", "SPI1")
  tbls.all <- lapply(tfs.oi, function(targetGene) build.model(targetGene, 2))
  names(tbls.all) <- tfs.oi

} # buildAll
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   buildAll()
