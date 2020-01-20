library(TrenaMultiScore)
library(TrenaProjectErythropoiesis)
library(TrenaProjectAD)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tmse")) {
   message(sprintf("--- creating instance of TrenaMultiScore"))
   tpe <- TrenaProjectErythropoiesis()
   tpad <- TrenaProjectAD()
   tmse <- TrenaMultiScore(tpe, "GATA2");
   tmsa <- TrenaMultiScore(tpad, "TREM2")
   }

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_getGeneHancerRegion()
   test_findOpenChromatin()
   test_findFimoTFBS()
   test_scoreMotifHitsForConservation()
   test_getTargetGeneInfo()
   test_addDistanceToTSS()
   test_scoreMotifHitsForGeneHancer()
   test_addGenicAnnotations()
   test_addChIP()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue("TrenaMultiScore" %in% is(tmsa))
   checkTrue("TrenaProjectAD" %in% is (getProject(tmsa)))

   checkTrue("TrenaMultiScore" %in% is(tmse))
   checkTrue("TrenaProjectErythropoiesis" %in% is (getProject(tmse)))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_getGeneHancerRegion <- function()
{
   message(sprintf("--- test_getGeneHancerRegion"))

   tbl.gata2 <- getGeneHancerRegion(tmse)
   checkTrue(is.data.frame(tbl.gata2))
   checkEquals(as.list(tbl.gata2), list(chrom="chr3", start=128073944, end=128621958))

   tbl.trem2 <- getGeneHancerRegion(tmsa)
   checkTrue(is.data.frame(tbl.trem2))
   checkEquals(as.list(tbl.trem2), list(chrom="chr6", start=41131001, end=41296400))

} # test_getGeneHancerRegion
#------------------------------------------------------------------------------------------------------------------------
test_findOpenChromatin <- function()
{
   message(sprintf("--- test_findOpenChromatin"))

   findOpenChromatin(tmse, "chr3", start=128470539, end=128502070)
   tbl.oc <- getOpenChromatin(tmse)
   dim(tbl.oc)
   checkTrue(nrow(tbl.oc) > 12)
   checkEquals(ncol(tbl.oc), 7)

   checkEquals(dim(getOpenChromatin(tmsa)), c(0,0))
   findOpenChromatin(tmsa, "chr3", start=128470539, end=128502070) # start=128478422, end=128494187)
   tbl.oc2 <- getOpenChromatin(tmsa)
   checkTrue(nrow(tbl.oc2) > 5000)

} # test_findOpenChromatin
#------------------------------------------------------------------------------------------------------------------------
test_findFimoTFBS <- function()
{
   message(sprintf("--- test_findFimoTFBS"))

   findOpenChromatin(tmse, "chr3", start=128483204, end=128483276)
   getOpenChromatin(tmse)

   findFimoTFBS(tmse)
   tbl.fimo <- getMultiScoreTable(tmse)
   checkEquals(dim(tbl.fimo), c(12, 9))
   checkTrue(all(c("SP2", "ZNF263") %in% tbl.fimo$tf))

   findFimoTFBS(tmse, fimo.threshold=1e-3)
   tbl.fimo <- getMultiScoreTable(tmse)

   checkTrue(nrow(tbl.fimo) > 90)
   checkTrue(all(c("SP2", "ZNF263", "CEBPB", "SP1") %in% tbl.fimo$tf)) # and many others

} # test_findFimoTFBS
#------------------------------------------------------------------------------------------------------------------------
test_scoreMotifHitsForConservation <- function()
{
   message(sprintf("--- test_scoreMotifHitsForConservation"))

   findOpenChromatin(tmse, "chr3", start=128481000, end=128489000)
   findFimoTFBS(tmse, fimo.threshold=1e-5)
   tbl.fimo <- getMultiScoreTable(tmse)
   checkEquals(dim(tbl.fimo), c(48, 9))

   scoreMotifHitsForConservation(tmse)

   tbl <- getMultiScoreTable(tmse)
   checkEquals(dim(tbl), c(48, 12))
   checkTrue(all(c("phast7", "phast30", "phast100") %in% colnames(tbl)))
   checkEqualsNumeric(mean(tbl$phast7), 0.57, tolerance=0.05)
   checkEqualsNumeric(mean(tbl$phast30), 0.66, tolerance=0.05)
   checkEqualsNumeric(mean(tbl$phast100), 0.77, tolerance=0.05)

} # test_scoreMotifHitsForConservation
#------------------------------------------------------------------------------------------------------------------------
test_getTargetGeneInfo <- function()
{
   message(sprintf("--- test_getTargetGeneInfo"))
   x <- getTargetGeneInfo(tmse)
   checkTrue(all(c("tss", "strand", "chrom", "start", "end") %in% names(x)))
   checkEquals(x$tss, 128493185)
   checkEquals(x$strand, -1)

} # test_getTargetGeneInfo
#------------------------------------------------------------------------------------------------------------------------
test_addDistanceToTSS <- function()
{
   message(sprintf("--- test_addDistanceToTSS"))

   tss <- getTargetGeneInfo(tmse)$tss
   shoulder <- 5000
   findOpenChromatin(tmse, "chr3", start=tss-(2*shoulder), end=tss+shoulder)

   findFimoTFBS(tmse, fimo.threshold=1e-5)
   tbl.fimo <- getMultiScoreTable(tmse)
   checkEquals(dim(tbl.fimo), c(123, 9))

   scoreMotifHitsForConservation(tmse)
   addDistanceToTSS(tmse)

   tbl <- getMultiScoreTable(tmse)
   checkEquals(head(sort(tbl$tss)), c(-4770, -4769, -4769, -4768, -4768, -4768))

} # test_addDistanceToTSS
#------------------------------------------------------------------------------------------------------------------------
test_scoreMotifHitsForGeneHancer <- function()
{
   message(sprintf("--- test_scoreMotifHitsForGeneHancer"))

   tss <- getTargetGeneInfo(tmse)$tss
   shoulder <- 100000
   findOpenChromatin(tmse, "chr3", start=tss-(2*shoulder), end=tss+shoulder)

   findFimoTFBS(tmse, fimo.threshold=1e-6)
   tbl.fimo <- getMultiScoreTable(tmse)
   # checkEquals(dim(tbl.fimo), c(5, 9))

   scoreMotifHitsForGeneHancer(tmse)

   tbl <- getMultiScoreTable(tmse)
   gh.values <- unique(tbl$gh)

     # when a motif is within the full range reported by genehancer, but does
     # not fall within a genehancer region, then it gets a zero score. check that
     # also check to see that multiple non-zero scores were picked up
   checkTrue(0 %in% gh.values)
   checkTrue(length(unique(gh.values)) == 5)

} # test_scoreMotifHitsForGeneHancer
#------------------------------------------------------------------------------------------------------------------------
test_addGenicAnnotations <- function()
{
   message(sprintf("--- test_addGenicAnnotations"))

   findOpenChromatin(tmse, "chr3", start=128481000, end=128489000)
   getOpenChromatin(tmse)
   findFimoTFBS(tmse, fimo.threshold=1e-3)
   addDistanceToTSS(tmse)

   scoreMotifHitsForConservation(tmse)
   addGenicAnnotations(tmse)

   tbl.fimo <- getMultiScoreTable(tmse)
   checkEquals(ncol(tbl.fimo), 15)
   checkTrue(all(c("annot.type", "annot.symbol") %in% colnames(tbl.fimo)))
   checkTrue(nrow(tbl.fimo) > 300)

} # test_addGenicAnnotations
#------------------------------------------------------------------------------------------------------------------------
test_addChIP <- function()
{
   message(sprintf("--- test_addChip"))

   findOpenChromatin(tmse, "chr3", start=128481000, end=128489000)
   tbl.oc <- getOpenChromatin(tmse)
   checkEquals(nrow(tbl.oc), 8)

   findFimoTFBS(tmse, fimo.threshold=1e-5)
   tbl.fimo <- getMultiScoreTable(tmse)
   checkEquals(dim(tbl.fimo), c(48, 9))
   tbl.fimo
   addChIP(tmse)
   tbl.fimo <- getMultiScoreTable(tmse)
   checkEquals(dim(tbl.fimo), c(48, 10))
   checkTrue(nrow(subset(tbl.fimo, chip)) > 15)

} # test_addChIP
#------------------------------------------------------------------------------------------------------------------------
test_erythropoeisis.hoxb4 <- function()
{
   message(sprintf("--- test_erythropoeisis.hoxb4"))

   tms.hoxb4 <- TrenaMultiScore(tpe, "HOXB4");
   getGeneHancerRegion(tms.hoxb4)
   findOpenChromatin(tms.hoxb4)
   findFimoTFBS(tms.hoxb4, fimo.threshold=1e-3)
   scoreMotifHitsForConservation(tms.hoxb4)
   scoreMotifHitsForGeneHancer(tms.hoxb4)
   addDistanceToTSS(tms.hoxb4)

   mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28")
   addGeneExpressionCorrelations(tms.hoxb4, mtx)
   addGenicAnnotations(tms.hoxb4)
   tbl <- getMultiScoreTable(tms.hoxb4)
   tbl.sub.neg <-  subset(tbl, p.value < 0.0001 & phast100 > 0.8 & cor < -0.4 & gh > 0)
   checkEquals(sort(unique(tbl.sub.neg$tf)), c("IRF1", "MXI1", "RXRG"))
   tbl.sub.pos <-
      subset(tbl, p.value < 0.00001 & phast100 > 0.8 & cor > 0.8 & gh > 0)
   checkEquals(sort(unique(tbl.sub.pos$tf)), c("MYC", "STAT1", "ZNF263"))

} # test_erythropoeisis.hoxb4
#------------------------------------------------------------------------------------------------------------------------
test_AD.trem2 <- function()
{
   message(sprintf("--- test_AD.trems"))

   tms.trem2 <- TrenaMultiScore(tpad, "TREM2");
   getGeneHancerRegion(tms.trem2)
   findOpenChromatin(tms.trem2)
   findFimoTFBS(tms.trem2, fimo.threshold=1e-3)
   scoreMotifHitsForConservation(tms.trem2)
   scoreMotifHitsForGeneHancer(tms.trem2)
   addDistanceToTSS(tms.trem2)
   addChIP(tms.trem2)

   mtx <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/mayo.tcx.16969x262.covariateCorrection.log+scale.RData"))
   addGeneExpressionCorrelations(tms.trem2, mtx)
   addGenicAnnotations(tms.trem2)
   tbl <- getMultiScoreTable(tms.trem2)

} # test_AD.trem2
#------------------------------------------------------------------------------------------------------------------------

if(!interactive())
   runTests()
