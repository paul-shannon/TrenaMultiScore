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
   test_getFimoTFBS()

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
test_getFimoTFBS <- function()
{
   message(sprintf("--- test_FimoTFBS"))

   findOpenChromatin(tmse, "chr3", start=128483204, end=128483276)
   getOpenChromatin(tmse)
   tbl.fimo <- getFimoTFBS(tmse)

} # test_getFimoTFBS
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
