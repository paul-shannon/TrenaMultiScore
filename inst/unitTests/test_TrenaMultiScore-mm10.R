library(TrenaMultiScore)
library(TrenaProjectMM10)
library(RUnit)
#----------------------------------------------------------------------------------------------------
if(!exists("tms.mm10")) {
   message(sprintf("--- creating instance of TrenaMultiScore, mm10"))
   tp <- TrenaProjectMM10()
   tms.mm10 <- TrenaMultiScore(tp, "Trem2");
   }
#----------------------------------------------------------------------------------------------------
runTests <- function()
{

} # test_getGeneHancerRegion
#----------------------------------------------------------------------------------------------------
test_getGeneHancerRegion <- function()
{
   message(sprintf("--- test_getGeneHancerRegion"))

   tbl.gata2 <- getGeneHancerRegion(tmse)
   checkTrue(is.data.frame(tbl.gata2))
   checkEquals(as.list(tbl.gata2), list(chrom="chr3", start=128073944, end=128621958, width=548015))

   tbl.trem2 <- getGeneHancerRegion(tmsa)
   checkTrue(is.data.frame(tbl.trem2))
   checkEquals(as.list(tbl.trem2), list(chrom="chr6", start=41131001, end=41296400, width=165400))

} # test_getGeneHancerRegion
#------------------------------------------------------------------------------------------------------------------------


