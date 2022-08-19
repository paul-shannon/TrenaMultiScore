library(TrenaMultiScore)
library(TrenaProjectErythropoiesis)
library(factoextra)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{

} # runTests
#------------------------------------------------------------------------------------------------------------------------
explore_addRBPtoTable <- function()
{
   tpe <- TrenaProjectErythropoiesis()
   tms <- TrenaMultiScore(tpe,  "GATA2")

   loc <- getGeneHancerRegion(tms)
   with(loc, 1 + end - start) # 548kb
   findOpenChromatin(tms)
   findFimoTFBS(tms, fimo.threshold=1e-3)   # 146k hits!
   "TBX15" %in% getMultiScoreTable(tms)$tf
   addChIP(tms)
   scoreMotifHitsForConservation(tms)
   scoreMotifHitsForGeneHancer(tms)
   addGenicAnnotations(tms)
   addDistanceToTSS(tms)

   mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28")
   addGeneExpressionCorrelations(tms, mtx)

   tbl.early <- getMultiScoreTable(tms)
   head(tbl.early)
   "TBX15" %in% tbl.early$tf
   subset(tbl.early, tf=="TBX15")


   tbl.rbp <- addRnaBindingProteins(tms)
   head(tbl.rbp)

   f <- function(tf){
     if(tf %in% rownames(mtx))
        return(cor(mtx[obj@targetGene,], mtx[tf,], method="spearman"))
     else return(NA)
     }






} # explore_addRBPtoTable
#------------------------------------------------------------------------------------------------------------------------
