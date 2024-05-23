library(TrenaMultiScore)
library(TrenaProjectErythropoiesis)
library(factoextra)  # Extract and Visualize the Results of Multivariate Data Analyses

tpe <- TrenaProjectErythropoiesis()
tms <- TrenaMultiScore(tpe,  "KLF1")

loc <- getGeneHancerRegion(tms)
with(loc, 1 + end - start) # 203kb
findOpenChromatin(tms)  # consults atac-seq
   # ~/github/TrenaProjectErythropoiesis/misc/multiScore/brandAtacSeqCuration/tbl.atac.fp.RData
   # ~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.atacMerged.RData

findFimoTFBS(tms, fimo.threshold=1e-3)   # 146k hits!
dim(getMultiScoreTable(tms))  # 3956 9
addChIP(tms)    # from ReMap 2018
scoreMotifHitsForConservation(tms)
scoreMotifHitsForGeneHancer(tms)
addGenicAnnotations(tms)
addDistanceToTSS(tms)

mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28")
addGeneExpressionCorrelations(tms, mtx)

tbl <- getMultiScoreTable(tms)
dim(tbl)  # 3956 17
head(tbl)
"TBX15" %in% tbl.early$tf
subset(tbl, tf=="TBX15")


as.data.frame(t(subset(tbl, chip & abs(cor) > 0.7 & fimo_pvalue < 1e-5)[4,]))

table(subset(tbl, chip & abs(cor) > 0.7 & fimo_pvalue < 1e-3 & abs(tss) < 1000)$tf)
  #  E2F4   SP2   SP4 TFDP1
  #  1     1     2     1

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4071162/
#  E2Fs are usually bound to DNA in a heterodimer with one of three
#  functional DNA-binding cofactors, TFDP1, 2, or 3 [16–20]; however, the
#  specific role of these coregulators in E2F-mediated cell cycle control
#  during terminal erythropoiesis has not been explored.

plot(mtx["KLF1",], type="b", col="blue", pch=16, axes=FALSE,
     xlab="day", ylab="mRNA", ylim=c(4,12), main="KLF1")
x = colnames(mtx)
x = sub("day", "", x)
x = sub("\\.r[12]", "", x)
axis(1, at=0:28, labels = c("", x))
axis(2, at=4:14, labels=as.character(4:14))

