## ----main, echo=TRUE------------------------------------------------------------------------------
suppressWarnings(
suppressMessages({
   library(TrenaMultiScore)
   library(TrenaProjectErythropoiesis)
   library(factoextra)
     # Extract and Visualize the Results of Multivariate Data Analyses

tpe <- TrenaProjectErythropoiesis()
tms <- TrenaMultiScore(tpe,  "KLF1")
}))


## ----enhancers, echo=TRUE-------------------------------------------------------------------------
loc <- getGeneHancerRegion(tms)
print(sprintf("genehancer genomic span: %d kb", with(loc, 1 + end - start)))


## ----atac, echo=TRUE------------------------------------------------------------------------------
findOpenChromatin(tms)  # consults atac-seq


## ----fimo, echo=TRUE------------------------------------------------------------------------------
findFimoTFBS(tms, fimo.threshold=1e-3)


## ----basicTable, echo=TRUE------------------------------------------------------------------------
tbl = getMultiScoreTable(tms)
rownames(tbl) = NULL
print(sprintf("rows: %d   cols: %d",   nrow(tbl), ncol(tbl)))


## ----comment='', echo=FALSE, results='asis'-------------------------------------------------------
 knitr::kable(head(tbl), caption = "First 6 rows", floating.environment="sidewaystable")


## ----chip, echo=TRUE------------------------------------------------------------------------------
addChIP(tms)
print(sprintf("ncol(tbl): %d", ncol(getMultiScoreTable(tms))))



## ----tss, echo=TRUE-------------------------------------------------------------------------------
addDistanceToTSS(tms)
print(sprintf("ncol(tbl): %d", ncol(getMultiScoreTable(tms))))


## ----expression, echo=TRUE------------------------------------------------------------------------
mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28")
addGeneExpressionCorrelations(tms, mtx)
print(sprintf("ncol(tbl): %d", ncol(getMultiScoreTable(tms))))


## ----subset1, echo=TRUE---------------------------------------------------------------------------
tbl <- getMultiScoreTable(tms)
tbl.strong <- subset(tbl, chip & abs(cor) > 0.7 & fimo_pvalue < 1e-3 & abs(tss) < 1000)
print(table(tbl.strong$tf))

plot(mtx["KLF1",], type="b", col="blue", pch=16, axes=FALSE,
      xlab="day", ylab="mRNA", ylim=c(4,14), main="KLF1")
dayNames = colnames(mtx)
dayNames = sub("day", "", dayNames)
dayNames = sub("\\.r[12]", "", dayNames)
axis(1, at=0:28, labels = c("", dayNames))
axis(2, at=4:14, labels=as.character(4:14))
lines(mtx["E2F4",], type="b", pch=16, col="red")
lines(mtx["TFDP1",], type="b", pch=16, col="darkgreen")
legend(11, 7, c("KLF1", "E2F4", "TFDP1"), c("blue", "red", "darkgreen"))
## ----KLF1 expression, echo=TRUE, out.width = '70%'------------------------------------------------
#plot(mtx["KLF1",], type="b", col="blue", pch=16, axes=FALSE,
#      xlab="day", ylab="mRNA", ylim=c(4,14), main="KLF1")
# x = colnames(mtx)
# x = sub("day", "", x)
# x = sub("\\.r[12]", "", x)
# axis(1, at=0:28, labels = c("", x))
# axis(2, at=4:14, labels=as.character(4:14))

