library(RUnit)
mtx.tmp <- get(load("~/github/TrenaProjectErythropoiesis/viz/srm.vs.mrna/shinyapps.io/srm.rna.averaged.clean.RData"))
targetGenes <- rownames(mtx.tmp)
length(targetGenes) # 103
#targetGene <- "CHD1"
#targetGene <- "BACH1"
#targetGene <- "CBP"
#targetGene <- "CEBPB"
targetGene <- targetGenes[12]
source("./tmsRunner.R")
trenaProject <- TrenaProjectErythropoiesis()
tbl.atac.merged <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.atacMerged.RData"))


x <- TMS$new(targetGene, trenaProject)
x$addGeneHancer()
tbl.gh <- x$getGeneHancer()
#checkEquals(dim(tbl.gh), c(25, 16))
tbl.ghRegion <- x$getGeneHancerRegion()
#checkEquals(tbl.ghRegion$start, min(tbl.gh$start) - 1000)
#checkEquals(tbl.ghRegion$end, max(tbl.gh$end) + 1000)
x$addOpenChromatin(promoter.only=TRUE)
tbl.oc <- x$getOpenChromatin()

# x$viewOpenChromatin()
#checkEquals(dim(tbl.oc), c(88, 7))
motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "hocomocov11-core"))
length(motifs)
x$addAndScoreFimoTFBS(fimo.threshold=1e-4, motifs=motifs)
tbl.tms <- x$getTfTable()
dim(tbl.tms)  # 165 16
length(unique(tbl.tms$tf))  # [1] 44
x$addRBP()
tbl.rbp <- x$getRbpTable()
dim(tbl.rbp)  # 7652  12
getExpressionMatrixNames(trenaProject)
mtx <- getExpressionMatrix(trenaProject, "brandLabDifferentiationTimeCourse-27171x28")
                           #"brandLabDifferentiationTimeCourse-27171x28-namesCorrected")
x$add.tf.mrna.correlations(mtx, featureName="cor.all")
tbl.tms <- x$getTfTable()
dim(tbl.tms)

fivenum(tbl.tms$cor.all)
x$add.rbp.mrna.correlations(mtx, featureName="cor.all")
tbl.rbp <- x$getRbpTable()
dim(tbl.rbp)
fivenum(tbl.rbp$cor.all)
dim(tbl.rbp)  # 7652   13

tfs <- unique(subset(tbl.tms, gh > 600 & abs(cor.all) > 0.5)$tf)
printf("candiate tfs: %d", length(tfs))
rbps <- unique(subset(tbl.rbp, celltype=="K562" & abs(cor.all) > 0.5)$gene)
printf("candidate rbps in K562 cells: %d", length(rbps))
#tbl.trena.tf <- x$build.trena.model(tfs, list(), mtx)
#dim(tbl.trena.tf)  # 8 7
tbl.trena.both <- x$build.trena.model(tfs, rbps, mtx)
head(tbl.trena.both[order(tbl.trena.both$spearmanCoeff, decreasing=TRUE),], n=20)
tbl.model <- x$build.linear.model(sort.by.column="spearmanCoeff")
# tbl.model <- x$build.linear.model(sort.by.column="betaLasso")
# tbl.model <- x$build.linear.model(sort.by.column="rfScore")
printf("--- modeling %s", targetGene)
print(subset(tbl.model, p.value <= 0.2))
printf("adjusted r-squared: %5.2f", x$get.lm.Rsquared())
