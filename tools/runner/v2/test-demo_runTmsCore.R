library(R6)
library(RUnit)
source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")

tbl.atac <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.atacMerged.RData"))
tbl.fimo <- get(load("tbl.fimo.BACH1.RData"))
dim(tbl.fimo)                 # 114021 9
tbl.fimo.promoter.only <- subset(tbl.fimo, chrom=="chr21" & start >= 29296789 & end <= 29304792)
dim(tbl.fimo.promoter.only)   #   6143 9
trenaProject <- TrenaProjectErythropoiesis()
targetGene <- "BACH1"
mtx.rna <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints/mtx-rna-pkfm-27170x14.RData"))
mtx.srm <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/harmonized-rna-srm-14-timepoints/mtx-srm-copyNumber-100x13.RData"))

#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_promoterOnlyFimoTable()
    test_inclusiveFimoTable()

} # runTests
#----------------------------------------------------------------------------------------------------
test_promoterOnlyFimoTable <- function()
{

   message(sprintf("--- test_promoterOnlyFimoTable"))

   tms <- TMS$new(trenaProject, targetGene, tbl.fimo.promoter.only, tbl.atac, quiet=FALSE)
   tms$scoreFimoTFBS()   # chip, conservation, genehancer, genic annotations, distance to tss
   tms$add.tf.mrna.correlations(mtx.rna, featureName="cor.all")

   tbl.tms <- tms$getTfTable()
   checkEquals(nrow(tbl.tms), nrow(tbl.fimo.promoter.only))

   dim(subset(tbl.tms, fimo_pvalue <= 1e-6 & chip))

   dim(subset(tbl.tms, fimo_pvalue < 1e-5 & abs(cor.all) > 0.8 & chip & gh > 600))

   tfs <- subset(tbl.tms, fimo_pvalue < 1e-5 & abs(cor.all) > 0.7 & gh > 1)$tf
   length(tfs)

   tms$addRBP()
   tms$add.rbp.mrna.correlations(mtx.rna, featureName="cor.all")   # added to tbl.rbp
   tbl.rbp <- tms$getRbpTable()
   dim(tbl.rbp)

   rbps <- unique(subset(tbl.rbp, abs(cor.all) > 0.3)$gene)
   printf("candidate rbps: %d", length(rbps))
   tbl.trena.tf <- tms$build.trena.model(tfs, list(), mtx.rna)

   checkEquals(unique(tbl.trena.tf$class), "tf")
   checkTrue(all(c("RXRA", "ZBTB7A", "SP2", "NFE2", "KLF1", "KLF16") %in% tbl.trena.tf$gene[1:10]))

   tbl.trena.both <- tms$build.trena.model(tfs, rbps, mtx.rna)

   checkTrue(all(c("tf", "rbp") %in% tbl.trena.both$class))
   checkTrue(all(c("STAU1", "RBM15", "MBNL2", "RXRA", "ALKBH5", "ZBTB7A") %in% tbl.trena.both$gene[1:10]))
   checkEquals("STAU1", tbl.trena.both$gene[1])

} # test_inclusiveFimoTable
#----------------------------------------------------------------------------------------------------
test_inclusiveFimoTable <- function()
{
   message(sprintf("--- test_inclusiveFimoTable"))

   tms <- TMS$new(trenaProject, targetGene, tbl.fimo, tbl.atac, quiet=FALSE)
   tms$scoreFimoTFBS()   # chip, conservation, genehancer, genic annotations, distance to tss
   tms$add.tf.mrna.correlations(mtx.rna, featureName="cor.all")

   tbl.tms <- tms$getTfTable()
   checkEquals(nrow(tbl.tms), nrow(tbl.fimo))

   dim(subset(tbl.tms, fimo_pvalue <= 1e-6 & chip))

   dim(subset(tbl.tms, fimo_pvalue < 1e-5 & abs(cor.all) > 0.8 & chip & gh > 600))

   tfs <- subset(tbl.tms, fimo_pvalue < 1e-5 & abs(cor.all) > 0.7 & gh > 1)$tf
   length(tfs)

   tms$addRBP()
   tms$add.rbp.mrna.correlations(mtx.rna, featureName="cor.all")   # added to tbl.rbp
   tbl.rbp <- tms$getRbpTable()
   dim(tbl.rbp)

   rbps <- unique(subset(tbl.rbp, abs(cor.all) > 0.3)$gene)
   printf("candidate rbps: %d", length(rbps))
   tbl.trena.tf <- tms$build.trena.model(tfs, list(), mtx.rna)

   checkEquals(unique(tbl.trena.tf$class), "tf")
   checkTrue(all(c("NR3C1", "ELF3", "MAFG", "BCL6", "KLF13", "IRF1") %in% tbl.trena.tf$gene[1:10]))

   tbl.trena.both <- tms$build.trena.model(tfs, rbps, mtx.rna)

   checkTrue(all(c("tf", "rbp") %in% tbl.trena.both$class))
   checkTrue(all(c("STAU1", "NR3C1", "ELF3", "RBM15", "MAFG", "BCL6") %in% tbl.trena.both$gene[1:10]))

} # test_inclusiveFimoTable
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()

# printf("candidate tfs: %d", length(tfs))
# rbps <- unique(subset(tbl.rbp, celltype=="K562" & abs(cor.all) > 0.5)$gene)
# printf("candidate rbps in K562 cells: %d", length(rbps))
# #tbl.trena.tf <- x$build.trena.model(tfs, list(), mtx.rna)
# #dim(tbl.trena.tf)  # 8 7
# tbl.trena.both <- x$build.trena.model(tfs, rbps, mtx.rna)
# dim(tbl.trena.both)  # 42 7
# linear.models <- list()
# linear.models[["spearman"]] <- x$build.linear.model(mtx.rna, sort.by.column="spearmanCoeff", candidate.regulator.max=15)
# linear.models[["pearson"]] <- x$build.linear.model(mtx.rna, sort.by.column="spearmanCoeff", candidate.regulator.max=15)
# linear.models[["betaLasso"]] <- x$build.linear.model(mtx.rna, sort.by.column="betaLasso", candidate.regulator.max=15)
# linear.models[["betaRidge"]] <- x$build.linear.model(mtx.rna, sort.by.column="betaRidge", candidate.regulator.max=15)
# linear.models[["rfScore"]] <- x$build.linear.model(mtx.rna, sort.by.column="rfScore", candidate.regulator.max=15)
# linear.models[["xgboost"]] <- x$build.linear.model(mtx.rna, sort.by.column="xgboost", candidate.regulator.max=15)
# printf("--- modeling %s", targetGene)
# #print(subset(tbl.model, p.value <= 0.2))
# #printf("adjusted r-squared: %5.2f", x$get.lm.Rsquared())
#
