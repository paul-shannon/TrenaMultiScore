library(TrenaMultiScore)
library(TrenaProjectErythropoiesis)
library(TrenaProjectAD)
library(RUnit)
library(factoextra)

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
# both haney and krumsiek report that GATA1 regulates TAL1 (aka SCL)
# can we see that?  brandLabDifferentiationTimeCourse-27171x28 shows 0.75 correlation
test_erythropoeisis.tal1 <- function()
{
   message(sprintf("--- test_erythropoeisis.tal1"))

   tms.tal1 <- TrenaMultiScore(tpe, "TAL1");
   getGeneHancerRegion(tms.tal1)   # 190kb
   findOpenChromatin(tms.tal1)
   findFimoTFBS(tms.tal1, fimo.threshold=1e-3)
   scoreMotifHitsForConservation(tms.tal1)
   scoreMotifHitsForGeneHancer(tms.tal1)
   addDistanceToTSS(tms.tal1)

   mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28")
   addGeneExpressionCorrelations(tms.tal1, mtx)
   addGenicAnnotations(tms.tal1)
   addChIP(tms.tal1)

   tbl <- getMultiScoreTable(tms.tal1)
   tbl$cor[which(is.na(tbl$cor))] <- 0
   tbl$motifScore <- round(-log10(tbl$p.value), 2)

   dim(tbl)
   cor.spread <- fivenum(abs(tbl$cor))
   cor.threshold <- cor.spread[4]
   cor.threshold <- 0.9 * cor.spread[5]
   tbl <- subset(tbl, abs(cor) >= cor.threshold)
   dim(tbl)

   tfoi <- unique(tbl$tf)
   length(tfoi)
   tblc <- data.frame(row.names=tfoi)
   tblc$cor <- as.numeric(lapply(tfoi, function(TF) mean(abs(subset(tbl, tf==TF)$cor))))
   tblc$gh <-  as.numeric(lapply(tfoi, function(TF) mean(subset(tbl, tf==TF)$gh)))
   tblc$count <- as.numeric(lapply(tfoi, function(TF) nrow(subset(tbl, tf==TF))))
   tblc$chip <- as.numeric(lapply(tfoi, function(TF) sum(subset(tbl, tf==TF)$chip)))
   tblc$fimo <- as.numeric(lapply(tfoi, function(TF) mean(subset(tbl, tf==TF)$motifScore)))
   tblc$phast <- as.numeric(lapply(tfoi, function(TF) mean(as.numeric(as.matrix(subset(tbl, tf==TF)[, c("phast7", "phast30", "phast100")])))))
   #tblc$tss   <- as.numeric(lapply(tfoi, function(TF) mean(abs(subset(tbl, tf==TF)$tss))))


   mtx <- as.matrix(tblc)

   pca <- prcomp(mtx, scale=TRUE)
   library(factoextra)
   fviz_eig(pca)

  fviz_pca_ind(pca,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
               )

  fviz_pca_var(pca,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
               )

   fviz_pca_biplot(pca, repel = TRUE,
                   col.var = "#2E9FDF", # Variables color
                   col.ind = "#696969"  # Individuals color
                   )


} # test_erythropoeisis.tal1
#------------------------------------------------------------------------------------------------------------------------
test_erythropoeisis.gata2 <- function()
{
   message(sprintf("--- test_erythropoeisis.gata2"))

   tms.gata2 <- TrenaMultiScore(tpe, "GATA2");
   getGeneHancerRegion(tms.gata2)   # 190kb
   findOpenChromatin(tms.gata2)
   findFimoTFBS(tms.gata2, fimo.threshold=1e-5)
   scoreMotifHitsForConservation(tms.gata2)
   scoreMotifHitsForGeneHancer(tms.gata2)
   addDistanceToTSS(tms.gata2)

   mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28")
   addGeneExpressionCorrelations(tms.gata2, mtx)
   addGenicAnnotations(tms.gata2)
   addChIP(tms.gata2)

   tbl <- getMultiScoreTable(tms.gata2)
   table(tbl$chip)
   tbl$cor[which(is.na(tbl$cor))] <- 0
   tbl$motifScore <- round(-log10(tbl$p.value), 2)

   dim(tbl)
   cor.spread <- fivenum(abs(tbl$cor))
   cor.threshold <- cor.spread[3]
   #cor.threshold <- 0.7 * cor.spread[5]
   tbl <- subset(tbl, abs(cor) >= cor.threshold)
   dim(tbl)

   tfoi <- unique(tbl$tf)
   length(tfoi)
   tblc <- data.frame(row.names=tfoi)
   #tblc$absCor <- as.numeric(lapply(tfoi, function(TF) mean(abs(subset(tbl, tf==TF)$cor))))
   #tblc$cor <- as.numeric(lapply(tfoi, function(TF) mean(subset(tbl, tf==TF)$cor)))
   tblc$gh <-  round(as.numeric(lapply(tfoi, function(TF) mean(subset(tbl, tf==TF)$gh))),0)
   tblc$count <- as.numeric(lapply(tfoi, function(TF) nrow(subset(tbl, tf==TF))))
   tblc$chip <- round(as.numeric(lapply(tfoi, function(TF) mean(subset(tbl, tf==TF)$chip))),0)
   tblc$fimo <- round(as.numeric(lapply(tfoi, function(TF) mean(subset(tbl, tf==TF)$motifScore))),2)
   tblc$phast <- round(as.numeric(lapply(tfoi,
           function(TF) mean(as.numeric(as.matrix(subset(tbl, tf==TF)[, c("phast7", "phast30", "phast100")]))))),2)
   #tblc$tss   <- as.numeric(lapply(tfoi, function(TF) mean(abs(subset(tbl, tf==TF)$tss))))

   mtx <- as.matrix(tblc)

   pca <- prcomp(mtx, scale=TRUE)
   library(factoextra)
   fviz_eig(pca)

   fviz_pca_ind(pca,
                col.ind = "cos2", # Color by the quality of representation
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = FALSE     # Avoid text overlapping
                )

   fviz_pca_var(pca,
                col.var = "contrib", # Color by contributions to the PC
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE     # Avoid text overlapping
                )

   fviz_pca_biplot(pca, repel = TRUE,
                   col.var = "red",     # #2E9FDF", # Variables color
                   col.ind = "#696969"  # Individuals color
                   )


} # test_erythropoeisis.gata2
#------------------------------------------------------------------------------------------------------------------------
test_AD.trem2 <- function()
{
   message(sprintf("--- test_AD.trems"))

   tms.trem2 <- TrenaMultiScore(tpad, "TREM2");
   getGeneHancerRegion(tms.trem2)
   findOpenChromatin(tms.trem2)
   findFimoTFBS(tms.trem2, fimo.threshold=1e-3)
   addChIP(tms.trem2)
   scoreMotifHitsForConservation(tms.trem2)
   scoreMotifHitsForGeneHancer(tms.trem2)
   addGenicAnnotations(tms.trem2)
   addDistanceToTSS(tms.trem2)
   mtx <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/mayo.tcx.16969x262.covariateCorrection.log+scale.RData"))
   addGeneExpressionCorrelations(tms.trem2, mtx)

   tbl.raw <- getMultiScoreTable(tms.trem2)
   dim(tbl.raw)
   tbl <- subset(tbl.raw, chip | p.value <= 1e-4)
   dim(tbl)
   length(unique(tbl$tf))

   table(tbl$chip)

   tbl$cor[which(is.na(tbl$cor))] <- 0
   tbl$motifScore <- round(-log10(tbl$p.value), 2)

   dim(tbl)
   hist(abs(tbl$cor))
   tbl <- subset(tbl, abs(cor) >= 0.1)
   dim(tbl)

   tfoi <- unique(tbl$tf)
   length(tfoi)
   tblc <- data.frame(row.names=tfoi)
   tblc$cor <- as.numeric(lapply(tfoi, function(TF) mean(abs(subset(tbl, tf==TF)$cor))))
   tblc$gh <-  as.numeric(lapply(tfoi, function(TF) mean(subset(tbl, tf==TF)$gh)))
   tblc$count <- as.numeric(lapply(tfoi, function(TF) nrow(subset(tbl, tf==TF))))
   #tblc$chip <- as.numeric(lapply(tfoi, function(TF) subset(tbl, tf==TF)$chip))
   tblc$chip <- as.integer(lapply(tfoi, function(TF) sum(subset(tbl, tf==TF)$chip)))
   tblc$fimo <- as.numeric(lapply(tfoi, function(TF) mean(subset(tbl, tf==TF)$motifScore)))
   tblc$phast <- as.numeric(lapply(tfoi, function(TF) mean(as.numeric(as.matrix(subset(tbl, tf==TF)[, c("phast7", "phast30", "phast100")])))))
   #tblc$tss   <- as.numeric(lapply(tfoi, function(TF) mean(abs(subset(tbl, tf==TF)$tss))))


   tblc$count[tblc$count > 10] <- 10
   tblc$chip[tblc$chip > 10] <- 10

   tblc.trimmed <- subset(tblc, cor > 0.32 | chip > 0 | fimo > 5)
   mtx <- as.matrix(tblc.trimmed[, c("chip", "count", "cor")])

   pca <- prcomp(mtx, scale=TRUE)
   fviz_eig(pca)

   fviz_pca_ind(pca,
                col.ind = "cos2", # Color by the quality of representation
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE     # Avoid text overlapping
                )

  fviz_pca_var(pca,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
               )

   fviz_pca_biplot(pca, repel = FALSE,
                   col.var = "RED", # Variables color
                   col.ind = "#696969"  # Individuals color
                   )

    # m2 <- mtx[rownames(subset(as.data.frame(pca$x), PC1 > 0 & PC2 > 0)),]
    # m2[order(m2[,"cor"], decreasing=TRUE),]
    #       chip count  cor
    # IKZF1    4    13 0.57
    # CEBPA   13    18 0.53
    # ESR2     0    23 0.51
    # ESR1    21    63 0.49
    # SPI1    67    78 0.47
    # FLI1    39    51 0.46
    # RUNX1   10    14 0.46
    # SP2      0    26 0.45
    # SP3      0    33 0.43
    # ERG     18    23 0.43
    # ELF1    18    27 0.41
    # ETS1     5    14 0.41
    # GABPA   17    22 0.38
    # PRDM1    3    20 0.38
    # RFX5     4    19 0.37
    # SP1     29   129 0.33
    # CREB1    5    15 0.32


} # test_AD.trem2
#------------------------------------------------------------------------------------------------------------------------

if(!interactive())
   runTests()
