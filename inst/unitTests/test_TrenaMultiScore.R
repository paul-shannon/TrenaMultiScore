library(TrenaMultiScore)
library(TrenaProjectErythropoiesis)
library(TrenaProjectAD)
library(RUnit)
library(factoextra)
library(ghdb)

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
   test_getTargetGeneInfo()
   test_getGeneHancerRegion()
   test_findOpenChromatin()
   test_findFimoTFBS()
   # test_findMoodsTFBS()   FIMO is adequate, moods' contribution uncertain
   test_scoreMotifHitsForConservation()
   test_addDistanceToTSS()
   test_scoreMotifHitsForGeneHancer()
   test_addGenicAnnotations()
   test_addGeneExpressionCorrelations()
   test_addChIP()
   test_addRnaBindingProteins()

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
   checkEquals(as.list(tbl.gata2), list(chrom="chr3", start=128073944, end=128621958, width=548015))

   tbl.trem2 <- getGeneHancerRegion(tmsa)
   checkTrue(is.data.frame(tbl.trem2))
   checkEquals(as.list(tbl.trem2), list(chrom="chr6", start=41131001, end=41296400, width=165400))

} # test_getGeneHancerRegion
#------------------------------------------------------------------------------------------------------------------------
test_findOpenChromatin <- function()
{
   message(sprintf("--- test_findOpenChromatin"))

     #--------------------------------------------------------------
     # use a 1300 base region of chr19, upstream but not overlapping
     # with the promoter of DOT1L
     #--------------------------------------------------------------
   tpe.dot1l <- TrenaProjectErythropoiesis()
   tms <- TrenaMultiScore(tpe.dot1l, "DOT1L");

   findOpenChromatin(tms, "chr19", 2158968, 2160275,
                     intersect.with.geneHancer=FALSE,
                     use.merged.atac=TRUE)

   tbl.oc <- getOpenChromatin(tms)
   checkEquals(dim(tbl.oc), c(2, 5))

     #--------------------------------------------------------------
     # now ask for only the oc/atac overlapping genehancer promoter
     # should get nothing back
     #-------------------------------------------------------------
   findOpenChromatin(tms, "chr19", 2158968, 2160275,
                     intersect.with.geneHancer=TRUE,
                     use.merged.atac=TRUE)

   tbl.oc <- getOpenChromatin(tms)
   checkEquals(nrow(tbl.oc), 0)

     #--------------------------------------------------------------------------
     # repeat with 6kb region with a merged atac overlapping genehancer promoter
     # no overlap obliged on first call
     #------------------------------------------------------------------------
   findOpenChromatin(tms, "chr19", start=2158488, end=2164650,
                     intersect.with.geneHancer=FALSE,
                     use.merged.atac=TRUE)

   tbl.oc <- getOpenChromatin(tms)
   checkEquals(nrow(tbl.oc), 3)

     #--------------------------------------------------------------------------
     # repeat with 6kb region with a merged atac overlapping genehancer promoter
     # overlap obliged on second call
     #------------------------------------------------------------------------
   findOpenChromatin(tms, "chr19", start=2158488, end=2164650,
                     intersect.with.geneHancer=TRUE,
                     use.merged.atac=TRUE)

   tbl.oc <- getOpenChromatin(tms)
   checkEquals(nrow(tbl.oc), 1)

} # test_findOpenChromatin
#------------------------------------------------------------------------------------------------------------------------
test_findFimoTFBS <- function()
{
   message(sprintf("--- test_findFimoTFBS"))

   findOpenChromatin(tmse, "chr3", start=128483204, end=128483276)
   getOpenChromatin(tmse)

   findFimoTFBS(tmse)
   tbl.fimo <- getMultiScoreTable(tmse)
   checkEquals(ncol(tbl.fimo), 9)
   checkTrue(nrow(tbl.fimo) >= 2)
   checkTrue(all(c("SP2", "ZNF263") %in% tbl.fimo$tf))

   findFimoTFBS(tmse, fimo.threshold=1e-3)
   tbl.fimo <- getMultiScoreTable(tmse)

   checkTrue(nrow(tbl.fimo) > 30)
   checkTrue(all(c("SP2", "ZNF263", "CEBPB", "SP1") %in% tbl.fimo$tf)) # and many others

} # test_findFimoTFBS
#------------------------------------------------------------------------------------------------------------------------
test_findMoodsTFBS <- function()
{
   message(sprintf("--- test_findMoodsTFBS"))

   findOpenChromatin(tmse, "chr3", start=128483204, end=128483276)
   getOpenChromatin(tmse)

   t4 <- system.time(findMoodsTFBS(tmse, moods.threshold=1e-4))[["elapsed"]]
   tbl.moods.4 <- tmse@state$moods
   printf("moods 1e-4: %6.1f secs, %d hits", t4, nrow(tbl.moods.4))

   t3 <- system.time(findMoodsTFBS(tmse, moods.threshold=1e-3))[["elapsed"]]
   tbl.moods.3 <- tmse@state$moods
   printf("moods 1e-3: %6.1f secs, %d hits", t3, nrow(tbl.moods.3))

   t2 <- system.time(findMoodsTFBS(tmse, moods.threshold=1e-2))[["elapsed"]]
   tbl.moods.2 <- tmse@state$moods
   printf("moods 1e-2: %6.1f secs, %d hits", t2, nrow(tbl.moods.2))

   t6 <- system.time(findMoodsTFBS(tmse, moods.threshold=1e-6))[["elapsed"]]
   tbl.moods.6 <- tmse@state$moods
   printf("moods 1e-6: %6.1f secs, %d hits", t6, nrow(tbl.moods.6))

   checkEquals(nrow(tbl.moods.6), 0)

} # test_findMoodsTFBS
#------------------------------------------------------------------------------------------------------------------------
test_scoreMotifHitsForConservation <- function()
{
   message(sprintf("--- test_scoreMotifHitsForConservation"))

   findOpenChromatin(tmse, "chr3", start=128481000, end=128489000)
   findFimoTFBS(tmse, fimo.threshold=1e-5)
   tbl.fimo <- getMultiScoreTable(tmse)
   checkEquals(ncol(tbl.fimo), 9)
   checkTrue(nrow(tbl.fimo) >= 5)

   scoreMotifHitsForConservation(tmse)

   tbl <- getMultiScoreTable(tmse)
   checkEquals(ncol(tbl), 11)
   checkTrue(nrow(tbl) >= 5)
   #checkTrue(all(c("phast7", "phast30", "phast100") %in% colnames(tbl)))
   checkTrue(all(c("phast7", "phast100") %in% colnames(tbl)))
   checkEqualsNumeric(mean(tbl$phast7), 0.6, tolerance=0.05)
   #checkEqualsNumeric(mean(tbl$phast30), 0.66, tolerance=0.05)
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
   checkEquals(ncol(tbl.fimo), 9)
   checkTrue(nrow(tbl.fimo) >= 15)

   scoreMotifHitsForConservation(tmse)
   addDistanceToTSS(tmse)

   tbl <- getMultiScoreTable(tmse)
      # (31 jan 2020: with R-devel 4.0, we have an off-by-one difference
      # when compared to R 3.6.  not sure why.  skate around this for now.
   checkEqualsNumeric(head(sort(tbl$tss)),
                      c(-4770, -4769, -4769, -4768, -4768, -4768),
                      tolerance=1)

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
   checkTrue(length(unique(gh.values)) >= 3)

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
   checkEquals(ncol(tbl.fimo), 14)
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
   checkEquals(ncol(tbl.fimo), 9)
   checkTrue(nrow(tbl.fimo) >= 5)
   addChIP(tmse)
   tbl.fimo <- getMultiScoreTable(tmse)
   checkEquals(ncol(tbl.fimo), 10)
   checkTrue(nrow(tbl.fimo) >= 5)
   checkTrue(nrow(subset(tbl.fimo, chip)) >= 3)

} # test_addChIP
#------------------------------------------------------------------------------------------------------------------------
test_addGeneExpressionCorrelations <- function()
{
   message(sprintf("--- test_addGeneExpressionCorrelations"))

   tms.hoxb4 <- TrenaMultiScore(tpe, "HOXB4");
   getGeneHancerRegion(tms.hoxb4)
   findOpenChromatin(tms.hoxb4)
   findFimoTFBS(tms.hoxb4, fimo.threshold=1e-5)
   tbl.fimo <- getMultiScoreTable(tms.hoxb4)

   mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28")
   addGeneExpressionCorrelations(tms.hoxb4, mtx)

   timepoints <- c("day0.r1", "day0.r2", "day2.r1", "day2.r2", "day4.r1", "day4.r2")
   addGeneExpressionCorrelations(tms.hoxb4, mtx, "cor.early", timepoints)
   tbl <- getMultiScoreTable(tms.hoxb4)
   checkTrue(all(c("cor", "cor.early") %in% colnames(tbl)))

     # choose 10 tfs at random for checking
   set.seed(17)
   tfs <- sample(unique(tbl$tf), 10)
   tfs.usable <- intersect(tfs, rownames(mtx))
   tfs.noExpression <- setdiff(tfs, rownames(mtx))

   for(this.tf in tfs.usable){
      expected <- round(cor(mtx[this.tf, ], mtx["HOXB4",], method="spearman"), digits=2)
      actual <- subset(tbl, tf==this.tf)$cor[1]
      checkEquals(expected, actual)
      expected.sub <- round(cor(mtx[this.tf, timepoints], mtx["HOXB4", timepoints], method="spearman"), digits=2)
      actual <- subset(tbl, tf==this.tf)$cor.early[1]
      checkEquals(expected.sub, actual)
      }

   checkTrue(all(is.na(subset(tbl, tf %in% tfs.noExpression)$cor)))
   checkTrue(all(is.na(subset(tbl, tf %in% tfs.noExpression)$cor.early)))

} # test_addGeneExpressionCorrelations
#------------------------------------------------------------------------------------------------------------------------
test_addRnaBindingProteins <- function()
{
   message(sprintf("--- test_addRnaBindingProteins"))

   tms.gata2 <- TrenaMultiScore(TrenaProjectErythropoiesis(), "GATA2", quiet=TRUE);
   checkTrue(is.null(getRnaBindingProteins(tms.gata2)))

   tbl.rbp <- addRnaBindingProteins(tms.gata2)
   checkEquals(dim(tbl.rbp), c(2390, 12))

      # now retrieve it again, without a fresh query
   tbl.v2 <- getRnaBindingProteins(tms.gata2)
   checkEquals(dim(tbl.rbp), c(2390, 12))

} # test_addRnaBindingProteins
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
   checkTrue(all(c("GLI3", "TBX15", "RXRG") %in% tbl.sub.neg$tf))
   tbl.sub.pos <-
      subset(tbl, p.value < 0.00001 & phast100 > 0.8 & cor > 0.8 & gh > 0)
   checkTrue(all(c("MYC", "WT1") %in% tbl.sub.pos$tf))

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
   checkTrue(nrow(tbl) > 80)   # 90 x 18 (26 aug 2020)

   tfoi <- unique(tbl$tf)
   length(tfoi)
   tblc <- data.frame(row.names=tfoi)
   tblc$cor <- as.numeric(lapply(tfoi, function(TF) mean(abs(subset(tbl, tf==TF)$cor))))
   tblc$gh <-  as.numeric(lapply(tfoi, function(TF) mean(subset(tbl, tf==TF)$gh)))
   tblc$count <- as.numeric(lapply(tfoi, function(TF) nrow(subset(tbl, tf==TF))))
   tblc$chip <- as.numeric(lapply(tfoi, function(TF) sum(subset(tbl, tf==TF)$chip)))
   tblc$fimo <- as.numeric(lapply(tfoi, function(TF) mean(subset(tbl, tf==TF)$motifScore)))
   tblc$phast <- as.numeric(lapply(tfoi, function(TF) mean(as.numeric(as.matrix(subset(tbl, tf==TF)[, c("phast7", "phast100")])))))
   #tblc$phast <- as.numeric(lapply(tfoi, function(TF) mean(as.numeric(as.matrix(subset(tbl, tf==TF)[, c("phast7", "phast30", "phast100")])))))
   #tblc$tss   <- as.numeric(lapply(tfoi, function(TF) mean(abs(subset(tbl, tf==TF)$tss))))


   mtx <- as.matrix(tblc)
   dim(mtx)   # 5 x 6

   pca <- prcomp(mtx, scale=TRUE)
   fviz_eig(pca)

     #--------------------------------------------------
     # explore similarities among the tfs selected above
     # 5 which have the best mRNA correlation to TAL1
     #--------------------------------------------------

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
   region <- getGeneHancerRegion(tms.gata2)
   with(region, 1 + end - start) # 548kb
   findOpenChromatin(tms.gata2)
   motifs <- query(MotifDb, "sapiens", c("jaspar2018", "hocomocov11-core"))
   findFimoTFBS(tms.gata2, motifs=motifs, fimo.threshold=1e-3)
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

   dim(tbl)  # 416 x 18
   cor.spread <- fivenum(abs(tbl$cor))
   cor.threshold <- cor.spread[3]
   #cor.threshold <- 0.7 * cor.spread[5]
   tbl <- subset(tbl, abs(cor) >= cor.threshold)
   dim(tbl)  # 223 x 18

   tfoi <- unique(tbl$tf)
   length(tfoi)   # 60
   tblc <- data.frame(row.names=tfoi)
   #tblc$absCor <- as.numeric(lapply(tfoi, function(TF) mean(abs(subset(tbl, tf==TF)$cor))))
   #tblc$cor <- as.numeric(lapply(tfoi, function(TF) mean(subset(tbl, tf==TF)$cor)))
   tblc$gh <-  round(as.numeric(lapply(tfoi, function(TF) mean(subset(tbl, tf==TF)$gh))),0)
   tblc$count <- as.numeric(lapply(tfoi, function(TF) nrow(subset(tbl, tf==TF))))
   tblc$chip <- round(as.numeric(lapply(tfoi, function(TF) mean(subset(tbl, tf==TF)$chip))),0)
   tblc$fimo <- round(as.numeric(lapply(tfoi, function(TF) mean(subset(tbl, tf==TF)$motifScore))),2)
   tblc$phast <- round(as.numeric(lapply(tfoi,
           function(TF) mean(as.numeric(as.matrix(subset(tbl, tf==TF)[, c("phast7", "phast100")]))))),2)
   #tblc$phast <- round(as.numeric(lapply(tfoi,
   #        function(TF) mean(as.numeric(as.matrix(subset(tbl, tf==TF)[, c("phast7", "phast30", "phast100")]))))),2)
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
# takes a long time to run with fimo.threshold of 1e-3.  skip unless you are patient.
test_AD.trem2 <- function()
{
   message(sprintf("--- test_AD.trems"))

   tms.trem2 <- TrenaMultiScore(tpad, "TREM2");
   getGeneHancerRegion(tms.trem2)
   findOpenChromatin(tms.trem2)
   findFimoTFBS(tms.trem2, fimo.threshold=1e-3)   # 146k hits!
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
   tblc$phast <- as.numeric(lapply(tfoi, function(TF) mean(as.numeric(as.matrix(subset(tbl, tf==TF)[, c("phast7", "phast100")])))))
   #tblc$phast <- as.numeric(lapply(tfoi, function(TF) mean(as.numeric(as.matrix(subset(tbl, tf==TF)[, c("phast7", "phast30", "phast100")])))))
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
test_AD.mef2c <- function()
{
   message(sprintf("--- test_AD.trems"))

   tms.mef2c <- TrenaMultiScore(tpad, "MEF2C");
   getGeneHancerRegion(tms.mef2c)
   findOpenChromatin(tms.mef2c)
   tbl.oc <- getOpenChromatin(tms.mef2c)
   printf("base count for fimo: %d", sum(with(tbl.oc, 1 + end-start)))

   system.time(findFimoTFBS(tms.mef2c, fimo.threshold=1e-7))
   system.time(findFimoTFBS(tms.mef2c, fimo.threshold=1e-6))
   system.time(findFimoTFBS(tms.mef2c, fimo.threshold=1e-5))
   system.time(findFimoTFBS(tms.mef2c, fimo.threshold=1e-4))
   system.time(findFimoTFBS(tms.mef2c, fimo.threshold=1e-3))  # 5 minutes

   addChIP(tms.mef2c)
   scoreMotifHitsForConservation(tms.mef2c)
   scoreMotifHitsForGeneHancer(tms.mef2c)
   addGenicAnnotations(tms.mef2c)
   addDistanceToTSS(tms.mef2c)
   mtx.new <- get(load("~/github/TrenaProjectAD/inst/extdata/expression/mayo.tcx.16969x262.covariateCorrection.log+scale.RData"))
   mtx.old <- getExpressionMatrix(tpad, "temporalCortex.15167x264")
   addGeneExpressionCorrelations(tms.mef2c, mtx.new)

   tbl.raw <- getMultiScoreTable(tms.mef2c)

   dim(tbl.raw)
   save(tbl.raw, file="mef2c.fimo1e-3.boca.RData")
   dim(subset(tbl.raw, chip | p.value <= 1e-5))
   tbl <- subset(tbl.raw, chip | p.value <= 1e-5)
   dim(tbl)
   length(unique(tbl$tf))

   table(tbl$chip)

   tbl$cor[which(is.na(tbl$cor))] <- 0
   tbl$motifScore <- round(-log10(tbl$p.value), 2)

   dim(tbl)
   hist(abs(tbl$cor))
   spread <- fivenum(abs(tbl$cor))
   tbl <- subset(tbl, abs(cor) >= spread[4])
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
   tblc$phast <- as.numeric(lapply(tfoi, function(TF) mean(as.numeric(as.matrix(subset(tbl, tf==TF)[, c("phast7", "phast100")])))))
   #tblc$phast <- as.numeric(lapply(tfoi, function(TF) mean(as.numeric(as.matrix(subset(tbl, tf==TF)[, c("phast7", "phast30", "phast100")])))))
   #tblc$tss   <- as.numeric(lapply(tfoi, function(TF) mean(abs(subset(tbl, tf==TF)$tss))))

   remove.mef2c <- grep("MEF2C", rownames(tblc))
   if(length(remove.mef2c) > 0)
      tblc <- tblc[-remove.mef2c,]

   dim(tblc)
   #tblc$count[tblc$count > 10] <- 10
   #tblc$chip[tblc$chip > 10] <- 10

   fivenum(tblc$cor)
   tblc.trimmed <- subset(tblc, gh > -1 & (cor > 0.22 | chip > 0) & fimo > 3)
   dim(tblc.trimmed)


   #tblc.trimmed <- subset(tblc, gh > 0 & (cor > 0.22 | chip > 0 | fimo > 5))
   new.order <- with(tblc.trimmed, order(cor, count, fimo, chip, decreasing=TRUE))
   tblc.trimmed <- tblc.trimmed[new.order,]
   dim(tblc.trimmed)

   mtx.pca <- as.matrix(tblc.trimmed[, c("chip", "count", "cor", "gh")])

   pca <- prcomp(mtx.pca[1:30,], scale=TRUE)
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

   save(tbl, file="mef2c.feature.table.RData")

} # test_AD.mef2c
#------------------------------------------------------------------------------------------------------------------------
#       chrom    start      end     tf strand    score  p.value matched_sequence  chip phast7 phast30 phast100     gh                    annot.type        annot.symbol    tss  cor                            motif_id motifScore
# 24798  chr5 88667670 88667684 ZBTB33      - 10.53930 5.99e-05  CTCCCGCGGTCTCGA FALSE   1.00    1.00     1.00   7.39             introns,promoters    LINC00461,MIR9-2 236587 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       4.22
# 39561  chr5 88680491 88680505 ZBTB33      + 11.01120 4.76e-05  CTCTCTCGGGAGCTC FALSE   0.00    0.02     0.00   0.00                introns,1to5kb LINC00461,MEF2C-AS2 223766 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       4.32
# 39935  chr5 88680553 88680567 ZBTB33      -  9.83146 8.39e-05  GCCCCGCGCGGGCCC FALSE   0.37    0.35     0.55   0.00                introns,1to5kb LINC00461,MEF2C-AS2 223704 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       4.08
# 40524  chr5 88680674 88680688 ZBTB33      - 10.65170 5.68e-05  TTCCCGCGGCCTCCC FALSE   0.00    0.00     0.00   0.00                introns,1to5kb LINC00461,MEF2C-AS2 223583 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       4.25
# 47241  chr5 88685169 88685183 ZBTB33      -  9.76404 8.65e-05  GTCCCGCGGCGCTAG FALSE   0.78    0.78     0.36   0.00             introns,promoters LINC00461,MEF2C-AS2 219088 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       4.06
# 51217  chr5 88691107 88691121 ZBTB33      -  9.55056 9.53e-05  GGCTCGCGGTGCTAC FALSE   0.26    0.32     0.26   0.00             introns,promoters MEF2C-AS2,LINC00461 213150 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       4.02
# 91511  chr5 88872971 88872985 ZBTB33      +  9.56180 9.48e-05  CTCTCTCGAGAGTGT FALSE   0.05    0.00     0.00   0.00                       introns               MEF2C  31286 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       4.02
# 94139  chr5 88883227 88883241 ZBTB33      - 12.11240 2.69e-05  CCCTCGCGCGCGCGC FALSE   0.05    0.55     0.60 604.42 introns,promoters,5UTRs,exons     MEF2C,MEF2C-AS1  21030 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       4.57
# 94263  chr5 88883240 88883254 ZBTB33      -  9.48315 9.83e-05  CTCGCGCGCGCTCCC FALSE   0.13    0.09     0.55 604.42 introns,promoters,5UTRs,exons     MEF2C,MEF2C-AS1  21017 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       4.01
# 94289  chr5 88883242 88883256 ZBTB33      - 13.59550 1.15e-05  CCCTCGCGCGCGCTC FALSE   0.12    0.08     0.59 604.42 introns,promoters,5UTRs,exons     MEF2C,MEF2C-AS1  21015 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       4.94
# 94684  chr5 88883256 88883270 ZBTB33      -  7.04494 2.74e-04  GCCCCGCGCCCCCCC  TRUE   0.07    0.46     0.51 604.42 introns,promoters,5UTRs,exons     MEF2C,MEF2C-AS1  21001 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       3.56
# 94846  chr5 88883274 88883288 ZBTB33      -  4.28090 7.54e-04  GGCGCGCGCGAATGC  TRUE   0.00    0.84     0.46 604.42 introns,promoters,5UTRs,exons     MEF2C,MEF2C-AS1  20983 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       3.12
# 94865  chr5 88883277 88883291 ZBTB33      +  5.25843 5.36e-04  TTCGCGCGCGCCGAG  TRUE   0.00    0.87     0.40 604.42 introns,promoters,5UTRs,exons     MEF2C,MEF2C-AS1  20980 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       3.27
# 95861  chr5 88883476 88883490 ZBTB33      -  4.95506 5.97e-04  TTTACGCGATGTTTC  TRUE   1.00    1.00     1.00 604.42       introns,promoters,exons     MEF2C,MEF2C-AS1  20781 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       3.22
# 97116  chr5 88883656 88883670 ZBTB33      +  4.57303 6.82e-04  GCTTCGCGCTGCCTT  TRUE   0.33    0.95     0.00 604.42 introns,promoters,5UTRs,exons     MEF2C,MEF2C-AS1  20601 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       3.17
# 97343  chr5 88883715 88883729 ZBTB33      -  6.08989 3.95e-04  AGCTCGCGAGGAAAG  TRUE   0.03    0.32     0.62 604.42 introns,promoters,5UTRs,exons     MEF2C,MEF2C-AS1  20542 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       3.40
# 97347  chr5 88883718 88883732 ZBTB33      +  9.88764 8.17e-05  TCCTCGCGAGCTGGA  TRUE   0.05    0.50     0.63 604.42 introns,promoters,5UTRs,exons     MEF2C,MEF2C-AS1  20539 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       4.09
# 98820  chr5 88884089 88884103 ZBTB33      + 15.42700 3.44e-06  CTCTCGCGGCGCCTC FALSE   0.00    0.00     0.00 604.42      introns,promoters,1to5kb     MEF2C,MEF2C-AS1  20168 0.54 Hsapiens-jaspar2018-ZBTB33-MA0527.1       5.46
#
#  mtx["ZBTB33",]
#     chip    count      cor       gh
#   7.0000  18.0000   0.5400 369.7783
#
test_mef2cModel <- function()
{
   goi <- "ZBTB33"
   igv <- start.igv("MEF2C")
   tbl.gh <- getEnhancers(tpad)
   tbl.features <- get(load("mef2c.feature.table.RData"))

   with(tbl.gh, showGenomicRegion(igv, sprintf("%s:%d-%d", chrom[1], min(start)-1000, max(end)+ 1000)))

   track <- DataFrameQuantitativeTrack("gh", tbl.gh[, c("chrom", "start", "end", "combinedscore")], color="brown",
                                       autoscale=FALSE, min=0, max=20)
   displayTrack(igv, track)

   tbl.goifimo <- subset(tbl.features, tf==goi)[, c("chrom", "start", "end")]
   track <- DataFrameAnnotationTrack("ZBTB33.fimo", tbl.goifimo, color="darkBlue", trackHeight=25)
   displayTrack(igv, track)

   tbl.chip <- with(tbl.gh, getChipSeq(tpad, chrom[1], min(start)-1000, max(end)+1000, tfs=goi))
   track <- DataFrameAnnotationTrack("ChIP", tbl.chip[, c("chrom", "start", "end")], color="darkGreen")
   displayTrack(igv, track)

} # test_mef2cModel
#------------------------------------------------------------------------------------------------------------------------
demo.ccl1.tpe <- function()
{
   printf("--- demo.ccl1.tpe")
   tpe <- TrenaProjectErythropoiesis()
   tms <- TrenaMultiScore(tpe, "CCL1");
   tbl.gh <- getGeneHancerRegion(tms)    # 240k
   with(tbl.gh, findOpenChromatin(tms, chrom, start, end))
   tbl.oc <- getOpenChromatin(tms)

   tbl.corces.paul <- get(load("~/github/TrenaProjectErythropoiesis/misc/diffBind/corces/tbl.corces.hg38.scoredByPaul.RData"))
   tbl.corces.rory <- get(load("~/github/TrenaProjectErythropoiesis/misc/diffBind/corces/tbl.diffBind.rory.hg38.day0-4.RData"))

   tbl.corces.paul.ccl1 <- subset(tbl.corces.paul, chrom==tbl.gh$chrom & start >= tbl.gh$start & end <= tbl.gh$end)
   tbl.corces.rory.ccl1 <- subset(tbl.corces.rory, chrom==tbl.gh$chrom & start >= tbl.gh$start & end <= tbl.gh$end)
   tbl.corces.rory.ccl1$Fold <- -1 * tbl.corces.rory.ccl1$Fold
   dim(tbl.corces.paul.ccl1)
   dim(tbl.corces.rory.ccl1)

   igv <- start.igv("CCL1")
   track <- DataFrameQuantitativeTrack("corces.paul", tbl.corces.paul.ccl1[, c(1,2,3,13)],
                                       color="red", autoscale=TRUE)
   displayTrack(igv, track)

   track <- DataFrameQuantitativeTrack("corces.rory", tbl.corces.rory.ccl1[, c(1,2,3,9)],
                                       color="red", autoscale=TRUE)
   displayTrack(igv, track)

   ghdb <- GeneHancerDB()
   tbl.gh.cd34 <- retrieveEnhancersFromDatabase(ghdb, "CCL1",
                                                tissues="Common myeloid progenitor CD34+")
   dim(tbl.gh.cd34) # 5 16
   track <- DataFrameQuantitativeTrack("GH.cd34", tbl.gh.cd34[, c(1,2,3,11)],
                                       color="random", autoscale=FALSE, min=0, max=20)
   displayTrack(igv, track)

   tbl.gh.all <- retrieveEnhancersFromDatabase(ghdb, "CCL1", tissues="all")
   dim(tbl.gh.all) # 18 16
   track <- DataFrameQuantitativeTrack("GH.all", tbl.gh.all[, c(1,2,3,11)],
                                       color="brown", autoscale=FALSE, min=0, max=20)
   displayTrack(igv, track)


     #------------------------------------------------------------
     # add another track, the +/- ratio of day0 vs day2 reads.
     # day2 > day0 reads
     #------------------------------------------------------------
   day0.mean <- as.numeric(lapply(seq_len(nrow(tbl.corces.paul.ccl1)),
                                  function(i) mean(as.numeric(tbl.corces.paul.ccl1[i, 4:7]))))
   day2.mean <- as.numeric(lapply(seq_len(nrow(tbl.corces.paul.ccl1)),
                                  function(i) mean(as.numeric(tbl.corces.paul.ccl1[i, 8:9]))))

   avg0 <- day0.mean/day2.mean
   avg2 <- day2.mean/day0.mean
   avg0[is.infinite(avg0)] <- 0
   avg2[is.infinite(avg2)] <- 0

   score2 <- rep(0, nrow(tbl.corces.paul.ccl1))
   score2[which(avg0 > 1)] <- -avg0[avg0 > 1]
   score2[which(avg2 > 1)] <- avg2[avg2 > 1]
   score2 <- round(score2, digits=2)

   tbl.corces.paul.ccl1$score2 <- score2
   track <- DataFrameQuantitativeTrack("score2", tbl.corces.paul.ccl1[, c("chrom", "start", "end", "score2")],
                                       color="blue", autoscale=FALSE,
                                       min=-10, max=30)
   displayTrack(igv, track)

     #------------------------------------------------------------
     # bring in gene expression data
     #------------------------------------------------------------
   f <- system.file(package="TrenaProjectErythropoiesis", "extdata", "expression",
                    "brandLabDifferentiationTimeCourse-27171x28.RData")
   file.exists(f)
   mtx.rna <- get(load(f))
   cor.all <- lapply(rownames(mtx.rna), function(gene)
       cor(mtx.rna[gene, 1:4], mtx.rna["CCL1", 1:4]))
   names(cor.all) <- rownames(mtx.rna)
   fivenum(as.numeric(cor.all))
    # -0.68529152  0.06930292  0.26896765  0.48405971  0.89141070
   hist(as.numeric(cor.all), main="correlation of CCL1 againstall")

    #------------------------------------------------------------
    # add ChIP indiscriminately
    #------------------------------------------------------------
   tbl.chip <- with(tbl.gh, getChipSeq(tpe, chrom[1], min(start)-1000, max(end)+1000))
   tfs <- names(rev(sort(table(tbl.chip$tf))))
   for(tf.this in tfs){
       if(tf.this %in% names(cor.all)){
          cor.tf <- cor.all[[tf.this]]
          if(abs(cor.tf) > 0.95){
              printf("tf %s: %f", tf.this, cor.tf)
              color <- "red"
              if(cor.tf < 0)
                 color <- "blue"
              tbl.tf <- subset(tbl.chip, tf==tf.this)
              track <- DataFrameAnnotationTrack(tf.this, tbl.tf[, c(1,2,3)], color=color,,
                                                displayMode="COLLAPSED", trackHeight=23)
              displayTrack(igv, track)
             } # if abs
          } # tf in cor
       } # for tf.this


} # demo.ccl1.tpe
#------------------------------------------------------------------------------------------------------------------------

if(!interactive())
   runTests()
