library(TrenaMultiScore)
library(TrenaProjectErythropoiesis)
library(RUnit)
library(factoextra)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tmse")) {
   message(sprintf("--- creating instance of TrenaMultiScore"))
   tpe <- TrenaProjectErythropoiesis()
   tmse <- TrenaMultiScore(tpe, "TAL1");
   }
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
   # cor.threshold <- 0.9 * cor.spread[5]
   tbl <- subset(tbl, abs(cor) >= cor.threshold)
   dim(tbl)


   tbl.best <- subset(tbl, motifScore > 3 & phast100 > 0.5 & chip & abs(tss) < 3000)
   tfoi <- unique(tbl.best$tf)
   length(tfoi)

   tbl.trimmed <- subset(tbl, tf %in% tfoi)

   tblc <- data.frame(row.names=tfoi)
   tblc$cor <- as.numeric(lapply(tfoi, function(TF) mean(abs(subset(tbl.trimmed, tf==TF)$cor))))
   tblc$gh <-  as.numeric(lapply(tfoi, function(TF) mean(subset(tbl.trimmed, tf==TF)$gh)))
   tblc$count <- as.numeric(lapply(tfoi, function(TF) nrow(subset(tbl.trimmed, tf==TF))))
   tblc$chip <- as.numeric(lapply(tfoi, function(TF) sum(subset(tbl.trimmed, tf==TF)$chip)))
   tblc$fimo <- as.numeric(lapply(tfoi, function(TF) mean(subset(tbl.trimmed, tf==TF)$motifScore)))
   tblc$phast <- as.numeric(lapply(tfoi, function(TF) mean(as.numeric(as.matrix(subset(tbl.trimmed, tf==TF)[, c("phast7", "phast30", "phast100")])))))

   new.order <- with(tblc, order(cor, count, fimo, chip, decreasing=TRUE))
   tblc <- tblc[new.order,]

      #         cor        gh count chip     fimo     phast
      # TAL1   1.00 206.53667     3    3 3.646667 0.3777778
      # KLF1   0.93 323.40526    19    8 3.263158 0.2329825
      # NFE2   0.93 365.89600     5    3 4.178000 0.4646667
      # E2F4   0.87 193.97538    26   16 3.424231 0.3635897
      # GATA1  0.87 265.92429     7    7 4.004286 0.5033333
      # TFDP1  0.86 185.06024    82   24 3.877927 0.2915041
      # ZNF143 0.82 605.99000     3    1 3.403333 0.3600000
      # SP1    0.80 221.66307   212  119 3.848019 0.3230503
      # ZBTB7A 0.79 207.65000    15   10 3.605333 0.3571111
      # KLF13  0.70 324.69105    19    6 3.681053 0.2884211
      # IRF9   0.64  92.58429     7    1 3.331429 0.2285714
      # GATA3  0.60 207.55500    12    5 3.818333 0.5211111
      # FLI1   0.58 168.11200    15   10 3.832000 0.3800000
      # ETS1   0.55  62.00909    11    8 3.549091 0.3063636
      # JUNB   0.51 605.99000     4    2 3.497500 0.7516667


   mtx.pca <- as.matrix(tblc[, c("cor", "count", "chip", "fimo", "phast")])
   dim(mtx.pca)
   pca <- prcomp(mtx.pca, scale=TRUE)
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
build.model <- function(targetGene)
{
   tms.tg <- TrenaMultiScore(tpe, targetGene);
   getGeneHancerRegion(tms.tg)   # 190kb
   findOpenChromatin(tms.tg)
   findFimoTFBS(tms.tg, fimo.threshold=1e-2)
   scoreMotifHitsForConservation(tms.tg)
   scoreMotifHitsForGeneHancer(tms.tg)
   addDistanceToTSS(tms.tg)

   mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28")
   addGeneExpressionCorrelations(tms.tg, mtx)
   addGenicAnnotations(tms.tg)
   addChIP(tms.tg)

   tbl <- getMultiScoreTable(tms.tg)
   tbl$cor[which(is.na(tbl$cor))] <- 0
   tbl$motifScore <- round(-log10(tbl$p.value), 2)

   dim(tbl)
   cor.spread <- fivenum(abs(tbl$cor))
   cor.threshold <- cor.spread[4]
   # cor.threshold <- 0.9 * cor.spread[5]
   tbl <- subset(tbl, abs(cor) >= cor.threshold)
   dim(tbl)


   tbl.best <- subset(tbl, motifScore > 3 & (phast100 > 0.5 | chip) & abs(tss) < 5000)
   tfoi <- unique(tbl.best$tf)
   length(tfoi)

   tbl.trimmed <- subset(tbl, tf %in% tfoi)
   browser()

   tblc <- data.frame(row.names=tfoi)
   tblc$cor <- as.numeric(lapply(tfoi, function(TF) subset(tbl.trimmed, tf==TF)$cor[1]))
   tblc$gh <-  as.numeric(lapply(tfoi, function(TF) mean(subset(tbl.trimmed, tf==TF)$gh)))
   tblc$count <- as.numeric(lapply(tfoi, function(TF) nrow(subset(tbl.trimmed, tf==TF))))
   tblc$chip <- as.numeric(lapply(tfoi, function(TF) sum(subset(tbl.trimmed, tf==TF)$chip)))
   tblc$fimo <- as.numeric(lapply(tfoi, function(TF) mean(subset(tbl.trimmed, tf==TF)$motifScore)))
   tblc$phastAvg <- as.numeric(lapply(tfoi,
           function(TF) mean(as.numeric(as.matrix(subset(tbl.trimmed, tf==TF)[, c("phast7", "phast30", "phast100")])))))
   tblc$phast100 <- as.numeric(lapply(tfoi,
           function(TF) mean(as.numeric(as.matrix(subset(tbl.trimmed, tf==TF)[, "phast100"])))))

   tblc$targetGene <- targetGene
   new.order <- with(tblc, order(cor, count, fimo, chip, decreasing=TRUE))

   tblc[new.order,]
} # build.model
#------------------------------------------------------------------------------------------------------------------------
getOtherTBX15targets <- function()
{
  if(!exists("haney.erythropoiesis.tfs"))
     source("~/github/regulatoryGenomePaper/demos/common.R")

  haney.tfs <- haney.erythropoiesis.tfs()
  mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28")
  cors <- lapply(tfs, function(tf) cor(mtx[tf,], mtx["TBX15",]))
  names(cors) <- tfs
  tbl.cors <- data.frame(cor=as.numeric(cors))
  rownames(tbl.cors) <- tfs
  tf.oi <- rownames(subset(tbl.cors, abs(cor) > 0.8))

     # subset(tbl.cors, abs(cor) > 0.8)
     # CBFA2T3 -0.8680747
     # GATA2   -0.8722489
     # HHEX    -0.8560287
     # HOXB4   -0.8254034
     # SMARCC1 -0.8271527

  tbls.all <- lapply(tf.oi, function(targetGene) build.model(targetGene))
  names(tbls.all) <- tf.oi

  model.hhex <- build.model("HHEX")
  model.smarcc1 <- build.model("SMARCC1")
  model.zbtb33 <- build.model("ZBTB33")

} # getOtherTBX15targets
#------------------------------------------------------------------------------------------------------------------------
