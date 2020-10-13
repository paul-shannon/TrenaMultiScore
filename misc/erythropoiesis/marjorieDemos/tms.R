library(TrenaMultiScore)
library(TrenaProjectErythropoiesis)
library(RUnit)
library(factoextra)
library(RSQLite)
library(ghdb)
library(httr)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("ghdb"))
    ghdb <- GeneHancerDB()
if(!exists("tbl.atac"))
   tbl.atac <- get(load("~/github/TrenaProjectErythropoiesis/misc/multiScore/brandAtacSeqCuration/tbl.atac.fp.RData"))
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tmse")) {
   message(sprintf("--- creating instance of TrenaMultiScore"))
   tpe <- TrenaProjectErythropoiesis()
   tmse <- TrenaMultiScore(tpe, "TAL1");
   }

# from marjorie (email 29 aug 2020) use HSC/MPP for Day0 and CMP for Day2
if(!exists("tbl.atac.corces")){
    tbl.corces.full <- get(load("~/github/TrenaProjectErythropoiesis/prep/import/buenrostro/GSE74912-atacSeq.RData"))
    hsc.cols <- grep("\\.HSC", colnames(tbl.corces.full), ignore.case=TRUE)
    mpp.cols <- grep("MPP", colnames(tbl.corces.full), ignore.case=TRUE)
    cmp.cols <- grep("CMP", colnames(tbl.corces.full), ignore.case=TRUE)
    tbl.corces <- tbl.corces.full[, c(1:3, hsc.cols, mpp.cols, cmp.cols)]
    dim(tbl.corces)
    colnames(tbl.corces) <- c("chrom", "start", "end", "day0.1", "day0.2", "day0.3", "day0.4", "day2.1", "day2.2")
    tbl.atac.corces <- tbl.corces
    tbl.atac.corces.day0.quiet <- subset(tbl.atac.corces, day0.1 == 0 & day0.2 == 0 & day0.3==0 & day0.4 == 0 & day2.1 > 0)
    tbl.atac.corces.day0.quiet <- subset(tbl.atac.corces, day0.1 == 0 & day0.2 == 0 & (day2.1 > 0 | day2.2 > 0))
    tbl.atac.corces.day0.quiet <- subset(tbl.atac.corces, (day0.1 + day0.2 + day0.3 + day0.4) < 3 & (day2.1 + day2.2 > 6))
    }

#------------------------------------------------------------------------------------------------------------------------
# genes with zero-to-high expression between day 0 and day 2, and open chromatin at day 2
findCandidates <- function()
{
   expected <- c("brandLabDifferentiationTimeCourse-16173x28", "brandLabDifferentiationTimeCourse-27171x28")
   checkTrue(all(expected %in% getExpressionMatrixNames(tpe)))
   mtx <- getExpressionMatrix(tpe, expected[1])

   cutoff.1 <- 0.5
   cutoff.2 <- 1

   goi.0 <- names(which(sapply(rownames(mtx), function(geneName) all(mtx[geneName, 1:2] < cutoff.1))))
   print(length(goi.0))
   goi.2 <- names(which(sapply(rownames(mtx), function(geneName) all(mtx[geneName, 3:4] > cutoff.2))))
   print(length(goi.2))

   length(intersect(goi.0, goi.2))
   goi <-intersect(goi.0, goi.2)
   length(goi)

   tissues <- "Common myeloid progenitor CD34+"
   tissues <- "all"
   tbls.gh <- lapply(goi, function(gene) retrieveEnhancersFromDatabase(ghdb, gene, tissues=tissues))

   tbls.ghp <- lapply(goi, function(gene){
       tbl.gh <- retrieveEnhancersFromDatabase(ghdb, gene, tissues=tissues)
       if(nrow(tbl.gh) == 0)
         return(data.frame())
       subset(tbl.gh, combinedscore > 500)
       })

   length(tbls.ghp)
   names(tbls.ghp) <- goi
   tbls.ghp <- tbls.ghp[names(which(lapply(tbls.ghp, nrow) > 0))]
   length(tbls.ghp)

   # displayTrack(igv, DataFrameQuantitativeTrack("GH.CCL1", tbls.gh[["CCL1"]][, c("chrom", "start", "end", "combinedscore")], color="random", autoscale=FALSE, min=0,max=50))

   ghAndATAC <- function(tbl.gh){
       tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.gh[1, c("chrom", "start", "end")]),
                                            GRanges(tbl.atac.corces.day0.quiet[,  c("chrom", "start", "end")]),
                                            type="any"))
       if(nrow(tbl.ov) > 0)
           return(tbl.atac.corces.day0.quiet[tbl.ov$subjectHits,])
       return(data.frame())
       } # ghAndATAC

   tbls.ghp.atacPattern <- lapply(tbls.ghp, ghAndATAC)
   tbls.ghp.atacPattern <- tbls.ghp.atacPattern[names(which(lapply(tbls.ghp.atacPattern, nrow) > 0))]
   length(tbls.ghp.atacPattern)

   tbl.ccl1 <- ghAndATAC(tbls.gh[["CCL1"]])
   displayTrack(igv, DataFrameQuantitativeTrack("CCL1-atac-HSC.r1", tbl.ccl1[, c(1:3, 4)], color="random", autoscale=FALSE, min=0, max=20))
   displayTrack(igv, DataFrameQuantitativeTrack("CCL1-atac-HSC.r2", tbl.ccl1[, c(1:3, 5)], color="random", autoscale=FALSE, min=0, max=20))
   displayTrack(igv, DataFrameQuantitativeTrack("CCL1-atac-MPP.r1", tbl.ccl1[, c(1:3, 6)], color="random", autoscale=FALSE, min=0, max=20))
   displayTrack(igv, DataFrameQuantitativeTrack("CCL1-atac-MPP.r2", tbl.ccl1[, c(1:3, 7)], color="random", autoscale=FALSE, min=0, max=20))
   displayTrack(igv, DataFrameQuantitativeTrack("CCL1-atac-CMP.r1", tbl.ccl1[, c(1:3, 8)], color="random", autoscale=FALSE, min=0, max=20))
   displayTrack(igv, DataFrameQuantitativeTrack("CCL1-atac-CMP.r2", tbl.ccl1[, c(1:3, 9)], color="random", autoscale=FALSE, min=0, max=20))

   genes.filtered <- names (tbls.ghp.atacPattern)
   uri <- sprintf("http://localhost:8000/goEnrich")
   body.jsonString <- sprintf('%s', toJSON(list(geneSymbols=genes.filtered)))

   r <- POST(uri, body=body.jsonString)

      #sprintf('{"geneSymbols": "%s"}', goi.string))
   tbl <- fromJSON(content(r)[[1]])
   dim(tbl)
   wdth(1000)
   head(tbl, n=10)


} # findCandidates
#------------------------------------------------------------------------------------------------------------------------
build.model <- function(targetGene, fimoThresholdAsNegativeExponent=5)
{
   printf("=========== building model for %s, fimoThreshold: %f", targetGene,
          fimoThresholdAsNegativeExponent)

   results.subDirectory <- sprintf("fimo%d", fimoThresholdAsNegativeExponent)
   filename <- sprintf("%s.RData", targetGene)
   if(!file.exists(results.subDirectory))
      dir.create(results.subDirectory)

   tms.tg <- TrenaMultiScore(tpe, targetGene);
   getGeneHancerRegion(tms.tg)

   findOpenChromatin(tms.tg)
   if(nrow(getOpenChromatin(tms.tg)) == 0){
      message(sprintf("no open chromatin for %s, bailing out, saving empty model", targetGene))
      tbl <- data.frame()
      save(tbl, file=file.path(results.subDirectory, filename))
      return(data.frame())
      }
   fimoThreshold <- 10^(-fimoThresholdAsNegativeExponent)
   findFimoTFBS(tms.tg, fimo.threshold=fimoThreshold)
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
   tbl$targetGene <- targetGene

   dim(tbl)
   filePath <- file.path(results.subDirectory, filename)
   save(tbl, file=filePath)
   message(sprintf("saving %d rows model for %s to %s", nrow(tbl), targetGene, filePath))
   invisible(tbl)

} # build.model
#------------------------------------------------------------------------------------------------------------------------
# goi <- function()
# {
#    tbl <- read.table("G2vsG3GenesUP.txt", sep="\t", as.is=TRUE, header=TRUE, nrow=100)
#    dim(tbl)
#    genes <- tbl$GeneName
#    additional.genes <- read.table(file="additionalGenes.txt", stringsAsFactors=FALSE)$V1
#    length(genes)
#    length(additional.genes)
#    all.goi <- sort(unique(c(genes, additional.genes)))
#    length(all.goi)
#    library(org.Hs.eg.db)
#    tbl.ref <- select(org.Hs.eg.db, key=all.goi, columns=c("SYMBOL", "ENTREZID"), keytype="SYMBOL")
#    successes <- tbl.ref[which(!is.na(tbl.ref$ENTREZID)),]$SYMBOL
#    failures <- tbl.ref[which(is.na(tbl.ref$ENTREZID)),]$SYMBOL
#    tbl.ref2 <- select(org.Hs.eg.db, key=failures, columns=c("SYMBOL", "ALIAS", "ENTREZID"), keytype="ALIAS")
#    successes.round2 <- tbl.ref2[which(!is.na(tbl.ref2$SYMBOL)), "SYMBOL"]
#    goi.final <- sort(unique(c(successes, successes.round2)))
#    length(goi.final)  # 112
#
#    return(goi.final)
#
#} # goi
#------------------------------------------------------------------------------------------------------------------------
buildAll <- function(goi, fimoThresholdAsNegativeExponent)
{
  if(!exists("haney.erythropoiesis.tfs"))
     source("~/github/regulatoryGenomePaper/demos/common.R")

  #tfs.oi <- c("GATA1", "GATA2", "FLI1", "SPI1")
  #tfs.oi <- goi()[24:112]
  #tfs.oi <- goi()
  printf("running tms on %d genes", length(goi))

  f <- function(targetGene){
      tryCatch({
         fimo.threshold <- fimoThresholdAsNegativeExponent
         results.file.path <- file.path(sprintf("./fimo%d", fimo.threshold),
                                        sprintf("%s.RData", targetGene))
         if(file.exists(results.file.path))
             printf("results file exists, skipping: %s", results.file.path)
         else
             build.model(targetGene, fimo.threshold)
         },
      error = function(e){
         print("tms error")
         system(sprintf("touch %s", results.file.path))
         print(targetGene)
         print(e)
         })
      } # f

  tbls.all <- lapply(goi, f)
  names(tbls.all) <- goi

} # buildAll
#------------------------------------------------------------------------------------------------------------------------
collectResults <- function(directory, outfile)
{
   full.path.to.directory <- sprintf("~/github/TrenaMultiScore/misc/erythropoiesis/marjorieDemos/%s", directory)
   rdata.files <- list.files(path=full.path.to.directory, pattern=".RData")
   printf("RData file count: %d", length(rdata.files))

   tbls <- lapply(head(rdata.files, n=-1), function(file){
       full.path <- file.path(directory, file)
       if(file.size(full.path) > 1000){
           tbl <- get(load(full.path))
           printf("%30s: %d rows, %d cols", file, nrow(tbl), ncol(tbl))
           return(tbl)
       }})

   tbl.all <- do.call(rbind, tbls)
   printf("about to save combined table (%d rows)  to %s", nrow(tbl.all), outfile)
   if(!is.null(outfile))
      save(tbl.all, file=outfile)
   invisible(tbl.all)

} # collectResults
#------------------------------------------------------------------------------------------------------------------------
to.sqlite <- function(tbl, sqlite.filename)
{
   db <- dbConnect(SQLite(), sqlite.filename, synchronous=NULL)

   system.time(dbWriteTable(db, name="tms", value=tbl, overwrite=TRUE))  # less than 30 seconds
   dbDisconnect(db)

} # to.sqlite
#------------------------------------------------------------------------------------------------------------------------
