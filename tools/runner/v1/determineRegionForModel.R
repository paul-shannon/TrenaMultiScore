source("tmsRunner.R")
mtx.tmp <- get(load("~/github/TrenaProjectErythropoiesis/viz/srm.vs.mrna/shinyapps.io/srm.rna.averaged.clean.RData"))
targetGenes <- rownames(mtx.tmp)
trenaProject <- TrenaProjectErythropoiesis()
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_determineRegionForModel()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
determineRegionForModel <- function(targetGene)
{
   coi <- c("tss", "strand", "chrom", "start", "end")
   geneRecognized <- setTargetGene(trenaProject, targetGene)
   if(!geneRecognized)
       return(data.frame(chrom=NA, start=NA, end=NA, width=0, stringsAsFactors=FALSE))

   tss.list <- as.list(getTranscriptsTable(trenaProject)[coi])
   tbl.tss <- data.frame(chrom=tss.list$chrom, start=tss.list$tss, end=tss.list$tss, stringsAsFactors=FALSE)
   tbl.fallback <- data.frame(chrom=tss.list$chrom, start=tss.list$tss-5000, end=tss.list$tss+5000, width=10001, stringsAsFactors=FALSE)

   ghdb <- GeneHancerDB()
   suppressWarnings(
       tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all") # "Common myeloid progenitor CD34+")
       )
   if(nrow(tbl.gh) == 0)
       return(tbl.fallback)

   tbl.gh <- subset(tbl.gh, combinedscore > 550)
   if(nrow(tbl.gh) == 0)
       return(tbl.fallback)

   tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.tss), GRanges(tbl.gh[, c("chrom", "start", "end")])))

   if(nrow(tbl.ov) == 0)
       return(tbl.fallback)

   tss.hit <- unique(tbl.ov$subjectHits)
   tbl.promoter <- tbl.gh[tss.hit, c("chrom", "start", "end")]
   tbl.promoter$width <- 1 + abs(tbl.promoter$end-tbl.promoter$start)
   tbl.promoter$gene <- targetGene
   tbl.promoter

} # determineRegionForModle
#------------------------------------------------------------------------------------------------------------------------
test_determineRegionForModel <- function()
{
   message(sprintf("--- test_determineRegionForModel"))

   goi <- c("BACH1", "BCL11A_XL_L", "CBP", "CEBPB", "CHD1")
   checkEquals(determineRegionForModel("BACH1")$width, 6335)
   checkEquals(determineRegionForModel("BCL11A_XL_L")$width, 0)
   checkEquals(determineRegionForModel("CBP")$width, 0)
   checkEquals(determineRegionForModel("CEBPB")$width, 10064)
   checkEquals(determineRegionForModel("CHD1")$width, 6481)
   checkEquals(determineRegionForModel("CHD3")$width, 7271)
   checkEquals(determineRegionForModel("CHD4")$width, 7408)
   checkEquals(determineRegionForModel("coREST")$width, 0)
   checkEquals(determineRegionForModel("CTBP2")$width, 17337)
   checkEquals(determineRegionForModel("CTCF")$width, 5845)
   checkEquals(determineRegionForModel("DNMT1")$width, 3372)
   checkEquals(determineRegionForModel("DOT1L")$width, 13653)
   checkEquals(determineRegionForModel("E12/E47")$width, 0)
   checkEquals(determineRegionForModel("E2F4")$width, 4079)
   failures <- c("ETO2", "GR", "HXB4", "MLL1", "MLL3", "MLL4", "NC2B", "PO2F1", "RPB1", "SET1B", "SETB1",
                 "SIR6", "SMCA4", "SMRC1", "SNF5", "SPT16", "STA5A", "T2FA", "TF2B", "TF3C2",
                 "UBF1", "UTX", "ZC11A")

   start <- 15

   for(targetGene in targetGenes[start:103]){
       width <- determineRegionForModel(targetGene)$width
       printf("%12s: %d", targetGene, width)
       if(targetGene %in% failures){
         checkEquals(width, 0)
       } else {
         checkTrue(width > 1000)
         }
       } # for targetGene

} # test_determineRegionForModel
#------------------------------------------------------------------------------------------------------------------------


