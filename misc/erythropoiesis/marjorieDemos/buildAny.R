library(httr)
library(jsonlite)
library(raster) # for cv
library(RUnit)
#-------------------------------------------------------------------------------------------------------------------
file <- "~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.corces.hg38.day0thru4.RData"
file.exists(file)
tbl.corces <- get(load(file))
source("tms.R")
tbl.diffbind <- get(load("~/github/TrenaProjectErythropoiesis/misc/diffBind/corces/tbl.diffBind.rory.hg38.day0-4.th10.RData"))
head(lapply(tbl.diffbind, class))
#-------------------------------------------------------------------------------------------------------------------
buildBindingTable <- function(goi, fimo.score.threshold, abs.cor.threshold, shoulder, tbl.oc=data.frame())
{
   printf("--- buildBindingTable: %s", goi)

   if(nrow(tbl.oc) == 0) {
      shoulder.startAndEnd <- shoulder
      uri <- sprintf("http://localhost:8000/geneLoc")
      body.jsonString <- sprintf('%s', toJSON(list(gene=goi, genome="hg38", shoulder=shoulder.startAndEnd)))
      r <- POST(uri, body=body.jsonString)
      loc <- fromJSON(httr::content(r)[[1]])
      if(is.na(loc$chrom)){
         printf("--- no genomic location for %s, returning empty table", goi)
        return(data.frame())
        }
      loc$width <- loc$end-loc$start
      printf("loc$width: %d", loc$width)
      if(nrow(tbl.oc) == 0){
        printf("--- no open chromatin  %s, returning empty table", goi)
        return(data.frame())
        } #
     } # if no (i.e., empty) tbl.oc provided

   browser()
   tbl <- build.model(goi, fimoThresholdAsNegativeExponent=fimo.score.threshold, tbl.openChromatin=tbl.oc)
   printf("--- fimo hits found: %d", nrow(tbl))
   tbl.sub <- subset(tbl, abs(cor) > abs.cor.threshold)
   printf("--- fimo hits with correlated gene expression, +/-: %d", nrow(tbl.sub))

    #------------------------------------------------------------
    # blend in the timepoint counts
    #------------------------------------------------------------


   tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.corces), GRanges(tbl.oc)))
   colnames(tbl.ov) <- c("corces", "oc")
   tbl.corces.sub <- tbl.corces[tbl.ov$corces, 4:9]
   #browser()
   tbl.ocWithCounts <- cbind(tbl.oc, tbl.corces.sub)
   tbl.ocWithCounts

   tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.sub), GRanges(tbl.ocWithCounts)))
   colnames(tbl.ov) <- c("tfs", "oc")
   tbl.oc.forMerge <- tbl.ocWithCounts[tbl.ov$oc, c("Fold","p.value","day0.1","day0.2","day0.3","day0.4","day2.1","day2.2")]
   #cbind(tbl.oc.forMerge, tbl.corces.sub)
   rownames(tbl.oc.forMerge) <- NULL
   colnames(tbl.oc.forMerge)[1:2] <- c("diffbind.score", "diffbind.pval")
   tbl.out <- cbind(tbl.sub, tbl.oc.forMerge)
   tbl.out

} # buildBindingTable
#------------------------------------------------------------------------------------------------------------------------
gois <- c("AFP",
          "AKAP4",
          "ATP6V1B1-AS1",
          "C10orf107",
          "C1orf194",
          "C8orf48",
          "CALR3",
          "CAMK2N2",
          "CCL1",
          "CLCA2",
          "CYP4F11",
          "DHRS2",
          "DRD2",
          "ENTHD1",
          "EPB41L4A-AS2",
          "HMMR-AS1",
          "KRT16",
          "LINC01602",
          "LUM",
          "MIP",
          "MIR5581",
          "MMP12",
          "NLRP11",
          "SLC25A48",
          "SVEP1",
          "TEX19",
          "TMEM26",
          "ZNF560")

#------------------------------------------------------------------------------------------------------------------------
runBuilds <- function()
{
    fimo <- 3
    abscor <- 0.25
    shoulder <- 10000
    tbls.out <- lapply(gois, function(goi) buildBindingTable(goi, fimo, abscor, shoulder))
    save(tbls.out, file=sprintf("tbls.out.%d.%f.%d.RData", fimo, abscor, shoulder))
    tbl <- do.call(rbind, tbls.out)
    mtx <- as.matrix(table(tbl$tf, tbl$targetGene))
    dim(mtx)
    list.cv <- sort(apply(mtx, 1, cv))

    subset(tbl, tf=="NR5A2" & abs(cor) > 0.3 & abs(diffbind.score) > 2)[, -c(8, 17)]
    subset(tbl, tf=="E2F3" & abs(cor) > 0.5 & abs(diffbind.score) > 2)[, -c(8, 17)]
    TF <- "ZBTB6"
    TF <- "TBX15"
    TF <- "ELF2"
    TF <- "CTCF"
    subset(tbl, tf==TF & abs(cor) > 0.2 & abs(diffbind.score) > 2)[, -c(8, 17)]
    length(unique(tbl$tf)) # [1] 460

} # runBuilds
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# Alpha fetoprotein (AFP) participates in the build up of hematopoietic cells in the early embryonic stage
# bmc diagnostic pathology, 2019:
# https://diagnosticpathology.biomedcentral.com/articles/10.1186/s13000-019-0858-5
# chr4   73331138   73556174   +    225k  (gene size + 100k both ways
test.afp <- function()
{
   loc <- list(chrom="chr4", start=73331138, end=73556174)
   fimo <- 3
   abscor <- 0.25
   tbl.oc <- subset(tbl.diffbind, chrom==loc$chrom & start >= loc$start & end <= loc$end)
   dim(tbl.oc) # 9 12
   printf("--- %d atac regions open in %dk bases for %s", nrow(tbl.oc), as.integer(loc$width/1000), goi)
   tbl.out <- buildBindingTable(goi, fimo, abscor, shoulder=0, tbl.oc)

   tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.out), GRanges(tbl.oc)))
   checkEquals(nrow(tbl.ov), nrow(tbl.out))
   table(tbl.ov$subjectHits)  #    1   2   3   4   5   6   7   8   9
                              #  546 569 441 485 454 513 407 544 493
   tbl.afp <- tbl.out[, coi.01]
   rownames(tbl.afp) <- NULL
   colnames(tbl.afp)[grep("targetGene", colnames(tbl.afp))] <- "gene"
   save(tbl.afp, file="tbl.afp.RData")

} # test.afp
#------------------------------------------------------------------------------------------------------------------------
coi.01 <- c("chrom", "start","end","tf","targetGene", "gh","tss","cor",
            "chip","motifScore","diffbind.score","diffbind.pval","day0.1",
            "day0.2","day0.3","day0.4","day2.1","day2.2", "strand","motif_id")


