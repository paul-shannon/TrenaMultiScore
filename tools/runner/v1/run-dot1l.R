source("tmsRunner.R")
targetGene <- "DOT1L"
chromosome <- "chr19"
start.generous <- 2000000
end.generous   <- 2500000


viz <- FALSE

if(viz){   # determine start.focused, end.focused
    library(ghdb)
    igv <- start.igv("DOT1L", "hg38")
    showGenomicRegion(igv, sprintf("%s:%d-%d", chromosome, start.generous, end.generous))
    ghdb <- GeneHancerDB()
    tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all") # "Common myeloid progenitor CD34+")
    track <- DataFrameQuantitativeTrack("gh", tbl.gh[, c("chrom", "start", "end", "combinedscore")], autoscale=TRUE, color="brown")
    displayTrack(igv, track)
    tbl.atac.merged <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.atacMerged.RData"))
    tbl.atac.sub <- subset(tbl.atac.merged, chrom==chromosome & start >= start.generous & end <= end.generous)
    dim(tbl.atac.sub)
    track <- DataFrameAnnotationTrack("atac", tbl.atac.sub, color="red")
    displayTrack(igv, track)
    }

start.focused <- 2160508
end.focused  <- 2174767

mtx.rna <- getExpressionMatrix(trenaProject, "brandLabDifferentiationTimeCourse-27171x28")
#mtx <- getExpressionMatrix(trenaProject, "brandLabDifferentiationTimeCourse-27171x28-namesCorrected")


dot1l.model <- runTMS(targetGene, chromosome, start.focused, end.focused, mtx.rna)


