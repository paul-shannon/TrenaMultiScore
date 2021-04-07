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

# first model, just atac in gh promoter: 14k
start.focused <- 2160508
end.focused  <- 2174767
sprintf("%s:%d-%d", chromosome, start.focused, end.focused)

# second model, 3 extra atac regions added: 20k
start.focused <- 2157837
end.focused <- 2177660

printf("width: %5.2f", (1 + end.focused - start.focused)/1000)

mtx.rna <- getExpressionMatrix(trenaProject, "brandLabDifferentiationTimeCourse-27171x28")
#mtx <- getExpressionMatrix(trenaProject, "brandLabDifferentiationTimeCourse-27171x28-namesCorrected")


x <- runTMS(targetGene, chromosome, start.focused, end.focused, mtx.rna)
lm.tables <- x$get.lm.tables()
lm.rsquareds <- x$get.lm.Rsquareds()

filename <- sprintf("%s-model-%s:%d-%d.RData", targetGene, chromosome,
                    start.focused, end.focused)
save(x, file=filename)

