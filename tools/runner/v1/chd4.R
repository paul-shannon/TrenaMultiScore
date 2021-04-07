source("tmsRunner.R")
targetGene <- "CHD4"
chromosome <- "chr12"
start.generous <-  6549433
end.generous   <-  6628028


viz <- FALSE

if(viz){   # determine start.focused, end.focused
    igv <- start.igv(targetGene, "hg38")
    showGenomicRegion(igv, sprintf("%s:%d-%d", chromosome, start.generous, end.generous))
    ghdb <- GeneHancerDB()
    tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all") # "Common myeloid progenitor CD34+")
    track <- DataFrameQuantitativeTrack("gh", tbl.gh[, c("chrom", "start", "end", "combinedscore")], autoscale=TRUE, color="brown")
    displayTrack(igv, track)
    tbl.atac.merged <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.atacMerged.RData"))
    loc <- getGenomicRegion(igv)
    tbl.atac.sub <- subset(tbl.atac.merged, chrom==chromosome & start >= loc$start & end <= loc$end)
    dim(tbl.atac.sub)
    track <- DataFrameAnnotationTrack("atac", tbl.atac.sub, color="red")
    displayTrack(igv, track)
    }

# first model, just atac in gh promoter: 12.5k
start.focused <- 6600262
end.focused   <- 6617901

# second model, many extra atac regions added: 71kb
start.focused <-   6520372
end.focused   <- 6627794

printf("width: %5.2f", (1 + end.focused - start.focused)/1000)


mtx.rna <- getExpressionMatrix(trenaProject, "brandLabDifferentiationTimeCourse-27171x28")
#mtx <- getExpressionMatrix(trenaProject, "brandLabDifferentiationTimeCourse-27171x28-namesCorrected")


x <- runTMS(targetGene, chromosome, start.focused, end.focused, mtx.rna)
lm.tables <- x$get.lm.tables()
lm.rsquareds <- x$get.lm.Rsquareds()

filename <- sprintf("%s-model-%s:%d-%d.RData", targetGene, chromosome,
                    start.focused, end.focused)
save(x, file=filename)

lm.tables <- x$get.lm.tables()
lapply(names(lm.tables), function(name) {
    printf("-------- %s", name)
    browser
    print(subset(lm.tables[[name]])) # , p.value < 0.25))
    })
print(x$get.lm.Rsquareds())
