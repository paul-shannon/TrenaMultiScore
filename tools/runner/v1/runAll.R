# https://statisticsbyjim.com/regression/overfitting-regression-models/
source("tmsRunner.R")
source("determineRegionForModel.R")

mtx.tmp <- get(load("~/github/TrenaProjectErythropoiesis/viz/srm.vs.mrna/shinyapps.io/srm.rna.averaged.clean.RData"))
targetGenes <- rownames(mtx.tmp)
mtx.rna <- getExpressionMatrix(trenaProject, "brandLabDifferentiationTimeCourse-27171x28")

for(targetGene in targetGenes[49:103]){
   tbl.region <- determineRegionForModel(targetGene)
   if(tbl.region$width == 0){
       printf("no promoter region for %s", targetGene)
       next;
       }
   chromosome <- tbl.region$chrom
   start <- tbl.region$start
   end <- tbl.region$end
   x <- tryCatch(
       runTMS(targetGene, tbl.region$chrom, tbl.region$start, tbl.region$end, mtx.rna),
       error=function(e){
           printf("failure with %s", targetGene)
           return(NA)
           }
       )
   if(!is.na(x)){
       lm.tables <- x$get.lm.tables()
       lm.rsquareds <- x$get.lm.Rsquareds()
       filename <- sprintf("%s-model-%s:%d-%d.RData", targetGene, chromosome, start, end)
       save(x, file=filename)
       lm.tables <- x$get.lm.tables()
       lapply(names(lm.tables), function(name) {
           printf("-------- %s", name)
           print(subset(lm.tables[[name]])) # , p.value < 0.25))
       })
       print(x$get.lm.Rsquareds())
       } # if proceed
   } # for targetGenes
