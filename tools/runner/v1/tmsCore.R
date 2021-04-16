library(R6)
library(TrenaProjectErythropoiesis)
library(TrenaMultiScore)
library(RPostgreSQL)
library(trena)
#----------------------------------------------------------------------------------------------------
TMS = R6Class("TMS",

    #--------------------------------------------------------------------------------
    private = list(targetGene=NULL,
                   tp=NULL,
                   tms=NULL,
                   tbl.gh=NULL,
                   ghRegion=NULL,
                   tbl.fimo=NULL,
                   tbl.rbp=NULL,
		   tbl.trena=NULL,
		   linear.model=NULL,
                   tbls.lm=NULL,
                   lm.adjustedRsquareds=NULL,
                   mtx.rna=NULL,
                   igv=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(

        initialize = function(targetGene, trenaProject){
           suppressWarnings(db.access.test <-
                                try(system("ping -c 1 khaleesi", intern=TRUE, ignore.stderr=TRUE)))
           if(length(db.access.test) == 0)
              stop("khaleesi database server unavailable")
           printf("initializing TMS('%s')", targetGene)
           private$targetGene <- targetGene
           private$tp <- trenaProject
           private$tms <- TrenaMultiScore(private$tp, targetGene)
           private$tbls.lm <- list()
           private$lm.adjustedRsquareds <- list()
           },
        addGeneHancer = function(){
           private$tbl.gh <- getEnhancers(private$tp, tissues="all", maxSize=20000)
           private$ghRegion <- getGeneHancerRegion(private$tms)
           },
        getGeneHancer = function(){
           private$tbl.gh
           },
        getGeneHancerRegion = function(){
           private$ghRegion
           },
        viewOpenChromatin = function(tbl.atac.merge=data.frame()){
            if(is.null(private$igv)){
               private$igv <- start.igv(private$targetGene, "hg38")
               }
           tbl.gh.tmp <- private$tbl.gh[, c("chrom", "start", "end", "combinedscore")]
           track <- DataFrameQuantitativeTrack("gh", tbl.gh.tmp, autoscale=TRUE, color="brown")
           displayTrack(private$igv, track)
           if(nrow(tbl.atac.merged) > 0){
              loc.chrom <- private$ghRegion$chrom[1]
              loc.start <- private$ghRegion$start[1]
              loc.end <- private$ghRegion$end[1]
              tbl.atac.gene <- subset(tbl.atac.merged, chrom==loc.chrom & start >= loc.start & end <= loc.end)
              track <- DataFrameAnnotationTrack("atac", tbl.atac.gene, color="red")
              displayTrack(private$igv, track)
              tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.gh.tmp), GRanges(tbl.atac.gene)))
              tbl.atac.gh <- tbl.atac.gene[unique(tbl.ov[,2]),]
              track <- DataFrameAnnotationTrack("atac+gh", tbl.atac.gh, color="darkblue")
              displayTrack(private$igv, track)
              }
           },
        addOpenChromatin = function(chrom, start, end,
                                    intersect.with.genehancer,
                                    promoter.only){
           if(promoter.only){
              tbl.gh.promoter <- subset(private$tbl.gh, combinedscore > 500)
              region.count <- nrow(tbl.gh.promoter)
              stopifnot(region.count > 0)
              widths <- with(tbl.gh.promoter, 1 + end - start)
              biggest <- which(widths == max(widths))
              message(sprintf("tbl.gh.promoter: %d rows, biggest: %d", region.count, max(widths)))
              chrom <- tbl.gh.promoter$chrom[biggest]
              start <- tbl.gh.promoter$start[biggest]
              end   <- tbl.gh.promoter$end[biggest]
              }
            invisible(findOpenChromatin(private$tms, chrom, start, end,
                                       intersect.with.genehancer,
                                       use.merged.atac=TRUE))
           },
        getOpenChromatin = function(){
           getOpenChromatin(private$tms)
           },
        addAndScoreFimoTFBS = function(fimo.threshold=1e-4, motifs){
           findFimoTFBS(private$tms, motifs=motifs, fimo.threshold=fimo.threshold)
           tbl.tmp <- getMultiScoreTable(private$tms)
           if(nrow(tbl.tmp) == 0){
               s <- sprintf("no motif hits at threshold %f", fimo.threshold)
               stop(s)
               }
           addChIP(private$tms)
           scoreMotifHitsForConservation(private$tms)
           scoreMotifHitsForGeneHancer(private$tms)
           addGenicAnnotations(private$tms)
           addDistanceToTSS(private$tms)
           },
        addChIP = function(){
           addChIP(private$tms)
           },
        getTfTable = function(){
           getMultiScoreTable(private$tms)
           },

        addRBP = function(){
           db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="genereg2021", host="khaleesi")
             #  AND celltype='K562'
           query <- sprintf("select * from rbp where target='%s'", private$targetGene)
           tbl.rbp <- dbGetQuery(db, query)
           colnames(tbl.rbp)[grep("endpos", colnames(tbl.rbp))] <- "end"
           private$tbl.rbp <- tbl.rbp
           dbDisconnect(db)
           },

        getRbpTable = function(){
            private$tbl.rbp
            },

        add.tf.mrna.correlations = function(mtx.mrna, featureName){
           private$mtx.rna <- mtx.mrna
           addGeneExpressionCorrelations(private$tms, mtx.mrna, featureName)
           },

        add.rbp.mrna.correlations = function(mtx.mrna, featureName){
           f <- function(rbp){
              if(rbp %in% rownames(mtx.mrna))
                 return(cor(mtx.mrna[private$targetGene,], mtx.mrna[rbp,], method="spearman"))
              else return(NA)
              }
           suppressWarnings(cor <- unlist(lapply(private$tbl.rbp$gene, f)))
           private$tbl.rbp[, featureName] <- cor
           },

	 build.trena.model = function(tfs, rbps, mtx, order.by.column="spearmanCoeff", decreasing.order=TRUE){
           candidates <- intersect(c(tfs, rbps), rownames(mtx))
           suppressWarnings({
              solver <- EnsembleSolver(mtx,
                            targetGene=private$targetGene,
                            candidateRegulators=candidates,
                            geneCutoff=1.0,
                            solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"))
              tbl.out <- run(solver)
              })
           new.order <- order(tbl.out[, order.by.column], decreasing=decreasing.order)
	   tbl.out <- tbl.out[new.order,]
	   private$tbl.trena <- tbl.out
           tbl.out
	   },

        get.trena.model = function(){
           private$tbl.trena
           },

         build.linear.model = function(mtx, sort.by.column="rfScore", candidate.regulator.max=20){
           trena.order <- order(private$tbl.trena[,sort.by.column], decreasing=TRUE)
           private$tbl.trena <- private$tbl.trena[trena.order,]
           top.tfs <- private$tbl.trena$gene
	   if(length(top.tfs) > candidate.regulator.max)
	     top.tfs <- head(top.tfs, n=candidate.regulator.max)
           formula <- as.formula(sprintf("%s ~ 0 + .", private$targetGene))
	   tbl.data <- as.data.frame(t(mtx))[, c(private$targetGene, top.tfs)]
           model <- lm(formula, data=tbl.data)
	   private$linear.model <- model
           tbl.lm <- as.data.frame(coefficients(summary(model)))
           colnames(tbl.lm) <- c("estimate", "stderr", "t.value", "p.value")
           new.order <- order(tbl.lm$p.value, decreasing=FALSE)
           tbl.lm <- tbl.lm[new.order,]
           class <- rep("tf", nrow(tbl.lm))
           rbp.indices <- which(rownames(tbl.lm) %in% private$tbl.rbp$gene)
           class[rbp.indices] <- "rbp"
           tbl.lm$class <- class
           tbl.tms <- self$getTfTable()
           tfCount <- unlist(lapply(rownames(tbl.lm), function(tf.name) nrow(subset(tbl.tms, tf==tf.name))))
           tbl.lm$chip  <- unlist(lapply(rownames(tbl.lm), function(tf.name) nrow(subset(tbl.tms, tf==tf.name & chip))))
           rbpCount <- rep(0, nrow(tbl.lm))
           rbp.indices <- which(tbl.lm$class == "rbp")
           if(length(rbp.indices) > 0){
              rbp.names <-  rownames(tbl.lm)[rbp.indices]
              counts <- unlist(lapply(rbp.names, function(name) nrow(subset(private$tbl.rbp, gene==name))))
              rbpCount[rbp.indices] <- counts
              }
           tbl.lm$count <- tfCount + rbpCount
           correlated.expression <- unlist(lapply(rownames(tbl.lm),
                function(gene) if(gene %in% rownames(mtx)) cor(mtx[gene, ], mtx[private$targetGene,]) else 0))
           correlated.expression <- round(correlated.expression, digits=2)
           tbl.lm$cor <- correlated.expression
           private$tbls.lm[[sort.by.column]] <- tbl.lm
           private$lm.adjustedRsquareds[[sort.by.column]] <- summary(model)$adj.r.squared
	   invisible(tbl.lm)
           },

        get.lm.tables = function(){
           private$tbls.lm
           },

        get.lm.Rsquareds = function(){
           private$lm.adjustedRsquareds
           },

        get.mtx.rna = function(){
           invisible(private$mtx.rna)
           }

        #------------------------------------------------------------
       ) # public
    ) # class

#--------------------------------------------------------------------------------
test.tmsCore <- function(){
    require("TrenaProjectErythropoiesis")
    trenaProject <- TrenaProjectErythropoiesis()
    targetGene <- "BACH1"
    x <- TMS$new(targetGene, trenaProject)
    x$addGeneHancer()
    tbl.gh <- x$getGeneHancer()
    max.score <- max(tbl.gh$combinedscore)
    tbl.gh.promoter <- subset(tbl.gh, combinedscore == max.score)[1,]
    chrom <- tbl.gh.promoter$chrom
    start <- tbl.gh.promoter$start
    end   <- tbl.gh.promoter$end
    x$addOpenChromatin(chrom, start, end,
                       intersect.with.genehancer=FALSE,
                       promoter.only=FALSE)
    tbl.oc <- x$getOpenChromatin()
    motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "hocomocov11-core"))
    length(motifs)
    x$addAndScoreFimoTFBS(fimo.threshold=1e-4, motifs=motifs)
    tbl.tms <- x$getTfTable()
    dim(tbl.tms)  # 165 16
    length(unique(tbl.tms$tf))  # [1] 44
    x$addRBP()
    tbl.rbp <- x$getRbpTable()
    dim(tbl.rbp)  # 7652  12
    mtx.rna <- getExpressionMatrix(trenaProject, "brandLabDifferentiationTimeCourse-27171x28")

    x$add.tf.mrna.correlations(mtx.rna, featureName="cor.all")
    tbl.tms <- x$getTfTable()
    dim(tbl.tms)

    fivenum(tbl.tms$cor.all)
    x$add.rbp.mrna.correlations(mtx.rna, featureName="cor.all")
    tbl.rbp <- x$getRbpTable()
    dim(tbl.rbp)
    fivenum(tbl.rbp$cor.all)
    dim(tbl.rbp)  # 7652   13

    tfs <- unique(subset(tbl.tms, abs(cor.all) > 0.5)$tf)
    printf("candidate tfs: %d", length(tfs))
    rbps <- unique(subset(tbl.rbp, celltype=="K562" & abs(cor.all) > 0.5)$gene)
    printf("candidate rbps in K562 cells: %d", length(rbps))
    #tbl.trena.tf <- x$build.trena.model(tfs, list(), mtx.rna)
    #dim(tbl.trena.tf)  # 8 7
    tbl.trena.both <- x$build.trena.model(tfs, rbps, mtx.rna)
    linear.models <- list()
    linear.models[["spearman"]] <- x$build.linear.model(mtx.rna, sort.by.column="spearmanCoeff", candidate.regulator.max=15)
    linear.models[["pearson"]] <- x$build.linear.model(mtx.rna, sort.by.column="spearmanCoeff", candidate.regulator.max=15)
    linear.models[["betaLasso"]] <- x$build.linear.model(mtx.rna, sort.by.column="betaLasso", candidate.regulator.max=15)
    linear.models[["betaRidge"]] <- x$build.linear.model(mtx.rna, sort.by.column="betaRidge", candidate.regulator.max=15)
    linear.models[["rfScore"]] <- x$build.linear.model(mtx.rna, sort.by.column="rfScore", candidate.regulator.max=15)
    linear.models[["xgboose"]] <- x$build.linear.model(mtx.rna, sort.by.column="xgboost", candidate.regulator.max=15)
    printf("--- modeling %s, successful", targetGene)

    } # test.tmsCore
