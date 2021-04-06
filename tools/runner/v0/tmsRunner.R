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
                   tbl.lm=NULL,
                   lm.adjustedRsquared=NULL,
                   igv=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(

        initialize = function(targetGene, trenaProject){
           suppressWarnings(db.access.test <-
                                try(system("/sbin/ping -c 1 khaleesi", intern=TRUE, ignore.stderr=TRUE)))
           if(length(db.access.test) == 0)
              stop("khaleesi database server unavailable")
           printf("initializing TMS('%s')", targetGene)
           private$targetGene <- targetGene
           private$tp <- trenaProject
           private$tms <- TrenaMultiScore(private$tp, targetGene)
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
              browser()
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
           addGeneExpressionCorrelations(private$tms, mtx.mrna, featureName)
           },

        add.rbp.mrna.correlations = function(mtx.mrna, featureName){
           f <- function(rbp){
              if(rbp %in% rownames(mtx))
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
	   private$tbl.trena <- tbl.out
           new.order <- order(tbl.out[, order.by.column], decreasing=decreasing.order)
	   tbl.out[new.order,]
	   },

        get.trena.model = function(){
           private$tbl.trena
           },

         build.linear.model = function(sort.by.column="rfScore", candidate.regulator.max=20){
           trena.order <- order(private$tbl.trena[,sort.by.column], decreasing=TRUE)
           private$tbl.trena <- private$tbl.trena[trena.order,]
           top.tfs <- private$tbl.trena$gene
	   if(length(top.tfs) > candidate.regulator.max)
	     top.tfs <- head(top.tfs, n=candidate.regulator.max)
           formula <- as.formula(sprintf("%s ~ 1 + .", private$targetGene))
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
           private$tbl.lm <- tbl.lm
           private$lm.adjustedRsquared <- summary(model)$adj.r.squared
	   invisible(tbl.lm)
           },

        get.lm.table = function(){
           private$tbl.lm
           },

        get.lm.Rsquared = function(){
           private$lm.adjustedRsquared
           }

        #------------------------------------------------------------
       ) # public
    ) # class

#--------------------------------------------------------------------------------
