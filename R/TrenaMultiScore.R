# 1) find tss, chrom, strand, min and max of genehancer regions (aka, "geneRegion")
# 2) load open chromatin in geneRegion: atac-seq or database footprint regions
# 3) find all TFBS in atac-seq regions
# 4) score all TFBS for conservation
# 5) determine distance from TSS for all TFBS
#------------------------------------------------------------------------------------------------------------------------
#' @import methods
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#' @import RPostgreSQL
#' @import ghdb
#'
#' @import GenomicRanges
#' @import GenomicScores
#' @import phastCons7way.UCSC.hg38
#' @import phastCons100way.UCSC.hg38
#' @import MotifDb
#' @import motifmatchr
#' @importFrom universalmotif convert_motifs
#' @importFrom TFBSTools PFMatrixList
#' @import annotatr
#'
#' @title TrenaMultiScore-class
#'
#' @name TrenaMultiScore-class
#' @rdname TrenaMultiScore-class
#' @aliases TrenaMultiScore
#' @exportClass TrenaMultiScore
#'

.TrenaMultiScore <- setClass("TrenaMultiScore",
                             representation=representation(
                                trenaProject="TrenaProject",
                                targetGene="character",
                                motifDb="MotifList",
                                state="environment"
                                )
                             )
#------------------------------------------------------------------------------------------------------------------------
setGeneric('getProject', signature='obj', function(obj) standardGeneric('getProject'))
setGeneric('getGeneHancerRegion', signature='obj', function(obj) standardGeneric('getGeneHancerRegion'))
setGeneric('findOpenChromatin', signature='obj', function(obj, chrom=NA, start=NA, end=NA,
                                                          intersect.with.geneHancer=FALSE,
                                                          use.merged.atac=FALSE)

              standardGeneric('findOpenChromatin'))
#setGeneric('explicitlySetOpenChromatin', signature='obj', function(obj, tbl) standardGeneric('explicitlySetOpenChromatin'))
setGeneric('getOpenChromatin', signature='obj', function(obj) standardGeneric('getOpenChromatin'))
setGeneric('findFimoTFBS', signature='obj', function(obj, motifs=NA, fimo.threshold=NA, chrom=NA, start=NA, end=NA, genome="hg38")
              standardGeneric('findFimoTFBS'))
setGeneric('findMoodsTFBS', signature='obj', function(obj, motifs=NA, moods.threshold=NA, chrom=NA, start=NA, end=NA)
              standardGeneric('findMoodsTFBS'))
setGeneric('scoreMotifHitsForConservation', signature='obj', function(obj) standardGeneric('scoreMotifHitsForConservation'))
setGeneric('getMultiScoreTable', signature='obj', function(obj) standardGeneric('getMultiScoreTable'))
setGeneric('getTargetGeneInfo', signature='obj', function(obj) standardGeneric('getTargetGeneInfo'))
setGeneric('addDistanceToTSS', signature='obj', function(obj) standardGeneric('addDistanceToTSS'))
setGeneric('scoreMotifHitsForGeneHancer', signature='obj', function(obj) standardGeneric('scoreMotifHitsForGeneHancer'))
setGeneric('addGeneExpressionCorrelations', signature='obj', function(obj, mtx, featureName="cor", colnames=list()) standardGeneric('addGeneExpressionCorrelations'))
setGeneric('addGenicAnnotations', signature='obj', function(obj) standardGeneric('addGenicAnnotations'))
setGeneric('addChIP', signature='obj', function(obj) standardGeneric('addChIP'))
setGeneric('addRnaBindingProteins', signature='obj', function(obj) standardGeneric('addRnaBindingProteins'))
setGeneric('getRnaBindingProteins', signature='obj', function(obj) standardGeneric('getRnaBindingProteins'))
setGeneric('scoreMotifHitsForOpenChromatin', signature='obj', function(obj) standardGeneric('scoreMotifHitsForOpenChromatin'))
#------------------------------------------------------------------------------------------------------------------------
#' Define an object of class TrenaMultiScore
#'
#' @description
#' Expression, variant and covariate data for the genes of interest (perhaps unbounded) for pre-term birth studies
#'
#' @rdname TrenaMultiScore-class
#' @param trenaProject a concrete subclass of TrenaProject
#' @param targetGene a character string
#' @param tbl.fimo a data.frame, empty or with fimo matches calculated separate
#' @param quiet logical
#'
#' @export
#'
#' @return An object of the TrenaMultiScore class
#'

TrenaMultiScore <- function(trenaProject, targetGene, tbl.fimo=data.frame(), tbl.oc=data.frame(), quiet=TRUE)
{
   setTargetGene(trenaProject, targetGene)
   state <- new.env(parent=emptyenv())
   state$quiet <- quiet
   state$openChromatin <- tbl.oc
   state$genehancer <- data.frame()
   state$fimo <- tbl.fimo

   .TrenaMultiScore(trenaProject=trenaProject, targetGene=targetGene, motifDb=MotifDb, state=state)

} # TrenaMultiScore, the constructor
#------------------------------------------------------------------------------------------------------------------------
#' return the TMS project
#'
#' @description
#' get a reference to the TrenaProject concrete subclass included in this object
#'
#' @rdname getProject
#' @return a TrenaProject instance
#'
#' @export
#'
setMethod('getProject', 'TrenaMultiScore',

     function(obj){
        obj@trenaProject
        }) # getProject

#------------------------------------------------------------------------------------------------------------------------
#' return the full extent of GeneHancer regions for the targetGene
#'
#' @description
#' GeneHancer often reports distal enhancers within ~1MB of the gene; return that region
#'
#' @rdname getGeneHancerRegion
#' @return a data.frame
#'
#' @export
#'
setMethod('getGeneHancerRegion', 'TrenaMultiScore',

     function(obj){
        ghdb <- GeneHancerDB()
        tbl.gh <- retrieveEnhancersFromDatabase(ghdb, obj@targetGene, tissues="all")
        obj@state$genehancer <- tbl.gh
        start <- min(tbl.gh$start) - 1000
        end <- max(tbl.gh$end) + 1000
        width <- 1 + end - start
        data.frame(chrom=tbl.gh$chrom[1],
                   start=start, end=end, width=width,
                   stringsAsFactors=FALSE)

        }) # getGeneHancerRegion

#------------------------------------------------------------------------------------------------------------------------
#' find open chromatin in the current gene region, or (if supplied) somewhere else
#'
#' @description
#' Open chromatin from DNase footprints or scATAC-seq is retrieved, saved for later use
#'
#' @rdname findOpenChromatin
#' @param obj a TrenaMultiScore object
#' @param chrom typically NA
#' @param chrom start, integer, chromosomal location, default NA
#' @param chrom end, integer, chromosomal location, default NA
#'
#' @export
#'
setMethod('findOpenChromatin', 'TrenaMultiScore',

     function(obj, chrom=NA, start=NA, end=NA, intersect.with.geneHancer= FALSE, use.merged.atac=FALSE){
         if(is.na(chrom)){
            tbl.x <- getGeneHancerRegion(obj)
            chrom <- tbl.x$chrom
            start <- tbl.x$start
            end   <- tbl.x$end
            }
        if("TrenaProjectErythropoiesis" %in% is(getProject(obj))){
           obj@state$openChromatin <- .queryBrandLabATACseq(chrom, start, end, use.merged.atac)
           if(!obj@state$quiet)
               message(sprintf("regions of Brand ATACseq open chromatin: %d", nrow(obj@state$openChromatin)))
           }
        else if("TrenaProjectAD" %in% is(getProject(obj))){
           #obj@state$openChromatin <- .queryBocaATACseq(chrom, start, end)
           obj@state$openChromatin <- .queryHintFootprintRegionsFromDatabase("brain_hint_16", chrom, start, end)
           if(!obj@state$quiet)
              message(sprintf("regions of open chromatin: %d", nrow(obj@state$openChromatin)))
           }
        else stop(sprintf("no support for open chromatin retrieval in %s", getProject(obj)@projectName))
        if(intersect.with.geneHancer){
            tbl.oc <- obj@state$openChromatin
            tbl.gh <- getGeneRegulatoryRegions(getProject(obj))
            tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.oc), GRanges(tbl.gh)))
            if(nrow(tbl.ov) == 0){
                tbl.filtered <- data.frame()
            } else {
               tbl.filtered <- tbl.oc[tbl.ov[,1],]
               }
            message(sprintf("intersecting oc with gh: from %d to %d regions", nrow(tbl.oc), nrow(tbl.filtered)))
            obj@state$openChromatin <- tbl.filtered
            }
        nrow(obj@state$openChromatin)
        })


#------------------------------------------------------------------------------------------------------------------------
#' retrieve currently assigned open chromatin
#'
#' @description
#' Return a data.frame with the open chromatin previously retrieved
#'
#' @rdname getOpenChromatin
#' @param obj a TrenaMultiScore object
#' @return data.frame
#'
#' @export
#'
setMethod('getOpenChromatin', 'TrenaMultiScore',

     function(obj){
         obj@state$openChromatin
         })


#------------------------------------------------------------------------------------------------------------------------
#' get scored FIMO motif binding sites
#'
#' @description
#' call FimoBatch to get data.frame of all matches above threshold
#'
#' @rdname getFimoTFBS
#'
#' @param motifs a list of MotifDb items, JASPAR2019 hsapiens by default
#' @param fimo.threshold 1e-4 by default
#' @param chrom geneHancerRegion chrom by default
#' @param start geneHancerRegion start by default
#' @param end geneHancerRegion end by default
#' @param genome character string, "hg38" by default
#'
#' @return a data.frame
#'
#' @export
#'
setMethod('findFimoTFBS', 'TrenaMultiScore',

    function(obj, motifs=list(), fimo.threshold=NA, chrom=NA, start=NA, end=NA, genome="hg38"){

           # fimo table may have been build separately, do not repeat.
       if(nrow(obj@state$fimo) > 1) return

       tbl.fp <- getOpenChromatin(obj)
       if(!obj@state$quiet)
            message(sprintf("--- findFimoTFBS in %d regions of open chromatin", nrow(tbl.fp)))

       if(nrow(tbl.fp) == 0)
          stop("TrenaMultiScore::getFimoTFBS error: no open chromatin regions previously identified.")

       if(length(motifs) == 0)
          motifs <- query(obj@motifDb, c("sapiens"), c("hocomocov11a-core", "jaspar2018"))

       if(is.na(fimo.threshold))
          fimo.threshold <- 1e-4
       tbl.gh <- getGeneHancerRegion(obj)
       if(is.na(chrom))
          chrom <- tbl.gh$chrom
       if(is.na(start))
          start <- tbl.gh$start
       if(is.na(end))
          end <- tbl.gh$end

       meme.file <- "tmp.meme"
       MotifDb::export(motifs, con=meme.file, format="meme")
       source("~/github/fimoService/batchMode/fimoBatchTools.R")
       if(!obj@state$quiet)
          message(sprintf("--- calling fimoBatch on %d regions, %e threshold", nrow(tbl.fp), fimo.threshold))
       tbl.fimo <- fimoBatch(tbl.fp, matchThreshold=fimo.threshold, genomeName=genome, pwmFile=meme.file)
       obj@state$fimo <- tbl.fimo
       if(!obj@state$quiet)
         message(sprintf("--- tbl.fimo stored with %d rows, %d columns", nrow(tbl.fimo), ncol(tbl.fimo)))
       }) # findFimoTFBS

#------------------------------------------------------------------------------------------------------------------------
#' get scored moods motif binding sites
#'
#' @description
#' call moods (wrapped by motifmatchr package) to get data.frame of all matches above threshold
#'
#' @rdname getMoodsTFBS
#'
#' @param motifs a list of MotifDb items, JASPAR2019 hsapiens by default
#' @param moods.threshold 1e-4 by default
#' @param chrom geneHancerRegion chrom by default
#' @param start geneHancerRegion start by default
#' @param end geneHancerRegion end by default
#'
#' @return a data.frame
#'
#' @export
#'
setMethod('findMoodsTFBS', 'TrenaMultiScore',

    function(obj, motifs=NA, moods.threshold=NA, chrom=NA, start=NA, end=NA){

       tbl.oc <- getOpenChromatin(obj)
       if(nrow(tbl.oc) == 0)
          stop("TrenaMultiScore::getMoodsTFBS error: no open chromatin regions previously identified.")

       if(is.na(motifs))
           motifs <- query(obj@motifDb, c("sapiens"), c("hocomoco", "jaspar2018"))

       if(is.na(moods.threshold))
          moods.threshold <- 1e-4

       gr.oc <- GRanges(tbl.oc)
       motifs.pfmatrix <- lapply(motifs, function(motif) convert_motifs(motif, "TFBSTools-PFMatrix"))
       motifs.pfmList <- do.call(PFMatrixList, motifs.pfmatrix)
       gr.match <- matchMotifs(motifs.pfmList, gr.oc, genome="hg38", out="positions", p.cutoff=moods.threshold)
       tbl.moods <- as.data.frame(gr.match)
       if(ncol(tbl.moods) == 8){  # the expected number
          if(nrow(tbl.moods) > 0){
             tfs <- mcols(MotifDb[tbl.moods$group_name])$geneSymbol
             tbl.moods$tf <- tfs
             }
          colnames(tbl.moods)[grep("seqnames", colnames(tbl.moods))] <- "chrom"
          colnames(tbl.moods)[grep("group_name", colnames(tbl.moods))] <- "motifName"
          colnames(tbl.moods)[grep("score", colnames(tbl.moods))] <- "moodsScore"
          tbl.moods$chrom <- as.character(tbl.moods$chrom)
          tbl.moods$strand <- as.character(tbl.moods$strand)
          }
       obj@state$moods <- tbl.moods
       sprintf("tbl.moods stored with %d rows, %d columns", nrow(tbl.moods), ncol(tbl.moods))
       }) # findMoodsTFBS

#------------------------------------------------------------------------------------------------------------------------
.queryBrandLabATACseq <- function(chrom.loc, start.loc, end.loc, use.merged.atac)
{
  tbl.atac <- get(load("~/github/TrenaProjectErythropoiesis/misc/multiScore/brandAtacSeqCuration/tbl.atac.fp.RData"))
  if(use.merged.atac){
      tbl.atac <- get(load("~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.atacMerged.RData"))
      }
  subset(tbl.atac, chrom==chrom.loc & start >= start.loc & end <= end.loc)

} # .queryBrandLabATACseq
#------------------------------------------------------------------------------------------------------------------------
.queryHintFootprintRegionsFromDatabase <- function(database.name, chrom.loc, start.loc, end.loc)
{
  suppressWarnings(
    db.access.test <- try(system("/sbin/ping -c 1 khaleesi", intern=TRUE, ignore.stderr=TRUE)))
  if(length(db.access.test) == 0)
     stop("khaleesi database server unavailable")

   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="brain_hint_16", host="khaleesi")
   query <- sprintf("select * from regions where chrom='%s' and start >= %d and endpos <= %d",
                    chrom.loc, start.loc, end.loc)
   tbl <- dbGetQuery(db, query)
   dbDisconnect(db)

   if(nrow(tbl) > 0){
      tbl <- tbl[, c("chrom", "start", "endpos")]
      colnames(tbl) <- c("chrom", "start", "end")
      }

   tbl

} # .queryHintFootprintRegionsFromDatabase
#------------------------------------------------------------------------------------------------------------------------
.queryBocaATACseq <- function(chrom.loc, start.loc, end.loc)
{
  f <- "~/github/TrenaProjectAD/misc/multiScore/boca/tbl.boca.RData"
  tbl.boca <- get(load(f))

  tbl.tmp <- subset(tbl.boca, chrom==chrom.loc & start >= start.loc & end <= end.loc)
  tbl.reduced <- as.data.frame(reduce(GRanges(tbl.tmp)), row.names=NULL)
  colnames(tbl.reduced)[1] <- "chrom"
  tbl.reduced$chrom <- as.character(tbl.reduced$chrom)

  tbl.reduced

} # .queryBocaATACseq
#------------------------------------------------------------------------------------------------------------------------
#' add open chromatin intersection info (none, partial, complete for the currently held fimo table
#'
#' @description
#'   adds "none", "partial", "complete" annotation to the "oc" column for each motif hit in the
#'   already built (or supplied) fimo table, using the current obj@state$openChromatin table
#'
#'
#' @rdname scoreMotifHitsForConservation
#'
#' @param obj a TrenaMultiScore object
#' @return a data.frame
#'
#' @export
#'
setMethod('scoreMotifHitsForOpenChromatin', 'TrenaMultiScore',

    function(obj){

      tbl.fimo <- obj@state$fimo
      tbl.oc <- obj@state$openChromatin
      stopifnot(nrow(tbl.oc) > 0)
      gr.fimo <- GRanges(seqnames=tbl.fimo$chrom, IRanges(tbl.fimo$start, tbl.fimo$end))
      gr.oc   <- GRanges(seqnames=tbl.oc$chrom, IRanges(tbl.oc$start, tbl.oc$end))
      tbl.ov  <- as.data.frame(findOverlaps(gr.oc, gr.fimo, type="any"))
      openChromatin <- rep(FALSE, nrow(tbl.fimo))
      if(nrow(tbl.ov) > 0){
         hits <- unique(tbl.ov$subjectHits)
         openChromatin[hits] <- TRUE
         }
      tbl.fimo$oc <- openChromatin
      rownames(tbl.fimo) <- NULL
      obj@state$fimo <- tbl.fimo
      }) # scoreMotifHitsForOpenChromatin


#------------------------------------------------------------------------------------------------------------------------
#' add conservation scores to the currently held fimo table
#'
#' @description
#'   adds several conservation scores, each an average, to each motif hit in the already build fimo table
#'
#' @rdname scoreMotifHitsForConservation
#'
#' @param obj a TrenaMultiScore object
#' @return a data.frame
#'
#' @export
#'
setMethod('scoreMotifHitsForConservation', 'TrenaMultiScore',

    function(obj){

      tbl.fimo <- obj@state$fimo
      if(nrow(tbl.fimo) == 0)
         stop("TrenaMultiScore::scoreMotifHitsForConservation error: no fimo hits yet identified.")

        # phastCons100way.UCSC.hg38
        # phastCons30way.UCSC.hg38

       # browser()
       tbl.cons <- as.data.frame(gscores(phastCons7way.UCSC.hg38,
                                         GRanges(tbl.fimo[, c("chrom", "start", "end")])), stringsAsFactors=FALSE)
       tbl.fimo$phast7 <- round(tbl.cons$default, digits=2)

       #if(grepl("khaleesi", Sys.info()[["nodename"]])){
       #   print(load("~/github/TrenaMultiScore/inst/extdata/phastCons30way.UCSC.hg38.downloaded.RData"))
       #   tbl.cons <- as.data.frame(gscores(phastCons30way.UCSC.hg38.downloaded,
       #                                     GRanges(tbl.fimo[, c("chrom", "start", "end")])), stringsAsFactors=FALSE)
       #   }
       #else{
       #   tbl.cons <- as.data.frame(gscores(phastCons30way.UCSC.hg38,
       #                                     GRanges(tbl.fimo[, c("chrom", "start", "end")])), stringsAsFactors=FALSE)
       #   }

       #tbl.fimo$phast30 <- round(tbl.cons$default, digits=2)

       tbl.cons <- as.data.frame(gscores(phastCons100way.UCSC.hg38,
                                         GRanges(tbl.fimo[, c("chrom", "start", "end")])), stringsAsFactors=FALSE)
       tbl.fimo$phast100 <- round(tbl.cons$default, digits=2)

       rownames(tbl.fimo) <- NULL
       obj@state$fimo <- tbl.fimo
       }) # scoreMotifHitsForConservation

#------------------------------------------------------------------------------------------------------------------------
#' returns the current state of scored motif matches
#'
#' @description
#'   can be called at any time
#'
#' @rdname getMultiScoreTable
#'
#' @param obj a TrenaMultiScore object
#'
#' @return a data.frame
#'
#' @export
#'
setMethod('getMultiScoreTable', 'TrenaMultiScore',

    function(obj){

      tbl.fimo <- obj@state$fimo
      if("score" %in% colnames(tbl.fimo))
          colnames(tbl.fimo)[grep("score", colnames(tbl.fimo))] <- "fimo_score"
      if("p.value" %in% colnames(tbl.fimo))
          colnames(tbl.fimo)[grep("p.value", colnames(tbl.fimo))] <- "fimo_pvalue"

      if("motif_id" %in% colnames(tbl.fimo)){ # can be quite long, move it to the end
         current.index <- grep("motif_id", colnames(tbl.fimo))
         colnames.in.preferred.order <- c(colnames(tbl.fimo)[-current.index], "motif_id")
         tbl.fimo <- tbl.fimo[,colnames.in.preferred.order]
         } # motif_id column

      invisible(unique(tbl.fimo))

      }) # getMultiScoreTable

#------------------------------------------------------------------------------------------------------------------------
#' append 'TSS' column to the accumulating table
#'
#' @description
#'   append 'TSS' column to the accumulating table, recording distance from the main
#'   transcript's TSS to the start of the motif hit, negative values for upstream locations
#'
#' @rdname addDistanceToTSS
#'
#' @param obj a TrenaMultiScore object
#'
#' @return None
#'
#' @export
#'
setMethod('addDistanceToTSS', 'TrenaMultiScore',

    function(obj){

      coi <- c("tss", "strand", "chrom", "start", "end")
      targetGene.info <- as.list(getTranscriptsTable(getProject(obj))[coi])
      tbl <- getMultiScoreTable(obj)
      tss <- NA_integer_
      if(length(targetGene.info) > 0)
         tss <- (targetGene.info$tss - tbl$start)  * (targetGene.info$strand) * -1
      tbl$tss <- tss
      obj@state$fimo <- tbl
      }) # addDistanceToTSS

#------------------------------------------------------------------------------------------------------------------------
#' get basic stats on the primary transcript
#'
#' @description
#'   get basic stats on the primary transcript

#' @rdname getTargetGeneInfo
#'
#' @param obj a TrenaMultiScore object
#'
#' @return None
#'
#' @export
#'
setMethod('getTargetGeneInfo', 'TrenaMultiScore',

    function(obj){

      coi <- c("geneSymbol", "tss", "strand", "chrom", "start", "end", "ensg")
      as.list(getTranscriptsTable(getProject(obj))[coi])
      }) # getTargetGeneInfo

#------------------------------------------------------------------------------------------------------------------------
#' attach a genehancer score to each motif hit, 0 if none
#'
#' @description
#'   adds several conservation scores, each an average, to each motif hit in the already build fimo table
#'
#' @rdname scoreMotifHitsForGeneHancer
#'
#' @param obj a TrenaMultiScore object
#' @return a data.frame
#'
#' @export
#'
setMethod('scoreMotifHitsForGeneHancer', 'TrenaMultiScore',

    function(obj){

      if(nrow(obj@state$genehancer) == 0)
         getGeneHancerRegion(obj)
      tbl.gh <- obj@state$genehancer
      if(!grepl("chr", tbl.gh$chrom[1]))
         tbl.gh$chrom <- paste0("chr", tbl.gh$chrom)
      tbl.fimo <- obj@state$fimo
      if(nrow(tbl.fimo) == 0)
         stop("TrenaMultiScore::scoreMotifHitsForGeneHancer error: no fimo hits yet identified.")

      tbl.fimo$gh <- rep(0, nrow(tbl.fimo))
      if(nrow(tbl.gh) > 0){
         tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.fimo), GRanges(tbl.gh)))
         colnames(tbl.ov) <- c("fimo", "gh")
         tbl.fimo$gh[tbl.ov$fimo] <- round(tbl.gh$combinedscore[tbl.ov$gh], digits=2)
         }

      rownames(tbl.fimo) <- NULL
      obj@state$fimo <- tbl.fimo
      }) # scoreMotifHitsForGeneHancer

#------------------------------------------------------------------------------------------------------------------------
#' add a tf/targetGene correlations score for those TFs also in the expression matrix
#'
#' @description
#' add a tf/targetGene correlations score for those TFs also in the expression matrix
#'
#' @rdname addGeneExpressionCorrelations
#'
#' @param obj a TrenaMultiScore object
#' @param mtx a numerical matrix, genes are rownames, samples are colnames
#' @param featureName a character string, default "cor", data.frame column where results are stored
#' @param colnames list of character strings, limits scope of correlation calcuation. default empty list: use all columns
#' @return None
#'
#' @export
#'
setMethod('addGeneExpressionCorrelations', 'TrenaMultiScore',

    function(obj, mtx, featureName="cor", colnames=list()){

      if(!obj@state$quiet){
         printf("TMS::addGeneExpressionCorrelations, '%s'", featureName)
         }

      if(length(colnames) > 0){
         stopifnot(all(colnames %in% colnames(mtx)))
         mtx <- mtx[, colnames]
         }

      tbl.fimo <- obj@state$fimo
      if(nrow(tbl.fimo) == 0)
         stop("TrenaMultiScore::addGeneExpressionCorrelations error: no fimo hits yet identified.")

      f <- function(tf){
         if(tf %in% rownames(mtx))
           return(cor(mtx[obj@targetGene,], mtx[tf,], method="pearson", use="pairwise.complete.obs"))
         else return(NA)
         }

      if(!obj@state$quiet){
         printf("about to run spearman cor on %d tfs from tbl.fimo", length(tbl.fimo$tf))
         }

      suppressWarnings({
          cor <- unlist(lapply(tbl.fimo$tf, f))
          })

      cor <- round(cor, digits=2)
      tbl.fimo[, featureName] <- cor
      obj@state$fimo <- tbl.fimo
      }) # addGeneExpressionCorrelations

#------------------------------------------------------------------------------------------------------------------------
#' add intron, exon, utr, promoter, etc annotations
#'
#' @description
#' add intron, exon, utr, promoter, etc annotations
#'
#' @rdname addGenicAnnotations
#'
#' @param obj a TrenaMultiScore object
#' @param mtx a numerical matrix, genes are rownames, samples are colnames
#' @return None
#'
#' @export
#'
setMethod('addGenicAnnotations', 'TrenaMultiScore',

    function(obj){

      tbl.fimo <- obj@state$fimo
      if(nrow(tbl.fimo) == 0)
         stop("TrenaMultiScore::addGenicAnnotations error: no fimo hits yet identified.")

      existing.colnames <- colnames(tbl.fimo)
      gr <- GRanges(tbl.fimo)
        # anno <- build_annotations(genome="hg38", annotations=c("hg38_basicgenes", "hg38_genes_intergenic"))
      anno <- get(load(system.file(package="TrenaMultiScore", "extdata", "genomeAnnotations.RData")))
      gr.annoResults <- annotate_regions(regions=gr, annotations=anno, ignore.strand=TRUE, quiet=FALSE)
      tbl.anno <- as.data.frame(gr.annoResults, row.names=NULL)
      coi <- c(existing.colnames, "annot.symbol", "annot.type")
      coi[1] <- "seqnames"  # bioc-speak for chromosome
      na.symbols <- which(is.na(tbl.anno$annot.symbol))
      if(length(na.symbols) > 0)
          tbl.anno$annot.symbol[na.symbols] <- "NA"
      tbl.aggregated <- aggregate(cbind(annot.type, annot.symbol) ~ seqnames+start+end+strand+tf,
                                  data=tbl.anno, function(x) paste(unique(x), collapse=","))
      tbl.aggregated$annot.type <- gsub("hg38_genes_", "", tbl.aggregated$annot.type)
      colnames(tbl.aggregated)[1] <- "chrom"
      tbl.aggregated$chrom <- as.character(tbl.aggregated$chrom)

      tbl.fimo <- merge(tbl.fimo, tbl.aggregated)
      obj@state$fimo <- tbl.fimo
      }) # addGenicAnnotations

#------------------------------------------------------------------------------------------------------------------------
#' add 1 if a motif hit intersects a ChIP-seq hit, otherwise 0.
#'
#' @description
#' add 1 if a motif hit intersects a ChIP-seq hit, otherwise 0.  full-with ChIP peaks are used.
#'
#' @rdname addChIP
#'
#' @param obj a TrenaMultiScore object
#' @param mtx a numerical matrix, genes are rownames, samples are colnames
#' @return None
#'
#' @export
#'
setMethod('addChIP', 'TrenaMultiScore',

    function(obj){

      tbl.fimo <- obj@state$fimo
      if(nrow(tbl.fimo) == 0)
         stop("TrenaMultiScore::addChIP error: no fimo hits yet identified.")

      has.ChIP <- rep(FALSE, nrow(tbl.fimo))
      chrom <- tbl.fimo$chrom[1]
      start <- min(tbl.fimo$start) - 1000
      end <- max(tbl.fimo$end) + 1000
      tbl.chip <- getChipSeq(getProject(obj), chrom, start, end)
      tbl.chip <- subset(tbl.chip, tf %in% tbl.fimo$tf)
      tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.fimo), GRanges(tbl.chip), type="any"))
      tbl.combined <- cbind(tbl.fimo[tbl.ov$queryHits,], tbl.chip[tbl.ov$subjectHits,])
      colnames(tbl.combined)[max(grep("tf", colnames(tbl.combined)))] <- "tf.chip"
      tbl.fimoWithChip <- subset(tbl.combined, tf==tf.chip)
          # some chip hits are reported in multiple cell lines, which we don't care about yet
          # remove them
      tbl.fimoWithChip <- tbl.fimoWithChip[-which(duplicated(tbl.fimoWithChip[, 1:4])),]
      rownames(tbl.fimoWithChip) <- NULL
      sig.chip <- with(tbl.fimoWithChip, sprintf("%s:%d-%d.%s", chrom, start, end, tf))
      sig.fimo <- with(tbl.fimo, sprintf("%s:%d-%d.%s", chrom, start, end, tf))
      chip.hits <- match(sig.chip, sig.fimo)
      if(length(chip.hits) > 0)
         has.ChIP[chip.hits] <- TRUE

      tbl.fimo$chip <- has.ChIP
      obj@state$fimo <- tbl.fimo
      }) # addChIP

#------------------------------------------------------------------------------------------------------------------------
#' query remote khaleesi (ingested POSTAR2) database for rna binding protein binding sites for the target gene
#'
#' @description
#' POSTAR2 curated many clip-seq (and related) data sources; we put them in a postgres database, queried here by gene
#'
#' @rdname addRnaBindingProteins
#'
#' @param obj a TrenaMultiScore object
#' @return data.frame
#'
#' @export
#'
setMethod('addRnaBindingProteins', 'TrenaMultiScore',

    function(obj){
      if(!obj@state$quiet)
        message(sprintf("querying rbp database for %s", obj@targetGene))
      geneRegDB <- dbConnect(PostgreSQL(), user= "trena", password="trena",
                             dbname="genereg2021", host="khaleesi")
      query <- sprintf("select * from rbp where target='%s'", obj@targetGene)
      obj@state$tbl.rbp <- dbGetQuery(geneRegDB, query)
      dbDisconnect(geneRegDB)
      if(!obj@state$quiet)
        message(sprintf("rbp binding sites found: %d", nrow(obj@state$tbl.rbp)))
      invisible(obj@state$tbl.rbp)
      }) # addRnaBindingProteins

#------------------------------------------------------------------------------------------------------------------------
#' return previously queried data.frame of rna binding protein binding sites for the target gene
#'
#' @description
#' POSTAR2 curated many clip-seq (and related) data sources; current value is returned here
#'
#' @rdname getRnaBindingProteins
#'
#' @param obj a TrenaMultiScore object
#' @return data.frame
#'
#' @export
#'
setMethod('getRnaBindingProteins', 'TrenaMultiScore',

    function(obj){
      invisible(obj@state$tbl.rbp)
      }) # getRnaBindingProteins

#------------------------------------------------------------------------------------------------------------------------


