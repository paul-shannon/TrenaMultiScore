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
#'
#' @import GenomicScores
#' @import phastCons7way.UCSC.hg38
#' @import phastCons100way.UCSC.hg38
#' @import phastCons30way.UCSC.hg38
#' @import MotifDb
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
setGeneric('findOpenChromatin', signature='obj', function(obj, chrom=NA, start=NA, end=NA)
              standardGeneric('findOpenChromatin'))
setGeneric('getOpenChromatin', signature='obj', function(obj) standardGeneric('getOpenChromatin'))
setGeneric('findFimoTFBS', signature='obj', function(obj, motifs=NA, fimo.threshold=NA, chrom=NA, start=NA, end=NA)
              standardGeneric('findFimoTFBS'))
setGeneric('scoreMotifHitsForConservation', signature='obj', function(obj) standardGeneric('scoreMotifHitsForConservation'))
setGeneric('getMultiScoreTable', signature='obj', function(obj) standardGeneric('getMultiScoreTable'))
setGeneric('getTargetGeneInfo', signature='obj', function(obj) standardGeneric('getTargetGeneInfo'))
setGeneric('addDistanceToTSS', signature='obj', function(obj) standardGeneric('addDistanceToTSS'))
#------------------------------------------------------------------------------------------------------------------------
#' Define an object of class TrenaMultiScore
#'
#' @description
#' Expression, variant and covariate data for the genes of interest (perhaps unbounded) for pre-term birth studies
#'
#' @rdname TrenaMultiScore-class
#' @param trenaProject a concrete subclass of TrenaProject
#'
#' @export
#'
#' @return An object of the TrenaMultiScore class
#'

TrenaMultiScore <- function(trenaProject, targetGene, quiet=TRUE)
{
   setTargetGene(trenaProject, targetGene)
   state <- new.env(parent=emptyenv())
   state$openChromatin <- data.frame()
   state$fimo <- data.frame()

   .TrenaMultiScore(trenaProject=trenaProject, targetGene=targetGene, motifDb=MotifDb,
                    state=state)

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
#' GeneHancer has distal enhancers withing ~1MB of the gene; return that region
#'
#' @rdname getGeneHancerRegion
#' @return a data.frame
#'
#' @export
#'
setMethod('getGeneHancerRegion', 'TrenaMultiScore',

     function(obj){
        tbl <- getEnhancers(getProject(obj))
        data.frame(chrom=tbl$chrom[1],
                   start=min(tbl$start) - 1000,
                   end=max(tbl$end) + 1000, stringsAsFactors=FALSE)

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

     function(obj, chrom=NA, start=NA, end=NA){
         if(is.na(chrom)){
            tbl.x <- getGeneHancerRegion(obj)
            chrom <- tbl.x$chrom
            start <- tbl.x$start
            end   <- tbl.x$end
            }
        if("TrenaProjectErythropoiesis" %in% is(getProject(obj))){
           obj@state$openChromatin <- .queryBrandLabATACseq(chrom, start, end)
           message(sprintf("regions of open chromatin: %d", nrow(obj@state$openChromatin)))
           }
        else if("TrenaProjectAD" %in% is(getProject(obj))){
           obj@state$openChromatin <- .queryHintFootprintRegionsFromDatabase("brain_hint_16", chrom, start, end)
           message(sprintf("regions of open chromatin: %d", nrow(obj@state$openChromatin)))
           }
        else stop(sprintf("no support for open chromatin retrieval in %s", getProject(obj)@projectName))
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
#'
#' @return a data.frame
#'
#' @export
#'
setMethod('findFimoTFBS', 'TrenaMultiScore',

    function(obj, motifs=NA, fimo.threshold=NA, chrom=NA, start=NA, end=NA){

       tbl.fp <- getOpenChromatin(obj)

       if(nrow(tbl.fp) == 0)
          stop("TrenaMultiScore::getFimoTFBS error: no open chromatin regions previously identified.")

       if(is.na(motifs))
           motifs <- query(obj@motifDb, c("sapiens", "jaspar2018"))

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
       export(motifs, con=meme.file, format="meme")
       source("~/github/fimoService/batchMode/fimoBatchTools.R")
       tbl.fimo <- fimoBatch(tbl.fp, matchThreshold=fimo.threshold, genomeName="hg38", pwmFile=meme.file)
       obj@state$fimo <- tbl.fimo
       sprintf("tbl.fimo stored with %d rows, %d columns", nrow(tbl.fimo), ncol(tbl.fimo))
       }) # findFimoTFBS

#------------------------------------------------------------------------------------------------------------------------
.queryBrandLabATACseq <- function(chrom.loc, start.loc, end.loc)
{

  tbl.atac <- get(load("~/github/TrenaProjectErythropoiesis/misc/multiScore/brandAtacSeqCuration/tbl.atac.fp.RData"))
  subset(tbl.atac, chrom==chrom.loc & start >= start.loc & end <= end.loc)

} # .queryBrandLabATACseq
#------------------------------------------------------------------------------------------------------------------------
.queryHintFootprintRegionsFromDatabase <- function(database.name, chrom.loc, start.loc, end.loc)
{
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


       tbl.cons <- as.data.frame(gscores(phastCons7way.UCSC.hg38,
                                         GRanges(tbl.fimo[, c("chrom", "start", "end")])), stringsAsFactors=FALSE)
       tbl.fimo$phast7 <- round(tbl.cons$default, digits=2)

       tbl.cons <- as.data.frame(gscores(phastCons30way.UCSC.hg38,
                                         GRanges(tbl.fimo[, c("chrom", "start", "end")])), stringsAsFactors=FALSE)
       tbl.fimo$phast30 <- round(tbl.cons$default, digits=2)

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

      invisible(obj@state$fimo)

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

      coi <- c("tss", "strand", "chrom", "start", "end")
      as.list(getTranscriptsTable(getProject(obj))[coi])
      }) # getTargetGeneInfo

#------------------------------------------------------------------------------------------------------------------------
