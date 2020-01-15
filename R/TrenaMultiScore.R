# 1) find tss, chrom, strand, min and max of genehancer regions (aka, "geneRegion")
# 2) load open chromatin in geneRegion: atac-seq or database footprint regions
# 3) find all TFBS in atac-seq regions
# 4) score all TFBS for conservation
# 5) determine distance from TSS for all TFBS
#------------------------------------------------------------------------------------------------------------------------
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
setGeneric('getFimoTFBS', signature='obj', function(obj, motifs=NA, fimo.threshold=NA, chrom=NA, start=NA, end=NA)
              standardGeneric('getFimoTFBS'))
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
           }
        else if("TrenaProjectAD" %in% is(getProject(obj))){
           obj@state$openChromatin <- .queryHintFootprintRegionsFromDatabase("brain_hint_16", chrom, start, end)
           }
         else stop(sprintf("no support for open chromatin retrieval in %s", getProject(ob)@projectName))
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
#' @rdname getGeneHancerRegion
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
setMethod('getFimoTFBS', 'TrenaMultiScore',

    function(obj, motifs=NA, fimo.threshold=NA, chrom=NA, start=NA, end=NA){

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

       tbl.fp <- getOpenChromatin(obj)
       meme.file <- "tmp.meme"
       export(motifs, con=meme.file, format="meme")
       source("~/github/fimoService/batchMode/fimoBatchTools.R")
       fimoBatch(tbl.fp, matchThreshold=fimo.threshold, genomeName="hg38", pwmFile=meme.file)
       }) # getFimoTFBS

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
