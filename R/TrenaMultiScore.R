#----------------------------------------------------------------------------------------------------
#' @import methods
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
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
                                targetGene="character"
                                )
                             )

#----------------------------------------------------------------------------------------------------
setGeneric('getProject', signature='obj', function(obj) standardGeneric('getProject'))
setGeneric('getGeneHancerRegion', signature='obj', function(obj) standardGeneric('getGeneHancerRegion'))
#----------------------------------------------------------------------------------------------------
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
   .TrenaMultiScore(trenaProject=trenaProject, targetGene=targetGene)

} # TrenaMultiScore, the constructor
#----------------------------------------------------------------------------------------------------
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

#----------------------------------------------------------------------------------------------------
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

#----------------------------------------------------------------------------------------------------
