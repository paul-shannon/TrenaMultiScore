library(R6)
library(RUnit)
#----------------------------------------------------------------------------------------------------
HeatMapMatrixFromTMSTable = R6Class("HeatMapMatrixFromTMSTable",

    #--------------------------------------------------------------------------------
    private = list(tmsTable=NULL,
                   quiet=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(

        initialize = function(tmsTable, quiet=TRUE){
           stopifnot(all(c("tf", "targetGene") %in% colnames(tmsTable)))
           private$tmsTable <- tmsTable
           private$quiet <- quiet
           },

        #------------------------------------------------------------
        calculate = function(){
           xtab <- table(private$tmsTable[, c("tf", "targetGene")])
           tbl.x <- as.data.frame(xtab)
           tbl.x$tf <- as.character(tbl.x$tf)
           tbl.x$targetGene <- as.character(tbl.x$targetGene)
           genes <- unique(as.character(tbl.x$targetGene))
           tfs <- unique(as.character(tbl.x$tf))
           mtx <-  matrix(0, nrow=length(genes), ncol=length(tfs), dimnames=list(genes, tfs))
           for(r in seq_len(nrow(tbl.x))){
              tf <- tbl.x$tf[r];
              gene <- tbl.x$targetGene[r];
              count <- tbl.x$Freq[r];
              mtx[gene, tf] <- count
              }

           deleters <- as.integer(which(rowSums(mtx) == 0))
           if(length(deleters) > 0){
              if(!private$quiet) message(sprintf("deleting %d/%d rows without binding sites",
                                                  length(deleters), nrow(mtx)))
              mtx <- mtx[-deleters,]
              }
           if(!private$quiet)
             message(sprintf("returning incidence matrix, %d genes, %d tfs", nrow(mtx), ncol(mtx)))

           invisible(mtx)
           } # calculate

       ) # public

    ) # class
#----------------------------------------------------------------------------------------------------
test_HeatMapMatrixFromTMSTable = R6Class("test_HeatMapMatrixFromTMSTable",

    #--------------------------------------------------------------------------------
    private = list(
       tmsTable = NULL,
       test.10 = function(){
          printf("--- test.10")
          x.10 <- HeatMapMatrixFromTMSTable$new(private$tmsTable[1:10,])
          mtx <- x.10$calculate()
          checkEquals(dim(mtx), c(1, 7))
          checkEquals(rowSums(mtx), c(AFP=10))
          },
       test.random.20 = function(){
          printf("--- test.random.20")
          set.seed(17)
          indices <- sample(seq_len(nrow(private$tmsTable)), 20)
          tbl.random <- private$tmsTable[indices,]
          x <- HeatMapMatrixFromTMSTable$new(tbl.random)
          mtx <- x$calculate()
          checkTrue(min(rowSums(mtx)) > 0)  # all genes returned should have at least one TF
          checkEquals(dim(mtx), c(12, 19))
          },
       test.1000 = function(){
          printf("--- test.random.20")
          set.seed(17)
          indices <- sample(seq_len(nrow(private$tmsTable)), 1000)
          tbl.random <- private$tmsTable[indices,]
          x <- HeatMapMatrixFromTMSTable$new(tbl.random)
          mtx <- x$calculate()
          checkEquals(length(unique(tbl.random$targetGene)), 24)
          row.counts <- rowSums(mtx)
          checkEquals(length(row.counts), 24)
          checkTrue(min(rowSums(mtx)) > 0)  # all genes returned should have at least one TF
          checkTrue(min(colSums(mtx)) > 0)  # all tfs returned should have at least one targetGene
          checkEquals(dim(mtx), c(24, 329))
          }
       ),

    #--------------------------------------------------------------------------------
    public = list(
       initialize = function(){
          dir <- "~/github/TrenaMultiScore/misc/erythropoiesis/marjorieDemos"
          file <- "tbl.3.0.250000.500000.RData"
          full.path <- file.path(dir, file)
          stopifnot(file.exists(full.path))
          private$tmsTable <- get(load(full.path))
          },
       runTests = function(){
          private$test.10()
          private$test.random.20()
          private$test.1000()
          xyx <- 99
          }
    ) # public
  ) # class

#----------------------------------------------------------------------------------------------------
# tester <- test_HeatMapMatrixFromTMSTable$new()
# tester$runTests()

