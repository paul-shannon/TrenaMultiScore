biocGet <- function(pkgs){
   library(BiocManager)
   BiocManager::install(pkgs)
   }

printf <- function(...) print(noquote(sprintf(...)))
code.pkgs <- c(
    'annotatr',
    'GenomicScores',
    'phastCons7way.UCSC.hg38',
    # 'phastCons30way.UCSC.hg38',  loaded from extdata, version/AH problem
    'phastCons100way.UCSC.hg38')

for(code.pkg in code.pkgs){
   suppressWarnings(
      needed <- !require(code.pkg, character.only=TRUE)
      )
   printf("%s needed? %s", code.pkg, needed)
   if(needed)
      biocGet(code.pkg)
   } # for
