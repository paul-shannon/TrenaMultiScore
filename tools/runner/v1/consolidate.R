library(RUnit)
library(org.Hs.eg.db)
library(httr)
library(jsonlite)


files = list.files(pattern=".RData")
length(files) # 75

tbls.spearman <- list()
tbls.lasso <- list()
tbls.ridge <- list()
tbls.rf    <- list()
tbls.boost  <- list()

for(filename in files){
    x <- get(load(filename))
    class(x)
    tbls <- x$get.lm.tables()
    length(tbls)
    names(tbls)
    geneName <- strsplit(filename, "-")[[1]][1]
    tbls[[1]]
    checkTrue(all(c("spearmanCoeff", "betaLasso", "betaRidge", "rfScore", "xgboost") %in% names(tbls)))
    tbl.spearman <- tbls[["spearmanCoeff"]]
    tbl.lasso <- tbls[["betaLasso"]]
    tbl.ridge <- tbls[["betaRidge"]]
    tbl.rf    <- tbls[["rfScore"]]
    tbl.boost  <- tbls[["xgboost"]]

    tbl.spearman$reg <- rownames(tbl.spearman)
    tbl.lasso$reg <- rownames(tbl.lasso)
    tbl.ridge$reg <- rownames(tbl.ridge)
    tbl.rf$reg <- rownames(tbl.rf)
    tbl.boost$reg <- rownames(tbl.boost)

    tbl.spearman$targetGene <- geneName
    tbl.lasso$targetGene <- geneName
    tbl.ridge$targetGene <- geneName
    tbl.rf$targetGene <- geneName
    tbl.boost$targetGene <- geneName

    tbls.spearman[[geneName]] <- tbl.spearman
    tbls.lasso[[geneName]] <- tbl.lasso
    tbls.ridge[[geneName]] <- tbl.ridge
    tbls.rf[[geneName]] <- tbl.rf
    tbls.boost[[geneName]] <- tbl.boost
}

tbl.spearman <- do.call(rbind, tbls.spearman)
dim(tbl.spearman) # 1120 9
rownames(tbl.spearman) <- NULL
tbl.spearman.sub <- subset(tbl.spearman, p.value <= 0.2)

tbl.freq <- as.data.frame(sort(table(tbl.spearman.sub$reg), decreasing=TRUE))
as.character(subset(as.data.frame(sort(table(tbl.spearman.sub$reg), decreasing=TRUE)), Freq > 3)$Var1)

# "GEMIN5" "RXRA"   "UCHL5"  "PATZ1"  "TCF3"   "U2AF2"

gemin5.targets <- subset(tbl.spearman.sub, reg=="GEMIN5")$targetGene
gemin5.targets

goi <- gemin5.targets
goi.string <- toJSON(goi)
uri <- sprintf("http://localhost:8000/goEnrich")
body.jsonString <- sprintf('%s', toJSON(list(geneSymbols=goi)))

r <- POST(uri, body=body.jsonString)

   #sprintf('{"geneSymbols": "%s"}', goi.string))
tbl.go <- fromJSON(content(r)[[1]])
head(tbl.go)

