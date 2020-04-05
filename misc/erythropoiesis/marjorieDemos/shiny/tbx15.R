directory <- "~/github/TrenaMultiScore/misc/erythropoiesis/marjorieDemos/fimo2"
rdata.files <- list.files(path=directory, pattern=".RData")

tbls <- lapply(head(rdata.files, n=-1), function(file){
    full.path <- file.path(directory, file)
    if(file.size(full.path) > 10000){
        tbl <- get(load(full.path))
        tbl.tbx15 <- subset(tbl, tf=="TBX15")
        return(tbl.tbx15)
        }})
tbl.tbx15 <- do.call(rbind, tbls)
dim(tbl.tbx15) # 40571 20
save(tbl.tbx15, file="tbx15.RData")

# table(subset(tbl.tbx15, gh > 500 & abs(cor) > 0.5 & abs(tss) < 500 & phast30 > 0.9)$targetGene)
#     BCL6    CPEB4    FOXJ3   MIR144  MIR4732    NATD1     OSR2 PDZK1IP1  RUNDC3A  SLC30A1     SPI1
#       31       10        2        4        4        3        6        6       39        3       15


