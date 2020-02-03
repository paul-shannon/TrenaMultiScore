files <- list.files(pattern=".RData")
tbls <- list()
for(file in files){
   tbl <- get(load(file))
   tbls[[file]] <- tbl
   }

tbl.fimo4 <- do.call(rbind, tbls)
rownames(tbl.fimo4) <- NULL
dim(tbl.fimo4)  # 94560

save(tbl.fimo4, file="tenGenesFimo4.RData")

