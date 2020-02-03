files <- list.files(pattern=".RData")
tbls <- list()
for(file in files){
   tbl <- get(load(file))
   tbls[[file]] <- tbl
   }

tbl.fimo3 <- do.call(rbind, tbls)
rownames(tbl.fimo3) <- NULL
dim(tbl.fimo3)  # fimo3: 741598  fimo4: 94560
save(tbl.fimo3, file="tenGenesFimo3.RData")

