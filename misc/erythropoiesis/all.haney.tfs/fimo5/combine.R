files <- dir(pattern="*.RData")
print(load("BPTF.RData"))
tbls <- list()
for(file in files){
   name <- sub(".RData", "", file)
   tbls[[name]] <- get(load(file))
   }
tbl <- do.call(rbind, tbls)
rownames(tbl) <- NULL
dim(tbl)
