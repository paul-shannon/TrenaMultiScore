library(RSQLite)
tbl <- get(load("tbl.105.RData"))
dim(tbl)  # 8874058      20
db.file.name <- "tms.brand105.fimo2.sqlite"
db <- dbConnect(SQLite(), db.file.name)

system.time(dbWriteTable(db, name="tms", value=tbl, overwrite=TRUE))  # less than 30 seconds
dbDisconnect(db)

# try it out
db <- dbConnect(SQLite(), db.file.name)
tbl.tbx15 <- dbGetQuery(db, "select * from tms where tf='TBX15'")
dim(tbl.tbx15)  # 45644
