library(shiny)
library(DT)
library(RSQLite)
#----------------------------------------------------------------------------------------------------
db.file.name <- "~/github/TrenaMultiScore/misc/erythropoiesis/marjorieDemos/collectResults/tms.brand105.fimo2.sqlite"
db <- dbConnect(SQLite(), db.file.name)
#----------------------------------------------------------------------------------------------------
#tbl <- get(load("tbx15.RData"))
#tbl <- get(load("tbl.105.RData"))
#tbl$absTSS <- abs(tbl$tss)
#rownames(tbl) <- NULL
coi <- c(
"tf",
"targetGene",
"cor",
"motifScore",
"tss",
"absTSS",
"chip",
#"score",
#"p.value",
"phast7",
"phast30",
"phast100",
"gh",
"annot.type",
"annot.symbol",
"motif_id",
"matched_sequence",
"chrom",
"start",
"end",
"strand")
#tbl <- tbl[, coi]

tfs <- c("any",
         sort(dbGetQuery(db, "select distinct tf from tms ")$tf))

targetGenes <- c("any",
                 sort(dbGetQuery(db, "select distinct targetGene from tms ")$targetGene))
length(targetGenes)
log.max.abs.tss <- 6 # round(log(max(abs(tbl$tss))) + 0.5)
max.gh.score <- 700 # round(max(tbl$gh) + 0.5)
#----------------------------------------------------------------------------------------------------
ui = fluidPage(
   titlePanel("Query trena multi-scored TFBS"),
   sidebarLayout(
      sidebarPanel(
         selectInput("tf", "Transcription Factor:", tfs, selectize=FALSE, selected="TBX15"),
         selectInput("targetGene", "Target Gene:", targetGenes, selectize=FALSE, selected="PDZK1IP1"),
         sliderInput("absCorrelation", "abs(cor):", min = 0, max = 1.0, value = c(0.0, 1.0)),
         sliderInput("absTSS", "log(abs(tss)):", min = 0, max = log.max.abs.tss,
                     value = c(0, log.max.abs.tss)),

         sliderInput("geneHancer", "GeneHancer combined score:", min = 0, max = max.gh.score,
                     value = c(0, max.gh.score)),
         sliderInput("phast30", "PhastCons30:", min = 0, max = 1, value = c(0, 1)),
         radioButtons("ChIP", "ChIP", choices = c("Yes", "No"), selected = "No",
                      inline = TRUE, width = NULL, choiceNames = NULL,  choiceValues = NULL),
         width=3),
      mainPanel(
         DTOutput("table")
         )
      ) # sidebarLayout
   ) # fluidPage

#----------------------------------------------------------------------------------------------------
server = function(session, input, output) {

   output$table <- DT::renderDataTable({
      queryElements <- c()
      if(input$tf != "any")
          queryElements <- c(queryElements, sprintf("tf='%s'", input$tf))
          #tbl <- subset(tbl, tf==input$tf)
      if(input$targetGene != "any")
         queryElements <- c(queryElements, sprintf("targetGene='%s'", input$targetGene))
         #tbl <- subset(tbl, targetGene == input$targetGene)
      min <- input$absCorrelation[1]
      max <- input$absCorrelation[2]
      queryElements <- c(queryElements, sprintf("abs(cor) >= %f AND abs(cor) <= %f", min, max))
      #tbl <- subset(tbl, abs(cor) >= min & abs(cor) <= max)
      min <- (input$absTSS[1])
      max <- (input$absTSS[2])
      #queryElements <- tbl <- subset(tbl, log10(abs(tss)) >= min & log10(abs(tss)) <= max)
      queryElements <- c(queryElements, sprintf("log10(abs(tss)) >= %f AND log10(abs(tss)) <= %f", min, max))

      #c hip <- input$ChIP == "Yes"
      # if(chip) tbl <- subset(tbl, chip)

      #min <- input$geneHancer[1]
      #max <- input$geneHancer[2]
      #tbl <- subset(tbl, gh >= min & gh <= max)

      #min <- input$phast30[1]
      #max <- input$phast30[2]
      #tbl <- subset(tbl, phast30 >= min & phast30 <= max)
      queryElementsCombined <- paste(queryElements, collapse=" AND ")
      queryString <- sprintf("select * from tms where %s", queryElementsCombined)
      printf("queryString: %s", queryString)
      tbl <- dbGetQuery(db, queryString)
      #browser()
      #tbl <- mtcars
      #queryString <- sprintf("select * from tms where %s", queryElementsCombined)
      #                       paste(queryElements, collapse=" AND "))
      #printf("queryString: %s", queryString)
      DT::datatable(tbl, options=list(pageLength=100, dom='<lfip<t>>'), class='nowrap display')
      })
   } # server

#----------------------------------------------------------------------------------------------------
runApp(shinyApp(ui=ui, server=server), port=8888)





#    output$table <- DT::renderDataTable({
#       queryElements <- c()
#       if(input$tf != "any")
#           queryElements <- paste(queryElements, sprintf("tf='%s'", input$tf))
#           #tbl <- subset(tbl, tf==input$tf)
#       if(input$targetGene != "any")
#          queryElements <- paste(queryElements, sprintf("targetGene='%s'", input$targetGene))
#          #tbl <- subset(tbl, targetGene == input$targetGene)
#       min <- input$absCorrelation[1]
#       max <- input$absCorrelation[2]
#       tbl <- subset(tbl, abs(cor) >= min & abs(cor) <= max)
#
#       min <- (input$absTSS[1])
#       max <- (input$absTSS[2])
#       tbl <- subset(tbl, log10(abs(tss)) >= min & log10(abs(tss)) <= max)
#
#       chip <- input$ChIP == "Yes"
#       if(chip) tbl <- subset(tbl, chip)
#
#       min <- input$geneHancer[1]
#       max <- input$geneHancer[2]
#       tbl <- subset(tbl, gh >= min & gh <= max)
#
#       min <- input$phast30[1]
#       max <- input$phast30[2]
#       tbl <- subset(tbl, phast30 >= min & phast30 <= max)
#
#       DT::datatable(tbl, options=list(pageLength=100, dom='<lfip<t>>'), class='nowrap display')
#       })
#
