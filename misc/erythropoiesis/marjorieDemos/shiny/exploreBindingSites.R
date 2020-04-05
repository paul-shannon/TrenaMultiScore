library(shiny)
library(DT)

tbl <- get(load("tbx15.RData"))
tbl$absTSS <- abs(tbl$tss)
rownames(tbl) <- NULL
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
tbl <- tbl[, coi]

tfs <- c("any", sort(unique(tbl$tf)))
targetGenes <- c("any",sort(unique(tbl$targetGene)))
length(targetGenes)
log.max.abs.tss <- round(log(max(abs(tbl$tss))) + 0.5)
max.gh.score <- round(max(tbl$gh) + 0.5)
#----------------------------------------------------------------------------------------------------
ui = fluidPage(
   titlePanel("Query trena multi-scored TFBS"),
   sidebarLayout(
      sidebarPanel(
         selectInput("tf", "Transcription Factor:", tfs, selectize=FALSE),
         selectInput("targetGene", "Target Gene:", targetGenes, selectize=FALSE),
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
         DTOutput("table") #, width=200)
         )
      ) # sidebarLayout
   ) # fluidPage

#----------------------------------------------------------------------------------------------------
server = function(session, input, output) {

   output$table <- DT::renderDataTable({
      if(input$tf != "any")
         tbl <- subset(tbl, tf==input$tf)
      if(input$targetGene != "any")
         tbl <- subset(tbl, targetGene == input$targetGene)
      min <- input$absCorrelation[1]
      max <- input$absCorrelation[2]
      tbl <- subset(tbl, abs(cor) >= min & abs(cor) <= max)

      min <- (input$absTSS[1])
      max <- (input$absTSS[2])
      printf("all abs(tss) > %f and < %f", min, max)
      tbl <- subset(tbl, log10(abs(tss)) >= min & log10(abs(tss)) <= max)

      chip <- input$ChIP == "Yes"
      if(chip) tbl <- subset(tbl, chip)

      min <- input$geneHancer[1]
      max <- input$geneHancer[2]
      tbl <- subset(tbl, gh >= min & gh <= max)

      min <- input$phast30[1]
      max <- input$phast30[2]
      tbl <- subset(tbl, phast30 >= min & phast30 <= max)

      tbl
      },
      class='nowrap display',
      options=list(pageLength=100)
      )

    #observeEvent(input$absCorrelation, ignoreInit=TRUE, {
    #    printf("%f - %f", min, max)
    #    })

   } # server

#----------------------------------------------------------------------------------------------------
runApp(shinyApp(ui=ui, server=server), port=8888)




