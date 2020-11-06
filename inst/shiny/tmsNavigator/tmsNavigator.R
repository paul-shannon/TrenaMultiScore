library(shiny)
library(R6)
library(dataTableWidget)
library(HeatmapWidget)
library(GOEnrichmentWidget)
library(shinyWidgets)   # for pickerInput
source("toHeatmapMatrix.R")
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
goi <- c("MIP", "AFP", "CAMK2N2", "CALR3", "TEX19", "CLCA2", "LINC01602", "CYP4F11",
         "MIR5581", "ENTD1", "LUM", "DRD2", "ZNF560", "MMP12", "TMEM26", "AKAP4")
#----------------------------------------------------------------------------------------------------
DemoApp = R6Class("DemoApp",

    #--------------------------------------------------------------------------------
    private = list(dtw = NULL,
                   heatmap = NULL,
                   goWidget = NULL,
                   tbl.tms = NULL,
                   tbl.sub = NULL,
                   tbl.currentSubset = NULL,
                      # calculated afresh on button click, from current filters:
                   heatmap.mtx = NULL,
                   eventMessage = NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(

        initialize = function(){
            printf("initializing demo")
            filename <- "~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.3.0.250000.500000.RData"
            private$tbl.tms <- get(load(filename))
            private$tbl.sub <- private$tbl.tms
            private$tbl.currentSubset <- private$tbl.tms
            rownames(private$tbl.tms) <- NULL
            private$dtw = dataTableWidget$new(id="dtw", private$tbl.tms,
                                              pageLength=15,
                                              lengthMenu=c(4, 10, 15, 20, 25, 30, 50),
                                              width="1600px", height="1000px")
            private$goWidget <- GOEnrichmentWidget$new(id="goWidget.1", title="GO Enrichment",
                                                      geneSymbols=goi)
            },

        #------------------------------------------------------------
        ui = function(){
           fluidPage(
              tabsetPanel(type="tabs",id="navigationTabSet", selected="openChromatinTab",
                 tabPanel(title="Introduction", value="introTab", includeHTML("intro.html")),
                 tabPanel(title="Open Chromatin", value="openChromatinTab",
                      sidebarLayout(
                        sidebarPanel(
                            pickerInput(inputId = "tfPicker",
                                        label = "Filter on TF",
                                        choices = sort(unique(private$tbl.tms$tf)),
                                        options = list(`actions-box`=TRUE, size=20,
                                                       `selected-text-format`="count > 6"),
                                        multiple = TRUE),
                            pickerInput(inputId = "targetGenePicker",
                                        label = "Filter on target gene",
                                        choices = sort(unique(private$tbl.tms$targetGene)),
                                        options = list(`actions-box`=TRUE, size=20,
                                                       `selected-text-format`="count > 6"),
                                        multiple = TRUE),
                            sliderInput("absDiffbind", "abs(diffbind):", min=0, max=4.0, value=c(0.0, 4.0), step=0.1),
                            sliderInput("absCorrelation", "abs(cor):", min=0, max=1.0, value=c(0.0, 1.0)),
                            sliderInput("absTSS", "log(abs(tss)):", min=0, max=10, value = c(0, 10), step=0.2),
                            sliderInput("motifScore", "motif score:", min=0, max=10, value=c(3, 10), step=0.5),
                            sliderInput("geneHancer", "GeneHancer combined score:", min = 0, max = 700,
                                        value = c(0, 700)),
                            radioButtons(
                                inputId="chipRadioButtons",
                                label="ChIP",
                                choices = c("yes", "no", "either"),
                                selected = "either",
                                inline = TRUE,
                                width = NULL,
                                choiceNames = NULL,
                                choiceValues = NULL),
                            actionButton("displayHeatmapButton", "Display TF/gene binding sites Heatmap"),
                            actionButton("sendGenesToGoWidgetButton", "Genes to GO Enrichment"),
                            width=2), # sidebarPanel
                        mainPanel(
                          private$dtw$ui(),
                          width=10
                          ) # mainPanel
                      ) # sidebarLayout for openChromatin
                      ), # tabPanel
                 tabPanel(title="Heatmap", value="heatmapTab",
                    wellPanel(
                        htmlOutput("heatmapClickReadout"),
                        style="background-color: white; margin-top: 10px; font-size: 18px; width=200px;"),
                    div(id="heatmapContainer", div(id="heatmapWrapper"))
                    ), # heatmap tabPanel
                 tabPanel(title="GO", value="GOtab",
                    private$goWidget$ui()
                    )
                 ) # tabsetPanel
            )},

        #------------------------------------------------------------
        server = function(input, output, session){

           private$dtw$server(input, output, session)
           private$goWidget$server(input, output, session)

           #---------------------------------------------------------------
           observe({   # any and all of the subsetting widgets
              printf("--- filtering triggered")
              tfs <- input$tfPicker

              if(is.null(tfs))
                private$tbl.sub <- private$tbl.tms
              else
                private$tbl.sub <- subset(private$tbl.tms, tf %in% tfs)

              targetGenes <- input$targetGenePicker

              if(!is.null(targetGenes))
                private$tbl.sub <- subset(private$tbl.sub, targetGene %in% targetGenes)

              chipChoice <- input$chipRadioButtons
              private$tbl.sub <- switch(chipChoice,
                                "yes"    = subset(private$tbl.sub, chip),
                                "no"     = subset(private$tbl.sub, !chip),
                                "either" = private$tbl.sub)

              absDiffbind <- input$absDiffbind
              private$tbl.sub <- subset(private$tbl.sub, abs(diffbind.score) >= absDiffbind[1] &
                                         abs(diffbind.score) <= absDiffbind[2])

              absCorrelation <- input$absCorrelation
              private$tbl.sub <- subset(private$tbl.sub, abs(cor) >= absCorrelation[1] & abs(cor) <= absCorrelation[2])

              absTSS <- input$absTSS
              lowerBound <- 10^(absTSS[1])
              upperBound <- 10^(absTSS[2])
              printf("upperBound: %f    lowerBound: %f", lowerBound, upperBound)
              private$tbl.sub <- subset(private$tbl.sub, abs(tss) >= lowerBound & abs(tss) <= upperBound)

              motifScores <- input$motifScore
              private$tbl.sub <- subset(private$tbl.sub, motifScore >= motifScores[1] & motifScore <= motifScores[2])

              geneHancer <- input$geneHancer
              private$tbl.sub <- subset(private$tbl.sub, gh >= geneHancer[1] & gh <= geneHancer[2])

              printf("--- observe, about to setTable of %d rows", nrow(private$tbl.sub))
              private$dtw$setTable(private$tbl.sub)

              #---------------------------------------------------------------------------------------
              observeEvent(input$displayHeatmapButton, ignoreInit=TRUE, {
                 tf.count <- length(unique(private$tbl.sub$tf))
                 targetGene.count <- length(unique(private$tbl.sub$targetGene))
                 printf("--- displayHeatmapButton, targetGenes %d, tfs %d", targetGene.count, tf.count)
                 x <- HeatMapMatrixFromTMSTable$new(private$tbl.sub)
                 private$heatmap.mtx <- x$calculate()
                 printf(" xtab matrix, %d x %d", nrow(private$heatmap.mtx), ncol(private$heatmap.mtx))
                 private$heatmap = HeatmapWidget$new(id="box1",
                                                     title="filtered",
                                                     private$heatmap.mtx,
                                                     width=800,
                                                     height=800)
                 output$heatmapClickReadout <- renderUI(HTML(""))
                 removeUI(selector="#heatmapWrapper", immediate=TRUE)
                 updateTabsetPanel(session, "navigationTabSet", selected="heatmapTab")
                 insertUI(selector="#heatmapContainer", where="beforeEnd",
                          div(id="heatmapWrapper"), immediate=TRUE)
                 insertUI(selector="#heatmapWrapper", where="beforeEnd", private$heatmap$ui(), immediate=TRUE)
                 private$heatmap$server(input, output, session)
                 })

              #---------------------------------------------------------------------------------------
              observeEvent(input$sendGenesToGoWidgetButton, ignoreInit=TRUE, {
                  goi <- unique(private$tbl.sub$targetGene) # [c(1:10, 1000:1010)])
                  count <- length(goi)
                  printf("genes for enrichment: %d", count)
                  printf("tbl.sub: %d, %d", nrow(private$tbl.sub), ncol(private$tbl.sub))
                  updateTabsetPanel(session, "navigationTabSet", selected="GOtab")
                  private$goWidget$setGenes(goi)
                 })

              #---------------------------------------------------------------------------------------
              observeEvent(input$heatmap_click, ignoreInit=TRUE, {
                  x <- input$heatmap_click
                  click.info <- private$heatmap$clickEventAsPosition(x)
                  if(length(click.info) == 0) return()
                  tbl.clickInfo <- as.data.frame(click.info)
                  row <- tbl.clickInfo$row_index[1]
                  col <- tbl.clickInfo$column_index[1]
                    # update the reactiveVal:
                  msg <- sprintf("%d %s binding sites for %s",
                                 private$heatmap.mtx[row, col],
                                 colnames(private$heatmap.mtx)[col],
                                 rownames(private$heatmap.mtx)[row])
                  output$heatmapClickReadout <- renderUI(HTML(msg))
                  })

              #---------------------------------------------------------------------------------------
              observeEvent(input$heatmap_brush, ignoreInit=TRUE, {
                 #printf("=== tmsNavigator sees heatmap brush")
                 #x <- input$heatmap_brush
                 #print(x)
                 })

               })
           } # server
       ) # public
    ) # class
#--------------------------------------------------------------------------------
app <- DemoApp$new()
x <- shinyApp(app$ui, app$server)
runApp(x, port=1156)

