library(shiny)
library(R6)
library(dataTableWidget)
library(HeatmapWidget)
library(msgBoxWidget)
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
                   tfGeneCountReadout = NULL,
                   heatmapWidget = NULL,
                   goWidget = NULL,
                   tbl.tms = NULL,
                   tbl.sub = NULL,
                   tbl.currentSubset = NULL,
                      # calculated afresh on button click, from current filters:
                   heatmap.mtx = NULL,
                   eventMessage = NULL,
                   input = NULL,
                   output = NULL,
                   session = NULL

                   ),

    #--------------------------------------------------------------------------------
    public = list(

        initialize = function(){
            printf("initializing demo")
            filename <- "~/github/TrenaProjectErythropoiesis/inst/extdata/genomicRegions/tbl.3.0.250000.500000.RData"
            tbl.tmp <- get(load(filename))
            diffbind.pval.score <- round(-log10(tbl.tmp$diffbind.pval), digits=3)
            tbl.tmp$diffbind.pScore <- diffbind.pval.score
            coi <- c("targetGene", "tf", "cor", "chip", "tss", "diffbind.pScore", "motifScore", "gh",
                     "day0.oc", "day2.oc", "day0.1", "day0.2", "day0.3", "day0.4", "day2.1", "day2.2")
            private$tbl.tms <- tbl.tmp[, coi]
            private$tbl.sub <- private$tbl.tms
            private$tbl.currentSubset <- private$tbl.tms
            rownames(private$tbl.tms) <- NULL
            private$tfGeneCountReadout <- msgBoxWidget$new(id="tfGeneCountReadout",
                                                           title="", boxWidth=200, boxHeight=24,
                                                           fontSize=16, backgroundColor="white")
            private$dtw = dataTableWidget$new(id="dtw", private$tbl.tms,
                                              pageLength=15,
                                              lengthMenu=c(4, 10, 15, 20, 25, 30, 50),
                                              #columnWidths=c(30, 30, 30),
                                              width="1600px", height="1000px")
            private$goWidget <- GOEnrichmentWidget$new(id="goWidget.1", title="GO Enrichment",
                                                      geneSymbols=goi)
            },

        #------------------------------------------------------------
        ui = function(){
           fluidPage(
              titlePanel("TrenaMultiScore:  Exploration of Early Hematopoiesis"),
              tabsetPanel(type="tabs",id="navigationTabSet", selected="openChromatinTab",
                 tabPanel(title="Introduction", value="introTab", includeHTML("intro.html")),
                 tabPanel(title="Open Chromatin", value="openChromatinTab",
                      sidebarLayout(
                        sidebarPanel(
                            private$tfGeneCountReadout$ui(),
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
                            radioButtons(
                                inputId="openChromatinPreference",
                                label="Chromatin Closed Early?",
                                choices = c("yes", "no"),
                                selected = "yes",
                                inline = TRUE,
                                width = NULL,
                                choiceNames = NULL,
                                choiceValues = NULL),
                            sliderInput("day0.oc.slider", "Day 0 open chromatin max reads:",
                                        min=0,
                                        max=20,
                                        value=2, step=1),
                            sliderInput("day2.oc.slider", "Day 2 open chromatin min reads:",
                                        min=0,
                                        max=20,
                                        value=20, step=1),
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
                            actionButton("displayHeatmapButton", "Selection to Heatmap"),
                            actionButton("sendGenesToGoWidgetButton", "Selected Genes to GO"),
                            actionButton("sendGenesAndTFsToGoWidgetButton", "Selected Genes+TFs to GO"),
                            width=2), # sidebarPanel
                        mainPanel(
                          div(private$dtw$ui(), style="margin:20px;"),
                          width=10
                          ) # mainPanel
                      ) # sidebarLayout for openChromatin
                      ), # tabPanel
                 tabPanel(title="Heatmap", value="heatmapTab",
                    wellPanel(
                        div(htmlOutput("heatmapClickReadout"),
                            style="background-color: gray; height: 80px; margin-top: 10px; font-size: 18px; width=200px;")
                        ),
                    div(id="heatmapContainer", div(id="heatmapWrapper"))
                    ), # heatmap tabPanel
                 tabPanel(title="GO", value="GOtab",
                    private$goWidget$ui()
                    )
                 ) # tabsetPanel
            )},

        #------------------------------------------------------------
        server = function(input, output, session){

           private$input <- input
           private$output <- output
           private$session <- session

           private$dtw$server(input, output, session)
           private$goWidget$server(input, output, session)
           private$tfGeneCountReadout$server(input, output, session)

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

                #------------------------------------------------------------
                #  interpret day0, day2 thresholds as min or max
                #  depending on users interest in closed/open, or open/closed
                #------------------------------------------------------------
              closedThenOpenSearch = input$openChromatinPreference == "yes"
              day0.oc.threshold <- input$day0.oc.slider
              day2.oc.threshold <- input$day2.oc.slider
              if(closedThenOpenSearch){
                private$tbl.sub <- subset(private$tbl.sub,
                                          day0.oc <= day0.oc.threshold & day2.oc >= day2.oc.threshold)
              }else{
                private$tbl.sub <- subset(private$tbl.sub,
                                          day0.oc >= day0.oc.threshold & day2.oc <= day2.oc.threshold)
                } # closedThenOpenSearch

              absDiffbind <- input$absDiffbind
              private$tbl.sub <- subset(private$tbl.sub, abs(diffbind.pScore) >= absDiffbind[1] &
                                         abs(diffbind.pScore) <= absDiffbind[2])

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

              private$tfGeneCountReadout$setText(sprintf("%4d targetGenes, %4d tfs",
                                                         length(unique(private$tbl.sub$targetGene)),
                                                         length(unique(private$tbl.sub$tf))))

              private$dtw$setTable(private$tbl.sub)
              }) # observe

           #---------------------------------------------------------------------------------------
             # change the label of the two open chromatin reads, min or max, reflecting user's choice
           observeEvent(input$openChromatinPreference, ignoreInit=TRUE, {
              day0.closedChromatin.day2.open <- input$openChromatinPreference == "yes"
              slider.label.day0 <- "Day 0 open chromatin max reads:"
              slider.label.day2 <- "Day 2 open chromatin min reads:"
              if(!day0.closedChromatin.day2.open){
                 slider.label.day0 <- "Day 0 open chromatin min reads:"
                 slider.label.day2 <- "Day 2 open chromatin max reads:"
                 }
              updateSliderInput(private$session, "day0.oc.slider", label = slider.label.day0)
              updateSliderInput(private$session, "day2.oc.slider", label = slider.label.day2)
              })

           #---------------------------------------------------------------------------------------
           observeEvent(input$displayHeatmapButton, ignoreInit=TRUE, {
              tf.count <- length(unique(private$tbl.sub$tf))
              targetGene.count <- length(unique(private$tbl.sub$targetGene))
              printf("--- displayHeatmapButton, targetGenes %d, tfs %d", targetGene.count, tf.count)
              x <- HeatMapMatrixFromTMSTable$new(private$tbl.sub)
              private$heatmap.mtx <- x$calculate()
              printf(" xtab matrix, %d x %d", nrow(private$heatmap.mtx), ncol(private$heatmap.mtx))
              mtx.xtab <- private$heatmap.mtx
              filename <- sprintf("mtx.xtab.%d.%d.RData", nrow(mtx.xtab), ncol(mtx.xtab))
              printf("saving mtx.xtab to %s", filename)
              save(mtx.xtab, file=filename)
              private$heatmapWidget = HeatmapWidget$new(id="box1",
                                                        title="Coordinate Regulation",
                                                        private$heatmap.mtx,
                                                        rowClusters=5,
                                                        colClusters=5,
                                                        rowTitle="Gene",
                                                        columnTitle="TF",
                                                        width="100%",
                                                        height=700)
              removeUI(selector="#heatmapWrapper", immediate=TRUE)
              updateTabsetPanel(session, "navigationTabSet", selected="heatmapTab")
              insertUI(selector="#heatmapContainer", where="beforeEnd",
                       div(id="heatmapWrapper"), immediate=TRUE)
              insertUI(selector="#heatmapWrapper", where="beforeEnd", private$heatmapWidget$ui(), immediate=TRUE)
              private$heatmapWidget$server(input, output, session)
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
           observeEvent(input$sendGenesAndTFsToGoWidgetButton, ignoreInit=TRUE, {

               goi <- c(unique(private$tbl.sub$targetGene),
                        unique(private$tbl.sub$tf))
               count <- length(goi)
               printf("genes for enrichment: %d", count)
               printf("tbl.sub: %d, %d", nrow(private$tbl.sub), ncol(private$tbl.sub))
               updateTabsetPanel(session, "navigationTabSet", selected="GOtab")
               private$goWidget$setGenes(goi)
              })

           #---------------------------------------------------------------------------------------
           #observeEvent(input$heatmap_click, ignoreInit=TRUE, {
           #    x <- input$heatmap_click
           #    click.info <- private$heatmapWidget$clickEventAsPosition(x)
           #    if(length(click.info) == 0) return()
           #    tbl.clickInfo <- as.data.frame(click.info)
           #    row <- tbl.clickInfo$row_index[1]
           #    col <- tbl.clickInfo$column_index[1]
           #      # update the reactiveVal:
           #    msg <- sprintf("%d %s binding sites for %s",
           #                   private$heatmap.mtx[row, col],
           #                   colnames(private$heatmap.mtx)[col],
           #                   rownames(private$heatmap.mtx)[row])
           #    output$heatmapClickReadout <- renderUI(HTML(msg))
           #    })

           #---------------------------------------------------------------------------------------
           #observeEvent(input$heatmap_brush, ignoreInit=TRUE, {
           #   #printf("=== tmsNavigator sees heatmap brush")
           #   #x <- input$heatmap_brush
           #   #print(x)
           #   })

            # })
           } # server
       ) # public
    ) # class
#--------------------------------------------------------------------------------
app <- DemoApp$new()
x <- shinyApp(app$ui, app$server)
runApp(x, port=1156)

