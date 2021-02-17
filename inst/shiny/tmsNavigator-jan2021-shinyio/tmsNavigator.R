library(shiny)
library(R6)
library(DataTableWidget)
library(HeatmapWidget)
library(MsgBoxWidget)
library(GOEnrichmentWidget)
library(shinyWidgets)   # for pickerInput
source("toHeatmapMatrix.R")
library(rsconnect)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
# goi <- c("MIP", "AFP", "CAMK2N2", "CALR3", "TEX19", "CLCA2", "LINC01602", "CYP4F11",
#         "MIR5581", "ENTD1", "LUM", "DRD2", "ZNF560", "MMP12", "TMEM26", "AKAP4")
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
            filename <- "tbl.all.ocCounts.RData"
            printf("--- loading %s", filename)
            tbl.tmp <- get(load(filename))
            printf("--- load complete")
            coi <- c("targetGene", "tf", "cor", "corEarly", "chip", "tss", "chrom", "start", "end" , "motifScore", "gh",
                     "day0.oc", "day2.oc", "itraq.rise", "phast7", "phast100", "annot.type", "annot.symbol")
                    # , "day0.1", "day0.2", "day0.3", "day0.4", "day2.1", "day2.2", "diffbind.pScore")
            coi <- intersect(coi, colnames(tbl.tmp))
            private$tbl.tms <- tbl.tmp[, coi]
            private$tbl.sub <- private$tbl.tms
            private$tbl.currentSubset <- private$tbl.tms
            rownames(private$tbl.tms) <- NULL
            private$dtw = dataTableWidget$new(id="dtw", private$tbl.tms,
                                              pageLength=15,
                                              lengthMenu=c(4, 10, 15, 20, 25, 30, 50),
                                              #columnWidths=c(30, 30, 30),
                                              width="1600px", height="1000px")
            private$goWidget <- GOEnrichmentWidget$new(id="goWidget.1", title="GO Enrichment",
                                                      geneSymbols=sort(unique(private$tbl.tms$targetGene)))
            },

        #------------------------------------------------------------
        ui = function(){
           fluidPage(
              titlePanel("TrenaMultiScore:  Exploration of Early Hematopoiesis (Binding Sites in Open Chromatin)"),
              tabsetPanel(type="tabs",id="navigationTabSet", selected="openChromatinTab",
                 tabPanel(title="Introduction", value="introTab", includeHTML("intro.html")),
                 tabPanel(title="Open Chromatin", value="openChromatinTab",
                      sidebarLayout(
                        sidebarPanel(
                            #private$tfGeneCountReadout$ui(),
                            h4(textOutput("tfGeneCountReadout")),
                            actionButton("displayHeatmapButton", "Selection to Heatmap"),
                            actionButton("sendGenesToGoWidgetButton", "Selected Genes to GO"),
                            pickerInput(
                                inputId = "targetGenePicker",
                                label = "Select a Target Gene",
                                choices = sort(unique(private$tbl.tms$targetGene)),
                                selected = sort(unique(private$tbl.tms$targetGene)),
                                multiple = TRUE,
                                options = pickerOptions(
                                   actionsBox = TRUE,
                                   liveSearch=TRUE,
                                   `selected-text-format`="count > 4"
                                   )
                            ),
                            pickerInput(
                                inputId = "tfPicker",
                                label = "Select candidate TFs",
                                choices = sort(unique(private$tbl.tms$tf)),
                                selected = sort(unique(private$tbl.tms$tf)),
                                multiple = TRUE,
                                options = pickerOptions(
                                   actionsBox = TRUE,
                                   liveSearch=TRUE,
                                   `selected-text-format`="count > 4"
                                   )
                               ),

                            radioButtons(
                                inputId="chipRadioButtons",
                                label="ChIP",
                                choices = c("yes", "no", "either"),
                                selected = "either",
                                inline = TRUE,
                                width = NULL,
                                choiceNames = NULL,
                                choiceValues = NULL),
                            sliderInput("geneHancer", "GeneHancer combined score:", min = 0, max = 700, value = c(0, 700)),
                            # radioButtons(
                            #     inputId="openChromatinPreference",
                            #     label="Chromatin Closed Early?",
                            #     choices = c("yes", "no"),
                            #     selected = "yes",
                            #     inline = TRUE,
                            #     width = NULL,
                            #     choiceNames = NULL,
                            #     choiceValues = NULL),
                            #sliderInput("day0.oc.slider", "Day 0 open chromatin max reads:",
                            #            min=0,
                            #            max=20,
                            #            value=1, step=0.1),
                            #sliderInput("day2.oc.slider", "Day 2 open chromatin min reads:",
                            #            min=0,
                            #            max=20,
                            #            value=4, step=0.1),
                            #sliderInput("absDiffbind", "abs(diffbind):", min=0, max=4.0, value=c(0.0, 4.0), step=0.1),
                            sliderInput("absCorrelationEarly", "abs(corEarly):", min=0, max=1.0, value=c(0.4, 1.0)),
                            sliderInput("absCorrelation", "abs(cor):", min=0, max=1.0, value=c(0.4, 1.0)),
                            sliderInput("absTSS", "log(abs(tss)):", min=0, max=10, value = c(0, 5), step=0.2),
                            sliderInput("motifScore", "motif score:", min=0, max=10, value=c(3, 10), step=0.5),
                            br(), br(),
                            actionButton("filterTableButton", "Filter by Above Settings"),
                            width=3), # sidebarPanel
                        mainPanel(
                          div(private$dtw$ui(), style="margin:20px;"),
                          width=9
                          ) # mainPanel
                      ) # sidebarLayout for openChromatin
                      ), # tabPanel
                 tabPanel(title="Heatmap", value="heatmapTab",
                    #wellPanel(
                    #    div(htmlOutput("heatmapClickReadout"),
                    #        style="background-color: gray; height: 80px; margin-top: 10px; font-size: 18px; width=200px;")
                    #    ),
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
           output$tfGeneCountReadout <- renderText(sprintf("%4d targetGene/s, %4d tf/s",
                                                           length(unique(private$tbl.tms$targetGene)),
                                                           length(unique(private$tbl.tms$tf))))


           #---------------------------------------------------------------
           #observe({   # any and all of the subsetting widgets
           observeEvent(input$filterTableButton, ignoreInit=TRUE, {

              printf("--- filtering triggered")

                  #-----------------------------------------------------
                  # always start filtering with the entire original table
                  #-----------------------------------------------------
              private$tbl.sub <- private$tbl.tms

              #rowCheck <- nrow(subset(private$tbl.sub, targetGene=="CHD4" & tf=="ZNF410"))
              #printf("rowCheck: %d", rowCheck)
              tfs <- input$tfPicker

              if(!is.null(tfs)){
                  printf("--- tfs is null")
                  private$tbl.sub <- subset(private$tbl.sub, tf %in% tfs)
                  }
              printf("--- tfs: %d", length(tfs))
              #print(tfs)

              #printf("rowCheck 1: %d", nrow(subset(private$tbl.sub, targetGene=="CHD4" & tf=="ZNF410")))

              targetGenes <- input$targetGenePicker
              if(!is.null(targetGenes)){
                private$tbl.sub <- subset(private$tbl.sub, targetGene %in% targetGenes)
                printf("targetGenes is null")
                }

              printf("--- targetGenes: %d", length(targetGenes))
              #print(targetGenes)


              #printf("rowCheck 2: %d", nrow(subset(private$tbl.sub, targetGene=="CHD4" & tf=="ZNF410")))
              chipChoice <- input$chipRadioButtons
              private$tbl.sub <- switch(chipChoice,
                                "yes"    = subset(private$tbl.sub, chip),
                                "no"     = subset(private$tbl.sub, !chip),
                                "either" = private$tbl.sub)

              #printf("rowCheck 3: %d", nrow(subset(private$tbl.sub, targetGene=="CHD4" & tf=="ZNF410")))
                #------------------------------------------------------------
                #  interpret day0, day2 thresholds as min or max
                #  depending on users interest in closed/open, or open/closed
                #------------------------------------------------------------
              closedThenOpenSearch = input$openChromatinPreference == "yes"
              #printf("rowCheck 4: %d", nrow(subset(private$tbl.sub, targetGene=="CHD4" & tf=="ZNF410")))
              #day0.oc.threshold <- input$day0.oc.slider
              #day2.oc.threshold <- input$day2.oc.slider
              #if(closedThenOpenSearch){
              #  private$tbl.sub <- subset(private$tbl.sub,
              #                            day0.oc <= day0.oc.threshold & day2.oc >= day2.oc.threshold)
              #}else{
              #  private$tbl.sub <- subset(private$tbl.sub,
              #                            day0.oc >= day0.oc.threshold & day2.oc <= day2.oc.threshold)
              #  } # closedThenOpenSearch

              #absDiffbind <- input$absDiffbind
              #private$tbl.sub <- subset(private$tbl.sub, abs(diffbind.pScore) >= absDiffbind[1] &
              #                           abs(diffbind.pScore) <= absDiffbind[2])

              absCorrelationEarly <- input$absCorrelationEarly
              private$tbl.sub <- subset(private$tbl.sub, abs(corEarly) >= absCorrelationEarly[1] &
                                                         abs(corEarly) <= absCorrelationEarly[2])
              absCorrelation <- input$absCorrelation
              private$tbl.sub <- subset(private$tbl.sub, abs(cor) >= absCorrelation[1] & abs(cor) <= absCorrelation[2])

              #printf("rowCheck 5: %d", nrow(subset(private$tbl.sub, targetGene=="CHD4" & tf=="ZNF410")))
              absTSS <- input$absTSS
              lowerBound <- 10^(absTSS[1])
              upperBound <- 10^(absTSS[2])
              #printf("lowerBound: %f    upperBound: %f", lowerBound, upperBound)
              private$tbl.sub <- subset(private$tbl.sub, abs(tss) >= lowerBound & abs(tss) <= upperBound)
              #printf("rowCheck 6: %d", nrow(subset(private$tbl.sub, targetGene=="CHD4" & tf=="ZNF410")))

              motifScores <- input$motifScore
              #printf("----------- motifScores")
              #print(motifScores)
              #printf("--- motifScores[1]: %f", motifScores[1])
              #printf("--- motifScores[2]: %f", motifScores[2])
              private$tbl.sub <- subset(private$tbl.sub, motifScore >= motifScores[1] & motifScore <= motifScores[2])

              #printf("rowCheck 7: %d", nrow(subset(private$tbl.sub, targetGene=="CHD4" & tf=="ZNF410")))

              geneHancer <- input$geneHancer
              private$tbl.sub <- subset(private$tbl.sub, gh >= geneHancer[1] & gh <= geneHancer[2])

              #printf("--- observe, about to setTable of %d rows", nrow(private$tbl.sub))

              output$tfGeneCountReadout <- renderText(sprintf("%4d targetGene/s, %4d tf/s",
                                                              length(unique(private$tbl.sub$targetGene)),
                                                              length(unique(private$tbl.sub$tf))))

              private$dtw$setTable(private$tbl.sub)
              tbl.sub <- private$tbl.sub
              #tmp.filename <- sprintf("tbl-%d.rows.RData", nrow(tbl.sub))
              #printf("=== saving tbl.sub to %s", tmp.filename)
              #save(tbl.sub, file=tmp.filename)
              }) # observeEvent: filterTableButton

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
              #updateSliderInput(private$session, "day0.oc.slider", label = slider.label.day0)
              #updateSliderInput(private$session, "day2.oc.slider", label = slider.label.day2)
              })

           #---------------------------------------------------------------------------------------
           observeEvent(input$displayHeatmapButton, ignoreInit=TRUE, {
              tf.count <- length(unique(private$tbl.sub$tf))
              targetGene.count <- length(unique(private$tbl.sub$targetGene))
              printf("--- displayHeatmapButton, targetGenes %d, tfs %d", targetGene.count, tf.count)
              x <- HeatMapMatrixFromTMSTable$new(private$tbl.sub)
              mtx.tmp <- x$calculate()
              printf(" xtab matrix, before, %d x %d", nrow(mtx.tmp), ncol(mtx.tmp))
              if(nrow(mtx.tmp) == 1)
                  mtx.tmp <- rbind(mtx.tmp, mtx.tmp)
              if(ncol(mtx.tmp) == 2){
                 printf("====== doubling columns")
                 mtx.tmp <- cbind(mtx.tmp, mtx.tmp)
                 }
              printf(" xtab matrix, after, %d x %d", nrow(mtx.tmp), ncol(mtx.tmp))
              private$heatmap.mtx <- mtx.tmp
              mtx.xtab <- private$heatmap.mtx
              #filename <- sprintf("mtx.xtab.%d.%d.RData", nrow(mtx.xtab), ncol(mtx.xtab))
              #printf("saving mtx.xtab to %s", filename)
              # save(mtx.xtab, file=filename)
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
deploy <-function()
{
   repos <- options("repos")[[1]]
   stopifnot(sort(names(repos)) == c("BioCann", "BioCsoft", "CRAN"))
   stopifnot(repos$BioCann=="https://bioconductor.org/packages/3.12/data/annotation")
   stopifnot(repos$BioCsoft=="https://bioconductor.org/packages/3.12/bioc")
   stopifnot(repos$CRAN=="https://cran.microsoft.com")

   require(devtools)
   Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
   install_github("PriceLab/BioShiny/MsgBoxWidget",      force=TRUE)
   install_github("PriceLab/BioShiny/HeatmapWidget",      force=TRUE)
   install_github("PriceLab/BioShiny/DataTableWidget",    force=TRUE)
   install_github("PriceLab/BioShiny/GOEnrichmentWidget", force=TRUE)
   install_github("dreamRs/shinyWidgets", force=TRUE)
   install_github("dreamRs/shinybusy", force=TRUE)
   #install_github("paul-shannon/igvShiny", force=TRUE)
   #install_github("PriceLab/BioShiny/igvWidget", force=TRUE)
   #install_github("PriceLab/BioShiny/GenomeTracksWidget", force=TRUE)


   require(rsconnect)
   deployApp(account="hoodlab",
             appName="tms-hematopoiesis",
             appTitle="TrenaMultiScore early Hematopoiesis",
             appFiles=c("tmsNavigator.R",
                        "toHeatmapMatrix.R",
                        "tbl.all.ocCounts.RData",
                        "intro.html"
                        ),
             appPrimaryDoc="tmsNavigator.R",
             forceUpdate=TRUE
             )

} # deploy
#------------------------------------------------------------------------------------------------------------------------
app <- DemoApp$new()
if(grepl("hagfish", Sys.info()[["nodename"]]) & !interactive()){
   runApp(shinyApp(app$ui(), app$server), port=1119)
   } else {
   shinyApp(app$ui(), app$server)
   }
