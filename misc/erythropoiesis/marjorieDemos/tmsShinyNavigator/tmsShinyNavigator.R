library(shiny)
library(R6)
library(dataTableWidget)
library(shinyWidgets)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
DemoApp = R6Class("DemoApp",

    #--------------------------------------------------------------------------------
    private = list(dtw = NULL,
                   tbl.tms = NULL,
                   tbl.currentSubset = NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(

        initialize = function(){
            printf("initializing demo")
            private$tbl.tms <- get(load("tbl.afp.srm.RData"))
            private$tbl.currentSubset <- private$tbl.tms
            rownames(private$tbl.tms) <- NULL
            private$dtw = dataTableWidget$new(id="dt", private$tbl.tms,
                                              pageLength=15,
                                              lengthMenu=c(4, 10, 15, 20, 25, 30, 50),
                                              width="1600px", height="1000px")
            },

        #------------------------------------------------------------
        ui = function(){
           fluidPage(
              titlePanel("Explore Regulation of Early Hematopoiesis Genes"),
              sidebarLayout(
                  sidebarPanel(
                   pickerInput(inputId = "tfPicker",
                                label = "Filter on TF",
                                choices = sort(unique(private$tbl.tms$tf)),
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
                      width=2), # sidebarPanel
               mainPanel(
                   private$dtw$ui(),
                   width=10
                 )
              ) # sidebarLayout
            )},

        #------------------------------------------------------------
        server = function(input, output, session){

           private$dtw$server(input, output, session)

           observe({
              printf("--- filtering triggered")
              tfs <- input$tfPicker

              if(is.null(tfs))
                tbl.sub <- private$tbl.tms
              else
                tbl.sub <- subset(private$tbl.tms, tf %in% tfs)

              chipChoice <- input$chipRadioButtons
              tbl.sub <- switch(chipChoice,
                                "yes"    = subset(tbl.sub, chip),
                                "no"     = subset(tbl.sub, !chip),
                                "either" = tbl.sub)

              absDiffbind <- input$absDiffbind
              tbl.sub <- subset(tbl.sub, abs(diffbind.score) > absDiffbind[1] &
                                         abs(diffbind.score) < absDiffbind[2])

              absCorrelation <- input$absCorrelation
              tbl.sub <- subset(tbl.sub, abs(cor) > absCorrelation[1] & abs(cor) < absCorrelation[2])

              absTSS <- input$absTSS
              lowerBound <- 10^(absTSS[1])
              upperBound <- 10^(absTSS[2])
              printf("upperBound: %f    lowerBound: %f", lowerBound, upperBound)
              tbl.sub <- subset(tbl.sub, abs(tss) >= lowerBound & abs(tss) <= upperBound)

              motifScores <- input$motifScore
              tbl.sub <- subset(tbl.sub, motifScore > motifScores[1] & motifScore < motifScores[2])

              geneHancer <- input$geneHancer
              tbl.sub <- subset(tbl.sub, gh > geneHancer[1] & gh < geneHancer[2])

              printf("--- observe, about to setTable of %d rows", nrow(tbl.sub))
              private$dtw$setTable(tbl.sub)
              })
           } # server
       ) # public
    ) # class
#--------------------------------------------------------------------------------
app <- DemoApp$new()
x <- shinyApp(app$ui, app$server)
runApp(x, port=1156)

