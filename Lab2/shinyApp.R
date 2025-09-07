############################################################
# Dylan Renard
# JHU Gene Expression Lab 2
# September 6th 2025
# Professor Brandon Higgs
#
# CDC15 Time Course Shiny App (Q8 on Lab 2)
#
# Description:
#   This application visualizes Spellman’s yeast cell cycle dataset
#   (CDC15 temperature-sensitive mutant arrest experiment).
#
# Features:
#   - Reads in spellman.txt and subsets the CDC15 arrays (23–46)
#   - Launches a Shiny web app with interactive controls:
#       * X Variable selector: pick any CDC15 time point
#       * Y Variable selector: pick any CDC15 time point
#       * Point Color selector: customize scatterplot aesthetics
#   - Generates scatter plots across all 6,178 yeast genes
#   - Drops missing values automatically (pairwise complete obs)
#
# Notes:
#   - Based on the Lecture 3 "iris dataset" Shiny demo,
#     but swaps iris columns for CDC15 experiment arrays
#   - Designed for exploratory correlation visualization
#   - No correlation statistics or p-values displayed
#
# Usage:
#   - Save this file as app.R in the same directory as spellman.txt
#   - Run in RStudio: click "Run App" OR use:
#         shiny::runApp("shinyApp.R")
#
# Source Data:
#   Spellman PT, et al. (1998). Comprehensive identification of
#   cell cycle–regulated genes of the yeast Saccharomyces cerevisiae
#   by microarray hybridization. Mol Biol Cell. 9:3273–3297.
############################################################


# ---- Load CDC15 data ----
# Assumes spellman.txt is in your working directory
spell <- read.table("spellman.txt",
                    header = TRUE, row.names = 1, check.names = FALSE,
                    sep = "", stringsAsFactors = FALSE)

# CDC15 experiment: arrays 23–46
dat <- as.matrix(spell[, 23:46])
storage.mode(dat) <- "double"

library(shiny)

# ---- UI ----
ui <- fluidPage(
  titlePanel("CDC15: Time Point vs Time Point"),
  sidebarLayout(
    sidebarPanel(
      selectInput('xcol', 'X Variable', dimnames(dat)[[2]], selected = dimnames(dat)[[2]][1]),
      selectInput('ycol', 'Y Variable', dimnames(dat)[[2]], selected = dimnames(dat)[[2]][2]),
      selectInput('ptcol', 'Point color',
                  choices = c("red", "blue", "black", "darkgreen", "orange"),
                  selected = "red")
    ),
    mainPanel(
      plotOutput('plot', height = "520px")
    )
  )
)

# ---- Server ----
server <- function(input, output) {

  selectedData <- reactive({
    d <- dat[, c(input$xcol, input$ycol), drop = FALSE]
    keep <- complete.cases(d[,1], d[,2])
    data.frame(x = as.numeric(d[keep, 1]),
               y = as.numeric(d[keep, 2]))
  })

  output$plot <- renderPlot({
    d <- selectedData()
    plot(d$x, d$y,
         pch = 16, cex = 0.6, col = input$ptcol,
         xlab = input$xcol, ylab = input$ycol,
         main = "All Genes: Time Point vs. Time Point")
    grid()
  })
}

# ---- Launch ----
shinyApp(ui = ui, server = server)