############################################################
# Dylan Renard
# Bioinformatics Lab 2 – Q8 (with tasteful extras)
# September 2025 • Professor Brandon Higgs
#
# CDC15 Time Course Shiny App
#
# Required feature:
#  - Choose any two CDC15 time points (arrays 23–46) and
#    plot all genes as a scatter plot (no stats shown).
#
# Optional extras (clearly labeled):
#  - Violin: distributions of the two selected time points
#  - Dendrogram: clusters of CDC15 time points (1 - Pearson r)
#  - Circos: chord diagram of strongest time-point correlations
#
# Usage:
#  - Save as app.R in the same folder as spellman.txt
#  - Run: shiny::runApp("app.R") or click "Run App" in RStudio
############################################################

# ---- Load data ----
spell <- read.table("spellman.txt",
                    header = TRUE, row.names = 1, check.names = FALSE,
                    sep = "", stringsAsFactors = FALSE)

# CDC15 experiment: arrays 23–46 (genes x timepoints)
dat <- as.matrix(spell[, 23:46])
storage.mode(dat) <- "double"

tp_choices <- colnames(dat)

# ---- Optional: precompute pieces for extras ----
# Dendrogram (time points): distance = 1 - Pearson correlation
cor_mat <- cor(dat, use = "pairwise.complete.obs", method = "pearson")
dist_tp <- as.dist(1 - cor_mat)
hc_tp   <- hclust(dist_tp, method = "average")

# ---- Shiny ----
library(shiny)

# helpers to avoid hard package dependency for extras
.has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

ui <- fluidPage(
  titlePanel("CDC15: Time Point vs Time Point"),
  sidebarLayout(
    sidebarPanel(
      selectInput('xcol', 'X Variable', tp_choices, selected = tp_choices[1]),
      selectInput('ycol', 'Y Variable', tp_choices, selected = tp_choices[2]),
      selectInput('ptcol', 'Point color',
                  choices = c("red","blue","black","darkgreen","orange"),
                  selected = "red"),
      hr(),
      # Circos controls (for extras tab)
      strong("Circos options (optional)"),
      sliderInput("corThresh", "Min |r| to draw edge:",
                  min = 0.0, max = 1.0, value = 0.85, step = 0.01),
      numericInput("maxEdges", "Max edges:", value = 80, min = 10, max = 300, step = 10),
      helpText("Circos uses time-point correlations; requires {circlize}."),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Scatter",
          br(),
          plotOutput('scatter', height = "520px")
        ),
        tabPanel("Violin (extra)",
          br(),
          uiOutput("violin_msg"),
          plotOutput('violin', height = "520px")
        ),
        tabPanel("Dendrogram (extra)",
          br(),
          plotOutput('dendro', height = "520px")
        ),
        tabPanel("Circos (extra)",
          br(),
          uiOutput("circos_msg"),
          plotOutput('circos', height = "600px")
        )
      ),
      width = 9
    )
  )
)

server <- function(input, output, session) {

  # ---- Core reactive: selected columns as complete pairs across genes ----
  selectedData <- reactive({
    d <- dat[, c(input$xcol, input$ycol), drop = FALSE]
    keep <- complete.cases(d[,1], d[,2])
    data.frame(x = as.numeric(d[keep, 1]),
               y = as.numeric(d[keep, 2]))
  })

  # ---- Required: Scatter ----
  output$scatter <- renderPlot({
    d <- selectedData()
    plot(d$x, d$y,
         pch = 16, cex = 0.6, col = input$ptcol,
         xlab = input$xcol, ylab = input$ycol,
         main = "All Genes: Time Point vs. Time Point")
    grid()
  })

  # ---- Extra: Violin (needs ggplot2) ----
  output$violin <- renderPlot({
    if (!.has_pkg("ggplot2")) return(invisible(NULL))
    d <- selectedData()
    df_long <- data.frame(
      value = c(d$x, d$y),
      tp    = factor(rep(c(input$xcol, input$ycol), each = nrow(d)),
                     levels = c(input$xcol, input$ycol))
    )
    ggplot2::ggplot(df_long, ggplot2::aes(x = tp, y = value)) +
      ggplot2::geom_violin(fill = input$ptcol, alpha = 0.6, trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
      ggplot2::labs(x = "Time point", y = "Expression (log2 ratio)",
                    title = "Distributions of Selected Time Points") +
      ggplot2::theme_minimal(base_size = 14)
  })

  output$violin_msg <- renderUI({
    if (.has_pkg("ggplot2")) return(NULL)
    HTML("<em>Install <code>ggplot2</code> to view the violin plot: <code>install.packages('ggplot2')</code></em>")
  })

  # ---- Extra: Dendrogram (precomputed) ----
  output$dendro <- renderPlot({
    plot(hc_tp, main = "CDC15 Time-Point Dendrogram (1 - Pearson r)",
         xlab = "", sub = "")
    grid()
  })

  # ---- Extra: Circos (needs circlize) ----
  output$circos <- renderPlot({
    if (!.has_pkg("circlize")) return(invisible(NULL))

    # Build edge list from correlation matrix with |r| >= threshold
    r <- cor_mat
    diag(r) <- NA
    thr <- input$corThresh
    # Take upper triangle indices
    upper <- which(upper.tri(r), arr.ind = TRUE)
    edges <- data.frame(
      from = colnames(r)[upper[,1]],
      to   = colnames(r)[upper[,2]],
      w    = r[upper]
    )
    # filter by threshold and limit to top edges by absolute weight
    edges <- edges[!is.na(edges$w) & abs(edges$w) >= thr, , drop = FALSE]
    if (nrow(edges) == 0) {
      plot.new(); title("No edges at current |r| threshold")
      return(invisible(NULL))
    }
    edges <- edges[order(-abs(edges$w)), , drop = FALSE]
    if (nrow(edges) > input$maxEdges) edges <- edges[seq_len(input$maxEdges), , drop = FALSE]

    # Draw chord diagram
    circlize::circos.clear()
    op <- par(mar = c(1,1,1,1))
    on.exit({ par(op); circlize::circos.clear() }, add = TRUE)
    circlize::chordDiagram(
      x = edges[, c("from","to","w")],
      directional = 0,
      grid.col = setNames(rep("grey30", length(tp_choices)), tp_choices),
      transparency = 0.3
    )
    title(main = sprintf("Circos: strongest correlations (|r| ≥ %.2f)", thr))
  })

  output$circos_msg <- renderUI({
    if (.has_pkg("circlize")) return(NULL)
    HTML("<em>Install <code>circlize</code> to view the circos plot: <code>install.packages('circlize')</code></em>")
  })
}

shinyApp(ui = ui, server = server)