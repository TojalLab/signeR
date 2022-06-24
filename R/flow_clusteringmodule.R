clustering_UI <- function(id) {
  ns <- NS(id)

  tabPanel(
    "Clustering",
    icon = icon("hubspot"),
    fluidRow(
      box(
        width = 12,
        title = span(strong("Unsupervised approaches"),
          style = "font-size:24px"
        ),
        p(
          "Finding sample groups is a usual task in exploratory data analysis.
          Clustering is an unsupervised approach and as such does not require
          additonal sample data (i.e. sample labels). signeRFlow offers
          hierarquical and fuzzy clustering algorithms to group samples
          according to their exposure profiles"
        )
      )
    ),
    fluidRow(
      box(
        title = p("Hierarquical Clustering"),
        width = 12, solidHeader = T, collapsible = T,
        status = "info", collapsed = F,
        fluidRow(
          box(
            width = 12,
            p(
              "signeRFlow generates a dendogram for each generated sample
              of the exposure matrix. Consensus results, i.e. branches that
              are recurrently found, are reported.
              Different distance metrics and clustering algorithms are
              available, please choose below.")
          )
        ),
        column(
          width = 4,
          fluidRow(
            box(
              width = 12, solidHeader = T, collapsible = F, background = "aqua",
              selectInput(
                inputId = ns("mdist"), label = "Method dist:",
                choices = c(
                  "euclidean", "maximum", "manhattan",
                  "canberra", "binary", "minkowski"
                ),
                selected = "euclidean", multiple = FALSE,
                size = 4, selectize = FALSE
              ),
              selectInput(
                inputId = ns("mhclust"), label = "Method hclust:",
                choices = c(
                  "ward.D", "ward.D2", "single", "complete",
                  "average", "mcquitty", "median", "centroid"
                ),
                selected = "average", multiple = FALSE,
                size = 4, selectize = FALSE
              )
            )
          )
        ),
        column(
          width = 8,
          fluidRow(
            box(
              width = 12, solidHeader = T,
              withSpinner(
                plotOutput(ns("hierarquical_plot")),
                color = "#0dc5c1"
              )
            )
          )
        )
      )
    ),
    fluidRow(
      box(
        title = p("Fuzzy Clustering"), width = 12, solidHeader = T,
        collapsible = T, status = "info", collapsed = F,
        fluidRow(
          box(
            width = 12,
            p("
              signeRFlow can apply the Fuzzy C-Means Clustering on each
              generated sample of the exposure matrix. Pertinence levels of
              samples to clusters are averaged over different runs of the
              algorithm. Means are considered  as the final pertinence levels
              and are shown in a heatmap."
            )
          )
        ),
        column(
          width = 4,
          fluidRow(
            box(
              width = 12, solidHeader = T, collapsible = F, background = "aqua",
              sliderInput(ns("liminf"), "Number of group from:", 0, 10, 0, 1),
              sliderInput(ns("limsup"), "Number of group to:", 0, 10, 0, 1),
              p("Set groups to 0 to let the algorithm to estimate."),
              actionButton(ns("startfuzzy"), label = "Run fuzzy", icon = NULL)
            )
          )
        ),
        column(
          width = 8,
          fluidRow(
            box(
              width = 12, solidHeader = T,
              conditionalPanel(
                condition = "input.startfuzzy",
                withSpinner(plotOutput(ns("fuzzy_plot")),
                  color = "#0dc5c1"
                )
              )
            )
          )
        )
      )
    )
  )
}

clustering <- function(input,
                       output,
                       session,
                       signatures) {
  ns <- session$ns

  sigs_obj <- reactive({
    req(signatures())
  })

  mdist <- reactive({
    req(input$mdist)
    return(input$mdist)
  })

  mhclust <- reactive({
    req(input$mhclust)
    return(input$mhclust)
  })

  output$hierarquical_plot <- renderPlot({
    method.dist <- mdist()
    method.hclust <- mhclust()
    sigs <- sigs_obj()
    # print(sigs$SignExposures@samples)
    if (is.null(sigs)) {
      return(NULL)
    }
    HClustExp(
      sigs$SignExposures,
      method.dist = method.dist,
      method.hclust = method.hclust
    )
  })

  limsup <- reactive({
    req(input$limsup)
    if (input$limsup) {
      if (is.null(input$liminf)) {
        input$liminf <- 1
      }
    }
    return(input$limsup)
  })

  fuzzing <- eventReactive(input$startfuzzy, {
    liminf <- input$liminf
    limsup <- limsup()
    sigs <- sigs_obj()
    # print(sigs$SignExposures@samples)
    if (is.null(sigs)) {
      return(NULL)
    }
    if (liminf == 0 & limsup == 0) {
      withProgress(
        message = "Running Fuzzy...",
        detail = "This operation may take a while...",
        value = 0,
        {
          FCE <- FuzzyClustExp(sigs$SignExposures)
        }
      )
    } else {
      clim <- c(liminf, limsup)
      # print(clim)
      # FCE <- FuzzyClustExp(sigs$SignExposures, Clim = clim)
      withProgress(
        message = "Running Fuzzy...",
        detail = "This operation may take a while...",
        value = 0,
        {
          FCE <- FuzzyClustExp(sigs$SignExposures, Clim = clim)
        }
      )
    }
  })

  observeEvent(input$startfuzzy, {
    fuzzing()
  })

  output$fuzzy_plot <- renderPlot({
    input$startfuzzy
    fce <- fuzzing()
    if (!is.null(fce)) {
      heatmap(fce$Meanfuzzy, Colv = NA)
    }
  })
}