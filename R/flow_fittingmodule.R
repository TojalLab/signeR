fitting_UI <- function(id) {
  ns <- NS(id)

  tagList(
    tabBox(
      id = "fitting", width = 12,
      tabPanel(
        "Fitting",
        icon = icon("brain"),
        fluidRow(
          box(
            width = 12,
            p(
              "Please ", strong("Upload your data"),
              " below. The counts of mutations are required, and should be
              organized in a matrix with 96 columns (corresponding to mutations
              types) and one line for each genome sample. Opptionally, a matrix
              with matching opportunities can be uploaded (",
              a(
                "see signeR documentation for details",
                href = "https://www.google.com/"
              ),
              "). signeR will estimate sample's exposures to known signatures,
              which can also be uploaded (a matrix with 96 lines, corresponding
              to mutations types, and one column for each signature).
              Oppitionally, a set of Cosmic signatures can be chosen.
              signeRFlow will estimate sample's exposure data and display them
              through interactive visualizations."
            )
          )
        ),
        fluidRow(
          box(
            title = "Upload data", width = 12, solidHeader = T,
            collapsible = T, status = "primary",
            messageBox(
              width = 12,
              "Upload a SNV matrix file with your own samples
              to use in signeR fitting module and previous known signatures
              (mandatories files).
              You can upload an opportunity file as well."
            ),
            fluidRow(
              box(
                width = 4, background = "orange",
                tags$head(tags$script(src = "message-handler.js")),
                  actionButton(ns("snvhelp"),
                    "SNV matrix help",
                    icon = icon("info-circle")
                  ),
                hr(),
                fileInput(ns("mutfile_fit"),
                  "SNV matrix*",
                  multiple = FALSE,
                  accept = c(
                    "text/csv", "text/comma-separated-values",
                    "text/plain", ".csv"
                  )
                )
              ),
              box(
                width = 4, background = "navy",
                tags$head(tags$script(src = "message-handler.js")),
                actionButton(ns("opphelp"),
                  " Opportunity help",
                  icon = icon("info-circle")
                ),
                hr(),
                fileInput(ns("oppfile_fit"),
                  "Opportunities",
                  multiple = FALSE,
                  accept = c(
                    "text/csv", "text/comma-separated-values",
                    "text/plain", ".csv"
                  )
                )
              ),
              box(
                width = 4, background = "orange",
                tags$head(tags$script(src = "message-handler.js")),
                actionButton(ns("knownsighelp"),
                  " Previous signatures help",
                  icon = icon("info-circle")
                ),
                hr(),
                fileInput(ns("kP"),
                  "Previous signatures*",
                  multiple = FALSE,
                  accept = c(
                    "text/csv", "text/comma-separated-values",
                    "text/plain", ".csv"
                  )
                )
              ),
            ),
          )
        ),
        fluidRow(
          box(
            width = 12, solidHeader = T, collapsible = F,
            fluidRow(
              box(
                width = 12,
                p(
                  strong(
                    "The algorithm computes a Bayesian approach to the fitting
                    of mutation counts to known mutational signatures, thus
                    estimating  exposures to mutational processes. If provided,
                    opportunities are used as weights for the factorization.
                    Further analysis parameters can be set below. Estimated
                    exposures can be visualized on different plots."
                  )
                )
              )
            ),
            tags$head(tags$script(src = "message-handler.js")),
            column(
              width = 4,
              fluidRow(
                box(
                  width = 12, solidHeader = T,
                  collapsible = F, background = "aqua",
                  fluidRow(
                    box(
                      width = 12, solidHeader = F,
                      collapsible = F, collapsed = F, background = "aqua",
                      div(
                        style = "display: inline-block;",
                        HTML("<b>Iterations</b>")
                      ),
                      div(
                        style = "display: inline-block;",
                        actionLink(ns("iterationhelp_fit"),
                          label = icon("info-circle")
                        )
                      ),
                      HTML("<br>"),
                      div(
                        style = "display: inline-block;width: 30%;",
                        numericInput(ns("em_iter_fit"), "EM", 10)
                      ),
                      div(
                        style = "display: inline-block;width: 30%;",
                        numericInput(ns("warm_fit"), "Warm-up", 10)
                      ),
                      div(
                        style = "display: inline-block;width: 30%;",
                        numericInput(ns("eval_fit"), "Final", 10)
                      )
                    )
                  ),
                  # hr(),
                  actionButton(
                    ns("startfitting"),
                    label = "Start Fitting analysis", icon = NULL
                  ),
                  uiOutput(ns("dwfitting"))
                ),
              )
            ),
            column(
              width = 8,
              fluidRow(
                box(
                  width = 12, solidHeader = T,
                    withSpinner(plotOutput(ns("bic_fitting")))
                )
              )
            )
          )
        ),
        fluidRow(
          # colocar uma minuatura dos gráficos ao invés de combobox
          box(
            title = "Plots", width = 12, solidHeader = T,
            collapsible = T, status = "primary", collapsed = F,
            column(
              width = 4,
              fluidRow(
                box(
                  width = 12, solidHeader = T,
                  collapsible = F, background = "aqua",
                  selectInput(
                    inputId = ns("whichplot_fit"), label = "Available Plots:",
                    choices = c(
                      "SignHeat", "ExposureBoxplot",
                      "ExposureBarplot", "ExposureBarplotRelative",
                      "ExposureHeat"
                    ),
                    selected = "SignHeat", multiple = FALSE,
                    size = 5, selectize = FALSE
                  )
                ),
                box(
                  width = 12, solidHeader = T,
                  collapsible = F, background = "olive",
                  uiOutput(ns("plot_help_fit"))
                )
              )
            ),
            column(
              width = 8,
              fluidRow(
                box(
                  width = 12, solidHeader = T,
                  withSpinner(
                    plotOutput(ns("choosen_plot_fit"),
                      width = "auto", height = "400px"
                    ),
                    color = "#0dc5c1"
                  )
                )
              )
            )
          )
        )
      ),
      clustering_UI(ns("clusteringmod")),
      covariate_UI(ns("covariatemod"))
    )
  )
}

fitting <- function(input,
                    output,
                    session,
                    width) {
  ns <- session$ns

  callModule(
    clustering,
    "clusteringmod",
    reactive(signatures_fitting())
  )

  callModule(
    covariate,
    "covariatemod",
    reactive(signatures_fitting())
  )

  mut_fit <- reactive({
    # req(input$mutfile_fit)
    if (is.null(input$mutfile_fit)) {
      showModal(modalDialog(
        title = "Oh no!",
        paste0("You have not uploaded a file, silly person!"),
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }
    df <- read.table(input$mutfile_fit$datapath, header=T,sep="\t",row.names=1,check.names=F)
    if (!validate_cnv(df)) {
      showModal(modalDialog(
        title = "Oh no!",
        paste0("You must upload a valid SNV matrix file."),
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }
    return(df)
  })

  kp <- reactive({
    # req(input$kP)
    if (is.null(input$kP)) {
      showModal(modalDialog(
        title = "Oh no!",
        paste0("You have not uploaded a file, silly person!"),
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }
    # read.table(input$kP$datapath, sep = "\t")
    df <- read.table(input$kP$datapath, header=T,sep="\t",row.names=1)
    if (!validate_knownsig(df)) {
      showModal(modalDialog(
        title = "Oh no!",
        paste0("You must upload a valid Previous signature file."),
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }
    return(df)
  })

  opp_fit <- reactive({
    # isso torno o input obrigatorio
    # req(input$oppfile_fit)
    if (is.null(input$oppfile_fit)) {
      return(NULL)
    }
    read.table(input$oppfile_fit$datapath)
  })

  observeEvent(input$iterationhelp_fit, {
    showModal(modalDialog(
      title = "Iterations parameters",
      includeMarkdown(
        system.file("extdata", "iterations.md", package = "signeR")
      ),
      size = "l", easyClose = TRUE
    ))
  })

  whichplot_fit <- reactive({
    req(input$whichplot_fit)
    if (is.null(input$whichplot_fit)) {
      return(NULL)
    }
    return(input$whichplot_fit)
  })

  output$plot_help_fit <- renderUI({
    req(input$whichplot_fit)
    if (!is.null(whichplot_fit())) {
      if (whichplot_fit() == "SignPlot") {
        p(
          strong("SignPlot help"), HTML("<br>"), HTML("<br>"),
          "Signatures barplot with error bars reflecting the sample percentiles
          0.05, 0.25, 0.75, and 0.95 for each entry."
        )
      } else if (whichplot_fit() == "SignHeat") {
        p(
          strong("SignHeat help"), HTML("<br>"), HTML("<br>"),
          "Heatmap showing the entries of each signature."
        )
      } else if (whichplot_fit() == "ExposureBoxplot") {
        p(
          strong("ExposureBoxplot help"), HTML("<br>"), HTML("<br>"),
          "The levels of exposure to each signature in all genome samples."
        )
      } else if (whichplot_fit() == "ExposureBarplot") {
        p(
          strong("ExposureBarplot help"), HTML("<br>"), HTML("<br>"),
          "Barplot showing the contributions of the signatures to genome
          samples mutation counts."
        )
      } else if (whichplot_fit() == "ExposureBarplotRelative") {
        p(
          strong("ExposureBarplot relative help"), HTML("<br>"), HTML("<br>"),
          "Barplot showing the relative contribution of signatures on each
          genome samples mutation counts."
        )
      } else if (whichplot_fit() == "ExposureHeat") {
        p(
          strong("ExposureHeat help"), HTML("<br>"), HTML("<br>"),
          "Heatmap showing the exposures for each genome sample.
          Samples are grouped according to their levels of exposure to the
          signatures, as can be seen in the dendrogram on the left."
        )
      }
    }
  })

  signatures_fitting <- eventReactive(input$startfitting, {
    req(mut_fit())
    mutation <- mut_fit()
    req(kp())
    knownsigs <- kp()
    # print(head(mutation))
    opportunity <- opp_fit()
    if (is.null(opportunity)) {
      opportunity <- NA
    }
    # print(opportunity)
    sigs_fit <- NULL
    # print(input$warm_fit)
    # print(input$eval_fit)
    # print(input$em_iter_fit)
    withProgress(
      message = "Running signeR fitting...",
      detail = "This operation may take a while...",
      value = 0,
      {
        sigs_fit <- signeR(
          M = mutation, Opport = opportunity,
          P = as.matrix(knownsigs), fixedP = TRUE, main_burn = input$warm_fit,
          main_eval = input$eval_fit, EM_eval = 50, EMit_lim = input$em_iter_fit
        )
      }
    )
  })

  observeEvent(input$startfitting, {
    signatures_fitting()
  })

  output$bic_fitting <- renderPlot({
    input$startfitting
    sigs_fit <- signatures_fitting()
    if (!is.null(sigs_fit)) {
      SignPlot(sigs_fit$SignExposures)
    }
  })

  output$choosen_plot_fit <- renderPlot({
    req(whichplot_fit())
    req(input$startfitting)
    if (input$startfitting) {
      input$startfitting
      sigs_fit <- signatures_fitting()
    }
    if (is.null(sigs_fit)) {
      return(NULL)
    }
    if (!is.null(whichplot_fit())) {
      # if (whichplot_fit() == "SignPlot") {
      #   withProgress(
      #     message = "Generating the signature barplot...",
      #     detail = "This operation may take a while...",
      #     value = 0,
      #     {
      #       SignPlot(sigs_fit$SignExposures)
      #     }
      #   )
      # } else 
      if (whichplot_fit() == "SignHeat") {
        withProgress(
          message = "Generating the signature heatmap...",
          detail = "This operation may take a while...",
          value = 0,
          {
            SignHeat(sigs_fit$SignExposures)
          }
        )
      } else if (whichplot_fit() == "ExposureBoxplot") {
        withProgress(
          message = "Generating the exposure boxplot...",
          detail = "This operation may take a while...",
          value = 0,
          {
            ExposureBoxplot(sigs_fit$SignExposures)
          }
        )
      } else if (whichplot_fit() == "ExposureBarplot") {
        withProgress(
          message = "Generating the exposure barplot...",
          detail = "This operation may take a while...",
          value = 0,
          {
            ExposureBarplot(sigs_fit$SignExposures)
          }
        )
      } else if (whichplot_fit() == "ExposureBarplotRelative") {
        withProgress(
          message = "Generating the relative exposure barplot...",
          detail = "This operation may take a while...",
          value = 0,
          {
            ExposureBarplot(sigs_fit$SignExposures, , relative = TRUE)
          }
        )
      } else if (whichplot_fit() == "ExposureHeat") {
        withProgress(
          message = "Generating the exposure heatmap...",
          detail = "This operation may take a while...",
          value = 0,
          {
            ExposureHeat(sigs_fit$SignExposures)
          }
        )
      }
      return(NULL)
    }
  })

  observeEvent(input$snvhelp, {
    showModal(modalDialog(
      title = "SNV matrix help",
      includeMarkdown(
        system.file("extdata", "snv_help.md", package = "signeR")
      ),
      size = "l", easyClose = TRUE
    ))
  })

  observeEvent(input$opphelp, {
    showModal(modalDialog(
      title = "Opportunity matrix help",
      includeMarkdown(
        system.file("extdata", "opp_help.md", package = "signeR")
      ),
      size = "l", easyClose = TRUE
    ))
  })

  observeEvent(input$knownsighelp, {
    showModal(modalDialog(
      title = "Previous signatures matrix help",
      includeMarkdown(
        system.file("extdata", "knownsigs_help.md", package = "signeR")
      ),
      size = "l", easyClose = TRUE
    ))
  })

  output$dwfitting <- renderUI({
    req(input$startfitting)
    if (input$startfitting) {
      input$startfitting
      sigs_fit <- signatures_fitting()
    }
    if (is.null(sigs_fit)) {
      return(NULL)
    }
    downloadButton(ns("btdwfitting"), "Download Rdata")
  })

  output$btdwfitting <- downloadHandler(
    filename = function() {
      paste("signeRFlow-fitting_", Sys.Date(), ".RData", sep = "")
    },
    content = function(file) {
      sigs_fit <- signatures_fitting()
      save("sigs_fit", file = file)
    }
  )

}
