tcgaexplorer_UI <- function(id) {
  ns <- NS(id)

  tagList(
    tabBox(
      id = "tcgaexplorer", width = 12,
      tabPanel(
        "Home",
        icon = icon("home", verify_fa=FALSE),
        fluidRow(
          box(
            width = 12, status = "success",
            messageBox(
              width = 12,
              includeMarkdown(
                system.file(
                  "extdata", "tcga_analysis.md", package = "signeR"
                )
              )
            )
          )
        ),
        fluidRow(
          box(
            width = 6,
            withSpinner(
              plotOutput(
                ns("signplot"),
                width = "auto",
                height = "400px"
              ),
              color = "#0dc5c1"
            )
          ),
          box(
            width = 6,
            withSpinner(
              plotOutput(
                ns("signheatplot"),
                width = "auto",
                height = "400px"
              ),
              color = "#0dc5c1"
            )
          )
        ),
        fluidRow(
          box(
            title = "Data summary", width = 12, solidHeader = T,
            collapsible = F, status = "primary",
            column(
              width = 9,
              DT::dataTableOutput(ns("tcga_clinical_data")),
            ),
            column(
              width = 3,
              uiOutput(ns("tcga_feature_class_table"))
            )
          )
        ),
        fluidRow(
          box(
            width = 12,
            collapsible = F, status = "primary",
            fluidRow(
              optionsBox(
                width = 3,
                column(
                  width = 12,
                  uiOutput(ns("feature_op")),
                  hr(),
                  selectInput(
                    inputId = ns("which_plot"), label = "Available Plots:",
                    choices = c(
                      "ExposureHeat", "ExposureBarplot", "ExposureBoxplot"
                    ),
                    selected = "ExposureHeat", multiple = FALSE,
                    size = 3, selectize = FALSE
                  )
                )
              ),
              plotBox(
                width = 9,
                column(
                  width = 12,
                  withSpinner(
                    plotOutput(
                      ns("plot_chosen"),
                      width = "auto",
                      height = "400px"
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
      tabPanel(
        "Covariate",
        icon = icon("question", verify_fa=FALSE),
        fluidRow(
          box(
            title = "Filtered data summary", width = 12, solidHeader = T,
            collapsible = F, status = "primary",
            messageBox(
              width = 12,
              uiOutput(ns("filtered_value"))
            ),
            p("
              Please choose a feature of interest from the ones available for
              TCGA samples."
            ),
            column(
              width = 9,
              DT::dataTableOutput(ns("filtered_tcga_clinical_data")),
            ),
            column(
              width = 3,
              uiOutput(ns("filtered_tcga_feature_class_table"))
            )
          )
        ),
        fluidRow(
          box(
            title = p("Plots"),
            width = 12, solidHeader = T,
            collapsible = F, status = "info",
            fluidRow(
              box(
                width = 12,
                p(
                  "After you select a feature, you will be able to
                  select one button below to show a plot."
                )
              )
            ),
            column(
              width = 2,
              fluidRow(
                box(
                  width = 12,
                  uiOutput(ns("radio_buttons"))
                )
              )
            ),
            column(
              width = 10,
              fluidRow(
                box(
                  width = 12, solidHeader = T,
                  uiOutput(ns("plot_options")),
                  withSpinner(
                    plotOutput(ns("covariate_plot"), height = "450px"),
                    color = "#0dc5c1"
                  )
                )
              )
            )
          )
        ),
      )
    )
  )
}

tcgaexplorer <- function(input,
                         output,
                         session,
                         analysis_type,
                         tcga_tumor,
                         signatures) {
  ns <- session$ns

  get_similarities_tcga <- function() {
    data <- tcga_similarities %>%
      filter(
        project == tcga_tumor()
      ) %>%
      select(-project) %>%
      column_to_rownames("sigs")
    return(data)
  }

  get_tcga_data <- function() {
    feature <- input$feature
    filter_op <- input$filter_op
    filter_value <- input$filter_value
    feature_row <- input$tcga_clinical_data_rows_selected

    data <- tcga_clinical %>% filter(
      project == tcga_tumor()
    )

    if (!is.null(feature) & !is.null(feature_row)) {
      print(paste0("filter data pela feature ", feature))

      col <- names(data[feature_row + 1])
      feature <- as.list(feature)

      data <- tcga_clinical %>%
        filter(
          project == tcga_tumor(),
          !!as.symbol(col) %in% feature
        )

      return(data)
    } else if (!is.null(filter_value) & !is.null(feature_row)) {
      col <- names(data[feature_row + 1])
      print(
        paste0(
          "filter data pela feature ", col, " ", filter_op, " ", filter_value
        )
      )
      data <- tcga_clinical %>%
        filter(
          project == tcga_tumor(),
          get(filter_op)(!!as.symbol(col), filter_value)
        )
    }

    return(data)
  }

  sigs_obj <- reactive({
    req(signatures())
    sigs <- signatures()
    feature <- input$feature
    filter_op <- input$filter_op
    filter_value <- input$filter_value
    feature_row <- input$tcga_clinical_data_rows_selected
    if (!is.null(feature) & !is.null(feature_row)) {
      print(paste0("filter sig pela feature ", feature))
      data <- tcga_clinical %>% filter(
        project == tcga_tumor()
      )

      col <- names(data[feature_row + 1])
      feature <- as.list(feature)

      data <- tcga_clinical %>%
        filter(
          project == tcga_tumor(),
          !!as.symbol(col) %in% feature
        )

      sigs$SignExposures <- Reorder_samples(
        sigs$SignExposures,
        ord = which(
          sigs$SignExposures@samples %in% data$barcode
        )
      )
    } else if (!is.null(filter_value) & !is.null(feature_row)) {
      data <- get_tcga_data()
      col <- names(data[feature_row + 1])
      print(
        paste0(
          "filter sig pela feature ", col, " ", filter_op, " ", filter_value
        )
      )

      sigs$SignExposures <- Reorder_samples(
        sigs$SignExposures,
        ord = which(
          sigs$SignExposures@samples %in% data$barcode
        )
      )
    }
    return(sigs)
  })

  callModule(
    clustering,
    "clusteringmod",
    reactive(sigs_obj())
  )

  tcga_clinical_data <- reactive({
    req(tcga_tumor())
    result <- tcga_clinical %>% filter(
      project == tcga_tumor()
    )
    if (is.data.frame(result)) {
      ff <- rownames(t(result))[-1]
      df <- data.frame()

      for (f in ff) {
        fn <- f
        s <- result %>%
          select(1, f) %>%
          gather("key", "value", 2) %>%
          filter(!is.na(value)) %>%
          nrow()
        fq <- round(s / length(unique(result[[1]])) * 100, 3)
        s_na <- result %>%
          select(1, f) %>%
          gather("key", "value", 2) %>%
          filter(is.na(value)) %>%
          nrow()
        fq_na <- round(s_na / length(unique(result[[1]])) * 100, 3)
        class <- ifelse(
          class(result[[f]]) == "character", "categoric", "numeric"
        )
        data <- data.frame(
          "feature" = fn, "class" = class, "count" = paste0(s, " (", fq, "%)"),
          "missing" = paste0(s_na, " (", fq_na, "%)")
        )
        df <- rbind(df, data)
      }

      return(df)
    } else {
      return(NA)
    }
  })

  filtered_tcga_clinical_data <- reactive({
    req(tcga_tumor())
    result <- get_tcga_data()
    if (is.data.frame(result)) {
      ff <- rownames(t(result))[-1]
      df <- data.frame()

      for (f in ff) {
        fn <- f
        s <- result %>%
          select(1, f) %>%
          gather("key", "value", 2) %>%
          filter(!is.na(value)) %>%
          nrow()
        fq <- round(s / length(unique(result[[1]])) * 100, 3)
        s_na <- result %>%
          select(1, f) %>%
          gather("key", "value", 2) %>%
          filter(is.na(value)) %>%
          nrow()
        fq_na <- round(s_na / length(unique(result[[1]])) * 100, 3)
        class <- ifelse(
          class(result[[f]]) == "character", "categoric", "numeric"
        )
        data <- data.frame(
          "feature" = fn, "class" = class, "count" = paste0(s, " (", fq, "%)"),
          "missing" = paste0(s_na, " (", fq_na, "%)")
        )
        df <- rbind(df, data)
      }

      return(df)
    } else {
      return(NA)
    }
  })

  feature_class <- reactive({
    data <- tcga_clinical %>% filter(
      project == tcga_tumor()
    )
    feature_row <- input$tcga_clinical_data_rows_selected

    if (!is.null(feature_row)) {
      col <- names(data[feature_row + 1])
      t <- data %>%
        select(col) %>%
        rownames_to_column() %>%
        select(-rowname) %>%
        with(class(get(col)))

      if (t == "character") {
        fqq <- data %>%
          select(col) %>%
          arrange(col) %>%
          with(unique(.))

        sss <- data %>%
          select(col) %>%
          group_by(fnn = get(col)) %>%
          filter(!is.na(fnn)) %>%
          summarise(n = n()) %>%
          mutate(freq = (n / sum(n)) * 100)

        df <- data.frame("groups" = NA, "n" = NA, "frequency" = NA)

        df %>%
          add_row(
            groups = sss$fnn,
            n = sss$n,
            frequency = paste0(round(sss$freq, 3), "%")
          ) %>%
          filter(!is.na(groups))
      } else if (t == "numeric") {
        df <- data %>%
          select(col) %>%
          summarise(
            min = min(get(col), na.rm = T),
            max = max(get(col), na.rm = T),
            mean = mean(get(col), na.rm = T),
            sd = sd(get(col), na.rm = T)
          )
      }
    } else {
      return(NULL)
    }
  })

  filtered_feature_class <- reactive({
    data <- get_tcga_data()
    feature_row <- input$filtered_tcga_clinical_data_rows_selected

    if (!is.null(feature_row)) {
      col <- names(data[feature_row + 1])
      t <- data %>%
        select(col) %>%
        rownames_to_column() %>%
        select(-rowname) %>%
        with(class(get(col)))

      if (t == "character") {
        fqq <- data %>%
          select(col) %>%
          arrange(col) %>%
          with(unique(.))

        sss <- data %>%
          select(col) %>%
          group_by(fnn = get(col)) %>%
          filter(!is.na(fnn)) %>%
          summarise(n = n()) %>%
          mutate(freq = (n / sum(n)) * 100)

        df <- data.frame("groups" = NA, "n" = NA, "frequency" = NA)

        df %>%
          add_row(
            groups = sss$fnn,
            n = sss$n,
            frequency = paste0(round(sss$freq, 3), "%")
          ) %>%
          filter(!is.na(groups))
      } else if (t == "numeric") {
        df <- data %>%
          select(col) %>%
          summarise(
            min = min(get(col), na.rm = T),
            max = max(get(col), na.rm = T),
            mean = mean(get(col), na.rm = T),
            sd = sd(get(col), na.rm = T)
          )
      }
    } else {
      return(NULL)
    }
  })

  output$tcga_clinical_data_selected <- renderPrint({
    input$tcga_clinical_data_rows_selected
  })

  output$filtered_tcga_clinical_data_selected <- renderPrint({
    input$filtered_tcga_clinical_data_rows_selected
  })

  output$tcga_feature_class_table <- renderTable({
    feature_class()
  })

  output$filtered_tcga_feature_class_table <- renderTable({
    filtered_feature_class()
  })

  observeEvent(input$tcga_clinical_data_rows_selected, {
    feature_class()
  })

  observeEvent(input$filtered_tcga_clinical_data_rows_selected, {
    filtered_feature_class()
    get_plots_choices()
  })

  output$tcga_feature_class_table <- renderTable({
    feature_class()
  })

  output$filtered_tcga_feature_class_table <- renderTable({
    filtered_feature_class()
  })

  output$tcga_clinical_data <- DT::renderDataTable(
    tcga_clinical_data(),
    server = FALSE, selection = list(mode = "single")
  )

  output$filtered_tcga_clinical_data <- DT::renderDataTable(
    filtered_tcga_clinical_data(),
    server = FALSE, selection = list(mode = "single")
  )

  output$feature_op <- renderUI({
    data <- tcga_clinical %>% filter(
      project == tcga_tumor()
    )
    feature_row <- input$tcga_clinical_data_rows_selected
    if (!is.null(feature_row)) {
      col <- names(data[feature_row + 1])
      t <- data %>%
        select(col) %>%
        rownames_to_column() %>%
        select(-rowname) %>%
        with(class(get(col)))

      if (t == "character") {
        fqq <- data %>%
          select(col) %>%
          arrange(col) %>%
          with(unique(.))

        pickerInput(
          inputId = ns("feature"),
          label = paste0("Options from ", col),
          choices = fqq[[col]],
          multiple = TRUE
        )
      } else {
        df <- data %>%
          select(col) %>%
          summarise(
            min = min(get(col), na.rm = T),
            max = max(get(col), na.rm = T),
            mean = mean(get(col), na.rm = T),
            sd = sd(get(col), na.rm = T)
          )

        fluidRow(
          width = 12,
          column(
            width = 6,
            pickerInput(
              inputId = ns("filter_op"),
              label = "Filter type:",
              choices = c(">", ">=", "<", "<="),
              multiple = FALSE
            )
          ),
          column(
            width = 6,
            sliderInput(
              ns("filter_value"),
              "Value:", df$min, df$max, df$mean, 1.0
            )
          )
        )
      }
    } else {
      HTML("You can select a feature to filter the dataset.</br>")
    }
  })

  output$signplot <- renderPlot({
    tcga_tumor()
    sigs <- sigs_obj()
    if (!is.null(sigs)) {
      withProgress(
        message = "Generating the signature barplot...",
        detail = "This operation may take a while...",
        value = 0,
        {
          SignPlot(sigs$SignExposures)
        }
      )
    }
  })

  output$signheatplot <- renderPlot({
    tcga_tumor()
    sigs <- sigs_obj()
    similarities <- get_similarities_tcga()
    if (!is.null(sigs)) {
      pheatmap(
        t(similarities),
        display_numbers = F,
        angle_col = 45, cluster_rows = F, cluster_cols = F,
        fontsize_row = 6
      )
    }
  })

  which_plot <- reactive({
    req(input$which_plot)
    if (is.null(input$which_plot)) {
      return(NULL)
    }
    return(input$which_plot)
  })

  output$plot_chosen <- renderPlot({
    tcga_tumor()
    sigs <- sigs_obj()
    if (!is.null(which_plot())) {
      if (which_plot() == "ExposureBoxplot") {
        withProgress(
          message = "Generating the exposure boxplot...",
          detail = "This operation may take a while...",
          value = 0,
          {
            ExposureBoxplot(sigs$SignExposures)
          }
        )
      } else if (which_plot() == "ExposureBarplot") {
        withProgress(
          message = "Generating the exposure barplot...",
          detail = "This operation may take a while...",
          value = 0,
          {
            ExposureBarplot(sigs$SignExposures)
          }
        )
      } else if (which_plot() == "ExposureHeat") {
        withProgress(
          message = "Generating the exposure heatmap...",
          detail = "This operation may take a while...",
          value = 0,
          {
            ExposureHeat(sigs$SignExposures)
          }
        )
      }
      return(NULL)
    }
  })

  diffexp_method <- reactive({
    req(input$diffexp_method)
    return(input$diffexp_method)
  })

  diffexp_quant <- reactive({
    req(input$diffexp_quant)
    return(input$diffexp_quant)
  })

  diffexp_cutoff <- reactive({
    req(input$diffexp_cutoff)
    return(input$diffexp_cutoff)
  })

  diffexp_padj <- reactive({
    req(input$diffexp_padj)
    return(input$diffexp_padj)
  })

  sclassif_method <- reactive({
    req(input$sclassif_method)
    return(input$sclassif_method)
  })

  survival_method <- reactive({
    req(input$survival_method)
    return(input$survival_method)
  })


  diffexpplot <- function() {
    output$covariate_plot <- renderPlot({
      data <- get_tcga_data()

      feature_row <- input$filtered_tcga_clinical_data_rows_selected
      if (!is.null(feature_row)) {
        col <- names(data[feature_row + 1])
        labels <- data[[col]]
        if (class(labels) == "character") {
          sigs <- sigs_obj()
          if (is.null(sigs)) {
            return(NULL)
          }
          if (!is.null(sigs)) {
            difexp_method <- diffexp_method()
            diffexp.quant <- diffexp_quant()
            diffexp.cutoff <- diffexp_cutoff()
            diffexp.padj <- diffexp_padj()

            DiffExp(
              sigs$SignExposures,
              labels = labels,
              quant = diffexp.quant, cutoff = diffexp.cutoff,
              p.adj = diffexp.padj
            )
          }
        }
      }
      return(NULL)
    })
  }

  sampleclassplot <- function() {
    output$covariate_plot <- renderPlot({
      data <- get_tcga_data()
      feature_row <- input$filtered_tcga_clinical_data_rows_selected
      if (!is.null(feature_row)) {
        col <- names(data[feature_row + 1])
        labels <- data[[col]]
        sclas_method <- sclassif_method()
        if (class(labels) == "character") {
          sigs <- sigs_obj()
          if (is.null(sigs)) {
            return(NULL)
          }
          if (!is.null(sigs)) {
            ExposureClassify(
              sigs$SignExposures,
              labels = labels,
              method = sclas_method
            )
          }
        }
      }
      return(NULL)
    })
  }

  correlationplot <- function() {
    output$covariate_plot <- renderPlot({
      data <- get_tcga_data()

      feature_row <- input$filtered_tcga_clinical_data_rows_selected
      if (!is.null(feature_row)) {
        col <- names(data[feature_row + 1])
        feature <- data[[col]]
        if (class(feature) == "numeric") {
          sigs <- sigs_obj()
          if (is.null(sigs)) {
            return(NULL)
          }
          if (!is.null(sigs)) {
            ExposureCorrelation(
              sigs$SignExposures,
              feature = feature
            )
          }
        }
      }
      return(NULL)
    })
  }

  linearregressionplot <- function() {
    output$covariate_plot <- renderPlot({
      data <- get_tcga_data()

      feature_row <- input$filtered_tcga_clinical_data_rows_selected
      if (!is.null(feature_row)) {
        col <- names(data[feature_row + 1])
        feature <- data[[col]]
        if (class(feature) == "numeric") {
          sigs <- sigs_obj()
          if (is.null(sigs)) {
            return(NULL)
          }
          if (!is.null(sigs)) {
            ExposureGLM(
              sigs$SignExposures,
              feature = feature
            )
          }
        }
      }
      return(NULL)
    })
  }

  survivalplot <- function() {
    output$covariate_plot <- renderPlot({
      data <- get_tcga_data()

      surv_method <- survival_method()
      if ("time" %in% names(data) && "status" %in% names(data)) {
        su <- as.matrix(data.frame(time = as.numeric(data$time), status = as.numeric(data$status)))
        sigs <- sigs_obj()
        if (is.null(sigs)) {
          return(NULL)
        }
        if (!is.null(sigs)) {
          ExposureSurvival(
            Exposures = sigs$SignExposures, surv = su,
            method = surv_method
          )
        }
      }
      return(NULL)
    })
  }

  coxplot <- function() {
    output$covariate_plot <- renderPlot({
      data <- get_tcga_data()

      if ("time" %in% names(data) && "status" %in% names(data)) {
        su <- as.matrix(data.frame(time = as.numeric(data$time), status = as.numeric(data$status)))
        sigs <- sigs_obj()
        if (is.null(sigs)) {
          return(NULL)
        }
        if (!is.null(sigs)) {
          ExposureSurvModel(
            Exposures = sigs$SignExposures, surv = su
          )
        }
      }
      return(NULL)
    })
  }

  diffexpui <- function() {
    req(input$plotid)
    output$plot_options <- renderUI({
      data <- get_tcga_data()

      dropdownButton(

        selectInput(
          inputId = ns("diffexp_method"), label = "Method:",
          choices = c("kruskal.test"),
          selected = "kruskal.test", multiple = FALSE,
        ) %>% shinyInput_label_embed(
          shiny::icon("info-circle", verify_fa=FALSE) %>%
            bs_embed_tooltip(title = "Method")
        ),
        selectInput(
          inputId = ns("diffexp_padj"), label = "P-value adjust:",
          choices = c("BH"),
          selected = "BH", multiple = FALSE,
        ) %>% shinyInput_label_embed(
          shiny::icon("info-circle", verify_fa=FALSE) %>%
            bs_embed_tooltip(title = "padj")
        ),
        numericInput(
          ns("diffexp_quant"), "P-value quantile", 0.5,
          min = 0, max = 1, step = 0.1
        ) %>% shinyInput_label_embed(
          shiny::icon("info-circle", verify_fa=FALSE) %>%
            bs_embed_tooltip(title = "quantile")
        ),
        numericInput(
          ns("diffexp_cutoff"), "P-value threshold", 0.5,
          min = 0, max = 1, step = 0.1
        ) %>% shinyInput_label_embed(
          shiny::icon("info-circle", verify_fa=FALSE) %>%
            bs_embed_tooltip(title = "threshold")
        ),
        circle = TRUE, status = "danger",
        icon = icon("gear", verify_fa=FALSE), width = "200px",
        tooltip = tooltipOptions(title = "Plot options")
      )
    })
  }

  sampleclassui <- function() {
    req(input$plotid)
    output$plot_options <- renderUI({
      data <- get_tcga_data()

      dropdownButton(

        selectInput(
          inputId = ns("sclassif_method"), label = "Method:",
          choices = c(
            "knn", "lvq", "logreg", "lda",
            "lasso", "nb", "svm", "rf", "ab"
          ),
          selected = "knn", multiple = FALSE
        ) %>% shinyInput_label_embed(
          shiny::icon("info-circle", verify_fa=FALSE) %>%
            bs_embed_tooltip(title = "Method")
        ),
        circle = TRUE, status = "danger",
        icon = icon("gear", verify_fa=FALSE), width = "200px",
        tooltip = tooltipOptions(title = "Plot options")
      )
    })
  }

  correlationui <- function() {
    req(input$plotid)
    output$plot_options <- renderUI({
      NULL
    })
  }

  linearregressionui <- function() {
    req(input$plotid)
    output$plot_options <- renderUI({
      NULL
    })
  }

  survivalui <- function() {
    req(input$plotid)
    output$plot_options <- renderUI({
      data <- get_tcga_data()

      dropdownButton(

        selectInput(
          inputId = ns("survival_method"), label = "Method:",
          choices = c("logrank", "cox"),
          selected = "logrank", multiple = FALSE,
        ) %>% shinyInput_label_embed(
          shiny::icon("info-circle", verify_fa=FALSE) %>%
            bs_embed_tooltip(title = "Method")
        ),
        circle = TRUE, status = "danger",
        icon = icon("gear", verify_fa=FALSE), width = "200px",
        tooltip = tooltipOptions(title = "Plot options")
      )
    })
  }

  coxui <- function() {
    req(input$plotid)
    output$plot_options <- renderUI({
      NULL
    })
  }

  observeEvent(input$plotid, {
    req(input$filtered_tcga_clinical_data_rows_selected)
    if (input$plotid == "de") {
      diffexpui()
      diffexpplot()
    } else if (input$plotid == "sc") {
      sampleclassui()
      sampleclassplot()
    } else if (input$plotid == "cor") {
      correlationui()
      correlationplot()
    } else if (input$plotid == "lr") {
      linearregressionui()
      linearregressionplot()
    } else if (input$plotid == "sv") {
      survivalui()
      survivalplot()
    } else if (input$plotid == "cx") {
      coxui()
      coxplot()
    }
  })

  # isso para que o spinner nao fique rodando
  output$covariate_plot <- renderPlot({
    return(NULL)
  })

  output$filtered_value <- renderUI({
    feature <- input$feature
    filter_op <- input$filter_op
    filter_value <- input$filter_value
    feature_row <- input$tcga_clinical_data_rows_selected
    data <- get_tcga_data()
    if (!is.null(feature_row) & !is.null(feature)) {
      col <- names(data[feature_row + 1])
      HTML(
        "You have filtered the <em>", tcga_tumor(),
        "</em> dataset by <strong>", col,
        "</strong> using the values: <span style='color:red'>",
        paste(feature, collapse = ","),
        "</span>."
      )
    } else if (!is.null(filter_value) & !is.null(feature_row)) {
      col <- names(data[feature_row + 1])
      HTML(
        "You have filtered the <em>", tcga_tumor(),
        "</em> dataset by <strong>", col,
        "</strong> using the values <span style='color:red'>",
        filter_op, " then ", filter_value, "</span>."
      )
    } else {
      HTML(
        "You are using the <span style='color:red'>complete</span>",
        "<strong>", tcga_tumor(), "</strong> dataset."
      )
    }
  })

  get_plots_choices <- function() {
    req(input$filtered_tcga_clinical_data_rows_selected)
    feature_row <- input$filtered_tcga_clinical_data_rows_selected
    data <- get_tcga_data()
    col <- names(data[feature_row + 1])
    t <- data %>%
      select(col) %>%
      rownames_to_column() %>%
      select(-rowname) %>%
      with(class(get(col)))
    if (col == "time" || col == "status") {
      return(c(
        `<i class='km-img'></i>` = "sv",
        `<i class='cox-img'></i>` = "cx"
      ))
    } else if (t == "character") {
      return(c(
        `<i class='bx-img'></i>` = "de",
        `<i class='sc-img'></i>` = "sc"
      ))
    } else if (t == "numeric") {
      return(c(
        `<i class='cor-img'></i>` = "cor",
        `<i class='lr-img'></i>` = "lr"
      ))
    }
  }

  output$radio_buttons <- renderUI({
    radioGroupButtons(
      inputId = ns("plotid"),
      label = NULL,
      choices = get_plots_choices(),
      selected = character(0),
      direction = "vertical"
    )
  })
}
