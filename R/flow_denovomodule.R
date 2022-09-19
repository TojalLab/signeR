denovo_UI <- function(id) {
  ns <- NS(id)

  tagList(
    tabBox(
      id = "denovo", width = 12,
      tabPanel(
        span(em("de novo")),
        icon = icon("brain"),
        fluidRow(
          box(
            width = 12,
            p(
              "Please ", strong("upload your data"),
              " below. The counts of mutations are required, as a VCF file or a 
              counts matrix. The counts of mutations should be
              organized in a matrix with 96 columns
              (corresponding to mutations types) and one line for each genome sample. 
              Opptionally, a matrix with matching opportunities can be
              uploaded, build with a BED file or you can use a already built genome opportunity matrix 
              (",
              a(
                "see signeR documentation for details",
                href = "https://bioconductor.org/packages/release/bioc/vignettes/signeR/inst/doc/signeR-vignette.html"
              ),
              ").
              signeRFlow will estimate signatures and sample's exposure data
              and display results through interactive visualizations."
            )
          )
        ),
        fluidRow(
          box(
            title = "Upload data", width = 12, solidHeader = T,
            collapsible = T, status = "primary",
            messageBox(
              width = 12,
              "Upload a VCF file or a SNV matrix file (mandatory) with your own samples
              to use in signeR de novo module.
              You can upload an opportunity file as well or use a already built genome opportunity.
              Also, you can upload a BED file to build an opportunity matrix."
            ),
            fluidRow(
              box(
                width = 6, background = "orange",
                tags$head(tags$script(src = "message-handler.js")),
                  actionButton(ns("snvhelp"),
                    "SNV matrix help",
                    icon = icon("info-circle")
                  ),
                hr(),
                prettyRadioButtons(
                  inputId = ns("genbuild"), label = "Genome build :", 
                  choiceNames = c("hg19/GRCh37", "hg38/GRCh38"),
                  choiceValues = c("hg19", "hg38"),
                  inline = TRUE, status = "primary", selected = "hg19",
                  fill=TRUE, 
                ),
                fileInput(ns("mutfile"),
                  "VCF file or SNV matrix*",
                  multiple = FALSE,
                  accept = c(
                    ".vcf",".vcf.gz","text/csv", "text/plain",
                    "text/comma-separated-values",".csv"
                  )
                )
              ),
              box(
                width = 6, background = "navy",
                tags$head(tags$script(src = "message-handler.js")),
                actionButton(ns("opphelp"),
                  " Opportunity help",
                  icon = icon("info-circle")
                ),
                hr(),
                uiOutput(ns("uigenopp")),
                fileInput(ns("oppfile"),
                  "Opportunities or Target file (BED)",
                  multiple = FALSE,
                  accept = c(
                    "text/csv", "text/comma-separated-values",
                    "text/plain", ".csv", ".bed"
                  )
                )
              )
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
                    "The algorithm computes a Bayesian approach to the
                    non-negative factorization (NMF) of the mutation counts
                    in a matrix product of mutational signatures and exposures
                    to mutational processes. If provided, opportunities are
                    used as weights for the factorization.  Further analysis
                    parameters can be set below. Results can be visualized on
                    different plots and found signatures can be compared to
                    the ones in Cosmic database. "
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
                  sliderInput(
                    ns("nsigs"), "Number of signatures (min and max):",
                    min = 1, max = 95, value = c(1, 10)
                  ),
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
                        actionLink(
                          ns("iterationhelp"),
                          label = icon("info-circle")
                        )
                      ),
                      HTML("<br>"),
                      div(
                        style = "display: inline-block;width: 30%;",
                        numericInput(ns("em_iter"), "EM", 10)
                      ),
                      div(
                        style = "display: inline-block;width: 30%;",
                        numericInput(ns("warm"), "Warm-up", 10)
                      ),
                      div(
                        style = "display: inline-block;width: 30%;",
                        numericInput(ns("eval"), "Final", 10)
                      )
                    )
                  ),
                  # hr(),
                  actionButton(
                    ns("startdenovo"),
                    label = "Start de novo analysis", icon = NULL
                  ),
                  uiOutput(ns("dwdenovo"))
                )
              )
            ),
            column(
              width = 8,
              fluidRow(
                box(
                  width = 12, solidHeader = T,
                    withSpinner(plotOutput(ns("bic_denovo")))
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
                    inputId = ns("whichplot"), label = "Available Plots:",
                    choices = c(
                      "SignPlot", "SignHeat", "ExposureBoxplot",
                      "ExposureBarplot", "ExposureBarplotRelative",
                      "ExposureHeat"
                    ),
                    selected = "SignPlot", multiple = FALSE,
                    size = 5, selectize = FALSE
                  )
                ),
                box(
                  width = 12, solidHeader = T,
                  collapsible = F, background = "olive",
                  uiOutput(ns("plot_help"))
                )
              )
            ),
            column(
              width = 8,
              fluidRow(
                box(
                  width = 12, solidHeader = T,
                  withSpinner(
                    plotOutput(
                      ns("choosen_plot"),
                      width = "auto",
                      height = "400px"
                    ),
                    color = "#0dc5c1"
                  )
                )
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "COSMIC Comparison", width = 12,
            solidHeader = T, collapsible = T, status = "info", collapsed = F,
            fluidRow(
              box(
                width = 12,
                p(
                  "The plot below ilustrates the cosine distance between
                  found signatures and Cosmic version 3.2."
                  )
              )
            ),
            fluidRow(
              box(
                width = 12, solidHeader = T,
                withSpinner(plotOutput(ns("comparison_plot"),width = "auto"),
                  color = "#0dc5c1"
                ),
                # uiOutput(ns("comparison_table"))
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

denovo <- function(input,
                   output,
                   session,
                   width) {
  ns <- session$ns

  callModule(
    clustering,
    "clusteringmod",
    reactive(signatures_denovo())
  )

  callModule(
    covariate,
    "covariatemod",
    reactive(signatures_denovo())
  )

  mut <- reactive({
    # req(input$mutfile)
    if (is.null(input$mutfile)) {
      showModal(modalDialog(
        title = "Oh no!",
        paste0("You have not uploaded a file, silly person!"),
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }
    
    #df <- read.table(input$mutfile$datapath, header=T,sep="\t",row.names=1,check.names=F)
    ext <- tools::file_ext(input$mutfile$datapath)
    if (ext == "vcf" || ext == "vcf.gz") {
      build <- input$genbuild
      if (build == "hg19"){
        if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19")) BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
        mygenome <- BSgenome.Hsapiens.UCSC.hg19
      } else {
        if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38")) BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
        mygenome <- BSgenome.Hsapiens.UCSC.hg38
      }

      vcfobj <- VariantAnnotation::readVcf(input$mutfile$datapath, build)

      df <- genCountMatrixFromVcf(mygenome, vcfobj)
      
    } else {
      df <- read.table(input$mutfile$datapath, header=T,sep="\t",row.names=1,check.names=F)
      if (!validate_cnv(df)) {
        showModal(modalDialog(
          title = "Oh no!",
          paste0("You must upload a valid SNV matrix file."),
          easyClose = TRUE,
          footer = NULL
        ))
        return(NULL)
      }
    }

    return(df)
  })

  opp <- reactive({
    # isso torno o input obrigatorio
    # req(input$oppfile)

    if (input$genopp == "yes" && is.null(input$oppfile)) {

      mutation <- mut()

      nsamples = 1
      if (!is.null(mutation)){
        nsamples = nrow(mutation)
      }

      withProgress(
        message = "Download genome opportunity...",
        detail = "This operation may take a while...",
        value = 0,
        {
            data <- download_opp_file(input$genbuild)
        }
      )

      opp <- as.matrix(read.table(data))
      opp <- opp[rep(1:nrow(opp), times=nsamples),]
      rownames(opp) <- rep(1:nrow(opp))

      return(opp)
    } else if(is.null(input$oppfile)) {

      return(NULL)
    } else {
      # read.table(input$oppfile$datapath)

      if(input$genopp == "yes") {
        showModal(modalDialog(
          title = "Opportunity conflict",
          paste0(
            "You have selected to use genome opportunity and uploaded a file.
            signeRFlow will use the uploaded file and ignore genome opportunity."
          ),
          easyClose = TRUE,
          footer = NULL
        ))
      }

      ext <- tools::file_ext(input$oppfile$datapath)
      if (ext == "bed") {
        build <- input$genbuild
        mutation <- mut()
        if (build == "hg19"){
          if (!require("BSgenome.Hsapiens.UCSC.hg19")) BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
          mygenome <- BSgenome.Hsapiens.UCSC.hg19
        } else {
          if (!require("BSgenome.Hsapiens.UCSC.hg38")) BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
          mygenome <- BSgenome.Hsapiens.UCSC.hg38
        }

        target_regions <- tryCatch(
          {
            rtracklayer::import(
              con=input$oppfile$datapath, format="bed", genome=build
            )
          },
          error=function(cond){
            print(cond)
          },
          warning=function(cond){
            print(cond)
          }
        )

        if (class(target_regions)[[1]] != "GRanges") {
          showModal(modalDialog(
            title = "BED error",
            paste0(
              "signerflow couldn't process your BED file.","\n",
              "Error message: ", target_regions
            ),
            easyClose = TRUE,
            footer = NULL
          ))
          return(NULL)
        }

        nsamples = 1
        if (!is.null(mutation)){
          nsamples = nrow(mutation)
        }

        opp <- genOpportunityFromGenome(
          mygenome,target_regions, nsamples=nsamples
        )

        return(opp)
      } else {
        opp <- read.table(input$oppfile$datapath)

        return(opp)
      }
    }
  })

  observeEvent(input$iterationhelp, {
    showModal(modalDialog(
      title = "Iterations parameters",
      includeMarkdown(
        system.file("extdata", "iterations.md", package = "signeR")
      ),
      size = "l", easyClose = TRUE
    ))
  })

  whichplot <- reactive({
    req(input$whichplot)
    if (is.null(input$whichplot)) {
      return(NULL)
    }
    return(input$whichplot)
  })

  output$plot_help <- renderUI({
    req(input$whichplot)
    if (!is.null(whichplot())) {
      if (whichplot() == "SignPlot") {
        p(
          strong("SignPlot help"), HTML("<br>"), HTML("<br>"),
          "Signatures barplot with error bars reflecting the sample percentiles
          0.05, 0.25, 0.75, and 0.95 for each entry."
        )
      } else if (whichplot() == "SignHeat") {
        p(
          strong("SignHeat help"), HTML("<br>"), HTML("<br>"),
          "Heatmap showing the entries of each signature."
        )
      } else if (whichplot() == "ExposureBoxplot") {
        p(
          strong("ExposureBoxplot help"), HTML("<br>"), HTML("<br>"),
          "The levels of exposure to each signature in all genome samples."
        )
      } else if (whichplot() == "ExposureBarplot") {
        p(
          strong("ExposureBarplot help"), HTML("<br>"), HTML("<br>"),
          "Barplot showing the contributions of the signatures to genome
          samples mutation counts."
        )
      } else if (whichplot() == "ExposureBarplotRelative") {
        p(
          strong("ExposureBarplot relative help"), HTML("<br>"), HTML("<br>"),
          "Barplot showing the relative contribution of signatures on each
          genome samples mutation counts."
        )
      } else if (whichplot() == "ExposureHeat") {
        p(
          strong("ExposureHeat help"), HTML("<br>"), HTML("<br>"),
          "Heatmap showing the exposures for each genome sample.
          Samples are grouped according to their levels of exposure to the
          signatures, as can be seen in the dendrogram on the left."
        )
      }
    }
  })

  signatures_denovo <- eventReactive(input$startdenovo, {
    req(mut())
    mutation <- mut()
    # print(head(mutation))
    opportunity <- opp()
    if (is.null(opportunity)) {
      opportunity <- NA
    }
    # print(opportunity)
    sigs <- NULL
    # print(input$nsigs)
    # print(input$warm)
    # print(input$eval)
    # print(input$em_iter)
    withProgress(
      message = "Running signeR...",
      detail = "This operation may take a while...",
      value = 0,
      {
        sigs <- signeR(
          M = mutation, Opport = opportunity,
          nlim = input$nsigs, main_burn = input$warm,
          main_eval = input$eval, EM_eval = 50, EMit_lim = input$em_iter
        )
      }
    )
  })

  observeEvent(input$startdenovo, {
    signatures_denovo()
  })

  output$bic_denovo <- renderPlot({
    input$startdenovo
    sigs <- signatures_denovo()
    if (!is.null(sigs)) {
      BICboxplot(sigs)
    }
  })

  output$choosen_plot <- renderPlot({
    req(whichplot())
    req(input$startdenovo)
    if (input$startdenovo) {
      input$startdenovo
      sigs <- signatures_denovo()
    }
    if (is.null(sigs)) {
      return(NULL)
    }
    if (!is.null(whichplot())) {
      if (whichplot() == "SignPlot") {
        withProgress(
          message = "Generating the signature barplot...",
          detail = "This operation may take a while...",
          value = 0,
          {
            SignPlot(sigs$SignExposures)
          }
        )
      } else if (whichplot() == "SignHeat") {
        withProgress(
          message = "Generating the signature heatmap...",
          detail = "This operation may take a while...",
          value = 0,
          {
            SignHeat(sigs$SignExposures)
          }
        )
      } else if (whichplot() == "ExposureBoxplot") {
        withProgress(
          message = "Generating the exposure boxplot...",
          detail = "This operation may take a while...",
          value = 0,
          {
            ExposureBoxplot(sigs$SignExposures)
          }
        )
      } else if (whichplot() == "ExposureBarplot") {
        withProgress(
          message = "Generating the exposure barplot...",
          detail = "This operation may take a while...",
          value = 0,
          {
            ExposureBarplot(sigs$SignExposures)
          }
        )
      } else if (whichplot() == "ExposureBarplotRelative") {
        withProgress(
          message = "Generating the relative exposure barplot...",
          detail = "This operation may take a while...",
          value = 0,
          {
            ExposureBarplot(sigs$SignExposures, relative = TRUE)
          }
        )
      } else if (whichplot() == "ExposureHeat") {
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

  # Cosmicsigs <- reactive({
  #   req(input$cosmic)
  #   if (is.null(input$cosmic)) {
  #     return(NULL)
  #   }
  #   Cosm <- read.table("Cosmic_SBS3.txt", sep = "\t", header = T)
  #   columns <- as.vector(input$cosmic)
  #   if ("All" %in% columns) {
  #     columns <- colnames(Cosm)
  #   } else {
  #     columns <- columns[columns %in% colnames(Cosm)]
  #   }
  #   return(Cosm[, columns])
  # })

  # output$comparison_table <- renderTable({
  #   # Depende de ter rodado análise denovo
  #   # de Cosmicsigs e de input$coscut.
  #   if (input$startdenovo) {
  #     input$startdenovo
  #     dist <- comparison()
  #     table <- data.frame(Found_signature = NA, Previous_signature = NA)
  #     for (i in 1:NROW(dist)) {
  #       for (j in 1:NCOL(dist)) {
  #         if (dist[i, j] > coscut) {
  #           table <- rbind(table, data.frame(
  #             Found_signature = rownames(dist)[i],
  #             Previous_signature = colnames(dist)[j]
  #           ))
  #         }
  #       }
  #     }
  #     return(table)
  #   }
  # })

  comparison <- function() {
      sigs <- signatures_denovo()
      # Normalizing signer signatures
      SE <- sigs[[6]]
      signer <- Median_sign(SE)
      rownames(signer) <- SE@mutations
      colnames(signer) <- paste0("S", 1:(sigs$Nsign))
      signer <- t(signer)
      signer <- signer / rowSums(signer)
      # Ordering cosmic table according to signeR order
      cosmic <- cosmic_data
      cosmic <- mutate(
        cosmic,
        t = paste0(
          # gsub(">", ".", Substitution.Type), ".", Trinucleotide
          Substitution.Type, ":", Trinucleotide
        )
      ) %>%
        dplyr::select(t, contains("SBS")) %>%
        column_to_rownames("t")
      cosmic <- cosmic[colnames(signer), ]
      # print(rownames(cosmic))
      # print(colnames(signer))
      stopifnot(rownames(cosmic) == colnames(signer))
      # Calculating similarities signer vs cosmic
      dist <- as.matrix(1 - proxy::dist(signer, t(cosmic), method = "cosine"))
      return(dist)
  }

  output$comparison_plot <- renderPlot({
    # Depende de ter rodado análise denovo
    # e de ler Cosmicsigs.
    if (input$startdenovo) {
      dist <- comparison()
      pheatmap(
        dist,
        display_numbers = F,
        main = "SigneR and COSMIC similarities",
        cluster_rows = F, cluster_cols = F
      )
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

  output$dwdenovo <- renderUI({
    req(input$startdenovo)
    if (input$startdenovo) {
      input$startdenovo
      sigs <- signatures_denovo()
    }
    if (is.null(sigs)) {
      return(NULL)
    }
    downloadButton(ns("btdwdenovo"), "Download Rdata")
  })

  output$uigenopp <- renderUI({
    req(input$genbuild)
    build = "hg19"
    if (input$genbuild == "hg38") {
      build = "hg38"
    }
    
    prettyRadioButtons(
      inputId = ns("genopp"), label = paste0("Use already built genome opportunity (",build,")?"), 
      choiceNames = c("Yes", "No"),
      choiceValues = c("yes", "no"),
      inline = TRUE, status = "primary", selected = "no",
      fill=TRUE, 
    )

  })

  output$btdwdenovo <- downloadHandler(
    filename = function() {
      paste("signeRFlow-denovo", Sys.Date(), ".RData", sep = "")
    },
    content = function(file) {
      sigs <- signatures_denovo()
      save("sigs", file = file)
    }
  )
}
