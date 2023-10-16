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
              " below. The counts of mutations are required, as a VCF, MAF or a 
              counts matrix file. The counts of mutations should be
              organized in a matrix with 96 columns (corresponding to mutations
              types) and one line for each sample. 
              Opptionally, a matrix with matching opportunities can be
              uploaded, build with a BED file or you can use a already built genome opportunity matrix (hg19 or hg38 only)
              (",
              a(
                "see signeR documentation for details",
                href = "https://bioconductor.org/packages/release/bioc/vignettes/signeR/inst/doc/signeR-vignette.html"
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
              "Upload a VCF, MAF or a SNV matrix file with your own samples
              to use in signeR fitting module and previous known signatures
              (mandatories files).
              You can upload an opportunity file as well or use a already built genome opportunity (hg19 or hg38 only).
              Also, you can upload a BED file to build an opportunity matrix."
            ),
            fluidRow(
              box(
                width = 4, background = "orange",
                tags$head(tags$script(src = "message-handler.js")),
                actionButton(ns("snvhelp"),
                  "SNV matrix help",
                  icon = icon("info-circle")
                ),
                actionButton(ns("genomehelp"),
                  "Genome installation help",
                  icon = icon("info-circle")
                ),
                hr(),
                uiOutput(ns("genomes_fit")),
                fileInput(ns("mutfile_fit"),
                  "VCF, MAF or SNV matrix*",
                  multiple = FALSE,
                  accept = c(
                    ".vcf",".vcf.gz","text/csv","text/plain",
                    "text/comma-separated-values",".csv",".maf",".maf.gz"
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
                uiOutput(ns("uigenopp_fit")),
                fileInput(ns("oppfile_fit"),
                  "Opportunities or Target file (BED)",
                  multiple = FALSE,
                  accept = c(
                    "text/csv", "text/comma-separated-values",
                    "text/plain", ".csv", ".bed"
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

  output$genomes_fit <- renderUI({

    genomes_available <- installed.genomes()

    if(length(genomes_available)==0){
      messageBox(
        title = "Warning:",
        solidHeader = TRUE,
        width = 12,
        paste0(
          "There is no genome installed. ",
          "If you need, install a genome using BSGenome (see help above)."
        )
      )
    }else{
      pickerInput(
        inputId = ns("genbuild_fit"),
        label = "Genome: ",
        choices = genomes_available,
        multiple = FALSE
      )
    }
  })

  mut_fit <- reactive({
    if (is.null(input$mutfile_fit)) {
      showModal(modalDialog(
        title = "Oh no!",
        paste0("You have not uploaded a file, silly person!"),
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }
    ext <- tools::file_ext(input$mutfile_fit$datapath)
    if (ext == "vcf" || ext == "vcf.gz") {
      req(input$genbuild_fit)

      mygenome <- getBSgenome(input$genbuild_fit)
      build <- unique(as.data.frame(GenomeInfoDb::seqinfo(mygenome))$genome)

      vcfobj <- VariantAnnotation::readVcf(input$mutfile_fit$datapath, build)

      df <- genCountMatrixFromVcf(mygenome, vcfobj)
      
    } else if (ext == "maf" || ext == "maf.gz") {
      req(input$genbuild_fit)

      mygenome <- getBSgenome(input$genbuild_fit)
      maf <- read_tsv(input$mutfile_fit$datapath)

      if (!validate_maf(maf)) {
        showModal(modalDialog(
          title = "Oh no!",
          paste0("You must upload a valid MAF file."),
          easyClose = TRUE,
          footer = NULL
        ))
        return(NULL)
      }

      df <- genCountMatrixFromMAF(mygenome, input$mutfile_fit$datapath)   

    } else {
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
    }

    return(df)
  })

  kp <- reactive({
    if (is.null(input$kP)) {
      showModal(modalDialog(
        title = "Oh no!",
        paste0("You have not uploaded a file, silly person!"),
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }
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

    if (input$genopp_fit == "yes" && is.null(input$oppfile_fit)) {

      mutation <- mut_fit()

      mygenome <- getBSgenome(input$genbuild_fit)
      build <- unique(as.data.frame(GenomeInfoDb::seqinfo(mygenome))$genome)

      nsamples = 1
      if (!is.null(mutation)){
        nsamples = nrow(mutation)
      }

      withProgress(
        message = "Download genome opportunity...",
        detail = "This operation may take a while...",
        value = 0,
        {
            data <- download_opp_file(build)
        }
      )

      opp <- as.matrix(read.table(data))
      opp <- opp[rep(1:nrow(opp), times=nsamples),]
      rownames(opp) <- rep(1:nrow(opp))

      return(opp)
    } else if(is.null(input$oppfile_fit)) {

      return(NULL)
    } else {

      if(input$genopp_fit == "yes") {
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

      ext <- tools::file_ext(input$oppfile_fit$datapath)
      if (ext == "bed") {
        
        mutation <- mut_fit()

        mygenome <- getBSgenome(input$genbuild_fit)
        build <- unique(as.data.frame(GenomeInfoDb::seqinfo(mygenome))$genome)

        target_regions <- tryCatch(
          {
            rtracklayer::import(
              con=input$oppfile_fit$datapath, format="bed", genome=build
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
          mygenome, target_regions, nsamples=nsamples
        )

        return(opp)
      } else {
        opp <- read.table(input$oppfile_fit$datapath)

        return(opp)
      }
    }
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
          "The levels of exposure to each signature in all samples."
        )
      } else if (whichplot_fit() == "ExposureBarplot") {
        p(
          strong("ExposureBarplot help"), HTML("<br>"), HTML("<br>"),
          "Barplot showing the contributions of the signatures to 
          samples mutation counts."
        )
      } else if (whichplot_fit() == "ExposureBarplotRelative") {
        p(
          strong("ExposureBarplot relative help"), HTML("<br>"), HTML("<br>"),
          "Barplot showing the relative contribution of signatures on each 
          samples mutation counts."
        )
      } else if (whichplot_fit() == "ExposureHeat") {
        p(
          strong("ExposureHeat help"), HTML("<br>"), HTML("<br>"),
          "Heatmap showing the exposures for each sample.
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
    opportunity <- opp_fit()
    if (is.null(opportunity)) {
      opportunity <- NA
    }
    sigs_fit <- NULL

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

  observeEvent(input$genomehelp, {
    showModal(modalDialog(
      title = "Genome installation help",
      includeMarkdown(
        system.file("extdata", "genome_help.md", package = "signeR")
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

  output$uigenopp_fit <- renderUI({
    req(input$genbuild_fit)

    mygenome <- getBSgenome(input$genbuild_fit)
    build <- unique(as.data.frame(GenomeInfoDb::seqinfo(mygenome))$genome)
    
    if (build %in% c('hg19', 'hg38')){
      prettyRadioButtons(
        inputId = ns("genopp_fit"), label = paste0("Use already built genome opportunity (",build,")?"), 
        choiceNames = c("Yes", "No"),
        choiceValues = c("yes", "no"),
        inline = TRUE, status = "primary", selected = "no",
        fill=TRUE, 
      )
    }

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
