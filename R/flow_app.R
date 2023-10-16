signeRFlow <- function() {
    if (!requireNamespace("shiny", quietly = TRUE)) {
        stop("signeRFlow requires 'shiny'. Please install it using
         install.packages('shiny')")
    }
  
    options(shiny.usecairo = F)
    # adjust limits for file update to 50MB
    options(shiny.maxRequestSize = 50 * 1024^2)

    headerTagList <- list(
        tags$style(
            type = "text/css",
            ".navbar .navbar-nav {float: right; font-size: 14px}
            .navbar .navbar-nav li a {font-size: 14px}
            .nav-tabs {font-size: 12px}"
        ),
        tags$style(
            "@import url(
                https://use.fontawesome.com/releases/v6.1.1/css/all.css);"
        ),
        tags$base(target = "_blank"),
        tags$script(
            '
            var dimension = [0, 0];
            $(document).on("shiny:connected", function(e) {
            dimension[0] = window.innerWidth;
            dimension[1] = window.innerHeight;
            Shiny.onInputChange("dimension", dimension);
            });
            $(window).resize(function(e) {
            dimension[0] = window.innerWidth;
            dimension[1] = window.innerHeight;
            Shiny.onInputChange("dimension", dimension);
            });
            '
        ),
        HTML(
            "
            <script>
            (function(i,s,o,g,r,a,m) {
            i['GoogleAnalyticsObject']=r;i[r]=i[r]||
            function() {
            (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();
            a=s.createElement(o), m=s.getElementsByTagName(o)[0];
            a.async=1;
            a.src=g;m.parentNode.insertBefore(a,m)
            })
            (window, document, 'script',
            '//www.google-analytics.com/analytics.js','ga');

            ga('create', '', 'auto');
            ga('send', 'pageview');

            </script>
            "
        ),
        includeCSS(
            system.file(
                "www", "custom.css",
                package = "signeR"
            )
        ),
        includeCSS(
            system.file(
                "www", "footer.css",
                package = "signeR"
            )
        ),
        includeCSS(
            system.file(
                "www", "bootstrapTable.min.css",
                package = "signeR"
            )
        )
    )

    footerTagList <- list(
        tags$footer(
            id = "acFooter",
            shiny::includeHTML(
                system.file(
                    "www", "footer.html",
                    package = "signeR"
                )
            )
        )
    )

    ui <- navbarPage(
        title = strong("signeRFlow"), selected = "Explore",
        tabPanel(
            "Explore", explorepage(),
            icon = icon("chart-simple", verify_fa = FALSE)
        ),
        tabPanel("About", aboutpage(), icon = icon("info-circle")),
        header = headerTagList,
        footer = footerTagList,
        collapsible = TRUE, inverse = TRUE,
        windowTitle = "signeRFlow",
    )

    server <- function(input, output, session) {
        ns <- session$ns

        callModule(
            tcgaexplorer,
            "tcgaexplorer",
            reactive(input$analysis_type),
            reactive(input$tcga_tumor),
            reactive(loadClinical()),
            reactive(loadSig())
        )

        callModule(
            denovo,
            "denovomod",
            reactive(width())
        )

        callModule(
            fitting,
            "fittingmod",
            reactive(width())
        )

        loadClinical <- function() {
            req(input$tcga_tumor)
            req(input$analysis_type)

            withProgress(
                message = "Download genome opportunity...",
                detail = "This operation may take a while...",
                value = 0,
                {
                    data <- download_clinical_file(input$tcga_tumor)
                }
            )
            clinical <- read_tsv(data)
            return(clinical)
        }

        loadSig <- function() {
            req(input$tcga_tumor)
            req(input$analysis_type)
            if (input$analysis_type == "de novo") {
                withProgress(
                    message = "Download TCGA dataset...",
                    detail = "This operation may take a while...",
                    value = 0,
                    {
                        data <- download_data_file("denovo", input$tcga_tumor)
                    }
                )
                print(data)
                load(data)
                sigs <- sig
                return(sigs)
            } else if (input$analysis_type == "fitting") {
                withProgress(
                    message = "Download TCGA dataset...",
                    detail = "This operation may take a while...",
                    value = 0,
                    {
                        data <- download_data_file("fitting", input$tcga_tumor)
                    }
                )
                print(data)
                load(data)
                sigs <- sig_test
                return(sigs)
            }
        }

        observeEvent(input$link_to_signer_denovo, {
            shinydashboard::updateTabItems(
                session, "explorertabs", "de_novo_tab"
            )
        })

        observeEvent(input$link_to_signer_fitting, {
            shinydashboard::updateTabItems(
                session, "explorertabs", "fitting_tab"
            )
        })

        observeEvent(input$link_to_tcga_explorer, {
            shinydashboard::updateTabItems(
                session, "explorertabs", "tcga_explorer_tab"
            )
        })
    }

    shinyApp(ui = ui, server = server)
}
