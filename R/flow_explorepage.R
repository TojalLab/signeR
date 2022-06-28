explorepage <- function() {
  dashboardPage(
    dashboardHeader(disable = TRUE),
    dashboardSidebar(
      sidebarMenu(
        id = "explorertabs",
        menuItem("Home",
          tabName = "dashboard",
          icon = icon("home")
        ),
        menuItem("signeR Analysis",
          icon = icon("chart-column", verify_fa = FALSE), startExpanded = TRUE,
          menuSubItem(
            span(em("de novo")),
            tabName = "de_novo_tab",
            icon = icon("cog")
          ),
          menuSubItem(
            "Fitting",
            tabName = "fitting_tab",
            icon = icon("cog")
          )
        ),
        menuItem("TCGA Explorer",
          icon = icon("database"), startExpanded = FALSE,
          tabName = "tcga_explorer_tab"
        ),
        span(
          h5(
            strong("TCGA Explorer Settings")
          ),
          style = "text-align: center"
        ),
        wellPanel(
          fluidRow(
            column(
              width = 12,
              p(
                "Show results according to analysis type",
                strong("(Analysis Type)"),
                "and by TCGA study",
                strong("(TCGA Study)"), "."
              )
            )
          ),
          selectInput(
            "analysis_type",
            "Analysis Type",
            choices = c("de novo", "fitting"),
            selected = "de novo"
          ),
          selectInput(
            "tcga_tumor",
            "TCGA Study",
            choices = tcga_tumors %>% arrange(projectID) %>% .$projectID,
          ),
          p(
            "This is a ", strong("global setting"), " used in TCGA Explorer
            module and can be changed at any time to update results."
          )
        )
      )
    ),
    dashboardBody(
      tabItems(
        tabItem(
          tabName = "dashboard",
          titleBox("signeRFlow"),
          textBox(
            width = 12,
            p(
              "Select a module to explore exposure data through
              interactive visualizations. You can analyze you own mutation
              data or explore the TCGA dataset."
            )
          ),
          sectionBox(
            title = "What's Inside",
            fluidRow(
              infoBox("Analysis Modules:", 3,
                width = 3, color = "olive",
                fill = TRUE, icon = icon("chart-simple", verify_fa = FALSE)
              ),
              infoBox("Analysis Plots:", 20,
                width = 3, color = "olive",
                fill = TRUE, icon = icon("square-poll-vertical", verify_fa = FALSE)
              ),
              infoBox("TCGA Cancers:", 33,
                width = 3, color = "teal",
                fill = TRUE, icon = icon("flask-vial", verify_fa = FALSE)
              ),
              infoBox("TCGA Samples:", "11,080",
                width = 3, color = "teal",
                fill = TRUE, icon = icon("users", verify_fa = FALSE)
              )
            )
          ),
          sectionBox(
            title = "Analysis Modules",
            messageBox(
              width = 12,
              p(
                "In case you want to analyze you own data, please click on
                signeR Analysis. Then you can infer signatures from your data (",
                span(em("de novo")), 
                ") or fit your data to known signatures (fitting).", br(),
                "Otherwise, click on TCGA Explorer to inquire exposure data
                from the public dataset."
              )
            ),
            fluidRow(
              imgLinkBox(
                width = 6,
                linkId = "link_to_signer_denovo",
                title = p("signeR ", em("de novo")),
                imgSrc = "assets/signeranalysis.jpg",
                boxText = "
                  This module provides access to signeR de novo
                  analysis to find signatures in your data,
                  estimating both signatures and related exposures.",
                linkText = "Open Module"
              ),
              imgLinkBox(
                width = 6,
                linkId = "link_to_signer_fitting",
                title = "signeR fitting",
                imgSrc = "assets/signeranalysis.jpg",
                boxText = "
                  This module provides access to signeR fitting
                  analysis to find exposures to known signatures in your data,
                  which can be uploaded or chosen from Cosmic database.
                  Exposures are estimated and can be explored.",
                linkText = "Open Module"
              )
            ),
            fluidRow(
              imgLinkBox(
                width = 6,
                linkId = "link_to_tcga_explorer",
                title = "TCGA Explorer",
                imgSrc = "assets/tcgaexplorer.jpg",
                boxText = "
                  This module provides access to the results of
                  signeR applications to TCGA datasets (33 cancer types).",
                linkText = "Open Module"
              )
            ),
          )
        ),
        tabItem(
          tabName = "de_novo_tab",
          denovo_UI("denovomod")
        ),
        tabItem(
          tabName = "fitting_tab",
          fitting_UI("fittingmod")
        ),
        tabItem(
          tabName = "tcga_explorer_tab",
          tcgaexplorer_UI("tcgaexplorer")
        )
      )
    )
  )
}
