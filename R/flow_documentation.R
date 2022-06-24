docspage <- function() {
  fluidPage(
    br(),
    titleBox("Documentation"),
    messageBox(
      width = 12,
      p(
        "The current version of the", strong("signeRFlow"),
        "was built in R using code hosted at",
        a("LINK TO REPO",href = "LINK TO REPO/"),
        ". The details below are copied directly from the",
        strong("README"), "for the code repository."
      )
    ),
    fluidRow(
      column(
        width = 12,
        column(
          width = 12,
          includeMarkdown(
            system.file("extdata", "docs.md", package = "signeR")
          )
        )
      )
    )
  )
}
