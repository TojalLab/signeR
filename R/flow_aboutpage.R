aboutpage <- function() {
  fluidPage(
    br(),
    titleBox("About"),
    fluidRow(
      column(
        width = 12,
        column(
          width = 12,
          includeMarkdown(
            system.file("www", "about.md", package = "signeR")
          )
        )
      )
    )
  )
}
