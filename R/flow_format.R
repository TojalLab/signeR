
titleBox <- function(title) {
  fluidRow(
    shinydashboard::box(width = 12, background = "yellow",
        span(strong(title),
             style = "font-size:24px")
    )
  )
}

textBox <- function(...) {
  shinydashboard::box(..., status = "success")
}

tableBox <- function(...) {
  shinydashboard::box(..., status = "warning")
}


sectionBox <- function(..., title) {
  fluidRow(
    shinydashboard::box(...,
        width = 12,
        title = title,
        solidHeader = TRUE, status = "warning", collapsible = TRUE
    )
  )
}

messageBox <- function(...) {
  shinydashboard::box(..., status = "danger", background = "green")
}

imgLinkBox <- function(..., linkId, title, imgSrc, boxText, linkText) {
  shinydashboard::box(
    ...,
    title = span(title, style = "font-size:15px"),
    solidHeader = TRUE, status = "primary",
    fluidRow(
      column(
        width = 4,
        shiny::img(src = imgSrc, width = "80%")
      ),
      column(
        width = 8,
        p(boxText),
        actionButton(inputId = linkId, label = linkText)
      )
    )
  )
}

optionsBox <- function(...) {
  shinydashboard::box(..., background = "navy")
}

tableBox <- function(...) {
  shinydashboard::box(..., status = "warning")
}

plotBox <- function(...) {
  shinydashboard::box(..., status = "warning")
}
