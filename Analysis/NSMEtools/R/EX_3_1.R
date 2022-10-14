library(shiny)

ui <- fluidPage(
  textInput("name", "What's your name?"),
  textOutput("greeting")
)
server1 <- function(input, output, server) {
  output$greeting <- renderText(paste0("Hello ", input$name))
}

server2 <- function(input, output, server) {
  greeting <- reactive(paste0("Hello ", input$name))
  output$greeting <- renderText(greeting())
}

server3 <- function(input, output, server) {
  output$greeting <- reactive(paste0("Hello ", input$name))
}

shinyApp(ui,server3)
