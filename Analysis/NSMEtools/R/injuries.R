library(shiny)
library(ggplot2)
library(vroom)
library(forcats)
library(dplyr)

# Load data, if needed.
download <- function(name) {
  url <- "https://github.com/hadley/mastering-shiny/blob/main/neiss/"
  download.file(paste0(url, name, "?raw=true"), paste0("neiss/", name), quiet = TRUE)
}
notLoaded=function(x) tryCatch(if (any(class(x) == 'tbl')) 0 else 0, error=function(e) 1) 

if ( notLoaded(injuries) ) {
  download("injuries.tsv.gz")
  injuries <- vroom::vroom("neiss/injuries.tsv.gz", show_col_types = FALSE)
}
if ( notLoaded(population) ) {
  download("population.tsv")
  population <- vroom::vroom("neiss/population.tsv", show_col_types = FALSE)
}
if ( notLoaded(products) ) {
  download("products.tsv")
  products <- vroom::vroom("neiss/products.tsv", show_col_types = FALSE)
}
# This is a named vector
prod_codes <- setNames(products$prod_code, products$title)

#<< count_top
count_top <- function(df, var, n = 5) {
  df %>%
    mutate({{ var }} := fct_lump(fct_infreq({{ var }}), n = n)) %>%
    group_by({{ var }}) %>%
    summarise(n = as.integer(sum(weight)))
}
#>>

ui <- fluidPage(
  fluidRow(
    column(6,selectInput("code", "Product", choices = prod_codes)
    )
  ),
  fluidRow(
    column(4, tableOutput("diag")),
    column(4, tableOutput("body_part")),
    column(4, tableOutput("location"))
  ),
  fluidRow(
    column(12, plotOutput("age_sex"))
  )
)

server <- function(input, output, session) {
  selected <- reactive(injuries %>% filter(prod_code == input$code))
  
  #<< tables
  output$diag <- renderTable(count_top(selected(), diag), width = "100%")
  output$body_part <- renderTable(count_top(selected(), body_part), width = "100%")
  output$location <- renderTable(count_top(selected(), location), width = "100%")
  #>>
  
  summary <- reactive({
    selected() %>%
      count(age, sex, wt = weight) %>%
      left_join(population, by = c("age", "sex")) %>%
      mutate(rate = n / population * 1e4)
  })
  
  output$age_sex <- renderPlot({
    summary() %>%
      ggplot(aes(age, n, colour = sex)) +
      geom_line() +
      labs(y = "Estimated number of injuries")
  }, res = 96)
}

shinyApp( ui, server )


