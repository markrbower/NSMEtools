download <- function(name) {
  url <- "https://github.com/hadley/mastering-shiny/blob/main/neiss/"
  download.file(paste0(url, name, "?raw=true"), paste0("neiss/", name), quiet = TRUE)
}
download("injuries.tsv.gz")
injuries <- vroom::vroom("neiss/injuries.tsv.gz", show_col_types = FALSE)
download("population.tsv")
population <- vroom::vroom("neiss/population.tsv", show_col_types = FALSE)
download("products.tsv")
products <- vroom::vroom("neiss/products.tsv", show_col_types = FALSE)
