options(repos=c(CRAN="https://cran.r-project.org"))
if (getRversion() < "3.5.0") {
  stop("The scripts in this repository require R >= 3.5.0")
}
install.packages("pak")
required_packages <- c(
  "ashr",
  "data.table",
  "glmnet",
  "glue",
  "ggplot2",
  "ggraph",
  "igraph",
  "magrittr",
  "ncvreg",
  "stringr"
)
pak::pkg_install(required_packages)
file.create(".r-installed")
