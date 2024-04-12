count_capital_letters <- function(string) {
  sum(gregexpr("[A-Z]", string, perl = TRUE)[[1]] > 0)
}