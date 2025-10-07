#00_functions.R

## 

library(tidyverse)
library(data.table)
library(kableExtra)
library(finalfit)
library(broom.mixed)

ExportRes <- function(modelname){
  #browser()
  a <- broom.mixed::tidy(modelname)
  b <- posterior_interval(modelname, 0.95) %>%
    as_tibble(rownames = "term")
  a <- a %>%
    inner_join(b)
  a <- a %>%
    mutate_at(vars(estimate, `2.5%`, `97.5%`), function(x) x %>%
                  exp() %>%
                  formatC(digits = 2, format = "f", flag = "0")) %>%
    mutate(res = paste0(estimate," (", `2.5%`, "-", `97.5%`, ")"))
  a %>%
    select(term,res)
}

