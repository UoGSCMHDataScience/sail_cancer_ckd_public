# 06_plots_and_figures.R

library(tidyverse)
library(tidycmprsk)
library(ggsurvfit)
library(ggpubr)

cohort_wimd  <- readRDS("S:/1214 - The Relationship Between Cancer and Kidney Disease/CKD_Cancer/Data/cohort_wimd.rds")
path <- "S:/1214 - The Relationship Between Cancer and Kidney Disease/CKD_Cancer/"

## density plot of kidney function men and women

cohort_wimd %>%
  ggplot(aes(x = eGFR, colour = sex)) +
  geom_density()


ss_advplot <- ss_adv_cancer %>%  
  filter(str_detect(term, "ckd_cat2")) %>%
  separate(term, into = c("term", "ckd_cat"), sep = 8) %>%
  separate(ckd_cat, into = c("egfr", "ckd_cat"), sep = 5) %>%
  view()
ss_adv_cancer %>% distinct(model)
ss_surv_cancer %>% distinct(model)


ss_survplot <- ss_surv_cancer %>%  
  filter(str_detect(term, "ckd_cat2")) %>%
  separate(term, into = c("term", "ckd_cat"), sep = 8) %>%
  separate(ckd_cat, into = c("egfr", "ckd_cat"), sep = 5) %>%
  view()
ss_survplot


