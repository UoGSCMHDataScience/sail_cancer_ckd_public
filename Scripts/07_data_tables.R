# 09_data_tables

library(tidyverse)
library(tableone)
library(dplyr)
library(stringr)

cohort_wimd  <- readRDS("S:/1214 - The Relationship Between Cancer and Kidney Disease/CKD_Cancer/Data/cohort_wimd.rds")
path <- "S:/1214 - The Relationship Between Cancer and Kidney Disease/CKD_Cancer/"

# create cohort omitting non-melanoma skin and haematological cancers
cohort_trim <- cohort_wimd %>%
  filter(str_detect(SITE_ICD10_O2, pattern = "C44|C76|C77|C78|C79|C80|C81|C82|C83|C84|C85|C86|C87|C88|C89|C90|C91|C92|C93|C94|C95|C96|C97", negate = TRUE))
cohort_colo <- cohort_wimd %>% 
  filter(str_detect(SITE_ICD10_O2, pattern = "C18|C19|C20", negate = FALSE))

# table1

base <- c("age", "eGFR", "ckd_cat", "SMOKING_STATUS", "mm_count_no_ckd",  "WIMD_2011_DECILE", "stage_simplified")
cat_base <- c("ckd_cat", "SMOKING_STATUS", "stage_simplified")

table1 <- CreateTableOne(vars = base, 
                         factorVars = cat_base, 
                         strata = "sex", 
                         includeNA = TRUE,
                         addOverall = TRUE,
                         data = cohort_colo)

write.csv(print(table1, 
        nonnormal = c("eGFR", "mm_count_no_ckd", "WIMD_2011_DECILE")), 
        "Outputs/table1.csv")

# table2

table2_overall <- cohort_trim %>%
  group_by(sex) %>%
  summarise(n = n(), 
            n_death = sum(death),
            pc_death = sum(death)/n()*100,
            time = median(ca_surv_time/365.25), 
            time_q1 = quantile(ca_surv_time/365.25, 0.25),
            time_q3 = quantile(ca_surv_time/365.25, 0.75)) %>%
  as_tibble() %>%
  mutate(surv = paste0(round(time, 1), " (", round(time_q1, 1), " - ", round(time_q3, 1), ")"),
         death = paste0(n_death, " (", round(pc_death, 1), ")")) %>%
  select(sex, n, death, surv)

table2_site <- cohort_trim %>%
  group_by(sex, ca_type) %>%
  summarise(n = n(), 
            n_death = sum(death),
            pc_death = sum(death)/n()*100,
            time = median(ca_surv_time/365.25), 
            time_q1 = quantile(ca_surv_time/365.25, 0.25),
            time_q3 = quantile(ca_surv_time/365.25, 0.75)) %>%
  as_tibble() %>%
  mutate(surv = paste0(round(time, 1), " (", round(time_q1, 1), " - ", round(time_q3, 1), ")"),
         death = paste0(n_death, " (", round(pc_death, 1), ")")) %>%
  select(sex, ca_type, n, death, surv)

table2 <- bind_rows(table2_overall, table2_site)

table2 <- table2 %>% 
  filter(n>400) %>% arrange(ca_type) %>% view()

write.csv(print(table2 %>% arrange(ca_type)),
          "Outputs/table2.csv")


table2_site_check <- cohort_trim %>%
  group_by(ca_type) %>%
  summarise(n = n(), 
            n_death = sum(death),
            pc_death = sum(death)/n()*100,
            time = median(ca_surv_time/365.25), 
            time_q1 = quantile(ca_surv_time/365.25, 0.25),
            time_q3 = quantile(ca_surv_time/365.25, 0.75)) %>%
  as_tibble() %>%
  mutate(surv = paste0(round(time, 1), " (", round(time_q1, 1), " - ", round(time_q3, 1), ")"),
         death = paste0(n_death, " (", round(pc_death, 1), ")")) %>%
  select(ca_type, n, death, surv)
table2_site_check

cohort_trim$ca_type %>% unique()



