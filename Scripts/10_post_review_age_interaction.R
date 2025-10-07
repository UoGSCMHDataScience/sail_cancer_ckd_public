## analysis_whole_group.R

library(tidyverse)
library(broom)
library(survival)
library(mgcv)
library(gratia)
library(oddsratio)
library(ggpubr)

cohort_trim  <- readRDS("S:/1214 - The Relationship Between Cancer and Kidney Disease/CKD_Cancer/Data/cohort_trim.rds")

cohort_trim <- cohort_trim %>%
  mutate(neg_age10 = -age/10) %>%
  mutate(high_egfr10 = if_else(eGFR > 90, eGFR/10, 8.2))

## sex-stratified analysis for cancer survival - per decrease in eGFR

m1 <- coxph(Surv(ca_surv_time, death) ~ neg_egfr_10 + neg_age10 + WIMD_2011_DECILE + SMOKING_STATUS + 
                        mm_count_no_ckd + ca_type + stage_simplified,
                      data = cohort_trim %>% filter(sex == 1 & eGFR <90))  %>% 
  tidy() %>%
  mutate(model = "m1") %>%
  mutate(or = round(exp(estimate), digits = 2),
         lower = round(exp(estimate - 1.96*std.error), digits = 2),
         upper = round(exp(estimate + 1.96*std.error), digits = 2)) %>%
  select(c(model, term, estimate, std.error, or, lower, upper, p.value)) %>%
  mutate(p.value = round(p.value, digits = 5)) %>%
  filter(str_detect(term, c("egfr|age10")))

m2 <- coxph(Surv(ca_surv_time, death) ~ neg_egfr_10 * neg_age10 + WIMD_2011_DECILE + SMOKING_STATUS + 
                        mm_count_no_ckd + ca_type + stage_simplified,
                      data = cohort_trim %>% filter(sex == 1 & eGFR <90))  %>% 
  tidy() %>%
  mutate(model = "m2") %>%
  mutate(or = round(exp(estimate), digits = 2),
         lower = round(exp(estimate - 1.96*std.error), digits = 2),
         upper = round(exp(estimate + 1.96*std.error), digits = 2)) %>%
  select(c(model, term, estimate, std.error, or, lower, upper, p.value)) %>%
  mutate(p.value = round(p.value, digits = 5)) %>%
  filter(str_detect(term, c("egfr|age10")))


f1 <- coxph(Surv(ca_surv_time, death) ~ neg_egfr_10 + neg_age10 + WIMD_2011_DECILE + SMOKING_STATUS + 
                      mm_count_no_ckd + ca_type + stage_simplified,
                    data = cohort_trim %>% filter(sex == 0 & eGFR <90)) %>% 
  tidy() %>%
  mutate(model = "f1") %>%
  mutate(or = round(exp(estimate), digits = 2),
         lower = round(exp(estimate - 1.96*std.error), digits = 2),
         upper = round(exp(estimate + 1.96*std.error), digits = 2)) %>%
  select(c(model, term, estimate, std.error, or, lower, upper, p.value)) %>%
  mutate(p.value = round(p.value, digits = 5)) %>%
  filter(str_detect(term, c("egfr|age10")))

f2 <- coxph(Surv(ca_surv_time, death) ~ neg_egfr_10 * neg_age10 + WIMD_2011_DECILE + SMOKING_STATUS + 
                          mm_count_no_ckd + ca_type + stage_simplified,
                        data = cohort_trim %>% filter(sex == 0 & eGFR <90)) %>% 
  tidy() %>%
  mutate(model = "f2") %>%
  mutate(or = round(exp(estimate), digits = 2),
         lower = round(exp(estimate - 1.96*std.error), digits = 2),
         upper = round(exp(estimate + 1.96*std.error), digits = 2)) %>%
  select(c(model, term, estimate, std.error, or, lower, upper, p.value)) %>%
  mutate(p.value = round(p.value, digits = 5)) %>%
  filter(str_detect(term, c("egfr|age10")))


## sex-stratified analysis for cancer survival - per increase in eGFR

m3 <- coxph(Surv(ca_surv_time, death) ~ high_egfr10 + neg_age10 + WIMD_2011_DECILE + SMOKING_STATUS + 
                      mm_count_no_ckd + ca_type + stage_simplified,
                    data = cohort_trim %>% filter(sex == 1 & eGFR >75)) %>% 
  tidy() %>%
  mutate(model = "m3") %>%
  mutate(or = round(exp(estimate), digits = 2),
         lower = round(exp(estimate - 1.96*std.error), digits = 2),
         upper = round(exp(estimate + 1.96*std.error), digits = 2)) %>%
  select(c(model, term, estimate, std.error, or, lower, upper, p.value)) %>%
  mutate(p.value = round(p.value, digits = 5)) %>%
  filter(str_detect(term, c("egfr|age10")))

m4 <- coxph(Surv(ca_surv_time, death) ~ high_egfr10 * neg_age10 + WIMD_2011_DECILE + SMOKING_STATUS + 
                          mm_count_no_ckd + ca_type + stage_simplified,
                        data = cohort_trim %>% filter(sex == 1 & eGFR >75)) %>% 
  tidy() %>%
  mutate(model = "m4") %>%
  mutate(or = round(exp(estimate), digits = 2),
         lower = round(exp(estimate - 1.96*std.error), digits = 2),
         upper = round(exp(estimate + 1.96*std.error), digits = 2)) %>%
  select(c(model, term, estimate, std.error, or, lower, upper, p.value)) %>%
  mutate(p.value = round(p.value, digits = 5)) %>%
  filter(str_detect(term, c("egfr|age10")))

anova(surv_age_m, surv_age_int_m)


f3 <- coxph(Surv(ca_surv_time, death) ~ high_egfr10 + neg_age10 + WIMD_2011_DECILE + SMOKING_STATUS + 
                      mm_count_no_ckd + ca_type + stage_simplified,
                    data = cohort_trim %>% filter(sex == 0 & eGFR >75)) %>% 
  tidy() %>%
  mutate(model = "f3") %>%
  mutate(or = round(exp(estimate), digits = 2),
         lower = round(exp(estimate - 1.96*std.error), digits = 2),
         upper = round(exp(estimate + 1.96*std.error), digits = 2)) %>%
  select(c(model, term, estimate, std.error, or, lower, upper, p.value)) %>%
  mutate(p.value = round(p.value, digits = 5)) %>%
  filter(str_detect(term, c("egfr|age10")))

f4 <- coxph(Surv(ca_surv_time, death) ~ high_egfr10 * neg_age10 + WIMD_2011_DECILE + SMOKING_STATUS + 
                          mm_count_no_ckd + ca_type + stage_simplified,
                        data = cohort_trim %>% filter(sex == 0 & eGFR >75)) %>% 
  tidy() %>%
  mutate(model = "f4") %>% 
  mutate(or = round(exp(estimate), digits = 2),
         lower = round(exp(estimate - 1.96*std.error), digits = 2),
         upper = round(exp(estimate + 1.96*std.error), digits = 2)) %>%
  select(c(model, term, estimate, std.error, or, lower, upper, p.value)) %>%
  mutate(p.value = round(p.value, digits = 5)) %>%
  filter(str_detect(term, c("egfr|age10"))) 

survmod_age_interactions <- bind_rows(m1, m2, m3, m4, f1, f2, f3, f4) 
write_csv(survmod_age_interactions, "Outputs/survmod_age_interactions.csv")

