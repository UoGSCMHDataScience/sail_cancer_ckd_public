## analysis_for_era.R

library(tidyverse)
library(broom)
library(survival)
library(mgcv)
library(gratia)
library(oddsratio)
library(ggpubr)

cohort_wimd  <- readRDS("S:/1214 - The Relationship Between Cancer and Kidney Disease/CKD_Cancer/Data/cohort_wimd.rds")
path <- "S:/1214 - The Relationship Between Cancer and Kidney Disease/CKD_Cancer/"

# create cohort omitting non-melanoma skin and haematological cancers
cohort_trim <- cohort_wimd %>%
  filter(str_detect(SITE_ICD10_O2, pattern = "C44|C76|C77|C78|C79|C80|C81|C82|C83|C84|C85|C86|C87|C88|C89|C90|C91|C92|C93|C94|C95|C96|C97", negate = TRUE))

## interaction terms for advanced cancer

abdo_adv <- 
  glm(cancer_stage_adv == 1 ~  ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "abdo")) %>% tidy() %>% mutate(model = "abdo")
dig_adv <- 
  glm(cancer_stage_adv == 1 ~  ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "dig")) %>% tidy()  %>% mutate(model = "dig")
#haem_adv <- 
#  glm(cancer_stage_adv == 1 ~  ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
#      family = "binomial",
#      data = cohort_trim %>% filter(ca_type == "haem")) %>% tidy()  %>% mutate(model = "haem")
headneck_adv <- 
  glm(cancer_stage_adv == 1 ~  ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "headneck")) %>% tidy() %>% mutate(model = "headneck")
kub_adv <- 
  glm(cancer_stage_adv == 1 ~  ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "kub")) %>% tidy() %>% mutate(model = "kub")
melanoma_adv <- 
  glm(cancer_stage_adv == 1 ~  ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "melanoma")) %>% tidy() %>% mutate(model = "melanoma")
#myeloma_adv <- 
#  glm(cancer_stage_adv == 1 ~  ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
#      family = "binomial",
#      data = cohort_trim %>% filter(ca_type == "myeloma")) %>% tidy() %>% mutate(model = "myeloma")
other_adv <- 
  glm(cancer_stage_adv == 1 ~  ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "other")) %>% tidy() %>% mutate(model = "other")
resp_adv <- 
  glm(cancer_stage_adv == 1 ~  ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "resp")) %>% tidy() %>% mutate(model = "resp")
#skin_adv <- 
#  glm(cancer_stage_adv == 1 ~  ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
#      family = "binomial",
#      data = cohort_trim %>% filter(ca_type == "skin")) %>% tidy() %>% mutate(model = "skin")
thyroid_adv <- 
  glm(cancer_stage_adv == 1 ~  ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "thyroid")) %>% tidy() %>% mutate(model = "thyroid")


ss_adv_cancer <- bind_rows(abdo_adv, dig_adv, headneck_adv, kub_adv, melanoma_adv,
                           other_adv, resp_adv, thyroid_adv)

ss_adv_cancer <- ss_adv_cancer %>%
  mutate(or = round(exp(estimate), digits = 2),
         lower = round(exp(estimate - 1.96*std.error), digits = 2),
         upper = round(exp(estimate + 1.96*std.error), digits = 2)) %>%
  select(c(model, term, estimate, std.error, or, lower, upper, p.value)) %>%
  mutate(p.value = round(p.value, digits = 5))

write_csv(ss_adv_cancer, "Outputs/ss_adv_cancer_interaction.csv")

## interaction terms for advanced cancer - continuous eGFR

abdo_adv_cont <- 
  glm(cancer_stage_adv == 1 ~  eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "abdo")) %>% tidy() %>% mutate(model = "abdo")
dig_adv_cont <- 
  glm(cancer_stage_adv == 1 ~  eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "dig")) %>% tidy()  %>% mutate(model = "dig")
#haem_adv_cont <- 
#  glm(cancer_stage_adv == 1 ~  eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
#      family = "binomial",
#      data = cohort_trim %>% filter(ca_type == "haem")) %>% tidy()  %>% mutate(model = "haem")
headneck_adv_cont <- 
  glm(cancer_stage_adv == 1 ~  eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "headneck")) %>% tidy() %>% mutate(model = "headneck")
kub_adv_cont <- 
  glm(cancer_stage_adv == 1 ~  eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "kub")) %>% tidy() %>% mutate(model = "kub")
melanoma_adv_cont <- 
  glm(cancer_stage_adv == 1 ~  eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "melanoma")) %>% tidy() %>% mutate(model = "melanoma")
#myeloma_adv_cont <- 
#  glm(cancer_stage_adv == 1 ~  eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
#      family = "binomial",
#      data = cohort_trim %>% filter(ca_type == "myeloma")) %>% tidy() %>% mutate(model = "myeloma")
other_adv_cont <- 
  glm(cancer_stage_adv == 1 ~  eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "other")) %>% tidy() %>% mutate(model = "other")
resp_adv_cont <- 
  glm(cancer_stage_adv == 1 ~  eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "resp")) %>% tidy() %>% mutate(model = "resp")
#skin_adv_cont <- 
#  glm(cancer_stage_adv == 1 ~  eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
#      family = "binomial",
#      data = cohort_trim %>% filter(ca_type == "skin")) %>% tidy() %>% mutate(model = "skin")
thyroid_adv_cont <- 
  glm(cancer_stage_adv == 1 ~  eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = cohort_trim %>% filter(ca_type == "thyroid")) %>% tidy() %>% mutate(model = "thyroid")


ss_adv_cancer_cont <- bind_rows(abdo_adv_cont, dig_adv_cont, headneck_adv_cont, kub_adv_cont, melanoma_adv_cont,
                           other_adv_cont, resp_adv_cont, thyroid_adv_cont)

ss_adv_cancer_cont <- ss_adv_cancer_cont %>%
  mutate(or = round(exp(estimate), digits = 2),
         lower = round(exp(estimate - 1.96*std.error), digits = 2),
         upper = round(exp(estimate + 1.96*std.error), digits = 2)) %>%
  select(c(model, term, estimate, std.error, or, lower, upper, p.value)) %>%
  mutate(p.value = round(p.value, digits = 5))

write_csv(ss_adv_cancer_cont, "Outputs/ss_adv_cancer_cont_interaction.csv")
  
## interaction terms for survival - ckd categories

abdo_surv <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = cohort_trim %>% filter(ca_type == "abdo")) %>% tidy() %>% mutate(model = "abdo")
dig_surv <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = cohort_trim %>% filter(ca_type == "dig")) %>% tidy()  %>% mutate(model = "dig")
#haem_surv <- 
#  coxph(Surv(ca_surv_time, death) ~ ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
#      data = cohort_trim %>% filter(ca_type == "haem")) %>% tidy()  %>% mutate(model = "haem")
headneck_surv <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = cohort_trim %>% filter(ca_type == "headneck")) %>% tidy() %>% mutate(model = "headneck")
kub_surv <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = cohort_trim %>% filter(ca_type == "kub")) %>% tidy() %>% mutate(model = "kub")
male_surv <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = cohort_trim %>% filter(ca_type == "male")) %>% tidy() %>% mutate(model = "male")
melanoma_surv <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = cohort_trim %>% filter(ca_type == "melanoma")) %>% tidy() %>% mutate(model = "melanoma")
#myeloma_surv <- 
#  coxph(Surv(ca_surv_time, death) ~ ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
#      data = cohort_trim %>% filter(ca_type == "myeloma")) %>% tidy() %>% mutate(model = "myeloma")
other_surv <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = cohort_trim %>% filter(ca_type == "other")) %>% tidy() %>% mutate(model = "other")
prostate_surv <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = cohort_trim %>% filter(ca_type == "prostate")) %>% tidy() %>% mutate(model = "prostate")
resp_surv <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = cohort_trim %>% filter(ca_type == "resp")) %>% tidy() %>% mutate(model = "resp")
#skin_surv <- 
#  coxph(Surv(ca_surv_time, death) ~ ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
#      data = cohort_trim %>% filter(ca_type == "skin")) %>% tidy() %>% mutate(model = "skin")
thyroid_surv <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = cohort_trim %>% filter(ca_type == "thyroid")) %>% tidy() %>% mutate(model = "thyroid")

ss_surv_cancer <- bind_rows(abdo_surv, dig_surv, headneck_surv, kub_surv, melanoma_surv,
                           other_surv, resp_surv, thyroid_surv)

ss_surv_cancer <- ss_surv_cancer  %>%
  mutate(or = round(exp(estimate), digits = 2),
         lower = round(exp(estimate - 1.96*std.error), digits = 2),
         upper = round(exp(estimate + 1.96*std.error), digits = 2)) %>%
  select(c(model, term, estimate, std.error, or, lower, upper, p.value)) %>%
  mutate(p.value = round(p.value, digits = 5))

write_csv(ss_surv_cancer, "Outputs/ss_surv_cancer_interaction.csv")

## interaction terms for survival - continuous egfr

abdo_surv_cont <- 
  coxph(Surv(ca_surv_time, death) ~ eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
        data = cohort_trim %>% filter(ca_type == "abdo")) %>% tidy() %>% mutate(model = "abdo")
dig_surv_cont <- 
  coxph(Surv(ca_surv_time, death) ~ eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
        data = cohort_trim %>% filter(ca_type == "dig")) %>% tidy()  %>% mutate(model = "dig")
#haem_surv_cont <- 
#  coxph(Surv(ca_surv_time, death) ~ eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
#        data = cohort_trim %>% filter(ca_type == "haem")) %>% tidy()  %>% mutate(model = "haem")
headneck_surv_cont <- 
  coxph(Surv(ca_surv_time, death) ~ eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
        data = cohort_trim %>% filter(ca_type == "headneck")) %>% tidy() %>% mutate(model = "headneck")
kub_surv_cont <- 
  coxph(Surv(ca_surv_time, death) ~ eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
        data = cohort_trim %>% filter(ca_type == "kub")) %>% tidy() %>% mutate(model = "kub")
male_surv_cont <- 
  coxph(Surv(ca_surv_time, death) ~ eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
        data = cohort_trim %>% filter(ca_type == "male")) %>% tidy() %>% mutate(model = "male")
melanoma_surv_cont <- 
  coxph(Surv(ca_surv_time, death) ~ eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
        data = cohort_trim %>% filter(ca_type == "melanoma")) %>% tidy() %>% mutate(model = "melanoma")
#myeloma_surv_cont <- 
#  coxph(Surv(ca_surv_time, death) ~ eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
#        data = cohort_trim %>% filter(ca_type == "myeloma")) %>% tidy() %>% mutate(model = "myeloma")
other_surv_cont <- 
  coxph(Surv(ca_surv_time, death) ~ eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
        data = cohort_trim %>% filter(ca_type == "other")) %>% tidy() %>% mutate(model = "other")
prostate_surv_cont <- 
  coxph(Surv(ca_surv_time, death) ~ eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
        data = cohort_trim %>% filter(ca_type == "prostate")) %>% tidy() %>% mutate(model = "prostate")
resp_surv_cont <- 
  coxph(Surv(ca_surv_time, death) ~ eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
        data = cohort_trim %>% filter(ca_type == "resp")) %>% tidy() %>% mutate(model = "resp")
#skin_surv_cont <- 
#  coxph(Surv(ca_surv_time, death) ~ eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
#        data = cohort_trim %>% filter(ca_type == "skin")) %>% tidy() %>% mutate(model = "skin")
thyroid_surv_cont <- 
  coxph(Surv(ca_surv_time, death) ~ eGFR*sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
        data = cohort_trim %>% filter(ca_type == "thyroid")) %>% tidy() %>% mutate(model = "thyroid")

ss_surv_cancer_cont <- bind_rows(abdo_surv_cont, dig_surv_cont, headneck_surv_cont, kub_surv_cont, melanoma_surv_cont,
                            other_surv_cont, resp_surv_cont, thyroid_surv_cont)

ss_surv_cancer_cont <- ss_surv_cancer_cont  %>%
  mutate(or = round(exp(estimate), digits = 2),
         lower = round(exp(estimate - 1.96*std.error), digits = 2),
         upper = round(exp(estimate + 1.96*std.error), digits = 2)) %>%
  select(c(model, term, estimate, std.error, or, lower, upper, p.value)) %>%
  mutate(p.value = round(p.value, digits = 5))

write_csv(ss_surv_cancer_cont, "Outputs/ss_surv_cancer_cont_interaction.csv")

ss_surv_cancer_cont %>% view()
