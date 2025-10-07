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
cohort_trim <- cohort_trim %>%
  mutate(sex2 = factor(case_when(
    sex == 0 ~ "Female",
    sex == 1 ~ "Male"),
    levels = c("Male", "Female")))

## compare men versus women for advanced cancer

mvf1 <- glm(cancer_stage_adv == 1 ~ sex + age, family = "binomial",
            data = cohort_trim) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "mvf1")

mvf2 <- glm(cancer_stage_adv == 1 ~ sex + age + WIMD_2011_DECILE, family = "binomial", 
            data = cohort_trim) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "mvf2")

mvf3 <- glm(cancer_stage_adv == 1 ~ age + sex * ckd_cat + SMOKING_STATUS + mm_count_no_ckd + WIMD_2011_DECILE, family = "binomial", 
            data = cohort_trim) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "mvf3")

mvf4 <- glm(cancer_stage_adv == 1 ~ age + sex * ckd_cat + SMOKING_STATUS + mm_count_no_ckd + WIMD_2011_DECILE + ca_type, family = "binomial", 
            data = cohort_trim) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "mvf4")

mvf_adv_cont <- glm(cancer_stage_adv == 1 ~ age + sex * eGFR + SMOKING_STATUS + mm_count_no_ckd + WIMD_2011_DECILE + ca_type, 
                    family = "binomial", data = cohort_trim) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "mvf_adv_cont")

mvf_advca <- bind_rows(mvf1, mvf2, mvf3, mvf4, mvf_adv_cont) 

write_csv(mvf_advca, "Outputs/advca_interaction.csv")

mvf_advca %>% filter(model == "mvf4") %>% view()

## compare men versus women for cancer survival

mvf5 <- coxph(Surv(ca_surv_time, death) ~ sex + age, data = cohort_trim) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "mvf5")

mvf6 <- coxph(Surv(ca_surv_time, death) ~ age + sex * ckd_cat, data = cohort_trim) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "mvf6")

mvf7 <- coxph(Surv(ca_surv_time, death) ~ age + sex * ckd_cat + SMOKING_STATUS + mm_count_no_ckd + WIMD_2011_DECILE, data = cohort_trim) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "mvf7")

mvf8 <- coxph(Surv(ca_surv_time, death) ~ age + sex * ckd_cat + SMOKING_STATUS + mm_count_no_ckd + WIMD_2011_DECILE + ca_type, data = cohort_trim) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "mvf8")

mvf9 <- coxph(Surv(ca_surv_time, death) ~ age + sex * ckd_cat + SMOKING_STATUS + mm_count_no_ckd + WIMD_2011_DECILE + ca_type + stage_simplified, data = cohort_trim) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "mvf9")

mvf_surv_cont <- coxph(Surv(ca_surv_time, death) ~ age + sex * eGFR + SMOKING_STATUS + mm_count_no_ckd + WIMD_2011_DECILE + ca_type + stage_simplified, 
                       data = cohort_trim) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "mvf_surv_cont")

mvf_surv_12m <- bind_rows(mvf5, mvf6, mvf7, mvf8, mvf9, mvf_surv_cont)

write_csv(mvf_surv_12m, "Outputs/surv_12m_interaction.csv")

mvf_surv_12m %>% filter(model == "mvf9") %>% view()

## test mod sensitivity analyses

cohort_trim <- cohort_trim %>%
  mutate(max_stage = if_else(stage_simplified == "Unknown" | stage_simplified == "Unable to stage", "4", stage_simplified), 
         min_stage = if_else(stage_simplified == "Unknown" | stage_simplified == "Unable to stage", "1", stage_simplified))

sens_max_stage <- coxph(Surv(ca_surv_time, death) ~ age + sex * ckd_cat + SMOKING_STATUS + mm_count_no_ckd + WIMD_2011_DECILE + ca_type + max_stage, 
                        data = cohort_trim) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "sens_max_stage")

sens_min_stage <- coxph(Surv(ca_surv_time, death) ~ age + sex * ckd_cat + SMOKING_STATUS + mm_count_no_ckd + WIMD_2011_DECILE + ca_type + min_stage, 
                        data = cohort_trim) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "sens_min_stage")

sens_stage <- bind_rows(sens_max_stage, sens_min_stage)

write_csv(sens_stage, "Outputs/sens_stage.csv")
