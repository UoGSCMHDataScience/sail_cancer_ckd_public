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

# create cohort omitting non-melanoma skin and 
cohort_trim <- cohort_wimd %>%
  filter(str_detect(SITE_ICD10_O2, pattern = "C44|C76|C77|C78|C79|C80|C81|C82|C83|C84|C85|C86|C87|C88|C89|C90|C91|C92|C93|C94|C95|C96|C97", negate = TRUE))

## sex-stratified analysis for advanced cancer

men <- cohort_trim %>% filter(sex == "1")
women <- cohort_trim %>% filter(sex == "0")

## men

abdo_adv_m <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = men %>% filter(ca_type == "abdo")) %>% tidy() %>% mutate(sex = "male", model = "abdo")
dig_adv_m <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = men %>% filter(ca_type == "dig")) %>% tidy()  %>% mutate(sex = "male", model = "dig")
#haem_adv_m <- 
#  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
#      family = "binomial",
#      data = men %>% filter(ca_type == "haem")) %>% tidy()  %>% mutate(sex = "male", model = "haem")
headneck_adv_m <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = men %>% filter(ca_type == "headneck")) %>% tidy() %>% mutate(sex = "male", model = "headneck")
kub_adv_m <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = men %>% filter(ca_type == "kub")) %>% tidy() %>% mutate(sex = "male", model = "kub")
male_adv_m <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = men %>% filter(ca_type == "male")) %>% tidy() %>% mutate(sex = "male", model = "male")
melanoma_adv_m <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = men %>% filter(ca_type == "melanoma")) %>% tidy() %>% mutate(sex = "male", model = "melanoma")
#myeloma_adv_m <- 
#  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
#     family = "binomial",
#      data = men %>% filter(ca_type == "myeloma")) %>% tidy() %>% mutate(sex = "male", model = "myeloma")
other_adv_m <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = men %>% filter(ca_type == "other")) %>% tidy() %>% mutate(sex = "male", model = "other")
prostate_adv_m <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = men %>% filter(ca_type == "prostate")) %>% tidy() %>% mutate(sex = "male", model = "prostate")
resp_adv_m <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = men %>% filter(ca_type == "resp")) %>% tidy() %>% mutate(sex = "male", model = "resp")
#skin_adv_m <- 
#  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
#      family = "binomial",
#      data = men %>% filter(ca_type == "skin")) %>% tidy() %>% mutate(sex = "male", model = "skin")
thyroid_adv_m <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = men %>% filter(ca_type == "thyroid")) %>% tidy() %>% mutate(sex = "male", model = "thyroid")


## women

abdo_adv_f <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = women %>% filter(ca_type == "abdo")) %>% tidy() %>% mutate(sex = "female", model = "abdo")
dig_adv_f <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = women %>% filter(ca_type == "dig")) %>% tidy() %>% mutate(sex = "female", model = "dig")
#haem_adv_f <- 
#  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
#      family = "binomial",
#      data = women %>% filter(ca_type == "haem")) %>% tidy() %>% mutate(sex = "female", model = "haem")
headneck_adv_f <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = women %>% filter(ca_type == "headneck")) %>% tidy() %>% mutate(sex = "female", model = "headneck")
kub_adv_f <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = women %>% filter(ca_type == "kub")) %>% tidy() %>% mutate(sex = "female", model = "kub")
female_adv_f <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = women %>% filter(ca_type == "female")) %>% tidy() %>% mutate(sex = "female", model = "female")
melanoma_adv_f <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = women %>% filter(ca_type == "melanoma")) %>% tidy() %>% mutate(sex = "female", model = "melanoma")
#myeloma_adv_f <- 
#  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
#      family = "binomial",
#      data = women %>% filter(ca_type == "myeloma")) %>% tidy() %>% mutate(sex = "female", model = "myeloma")
other_adv_f <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = women %>% filter(ca_type == "other")) %>% tidy() %>% mutate(sex = "female", model = "other")
breast_adv_f <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = women %>% filter(ca_type == "breast")) %>% tidy() %>% mutate(sex = "female", model = "breast")
resp_adv_f <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = women %>% filter(ca_type == "resp")) %>% tidy() %>% mutate(sex = "female", model = "resp")
#skin_adv_f <- 
#  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
#      family = "binomial",
#      data = women %>% filter(ca_type == "skin")) %>% tidy() %>% mutate(sex = "female", model = "skin")
thyroid_adv_f <- 
  glm(cancer_stage_adv == 1 ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd, 
      family = "binomial",
      data = women %>% filter(ca_type == "thyroid")) %>% tidy() %>% mutate(sex = "female", model = "thyroid")

ss_adv_cancer <- bind_rows(abdo_adv_m, dig_adv_m, headneck_adv_m, kub_adv_m, male_adv_m, melanoma_adv_m,
                           other_adv_m, prostate_adv_m, resp_adv_m, thyroid_adv_m,
                           abdo_adv_f, dig_adv_f, headneck_adv_f, kub_adv_f, female_adv_f, melanoma_adv_f,
                           other_adv_f, breast_adv_f, resp_adv_f, thyroid_adv_f)

write_csv(ss_adv_cancer, "Outputs/ss_adv_cancer.csv")
  
## sex-stratified analysis for survival

## men

abdo_surv_m <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = men %>% filter(ca_type == "abdo")) %>% tidy() %>% mutate(sex = "male", model = "abdo")
dig_surv_m <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = men %>% filter(ca_type == "dig")) %>% tidy()  %>% mutate(sex = "male", model = "dig")
#haem_surv_m <- 
#  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
#      data = men %>% filter(ca_type == "haem")) %>% tidy()  %>% mutate(sex = "male", model = "haem")
headneck_surv_m <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = men %>% filter(ca_type == "headneck")) %>% tidy() %>% mutate(sex = "male", model = "headneck")
kub_surv_m <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = men %>% filter(ca_type == "kub")) %>% tidy() %>% mutate(sex = "male", model = "kub")
male_surv_m <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = men %>% filter(ca_type == "male")) %>% tidy() %>% mutate(sex = "male", model = "male")
melanoma_surv_m <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = men %>% filter(ca_type == "melanoma")) %>% tidy() %>% mutate(sex = "male", model = "melanoma")
#myeloma_surv_m <- 
#  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
#      data = men %>% filter(ca_type == "myeloma")) %>% tidy() %>% mutate(sex = "male", model = "myeloma")
other_surv_m <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = men %>% filter(ca_type == "other")) %>% tidy() %>% mutate(sex = "male", model = "other")
prostate_surv_m <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = men %>% filter(ca_type == "prostate")) %>% tidy() %>% mutate(sex = "male", model = "prostate")
resp_surv_m <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = men %>% filter(ca_type == "resp")) %>% tidy() %>% mutate(sex = "male", model = "resp")
#skin_surv_m <- 
#  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
#      data = men %>% filter(ca_type == "skin")) %>% tidy() %>% mutate(sex = "male", model = "skin")
thyroid_surv_m <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = men %>% filter(ca_type == "thyroid")) %>% tidy() %>% mutate(sex = "male", model = "thyroid")


## women

abdo_surv_f <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = women %>% filter(ca_type == "abdo")) %>% tidy() %>% mutate(sex = "female", model = "abdo")
dig_surv_f <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = women %>% filter(ca_type == "dig")) %>% tidy() %>% mutate(sex = "female", model = "dig")
#haem_surv_f <- 
#  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
#      data = women %>% filter(ca_type == "haem")) %>% tidy() %>% mutate(sex = "female", model = "haem")
headneck_surv_f <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = women %>% filter(ca_type == "headneck")) %>% tidy() %>% mutate(sex = "female", model = "headneck")
kub_surv_f <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = women %>% filter(ca_type == "kub")) %>% tidy() %>% mutate(sex = "female", model = "kub")
female_surv_f <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = women %>% filter(ca_type == "female")) %>% tidy() %>% mutate(sex = "female", model = "female")
melanoma_surv_f <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = women %>% filter(ca_type == "melanoma")) %>% tidy() %>% mutate(sex = "female", model = "melanoma")
#myeloma_surv_f <- 
#  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
#      data = women %>% filter(ca_type == "myeloma")) %>% tidy() %>% mutate(sex = "female", model = "myeloma")
other_surv_f <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = women %>% filter(ca_type == "other")) %>% tidy() %>% mutate(sex = "female", model = "other")
breast_surv_f <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = women %>% filter(ca_type == "breast")) %>% tidy() %>% mutate(sex = "female", model = "breast")
resp_surv_f <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = women %>% filter(ca_type == "resp")) %>% tidy() %>% mutate(sex = "female", model = "resp")
#skin_surv_f <- 
#  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
#      data = women %>% filter(ca_type == "skin")) %>% tidy() %>% mutate(sex = "female", model = "skin")
thyroid_surv_f <- 
  coxph(Surv(ca_surv_time, death) ~ ckd_cat2 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + stage_simplified, 
      data = women %>% filter(ca_type == "thyroid")) %>% tidy() %>% mutate(sex = "female", model = "thyroid")

ss_surv_cancer <- bind_rows(abdo_surv_m, dig_surv_m, headneck_surv_m, kub_surv_m, male_surv_m, melanoma_surv_m,
                           other_surv_m, prostate_surv_m, resp_surv_m, thyroid_surv_m,
                           abdo_surv_f, dig_surv_f, headneck_surv_f, kub_surv_f, female_surv_f, melanoma_surv_f,
                           other_surv_f, breast_surv_f, resp_surv_f, thyroid_surv_f)

write_csv(ss_surv_cancer, "Outputs/ss_surv_cancer.csv")
