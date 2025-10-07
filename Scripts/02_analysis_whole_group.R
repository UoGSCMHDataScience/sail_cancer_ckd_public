## analysis_whole_group.R

library(tidyverse)
library(broom)
library(survival)
library(mgcv)
library(gratia)
library(oddsratio)
library(ggpubr)

cohort_wimd  <- readRDS("S:/1214 - The Relationship Between Cancer and Kidney Disease/CKD_Cancer/Data/cohort_wimd.rds")
path <- "S:/1214 - The Relationship Between Cancer and Kidney Disease/CKD_Cancer/"

# create cohort omitting non-melanoma skin, haematological, unspecified and multiple site cancers
cohort_trim <- cohort_wimd %>%
  filter(str_detect(SITE_ICD10_O2, pattern = "C44|C76|C77|C78|C79|C80|C81|C82|C83|C84|C85|C86|C87|C88|C89|C90|C91|C92|C93|C94|C95|C96|C97", negate = TRUE))

## create new CKD categories

#cohort_wimd <- cohort_wimd %>%
#  mutate(ckd_cat = factor(case_when(
#    eGFR <= 30 ~ "eGFR <30", 
#    eGFR >30 & eGFR <= 45 ~ "eGFR 30 - <45",
#    eGFR >45 & eGFR <= 60 ~ "eGFR 45 - <60", 
#    eGFR >60 & eGFR <= 75 ~ "eGFR 60 - <75",
#    eGFR >75 & eGFR <= 90 ~ "eGFR 75 - <90", 
#    eGFR >90 & eGFR <= 105 ~ "eGFR 90 - <105",
#    eGFR >105 & eGFR <=120 ~ "eGFR 105 - <120",
#    eGFR >120 ~ "eGFR >120"),
#    levels = c("eGFR 75 - <90", "eGFR >120", "eGFR 105 - <120", "eGFR 90 - <105", "eGFR 60 - <75", "eGFR 45 - <60", "eGFR 30 - <45", "eGFR <30")
#  ))

## comparison of male versus female

mvf <- glm(cancer_stage_adv == 1 ~ sex + age,
             family = binomial, 
             data = cohort_trim) %>% tidy(exponentiate = TRUE, conf.int = TRUE)
mvf1 <- glm(cancer_stage_adv == 1 ~ sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + ca_type,
           family = binomial, 
           data = cohort_trim) %>% tidy(exponentiate = TRUE, conf.int = TRUE)

mvf_s <- coxph(Surv(ca_surv_time, death) ~ sex + age,
           data = cohort_trim) %>% tidy(exponentiate = TRUE, conf.int = TRUE)
mvf1_s <- coxph(Surv(ca_surv_time, death) ~ sex + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + ca_type,
            data = cohort_trim) %>% tidy(exponentiate = TRUE, conf.int = TRUE)

mvf <- bind_rows(mvf %>% mutate(model = "mvf_adv_age"),
                 mvf1 %>% mutate(model = "mvf_adv_full"),
                 mvf_s %>% mutate(model = "mvf_surv_age"),
                 mvf1_s %>% mutate(model = "mvf_surv_full"))

write_csv(mvf, "Outputs/mvf.csv")


## sex-stratified analysis for advanced cancer

cohort_nst <- cohort_trim %>%
  group_by(sex) %>%
  nest() %>%
  ungroup()

modlist <- list(
  mod1 = cancer_stage_adv == 1 ~ ckd_cat + age,
  mod2 = cancer_stage_adv == 1 ~ ckd_cat + age + WIMD_2011_DECILE,
  mod3 = cancer_stage_adv == 1 ~ ckd_cat + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd,
  mod4 = cancer_stage_adv == 1 ~ ckd_cat + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + ca_type,
  mod4a = cancer_stage_adv == 1 ~ neg_egfr_10 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + ca_type
)

advca_modelres <- map(modlist, function(form) {
  map(cohort_nst$data, function(data_choose) glm(formula = form, family = "binomial", data = data_choose))
})

advca_modelres <- do.call(c, advca_modelres)
advca_modelres <- tibble(tosep = names(advca_modelres), res = advca_modelres)
advca_modelres$advca_smry <- map(advca_modelres$res, tidy)

advca_unnst <- advca_modelres %>%
  select(tosep, advca_smry) %>%
  unnest(cols = c(advca_smry))

advca_modelcofs <- advca_unnst %>%
  separate(tosep, into = c("model", "sex"),
           sep = 4) %>%
  mutate(sex_name = if_else(sex == 1, "Male", "Female")) %>%
  mutate(or = round(exp(estimate), digits = 2),
         lower = round(exp(estimate - 1.96*std.error), digits = 2),
         upper = round(exp(estimate + 1.96*std.error), digits = 2)) %>%
  mutate(p.value = round(p.value, digits = 5))

advca <- advca_modelcofs %>%
  filter(str_detect(term, "ckd_cat")) %>%
  separate(term, into = c("term", "ckd_cat"), sep = 7) %>%
  add_row(tibble(model = "mod1", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Male", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod1", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Female", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod2", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Male", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod2", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Female", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod3", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Male", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod3", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Female", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod4", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Male", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod4", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Female", or = 1, lower = 1, upper = 1, p.value = 1)) 

advca$ckd_cat <- factor(advca$ckd_cat, 
                        levels = c("eGFR <30", "eGFR 30 - <45", "eGFR 45 - <60", "eGFR 60 - <75", "eGFR 75 - <90", "eGFR 90 - <105", "eGFR 105 - <120", "eGFR >120"),
                        labels = c("<30", "30 - <45", "45 - <60", "60 - <75", "75 - <90", "90 - <105", "105 - <120", ">120"))

write_csv(advca, "Outputs/advanced_cancer_coefficients_sex_stratified.csv")

advplot <- advca %>%  
  filter(model == "mod4") %>%
  ggplot(aes(x = ckd_cat, y = or, ymin = lower, ymax = upper, colour = sex_name)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) + 
  scale_y_continuous(trans = "log10", limits = c(0.8, 2)) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), size = 0.3, width = 0.3, cex = 0.2, pch = 20, position = position_dodge(width = 0.4)) + 
  labs(x = "eGFRcr prior to cancer diagnosis", y = "Adjusted Odds Ratio (95% confidence interval)", title = "Diagnosis of advanced cancer", colour = "Sex") +
  theme_bw() + 
  theme(text = element_text(size = 20))
advplot  

## sex-stratified analysis for cancer survival

cohort_nst <- cohort_trim %>%
  group_by(sex) %>%
  nest() %>%
  ungroup()

modlist2 <- list(
  mod1 = Surv(ca_surv_time, death) ~ ckd_cat + age,
  mod2 = Surv(ca_surv_time, death) ~ ckd_cat + age + WIMD_2011_DECILE,
  mod3 = Surv(ca_surv_time, death) ~ ckd_cat + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd,
  mod4 = Surv(ca_surv_time, death) ~ ckd_cat + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + ca_type,
  mod5 = Surv(ca_surv_time, death) ~ ckd_cat + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + ca_type + stage_simplified,
  mod5a = Surv(ca_surv_time, death) ~ neg_egfr_10 + age + WIMD_2011_DECILE + SMOKING_STATUS + mm_count_no_ckd + ca_type + stage_simplified
)

survca_modelres <- map(modlist2, function(form) {
  map(cohort_nst$data, function(data_choose) coxph(formula = form, data = data_choose))
})

survca_modelres <- do.call(c, survca_modelres)
survca_modelres <- tibble(tosep = names(survca_modelres), res = survca_modelres)
survca_modelres$survca_smry <- map(survca_modelres$res, tidy)

survca_unnst <- survca_modelres %>%
  select(tosep, survca_smry) %>%
  unnest(cols = c(survca_smry))

survca_modelcofs <- survca_unnst %>%
  separate(tosep, into = c("model", "sex"),
           sep = 4) %>%
  mutate(sex_name = if_else(sex == 1, "Male", "Female")) %>%
  mutate(or = round(exp(estimate), digits = 2),
         lower = round(exp(estimate - 1.96*std.error), digits = 2),
         upper = round(exp(estimate + 1.96*std.error), digits = 2)) %>%
  mutate(p.value = round(p.value, digits = 5))

survca <- survca_modelcofs %>%
  filter(str_detect(term, "ckd_cat")) %>%
  separate(term, into = c("term", "ckd_cat"), sep = 7) %>%
  add_row(tibble(model = "mod1", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Male", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod1", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Female", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod2", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Male", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod2", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Female", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod3", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Male", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod3", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Female", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod4", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Male", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod4", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Female", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod5", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Male", or = 1, lower = 1, upper = 1, p.value = 1)) %>%
  add_row(tibble(model = "mod5", term = "ckd_cat", ckd_cat = "eGFR 75 - <90", sex_name = "Female", or = 1, lower = 1, upper = 1, p.value = 1)) 

survca$ckd_cat <- factor(survca$ckd_cat, 
                        levels = c("eGFR <30", "eGFR 30 - <45", "eGFR 45 - <60", "eGFR 60 - <75", "eGFR 75 - <90", "eGFR 90 - <105", "eGFR 105 - <120", "eGFR >120"),
                        labels = c("<30", "30 - <45", "45 - <60", "60 - <75", "75 - <90", "90 - <105", "105 - <120", ">120"))

write_csv(survca, "Outputs/survival_coefficients_sex_stratified.csv")

survplot <- survca %>%  
  filter(model == "mod5") %>%
  ggplot(aes(x = ckd_cat, y = or, ymin = lower, ymax = upper, colour = sex_name)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) + 
  scale_y_continuous(trans = "log10", limits = c(0.8, 4.5)) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), size = 0.3, width = 0.3, cex = 0.2, pch = 20, position = position_dodge(width = 0.4)) + 
  labs(x = "eGFRcr prior to cancer diagnosis", y = "Adjusted Hazard Ratio (95% confidence interval)", title = "Death after cancer diagnosis", colour = "Sex") +
  theme_bw() + 
  theme(text = element_text(size = 20))
survplot

ggarrange(advplot, survplot, 
          common.legend = TRUE)


## standardised mortality ratios

cohort_trim$age %>% summary()
cohort_trim <- cohort_trim %>%
  mutate(age_sex_cat = factor(case_when(
    sex == 1 & age < 60 ~ "Male <60",
    sex == 1 & age >= 60 & age < 70 ~ "Male 60-<70",
    sex == 1 & age >= 70 & age < 80 ~ "Male 70-<80",
    sex == 1 & age >= 80 ~ "Male >=80",
    sex == 0 & age < 60 ~ "Female <60",
    sex == 0 & age >= 60 & age < 70 ~ "Female 60-<70",
    sex == 0 & age >= 70 & age < 80 ~ "Female 70-<80",
    sex == 0 & age >= 80 ~ "Female >=80"),
    levels = c("Male <60", "Male 60-<70", "Male 70-<80", "Male >=80", "Female <60", "Female 60-<70", "Female 70-<80", "Female >=80")))

cohort_trim$age_sex_cat %>% as.factor() %>% summary()
group_tbl <-
  cohort_trim %>%
  group_by(age_sex_cat, ckd_cat2, death) %>%
  tally() %>%
  spread(death, n)

group_tbl <- group_tbl %>%
  mutate(death_rate = `1`/(`1`+`0`)) %>%
  select(age_sex_cat, ckd_cat2, death_rate) %>%
  mutate(death_rate = round(death_rate, 2)) 

group_tbl <- group_tbl %>%
  pivot_wider(names_from = age_sex_cat, values_from = death_rate)
  

write_csv(group_tbl, "Outputs/group_tbl.csv")
