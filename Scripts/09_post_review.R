# 09_post_review

library(tidyverse)
library(lme4)
library(nephro)
library(RODBC)
library(stringr)
library(tcltk)
library(knitr)
library(tableone)
library(ggpubr)
library(broom)
library(survival)

#cohort  <- readRDS("S:/1214 - The Relationship Between Cancer and Kidney Disease/CKD_Cancer/Data/cohort.rds")
#cohort_trim  <- readRDS("S:/1214 - The Relationship Between Cancer and Kidney Disease/CKD_Cancer/Data/cohort_trim.rds")
cohort_add  <- readRDS("Scratch_data/cohort_add.rds")

## login and extract ----

getLogin <- function(userName=''){
  require(tcltk);
  wnd <- tktoplevel();
  tclVar(userName) -> user;
  tclVar("") -> passVar;
  
  tkgrid(tklabel(wnd, text="Username:"));
  tkgrid(tkentry(wnd, textvariable=user) -> passBox);
  
  tkgrid(tklabel(wnd, text="Password:"));
  tkgrid(tkentry(wnd, textvariable=passVar, show="*") -> passBox);
  
  tkbind(passBox, "<Return>",function()
    tkdestroy(wnd));
  
  tkwait.window(wnd);
  password <- tclvalue(passVar);
  userName <- tclvalue(user);
  return(c(userName, password));
}

login <- getLogin('')
channel <- odbcConnect('PR_SAIL',login[1],login[2])
rm(login)

## pull extra data - sodium, urea, uric acid, albumin ----

xmn <- sqlTables(channel)
xmn_us <- xmn[ str_detect(xmn$TABLE_SCHEM, "0830"), ]
rm(xmn)
xmn_us

included <- cohort_trim %>%
  select(ALF_PE, closest_date)

alb = sqlQuery(channel, 
               "SELECT P.ALF_PE, EVENT_DT as alb_date, EVENT_VAL as alb
                      FROM SAIL1214V.WLGP_GP_EVENT_CLEANSED_20201001 E
                      LEFT JOIN SAIL1214V.WLGP_PATIENT_ALF_CLEANSED_20201001 P
                      ON E.PRAC_CD_PE = P.PRAC_CD_PE
                      AND E.LOCAL_NUM_PE = P.LOCAL_NUM_PE
                      WHERE EVENT_CD IN ('44M4.', '44MI.');")

alb2 <- alb %>%
  right_join(included) %>%
  filter(!is.na(closest_date)) %>%
  mutate(timediff = (closest_date - ALB_DATE)) %>%
  filter(timediff >= 0) %>%
  group_by(ALF_PE) %>%
  slice_min(timediff, n=1) %>%
  ungroup() %>%
  distinct(ALF_PE, .keep_all = TRUE)

urate = sqlQuery(channel, 
               "SELECT P.ALF_PE, EVENT_DT as urate_date, EVENT_VAL as urate
                      FROM SAIL1214V.WLGP_GP_EVENT_CLEANSED_20201001 E
                      LEFT JOIN SAIL1214V.WLGP_PATIENT_ALF_CLEANSED_20201001 P
                      ON E.PRAC_CD_PE = P.PRAC_CD_PE
                      AND E.LOCAL_NUM_PE = P.LOCAL_NUM_PE
                      WHERE EVENT_CD IN ('44K..', '44K5.', '44K6.');")

urate2 <- urate %>%
  right_join(included) %>%
  filter(!is.na(closest_date)) %>%
  mutate(timediff = (closest_date - URATE_DATE)) %>%
  filter(timediff >= 0) %>%
  group_by(ALF_PE) %>%
  slice_min(timediff, n=1) %>%
  ungroup() %>% 
  distinct(ALF_PE, .keep_all = TRUE)

sod = sqlQuery(channel, 
               "SELECT P.ALF_PE, EVENT_DT as sod_date, EVENT_VAL as sod
                      FROM SAIL1214V.WLGP_GP_EVENT_CLEANSED_20201001 E
                      LEFT JOIN SAIL1214V.WLGP_PATIENT_ALF_CLEANSED_20201001 P
                      ON E.PRAC_CD_PE = P.PRAC_CD_PE
                      AND E.LOCAL_NUM_PE = P.LOCAL_NUM_PE
                      WHERE EVENT_CD IN ('44h1.', '44h6.', '44I5.', '4Q43.');")

sod2 <- sod %>%
  right_join(included) %>%
  filter(!is.na(closest_date)) %>%
  mutate(timediff = (closest_date - SOD_DATE)) %>%
  filter(timediff >= 0) %>%
  group_by(ALF_PE) %>%
  slice_min(timediff, n=1) %>%
  ungroup() %>%
  distinct(ALF_PE, .keep_all = TRUE)

urea = sqlQuery(channel, 
              "SELECT P.ALF_PE, EVENT_DT as urea_date, EVENT_VAL as urea
                      FROM SAIL1214V.WLGP_GP_EVENT_CLEANSED_20201001 E
                      LEFT JOIN SAIL1214V.WLGP_PATIENT_ALF_CLEANSED_20201001 P
                      ON E.PRAC_CD_PE = P.PRAC_CD_PE
                      AND E.LOCAL_NUM_PE = P.LOCAL_NUM_PE
                      WHERE EVENT_CD IN ('44120', '44121', '44J1.', '44J2.', '44JZ.', '44J..', '44J8.', '44J9.', '44JA.', '44JB.');")
  
urea2 <- urea %>%
  filter(!is.na(UREA) & UREA < 100) %>%
  filter(UREA_DATE > as.Date("2000-01-01")) %>%
  right_join(included) %>%
  mutate(timediff = closest_date - UREA_DATE) %>%
  filter(timediff >= 0) %>%
  group_by(ALF_PE) %>%
  slice_min(timediff, n=1) %>%
  ungroup() %>%
  distinct(ALF_PE, .keep_all = TRUE)

rm(alb, urate, na, urea, included, xmn_us)

saveRDS(alb2, "Scratch_data/alb.Rds")
saveRDS(urea2, "Scratch_data/urea.Rds")
saveRDS(sod2, "Scratch_data/sod.Rds")
saveRDS(urate2, "Scratch_data/urate.Rds")

## bind together ----

cohort_trim %>% group_by(ALF_PE) %>% n_distinct()

cohort_add <- cohort_trim %>%
  left_join(alb2 %>% select(-c(closest_date, timediff)), join_by = "ALF_PE") %>%
  left_join(urate2 %>% select(-c(closest_date, timediff)), join_by = "ALF_PE") %>%
  left_join(urea2 %>% select(-c(closest_date, timediff)), join_by = "ALF_PE") %>%
  left_join(sod2 %>% select(-c(closest_date, timediff)), join_by = "ALF_PE") %>%
  group_by(ALF_PE) %>%
  slice_head(n=1)

cohort_add <- cohort_add %>%
  mutate(high_gfr = if_else(eGFR >100, "eGFR >=100", "eGFR <100")) %>%
  mutate(sex_name = if_else(sex == 0, "Female", "Male"),
         stage_name = factor(case_when(
           stage_simplified == 1 | stage_simplified == 2 ~ "Not advanced", 
           stage_simplified == 3 | stage_simplified == 4 ~ "Advanced",
           TRUE ~ "Unknown"),
           levels = c("Not advanced", "Advanced", "Unknown"))) %>%
  mutate(CREAT = sCr / 0.01131222)

saveRDS(cohort_add, "Scratch_data/cohort_add.Rds")

# data tables ----

## what data are available?----

missing_data <- cohort_add %>%
  group_by(high_gfr, sex_name) %>%
  summarise(N = n(), 
            urea = sum(is.na(UREA))/n()*100, 
            sod = sum(is.na(SOD))/n()*100, 
            alb = sum(is.na(ALB))/n()*100, 
            urate = sum(is.na(URATE)/n()*100))

write_csv(missing_data, "Outputs/missing_biochemistry_data.csv")

## is there a significant difference between eGFR >=100 and eGFR <100 ? ----

base <- c("age", "UREA", "CREAT", "ALB", "SOD")

table_high_gfr_adv <- CreateTableOne(vars = base, 
                         #factorVars = cat_base, 
                         strata = c("high_gfr"), 
                         includeNA = TRUE,
                         addOverall = TRUE,
                         data = cohort_add %>% filter(cancer_stage_adv == 1))  
print(table_high_gfr_adv, 
      nonnormal = c("age", "UREA", "CREAT", "ALB", "SOD"))


write.csv(
  print(table_high_gfr_adv),
        "Outputs/table_high_gfr_adv.csv")

table_high_gfr_not_adv <- CreateTableOne(vars = base, 
                                     #factorVars = cat_base, 
                                     strata = c("high_gfr"), 
                                     includeNA = TRUE,
                                     addOverall = TRUE,
                                     data = cohort_add %>% filter(cancer_stage_adv == 0))  
print(table_high_gfr_not_adv, 
      nonnormal = c("age", "UREA", "CREAT", "ALB", "SOD"))
write.csv(
  print(table_high_gfr_not_adv),
        "Outputs/table_high_gfr_not_adv.csv")

table_high_gfr_unknown <- CreateTableOne(vars = base, 
                                         factorVars = cat_base, 
                                         strata = c("high_gfr"), 
                                         includeNA = TRUE,
                                         addOverall = TRUE,
                                         data = cohort_add %>% filter(is.na(cancer_stage_adv)))  
table_high_gfr_unknown
write.csv(
  print(table_high_gfr_unknown),
        "Outputs/table_high_gfr_unknown.csv")

## draw density plots for age, urea, creatinine, sodium, albumin ----

age_plot <- cohort_add %>%
  ggplot(aes(x=age, colour = high_gfr)) +
  geom_density() +
  facet_grid(sex_name ~ stage_name) +
  xlim(20,100) +
  labs(x = "Age (years)", y = "Density", tag = "A", colour = "") +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw() +
  theme(text = element_text(size = 14))

urea_plot <- cohort_add %>%
  ggplot(aes(x=UREA, colour = high_gfr)) +
  geom_density() +
  facet_grid(sex_name ~ stage_name) +
  xlim(0,15) +
  labs(x = "Urea (mmol/L)", y = "Density", tag = "B", colour = "") +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw() +
  theme(text = element_text(size = 14))

creat_plot <- cohort_add %>%
  ggplot(aes(x=CREAT, colour = high_gfr)) +
  geom_density() +
  facet_grid(sex_name ~ stage_name) +
  xlim(20,200) +
  labs(x = "Creatinine (micromol/L)", y = "Density", tag = "C", colour = "") +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw() +
  theme(text = element_text(size = 14))

sod_plot <- cohort_add %>%
  ggplot(aes(x=SOD, colour = high_gfr)) +
  geom_density() +
  facet_grid(sex_name ~ stage_name) +
  xlim(130,150) +
  labs(x = "Sodium (mmol/L)", y = "Density", tag = "D", colour = "") +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw() +
  theme(text = element_text(size = 14))

alb_plot <- cohort_add %>%
  ggplot(aes(x=ALB, colour = high_gfr)) +
  geom_density() +
  facet_grid(sex_name ~ stage_name) +
  xlim(20,60) +
  labs(x = "Albumin (g/L)", y = "Density", tag = "E", colour = "") +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw() +
  theme(text = element_text(size = 14))

biochem_plot <- ggarrange(age_plot, urea_plot, creat_plot, sod_plot, alb_plot, 
          ncol = 1,
          common.legend = TRUE)
biochem_plot

ggsave(biochem_plot, 
       file = "Outputs/biochem_plot.tiff", 
       height = 20, 
       width = 16)

## check attentuation of results by adding SIADH and frailty markers ----

## test for interactions with sex

glm(cancer_stage_adv == 1 ~ age + sex * ckd_cat + SMOKING_STATUS + mm_count_no_ckd + WIMD_2011_DECILE + ca_type + SOD + UREA + ALB, family = "binomial", data = cohort_add) %>% 
  tidy() %>%
  mutate(model = "expl_adv_int") %>% 
  write_csv("Outputs/expl_adv_int.csv")

coxph(Surv(ca_surv_time, death) ~ age + sex * ckd_cat + SMOKING_STATUS + mm_count_no_ckd + WIMD_2011_DECILE + ca_type + stage_simplified  + SOD + UREA + ALB, 
        data = cohort_add) %>%
  tidy() %>%
  mutate(model = "expl_surv_int") %>% 
  write_csv("Outputs/expl_surv_int.csv")

## run nested models  

cohort_nst <- cohort_add %>%
  group_by(sex_name) %>%
  nest() %>%
  ungroup()

expl_adv <- map(cohort_nst$data, 
             ~ glm(cancer_stage_adv == 1 ~ age + ckd_cat + SMOKING_STATUS + mm_count_no_ckd + WIMD_2011_DECILE + ca_type
                 + SOD + UREA + ALB, family = "binomial", 
             data = .x)) 

expl_adv <- expl_adv %>% map(tidy) 

expl_adv_m <- expl_adv[[1]] %>% mutate(or = exp(estimate),
                              lower = exp(estimate - 1.96*std.error),
                              upper = exp(estimate + 1.96*std.error)) %>%
  mutate(sex_name = "Male")

expl_adv_f <- expl_adv[[2]] %>% mutate(or = exp(estimate),
                                       lower = exp(estimate - 1.96*std.error),
                                       upper = exp(estimate + 1.96*std.error)) %>%
  mutate(sex_name = "Female")

expl_adv <- bind_rows(expl_adv_m, expl_adv_f)
expl_adv <- expl_adv %>%  
  filter(str_detect(term, "ckd_cateGFR")) %>%
  separate(term, into = c("term", "ckd_cat"), sep = 7) %>%
  separate(ckd_cat, into = c("egfr", "ckd_cat"), sep = 5) 
expl_adv <- expl_adv %>%
  add_row(tibble(ckd_cat = "75 - <90", sex_name = "Male", or = 1, lower = 1, upper = 1)) %>%
  add_row(tibble(ckd_cat = "75 - <90", sex_name = "Female", or = 1, lower = 1, upper = 1))
expl_adv$ckd_cat <- factor(expl_adv$ckd_cat, 
                           levels = c("<30", "30 - <45", "45 - <60", "60 - <75", "75 - <90", "90 - <105", "105 - <120", ">120"))


expl_advplot <- expl_adv %>%
  ggplot(aes(x = ckd_cat, y = or, ymin - lower, ymax = upper, colour = sex_name)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  scale_y_continuous(trans = "log10", limits = c(0.5, 5)) +
  geom_hline(yintercept = 1, lty = 2) + 
  geom_errorbar(aes(ymin = if_else(lower <= 0.5, 0.5, lower), ymax = if_else(upper >= 5, 5, upper)), linewidth = 0.8, 
                width = 0.3, cex = 0.3, pch = 20, position = position_dodge(width = 0.4)) + 
  labs(x = "eGFRcr prior to cancer diagnosis", y = "Adjusted Odds Ratio (95% confidence interval)", title = "Advanced cancer at presentation", colour = "Sex") +
  theme_bw() + 
  theme(text = element_text(size = 20))


## survival

expl_surv <- map(cohort_nst$data, 
                 ~ coxph(Surv(ca_surv_time, death) ~ age + ckd_cat + SMOKING_STATUS + mm_count_no_ckd + WIMD_2011_DECILE + ca_type + stage_simplified  + SOD + UREA + ALB, 
               data = .x)) 

expl_surv <- expl_surv %>% map(tidy) 

expl_surv_m <- expl_surv[[1]] %>% mutate(or = exp(estimate),
                                       lower = exp(estimate - 1.96*std.error),
                                       upper = exp(estimate + 1.96*std.error)) %>%
  mutate(sex_name = "Male")

expl_surv_f <- expl_surv[[2]] %>% mutate(or = exp(estimate),
                                       lower = exp(estimate - 1.96*std.error),
                                       upper = exp(estimate + 1.96*std.error)) %>%
  mutate(sex_name = "Female")

expl_surv <- bind_rows(expl_surv_m, expl_surv_f)
expl_surv <- expl_surv %>%  
  filter(str_detect(term, "ckd_cateGFR")) %>%
  separate(term, into = c("term", "ckd_cat"), sep = 7) %>%
  separate(ckd_cat, into = c("egfr", "ckd_cat"), sep = 5) 
expl_surv <- expl_surv %>%
  add_row(tibble(ckd_cat = "75 - <90", sex_name = "Male", or = 1, lower = 1, upper = 1)) %>%
  add_row(tibble(ckd_cat = "75 - <90", sex_name = "Female", or = 1, lower = 1, upper = 1))
expl_surv$ckd_cat <- factor(expl_surv$ckd_cat, 
                            levels = c("<30", "30 - <45", "45 - <60", "60 - <75", "75 - <90", "90 - <105", "105 - <120", ">120"))

expl_survplot <- expl_surv %>%
  ggplot(aes(x = ckd_cat, y = or, ymin - lower, ymax = upper, colour = sex_name)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  scale_y_continuous(trans = "log10", limits = c(0.5, 5)) +
  geom_hline(yintercept = 1, lty = 2) + 
  geom_errorbar(aes(ymin = if_else(lower <= 0.5, 0.5, lower), ymax = if_else(upper >= 5, 5, upper)), linewidth = 0.8, 
                width = 0.3, cex = 0.3, pch = 20, position = position_dodge(width = 0.4)) +
  labs(x = "eGFRcr prior to cancer diagnosis", y = "Adjusted Hazard Ratio (95% confidence interval)", title = "All-cause mortality", colour = "Sex") +
  theme_bw() + 
  theme(text = element_text(size = 20)) + 
  annotate("text", x = "<30", y = 0.5, label = "*", size = 15) + 
  annotate("text", x = "30 - <45", y = 0.5, label = "*", size = 15) 

## save model outputs

write_csv(expl_adv, "Outputs/expl_adv.csv")
write_csv(expl_surv, "Outputs/expl_surv.csv")
