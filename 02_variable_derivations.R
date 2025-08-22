# data variable derivations and labeling

# change in ukb

library(dplyr)
library(tidyr)
library(forcats)
library(gtsummary)
library(ukbtools)
library(stringr)

source("labellist.R")

system('dx download pheno.csv') # downloaded phenotype data
system('dx download mets.csv')
system("dx download first_occurrences_processed.csv")

pheno <- read_csv("pheno.csv", name_repair = "universal_quiet")

mets <- read_csv("mets.csv", name_repair = "universal_quiet")

load("analysis_variable_keys.Rdata")


med_dat <- pheno %>%
  dplyr::select(eid, all_of(med_key$col.name))

med_dat2 <- med_dat %>%
  mutate(instance1 = str_replace( paste0(medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_f6153_0_0, medication_for_cholesterol_blood_pressure_or_diabetes_f6177_0_0), "NA", ""), 
         instance2 = str_replace( paste0(medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_f6153_0_1, medication_for_cholesterol_blood_pressure_or_diabetes_f6177_0_1), "NA", ""),
         instance3 = str_replace(paste0(medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_f6153_0_2, medication_for_cholesterol_blood_pressure_or_diabetes_f6177_0_2), "NA", ""),
         instance4 = paste(medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_f6153_0_3)) %>%
  dplyr::select(-c(medication_for_cholesterol_blood_pressure_diabetes_or_take_exogenous_hormones_f6153_0_0:medication_for_cholesterol_blood_pressure_or_diabetes_f6177_0_2))


med_dat3 <- med_dat2 %>%
  pivot_longer(-eid) %>%
  #filter(value != "NA") %>%
  dplyr::select(-name) %>%
  unique() %>%
  mutate(has = 1) %>%
  pivot_wider(names_from = value, values_from = has, values_fill = 0) %>%
  dplyr::select(-c("NA", "None of the above", "Do not know", "Prefer not to answer"))



dat <- dat %>%
  dplyr::select(-c(all_of(med_key$col.name))) %>%
  left_join(med_dat3)

dat <- pheno %>%
  mutate(SBP = ifelse(!is.na(Systolic.blood.pressure..automated.reading...Instance.0...Array.0), 
                      Systolic.blood.pressure..automated.reading...Instance.0...Array.0, 
                      Systolic.blood.pressure..manual.reading...Instance.0...Array.0))

dat <- dat %>%
  #filter(metabolite_status == "Phase 1") %>%
  # make age category
  mutate(age_cat = ifelse(
    Age.when.attended.assessment.centre...Instance.0 <50, "<50",
    ifelse(Age.when.attended.assessment.centre...Instance.0 >=50 & 
             Age.when.attended.assessment.centre...Instance.0 <60, "50-59",
           ifelse(Age.when.attended.assessment.centre...Instance.0 >=60, "60+", NA))),
    menopause_cat = ifelse(
      Age.when.attended.assessment.centre...Instance.0 <55 & 
        Had.menopause...Instance.0 %in% c("Prefer not to answer", "Not sure - other reason", "Not sure - had a hysterectomy", NA), "Undetermined",
      ifelse(Age.when.attended.assessment.centre...Instance.0 >=55 & 
               Had.menopause...Instance.0 %in% c("Prefer not to answer", "Not sure - other reason", "Not sure - had a hysterectomy", "Yes", NA), "Yes",
             ifelse(Sex == "Male" | Had.menopause...Instance.0 == "No", "No",
                    ifelse(Had.menopause...Instance.0 == "Yes", "Yes",
                           Had.menopause...Instance.0)))))

# now recode race

dat <- dat %>%
  mutate(
    race = ifelse(
      Ethnic.background...Instance.0 %in% c("White", "British", "Irish", "Any other white background"), "White",
      ifelse(
        Ethnic.background...Instance.0 %in% c("Black or Black British", "Caribbean", "African", "Any other Black background"), "Black",
        ifelse(Ethnic.background...Instance.0 %in% c("Asian or Asian British", "Chinese", "Any other Asian background"), "Asian",
               ifelse(Ethnic.background...Instance.0 %in% c("Indian", "Pakistani", "Bangladeshi"), "Southeast Asian",
                      ifelse(Ethnic.background...Instance.0 %in% c("Mixed", "White and Black Caribbean", "White and Black African", "White and Asian", "Any other mixed background"), "Mixed",
                             ifelse((Ethnic.background...Instance.0 %in% c("Prefer not to answer", "Do not know"))|is.na(Ethnic.background...Instance.0), "Unknown", 
                                    ifelse(
                                      Ethnic.background...Instance.0 %in% c("Other ethnic group"), "Other", NA)))))
      )
    )
  )

dat <- dat %>%
  mutate(race = forcats::fct_relevel(as.factor(race), "White", "Black",
                                     "Southeast Asian", "Asian", "Other", "Mixed", "Unknown"))

dat <- dat %>%
  mutate(bmi = as.numeric(scale(Body.mass.index..BMI....Instance.0, center = FALSE))) 


# now pare down any vars we don't need 

dat_metab_sub <- dat %>%
  #filter(withdraw == FALSE) %>% 
  mutate(
    metabolite_status = ifelse(
      Sample.Measured.Date.and.Time...Instance.0 < "2020-04-16", 
      "Phase 1", "Phase 2"
    ),
    SMOKING = case_when(
      (Current.tobacco.smoking...Instance.0 == "No" | is.na(Current.tobacco.smoking...Instance.0)) 
      & (past_tobacco_smoking_f1249_0_0 %in% c("Smoked on most or all days", "Smoked occasionally")) ~ "Previously smoked", 
      is.na(Current.tobacco.smoking...Instance.0) ~ "Unknown",
      .default = Current.tobacco.smoking...Instance.0
    )) %>%
  group_by(metabolite_status, age_cat) %>%
  mutate(
    ##TODO need to do this imputation within age strata and phase separately
    BMI = case_when(
      is.na(bmi) ~ mean(bmi, na.rm = TRUE),
      .default = bmi
    ),
    SBP = case_when(
      is.na(SBP) ~ mean(SBP, na.rm = TRUE),
      .default = SBP
    ),
    whr = waist_circumference_f48_0_0/hip_circumference_f49_0_0,
    WHR = case_when(
      is.na(whr) ~ mean(whr, na.rm = TRUE),
      .default = whr
    ),
    HDL = case_when(
      is.na(hdl_cholesterol_f30760_0_0) ~ mean(hdl_cholesterol_f30760_0_0, na.rm = TRUE),
      .default = hdl_cholesterol_f30760_0_0
    ),
    HBA1C = case_when(
      is.na(glycated_haemoglobin_hba1c_f30750_0_0) ~ mean(glycated_haemoglobin_hba1c_f30750_0_0, na.rm = TRUE),
      .default = glycated_haemoglobin_hba1c_f30750_0_0
    ),
    TC = case_when(
      is.na(cholesterol_f30690_0_0) ~ mean(cholesterol_f30690_0_0, na.rm = TRUE),
      .default = cholesterol_f30690_0_0
    ),
    menopause_cat = case_when(Sex == "Male" ~ "No",
                              .default = menopause_cat)
  ) %>%
  ungroup() %>%
  mutate(
    age_cat = fct_relevel(as.factor(age_cat), "<50", "50-59", "60+"),
    age = Age.when.attended.assessment.centre...Instance.0
  ) 

dat_metab_sub2 <- dat_metab_sub %>%
  left_join(final_fo %>% select(Participant.ID, prevalent_cad, prevalent_diabetes_t2, 
                                prevalent_hypertension, incident_cad))

var_label(dat_sub_metab2) <- outcomedatlist_fo
var_label(dat_sub_metab2) <- outcomedatlist
var_label(dat_sub_metab2) <- covardatlist
var_label(dat_sub_metab2) <- outcomedatlist_inc
var_label(dat_sub_metab2) <- outcomedatlist_prev

### CHECK dat_metabsub2 and see if dimensions are right 
save(dat_metab_sub, file = "processed_pheno.Rdata")

tabdat <- dat_metab_sub %>%
  filter(metabolite_status == "Phase 1") %>%
  dplyr::select(c(age_cat, Age.when.attended.assessment.centre...Instance.0, 
                  race, Sex, 
                   BMI, WHR, SBP, HDL, TC, HBA1C, SMOKING,
                  `Blood pressure medication`:Insulin,
                  prevalent_hypertension, prevalent_cad, prevalent_diabetes_t2)) %>%
  group_by(age_cat) %>%
  nest() %>%
  arrange(age_cat) %>%
  mutate(
    demo_tab = map(data, ~tbl_summary(., by = Sex))
  )

tbl_merge(tabdat$demo_tab, tab_spanner = as.character(tabdat$age_cat))

tabdat2 <- dat_metab_sub %>%
  filter(metabolite_status == "Phase 2") %>%
  dplyr::select(c(age_cat, Age.when.attended.assessment.centre...Instance.0, 
                  race, Sex, 
                  BMI, WHR, SBP, HDL, TC, HBA1C, SMOKING,
                  `Blood pressure medication`:Insulin,
                  prevalent_hypertension, prevalent_cad, prevalent_diabetes_t2)) %>%
  group_by(age_cat) %>%
  nest() %>%
  arrange(age_cat) %>%
  mutate(
    demo_tab = map(data, ~tbl_summary(., by = Sex))
  )
  
tbl_merge(tabdat2$demo_tab, tab_spanner = as.character(tabdat2$age_cat))

