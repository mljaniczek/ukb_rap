# data variable derivations and labeling

# change in ukb

library(dplyr)
library(tidyr)
library(forcats)
library(gtsummary)
library(ukbtools)
library(stringr)

load("data_pheno.Rdata")

load("data_metabolites.Rdata")

load("analysis_variable_keys.Rdata")


med_dat <- dat %>%
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

dat <- dat %>%
  mutate(SBP = ifelse(!is.na(systolic_blood_pressure_automated_reading_f4080_0_0), 
                      systolic_blood_pressure_automated_reading_f4080_0_0, 
                      systolic_blood_pressure_manual_reading_f93_0_0))

dat <- dat %>%
  #filter(metabolite_status == "Phase 1") %>%
  # make age category
  mutate(age_cat = ifelse(
    age_when_attended_assessment_centre_f21003_0_0 <50, "<50",
    ifelse(age_when_attended_assessment_centre_f21003_0_0 >=50 & 
             age_when_attended_assessment_centre_f21003_0_0 <60, "50-59",
           ifelse(age_when_attended_assessment_centre_f21003_0_0 >=60, "60+", NA))),
    menopause_cat = ifelse(
      age_when_attended_assessment_centre_f21003_0_0 <55 & 
        had_menopause_f2724_0_0 %in% c("Prefer not to answer", "Not sure - other reason", "Not sure - had a hysterectomy", NA), "Undetermined",
      ifelse(age_when_attended_assessment_centre_f21003_0_0 >=55 & 
               had_menopause_f2724_0_0 %in% c("Prefer not to answer", "Not sure - other reason", "Not sure - had a hysterectomy", "Yes", NA), "Yes",
             ifelse(sex_f31_0_0 == "Male" | had_menopause_f2724_0_0 == "No", "No",
                    ifelse(had_menopause_f2724_0_0 == "Yes", "Yes",
                           had_menopause_f2724_0_0)))))

# now recode race

dat <- dat %>%
  mutate(
    race = ifelse(
      ethnic_background_f21000_0_0 %in% c("White", "British", "Irish", "Any other white background"), "White",
      ifelse(
        ethnic_background_f21000_0_0 %in% c("Black or Black British", "Caribbean", "African", "Any other Black background"), "Black",
        ifelse(ethnic_background_f21000_0_0 %in% c("Asian or Asian British", "Chinese", "Any other Asian background"), "Asian",
               ifelse(ethnic_background_f21000_0_0 %in% c("Indian", "Pakistani", "Bangladeshi"), "Southeast Asian",
                      ifelse(ethnic_background_f21000_0_0 %in% c("Mixed", "White and Black Caribbean", "White and Black African", "White and Asian", "Any other mixed background"), "Mixed",
                             ifelse((ethnic_background_f21000_0_0 %in% c("Prefer not to answer", "Do not know"))|is.na(ethnic_background_f21000_0_0), "Unknown", 
                                    ifelse(
                                      ethnic_background_f21000_0_0 %in% c("Other ethnic group"), "Other", NA)))))
      )
    )
  )

dat <- dat %>%
  mutate(race = forcats::fct_relevel(as.factor(race), "White", "Black",
                                     "Southeast Asian", "Asian", "Other", "Mixed", "Unknown"))

dat <- dat %>%
  mutate(bmi = as.numeric(scale(body_mass_index_bmi_f21001_0_0, center = FALSE))) 


# now pare down any vars we don't need 


dat_metab <- dat %>% 
  filter(has_metabolites_phase2 == TRUE) 

metabs_sub <- metabs %>%
  filter(eid %in% dat_metab$eid)


dat_metab %>%
  dplyr::select(-eid) %>%
  tbl_summary(by = metabolite_status)


dat_metab_sub <- dat_metab %>%
  filter(withdraw == FALSE) %>%
  mutate(
    SMOKING = case_when(
      (current_tobacco_smoking_f1239_0_0 == "No" | is.na(current_tobacco_smoking_f1239_0_0)) 
      & (past_tobacco_smoking_f1249_0_0 %in% c("Smoked on most or all days", "Smoked occasionally")) ~ "Previously smoked", 
      is.na(current_tobacco_smoking_f1239_0_0) ~ "Unknown",
      .default = current_tobacco_smoking_f1239_0_0
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
    menopause_cat = case_when(sex_f31_0_0 == "Male" ~ "No",
                              .default = menopause_cat)
  ) %>%
  ungroup() %>%
  mutate(
    age_cat = fct_relevel(as.factor(age_cat), "<50", "50-59", "60+")
  )


save(dat_metab_sub, file = "processed_pheno.Rdata")

tabdat <- dat_metab_sub %>%
  filter(metabolite_status == "Phase 1") %>%
  dplyr::select(c(age_cat, age_when_attended_assessment_centre_f21003_0_0, 
                  race, sex_f31_0_0, 
                   BMI, WHR, SBP, HDL, TC, HBA1C, SMOKING,
                  `Blood pressure medication`:Insulin,
                  prevalent_hypertension, prevalent_cad, prevalent_diabetes_t2)) %>%
  group_by(age_cat) %>%
  nest() %>%
  arrange(age_cat) %>%
  mutate(
    demo_tab = map(data, ~tbl_summary(., by = sex_f31_0_0))
  )

tbl_merge(tabdat$demo_tab, tab_spanner = as.character(tabdat$age_cat))

tabdat2 <- dat_metab_sub %>%
  filter(metabolite_status == "Phase 2") %>%
  dplyr::select(c(age_cat, age_when_attended_assessment_centre_f21003_0_0, 
                  race, sex_f31_0_0, 
                  BMI, WHR, SBP, HDL, TC, HBA1C, SMOKING,
                  `Blood pressure medication`:Insulin,
                  prevalent_hypertension, prevalent_cad, prevalent_diabetes_t2)) %>%
  group_by(age_cat) %>%
  nest() %>%
  arrange(age_cat) %>%
  mutate(
    demo_tab = map(data, ~tbl_summary(., by = sex_f31_0_0))
  )
  
tbl_merge(tabdat2$demo_tab, tab_spanner = as.character(tabdat2$age_cat))

var_label(dat_sub_metab2) <- outcomedatlist_fo
var_label(dat_sub_metab2) <- outcomedatlist
var_label(dat_sub_metab2) <- covardatlist
var_label(dat_sub_metab2) <- outcomedatlist_inc
var_label(dat_sub_metab2) <- outcomedatlist_prev