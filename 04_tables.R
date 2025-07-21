tab_vars <- c("sex_f31_0_0", "prevalent_diabetes_t2", "prevalent_hypertension",)

tabdat <- dat_metab_sub %>%
  select(eid, all_of(char_key$col.name), metabolite_status)

final3 <- final2 %>%
  right_join(tabdat %>% select(eid, sex_f31_0_0, metabolite_status))

final3[is.na(final3)] <- 0
var_label(final3) <- outcomedatlist

metab_tab_icd10_phase1 <- final3 %>%
  filter(metabolite_status == "Phase 1") %>%
  select(sex_f31_0_0, all_of(names(unlist(outcomedatlist)))) %>%
  tbl_summary(by = sex_f31_0_0, missing = "no") %>%
  add_overall()

metab_tab_icd10_phase2 <- final3 %>%
  filter(metabolite_status == "Phase 2") %>%
  select(sex_f31_0_0, all_of(names(unlist(outcomedatlist)))) %>%
  tbl_summary(by = sex_f31_0_0, missing = "no") %>%
  add_overall()

all_tab_icd10 <- final3 %>%
  select(sex_f31_0_0, all_of(names(unlist(outcomedatlist)))) %>%
  tbl_summary(by = sex_f31_0_0, missing = "no") %>%
  add_overall()

tbl_merge_icd10 <-
  tbl_merge(tbls = list(all_tab_icd10, metab_tab_icd10_phase1, metab_tab_icd10_phase2),
            tab_spanner = c("All Samples", "Phase 1", "Phase 2"))

library(gt)

gtsave(as_gt(tbl_merge_icd10), file = "dx_prevalence_table_080123.rtf")

save(tbl_merge_icd10, file = "dx_prevalence_table_icd10_080123.Rdata")

sum(tabdat$has_metabolites)

dontwant <- c("weight_method_f21_0_0", "treatmentmedication_code_f20003_0_0", "forced_expiratory_volume_in_1second_fev1_f3063_0_0",
              "forced_expiratory_volume_in_1second_fev1_pilot_f10695_0_1",
              "forced_vital_capacity_fvc_f3062_0_0", "forced_vital_capacity_fvc_pilot_f10694_0_1")
deal_later <- c("past_tobacco_smoking_f1249_0_0", 
                "had_menopause_f2724_0_0",
                "ever_used_hormonereplacement_therapy_hrt_f2814_0_0"
                ,                               "ever_had_hysterectomy_womb_removed_f3591_0_0",
                "smokingsmokers_in_household_f1259_0_0", "ethnic_background_f21000_0_0")

var_label(tabdat) <- covardatlist
metab_tab_phase1 <- tabdat %>%
  filter(metabolite_status== "Phase 1") %>%
  #select(sex:had_menopause, pulse_rate_automated_reading_f102_0_0:pulse_rate_during_bloodpressure_measurement_f95_0_0) %>%
  #select(!all_of(dontwant), !all_of(deal_later)) %>%
  select(sex_f31_0_0, all_of(names(unlist(covardatlist)))) %>%
  tbl_summary(by = sex_f31_0_0, missing = "no") %>%
  bold_labels() %>%
  add_n() %>%
  add_overall()

metab_tab_phase2 <- tabdat %>%
  filter(metabolite_status == "Phase 2") %>%
  select(sex_f31_0_0, all_of(names(unlist(covardatlist)))) %>%
  tbl_summary(by = sex_f31_0_0, missing = "no") %>%
  bold_labels() %>%
  add_n() %>%
  add_overall()

all_tab <- tabdat %>%
  #filter(has_metabolites == 1) %>%
  #select(sex:had_menopause, pulse_rate_automated_reading_f102_0_0:pulse_rate_during_bloodpressure_measurement_f95_0_0) %>%
  #select(!all_of(dontwant), !all_of(deal_later)) %>%
  select(sex_f31_0_0, all_of(names(unlist(covardatlist)))) %>%
  tbl_summary(by =sex_f31_0_0, missing = "no") %>%
  bold_labels() %>%
  add_n() %>%
  add_overall()


tbl_merge1 <- 
  tbl_merge(tbls = list(all_tab, metab_tab_phase1, metab_tab_phase2),
            tab_spanner = c("All Samples", "Phase 1", "Phase 2"))# %>%
#modify_spanning_header(all_stat_cols() ~ "**Baseline Measurements**")

library(gt)

gtsave(as_gt(tbl_merge1), file = "baseline_table_081023.rtf")


save(tbl_merge1, file = "testtab2_080123.Rdata")





tab_disease <- dat_sub_metab2 %>%
  select(sex_f31_0_0, metabolite_status, all_of(names(unlist(outcomedatlist))), all_of(names(unlist(outcomedatlist_fo))))

metab_tab_phase1_sex <- tab_disease %>%
  filter(metabolite_status == "Phase 1") %>%
  select(-metabolite_status) %>%
  tbl_summary(by = sex_f31_0_0, missing = "no") %>%
  bold_labels() %>%
  add_n() %>%
  add_overall()

all_tab_phase1_outcome <- tab_disease %>%
  select(-metabolite_status) %>%
  tbl_summary(by = sex_f31_0_0, missing = "no") %>%
  bold_labels() %>%
  add_n() %>%
  add_overall()