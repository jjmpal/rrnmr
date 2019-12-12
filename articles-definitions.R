FR_definitions <- function(dset, variables) {
    dset <- dplyr::mutate(dset, age = case_when(VUOSI == 2017 ~ IKA_POIMINTA,
                                        TRUE ~ IKA),
                  female = as.factor(case_when(SUKUP == 2 ~ 1,
                                               SUKUP == 1 ~ 0)),
                  smoker = as.factor(case_when(TUPI3 == 4 ~ 1,
                                               TUPI3 %in% c(1, 2, 3) ~ 0)),
                  diabetes = as.factor(case_when(Q34 == 3 ~ 1,
                                                 Q34 %in% c(1, 2) ~ 0,
                                                 FR07_38 %in% c(1, 2, 6) ~ 0,
                                                 FR07_38 %in% c(3, 4, 5) ~ 1,
                                                 FR12_20A == 1 | FR12_20B == 1 | FR12_20F == 1 ~ 0,
                                                 FR12_20C == 1 | FR12_20D == 1 | FR12_20E == 1 ~ 1,
                                                 FT17_L1_7_1 == 1 | FT17_L1_7_2 == 1 | FT17_L1_7_6 == 1 ~ 0,
                                                 FT17_L1_7_3 == 1 | FT17_L1_7_4 == 1 | FT17_L1_7_5 == 1 ~ 1)),
                  bmi = case_when(VUOSI == 2017 ~ FT17_TT2_12_BMI,
                                  TRUE ~ BMI),
                  leisure = as.factor(case_when(VUOSI == 2017 ~ FT17_L1_55,
                                                TRUE ~ Q57)),
                  BP_EVERHIGH = case_when(VUOSI == 2017 ~ FT17_L1_15,
                                          TRUE ~ Q27),
                  BP_TREATED = case_when(VUOSI == 2017 ~ FT17_L1_16,
                                         TRUE ~ Q29),
                  BP_TREATLAST = case_when(VUOSI == 2017 ~ FT17_L1_17,
                                           TRUE ~ Q30),
                  bptreat = as.factor(case_when(BP_EVERHIGH == 1 ~ 0,
                                                BP_TREATED == 1 ~ 0,
                                                BP_TREATLAST %in% seq(2,6) ~ 0,
                                                BP_TREATLAST == 1 ~ 1)),
                  KOL_MEASURED = case_when(VUOSI == 2017 ~ FT17_L1_10,
                                       TRUE ~ Q23),
                  KOL_TREATED  = case_when(VUOSI == 2017 ~ FT17_L1_13,
                                           TRUE ~ K34),
                  koltreat = as.factor(case_when(KOL_MEASURED %in% c(5,6) ~ 0,
                                                 KOL_TREATED == 2 ~ 1,
                                                 KOL_TREATED == 1 ~ 0)),
                  sys = case_when(VUOSI == 2017 ~ (FT17_TT1_SYS1 + FT17_TT1_SYS2)/2,
                                  TRUE ~ (SYS1 + SYS2)/2),
                  dias = case_when(VUOSI == 2017 ~ (FT17_TT1_DIAS1 + FT17_TT1_DIAS2)/2,
                                   TRUE ~ (DIAS1 + DIAS2)/2),
                  htn = factor(ifelse(sys >= 140 | dias >= 90 | bptreat == 1, 1, 0)),
                  htn3 = factor(case_when(bptreat == 1 | sys >= 160 | dias >= 100 ~ "Stage 2",
                                          sys >= 140 | dias >= 90 ~ "Stage 1",
                                          sys <140 & dias < 90 ~ "Normotensive")),
                  cohort = as.character(VUOSI),
                  sampleid = sprintf("%s_%s", cohort, HAVTUN)) %>%
        dplyr::select(variables)
}

DILGOM_definitions <- function(dset.full, variables) {
    stopifnot(!missing(dset.full))
    filter(dset.full, VUOSI == 2007) %>%
        dplyr::mutate(age = IKA + 5,
                      female = as.factor(case_when(SUKUP == 2 ~ 1,
                                                   SUKUP == 1 ~ 0)),
                      smoker = as.factor(case_when(Q67_2014 == 1 ~ 0,
                                                   Q69_2014 %in% c(2, 3) ~ 0,
                                                   Q69_2014 == 1 ~ 1)),
                      diabetes = as.factor(case_when(DLGM14_K18_1 == 1 |
                                                     DLGM14_K18_2 == 1 |
                                                     DLGM14_K18_6 == 1 ~ 0,
                                                     DLGM14_K18_3 == 1 |
                                                     DLGM14_K18_4 == 1 |
                                                     DLGM14_K18_5 == 1 ~ 1)),
                      bmi = B_BMI_BIOL_2014,
                      leisure = as.factor(Q57_2014),
                      bptreat = as.factor(case_when(Q27_2014 == 1 ~ 0,
                                                    Q29_2014 == 1 ~ 0,
                                                    Q30_2014 %in% seq(2,6) ~ 0,
                                                    Q30_2014 == 1 ~ 1)),
                      koltreat = as.factor(case_when(Q24_2014 == 1 ~ 0,
                                                     K34_2014 == 2 ~ 1,
                                                     K34_2014 == 1 ~ 0)),
                      sys = (DLGM14_VP_OMRON1_SYS + DLGM14_VP_OMRON2_SYS)/2,
                      dias = (DLGM14_VP_OMRON1_DIAS + DLGM14_VP_OMRON2_DIAS)/2,
                      htn = factor(ifelse(sys >= 140 | dias >= 90 | bptreat == 1, 1, 0)),
                      htn3 = factor(case_when(bptreat == 1 | sys >= 160 | dias >= 100 ~ "Stage 2",
                                              sys >= 140 | dias >= 90 ~ "Stage 1",
                                              TRUE ~ "Normotensive")),
                      VUOSI = 2014,
                      cohort = as.character(VUOSI),
                      sampleid = sprintf("%s_%s", cohort, HAVTUN)) %>%
        dplyr::select(variables)
}

HEALTH_definitions <- function(pheno, phenoextra, metabolite, variables, cohort) {
    stopifnot(!missing(cohort))
    dset.pheno <- read_sas(pheno) %>%
        left_join(., read_sas(phenoextra) %>% dplyr::select(-SP2), by = "_RANDOM_ID_") %>%
        dplyr::mutate(age = IKA2,
               female = as.factor(case_when(SP2 == 2 ~ 1,
                                            SP2 == 1 ~ 0)),
               smoker = as.factor(case_when(FB01 == 0 ~ 0,
                                            FB05 == 1 ~ 1,
                                            FB05 == 2 ~ 0,
                                            FB05 == 3 ~ 0)),
               diabetes = as.factor(case_when(BA26 == 1 ~ 1,
                                              BA26 == 0 ~ 0)),
               bmi = BMII_BMI,
               leisure = as.factor(KYS1_K27),
               bptreat = as.factor(case_when(BA13 == 0 ~ 0,
                                             BA13C == 1 ~ 1,
                                             BA13C == 0 ~ 0)),
               koltreat = as.factor(ATC_C10A),
               sys = (MIT1_SYSTBP1 + MIT1_SYSTBP2)/2,
               dias = (MIT1_DIASTBP1 + MIT1_DIASTBP2)/2,
               htn = factor(ifelse(sys >= 140 | dias >= 90 | bptreat == 1, 1, 0)),
               htn3 = factor(case_when(bptreat == 1 | sys >= 160 | dias >= 100 ~ "Stage 2",
                                       sys >= 140 | dias >= 90 ~ "Stage 1",
                                       sys <140 & dias < 90 ~ "Normotensive")),
               cohort = as.character(cohort),
               sampleid = dropattributes(`_RANDOM_ID_`)) %>%
        dplyr::select(variables, -`_RANDOM_ID_`) 

    dset.nmr <- my_read_nmr(metabolite) %>%
        rename(sampleid = `_RANDOM_ID_`) %>%
      
    inner_join(dset.pheno, dset.nmr, by = c("sampleid" = "sampleid"))
}

