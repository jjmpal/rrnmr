my_read_sas <- function(file) {
    read_sas(file) %>%
        as_tibble %>%
        rename_all(toupper)
}

my_read_nmr <- function(file, exclude = c("GLUCONLAC",
                                          "GLUCONLAC2",
                                          "HAVTUN5",
                                          "HIETOH",
                                          "HIGHLAC",
                                          "HIGHLAC2",
                                          "HIGHPYR",
                                          "HIGHPYR2",
                                          "ISOPROPANOL",
                                          "LOGLTMINHIGLTMAT",
                                          "LOWGLUC",
                                          "LOWGLUC2",
                                          "LOWPROT",
                                          "HIGH_PYRUVATE",
                                          "ISOPROPYL_ALCOHOL",
                                          "LOW_PROTEIN_CONTENT",
                                          "HIGH_ETANOL")) {
    my_read_sas(file) %>%
        rename_all(toupper) %>%
        rename_all(~metabolitemapping(.)) %>%
        mutate_at(vars(starts_with("NMR_")), ~tryCatch(suppressWarnings(as.numeric(.)))) %>%
        mutate_at(vars(starts_with("NMR_")), ~ifelse(. > 0, ., NA)) %>%
        select(., -one_of(exclude[exclude %in% colnames(.)]))
}


remove_reexamined <- function(dset, levels = c(2017, 2012, 2007, 2002, 1997)) {
    dplyr::group_by(dset, FR_ID) %>%
        dplyr::arrange(factor(VUOSI, levels = levels)) %>%
        dplyr::mutate(remove = case_when(VUOSI == 2017 ~ 0,
                                         row_number() == 1 ~ 0,
                                         TRUE ~ 1)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(remove != 1) %>%
        dplyr::select(-remove)
}

mydescribe <- function(df) {
    tryCatch(suppressWarnings(select(df, -one_of("sampleid", "cohort", "htn3")))) %>%
        gather(key, value) %>%
        group_by(key) %>%
        summarise(Mean = mean(as.numeric(value), na.rm = TRUE),
                           SD = sd(as.numeric(value), na.rm = TRUE),
                           N = n(),
                           NAs = sum(is.na(value)),
                           ratio = sum(is.na(value))/n())
}

listna <- function(df) {
    lapply(df, function(x) is.na(x)) %>%
        Reduce("|", .)
}

missingbycohort <- function(dset, var.fr) {
    var.nmr <- colnames(dset) %>% mygrep(word = "NMR_")
    dset %>%
        dplyr::mutate(NA.fr = listna(.[var.fr]),
               NA.nmr = listna(.[var.nmr]),
               NA.all = NA.fr | NA.nmr) %>%
        group_by(cohort) %>%
        summarize(N.raw = n(),
                  missing.meta = sum(NA.fr),
                  missing.nmr = sum(NA.nmr),
                  N = N.raw-sum(NA.all))
}

getnamesbyclass <- function(groups) {
    ggforestplot::df_NG_biomarker_metadata %>%
        mutate(id = paste0("NMR_", machine_readable_name)) %>%
        select(id, name, group) %>%
        filter(group %in% groups) %>%
        pull(id)
}

df2l <- function(df) {
    list <- df %>% pull(machine_readable_name)
    names(list) <- df %>% pull(alt)
    return(list)
}

excludemetabolites <- function(bygroup, byname, byexpression) {
    ret.group <- getnamesbyclass(bygroup)
    ret.expr <- ggforestplot::df_NG_biomarker_metadata %>%
        filter(grepl(byexpression, machine_readable_name)) %>%
        pull(machine_readable_name)
    c(ret.group, paste0("NMR_", c(byname, ret.expr))) %>% unique
}


metabolitemapping <- function(..., prefix = "NMR_") {
    list <- ggforestplot::df_NG_biomarker_metadata %>%
        mutate(alt = purrr::map2(alternative_names,
                                 machine_readable_name,
                                 ~toupper(c(.x,
                                            gsub("[-/]", "_", c(.x,
                                                                gsub("%", "PROS", .x),
                                                                gsub("%", "PERCENT", .x))),
                                            .y)))) %>%
        select(machine_readable_name, alt) %>%
        unnest(alt) %>%
        distinct %>%
        df2l
    if (missing(...)) {
        paste0(prefix, unique(unname(unlist(list))))
    } else {
        ifelse(... %in% names(list), paste0(prefix, list[...]), ...)
    }
}

mycharacteristics <- function(dset.crosssectional, dset.longitudinal) {
    characteristics <- rbind(dset.crosssectional %>% mutate(type = "crosssectional"),
                             dset.longitudinal %>% mutate(type = "longitudinal")) %>%
        select(type, starts_with("NMR_")) %>%
        rename_at(vars(starts_with("NMR_")), ~bionames(.)) %>%
        tableone::CreateTableOne(data = .,
                                 strata = "type")
    print(characteristics,
                      exact = "stage",
                      quote = FALSE,
                      noSpaces = TRUE,
                      printToggle = FALSE,
                      digits = 2,
                      pDigits = 3,
                      contDigits = 3)  %>%
        as.data.frame %>%
        tibble::rownames_to_column(var = "rowname") %>%
        format(justify = "left", trim = TRUE) %>%
        select(-one_of("test")) %>%
        filter(grepl("mean", rowname)) %>%
        mutate(rowname = gsub(" \\(mean \\(SD\\)\\)", "", rowname)) %>%
        arrange(rowname)
}

followup.pairs <- function(dset) {
    dset <- drop_na(dset) 

    health.baseline <- dset %>% filter(cohort == "T2000") %>% pull(sampleid)
    health.followup <- dset %>% filter(cohort == "T2011") %>% pull(sampleid)

    dilgom.baseline <- dset %>% filter(cohort == "F2007") %>% pull(sampleid)
    dilgom.followup <- dset %>% filter(cohort == "F2014") %>% pull(sampleid)

    data.frame(cohort = c("HEALTH", "DILGOM"),
               valid = c(health.baseline[health.baseline %in% health.followup] %>% length,
                         dilgom.baseline[dilgom.baseline %in% dilgom.followup] %>% length),
               missing = c(health.baseline[!health.baseline %in% health.followup] %>% length,
                            dilgom.baseline[!dilgom.baseline %in% dilgom.followup] %>% length))
}

four_iqr_outlier <- function(list) {
    (list > median(list) + IQR(list) * 4) | (list < median(list) - IQR(list) * 4)
}

nmr_normalization <- function(df) {
    df %>%
        mutate_at(vars(starts_with("NMR_")), ~log(.)) %>%
        filter_at(vars(starts_with("NMR_")), all_vars(!four_iqr_outlier(.))) %>%
        mutate_at(vars(starts_with("NMR_")), ~scale(.))
}
