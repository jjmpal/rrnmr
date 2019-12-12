c2l <- function(...) {
    l <- as.list(c(...))
    names(l) <- c(...)
    l
}

pub.p <- function(p) {
  p <- as.numeric(p)
  ifelse(p < 0.01, ifelse(p<0.001, "<0.001", sprintf("%.3f", p)), sprintf("%.2f", p))
}

getfactorvariables <- function(df, vars) {
    classes <- lapply(c2l(vars), function(x) class(df[[x]]))
    classes[classes == "factor"] %>% names
}

mygrep <- function(..., word, ignorecase = TRUE, complement = FALSE) {
    c(...)[xor(grepl(word, c(...), ignore.case = ignorecase), (complement == TRUE))]
}

betacip <- function(df, exclude = c("estimate", "conf.low", "conf.high", "qval", "pval",
                                    "statistic", "std.error")) {
    dplyr::mutate(df,
                  mean_ci = sprintf("%.2f (%.2f–%.2f)",
                                    estimate,
                                    conf.low,
                                    conf.high),
                  p.value = pub.p(qval)) %>%
        dplyr::select(-one_of(exclude))
}

standardnames <- function(x, prefix = "NMR_") {
    toupper(paste0(prefix, gsub("[/-]", "_", x)))
}

bionames <- function(metabolites) {
    md <- as.data.frame(ggforestplot::df_NG_biomarker_metadata) 
    metabolites %>%
        purrr::map_chr(function(id) {
            id <- gsub("NMR_", "", id)
            names <- md %>% filter(machine_readable_name == id)
            return(descriptivenames(names$abbreviation[1], names$description[1]))
        })
}

biogroups <- function(metabolites) {
    md <- as.data.frame(ggforestplot::df_NG_biomarker_metadata) 
    metabolites %>%
        purrr::map_chr(function(id) {
            id <- gsub("NMR_", "", id)
            names <- md %>% filter(machine_readable_name == id)
            return(names$group[1])
        })
}


metanames <- function(x) {
    case_when(x == "sys" ~ "Systolic BP",
              x == "dias" ~ "Diastolic BP",
              x == "age" ~ "Age",
              x == "bmi" ~ "BMI",
              x == "leisure" ~ "Exercise",
              x == "female" ~ "Female",
              TRUE ~ x)
}

descriptivenames <- function(abbr, desc, length = 20) {
    ifelse(nchar(desc) < length, desc, abbr)
}

myscale <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

mylogscale <- function(x) {
    minimumpositive <- min(x[x > 0], na.rm = TRUE)
    ifelse(x > 0,  x, minimumpositive/2) %>%
        log %>%
        myscale
}

getmetabolites <- function(dset) {
     dset %>% dplyr::select(starts_with("NMR_")) %>% colnames %>% sort
}

permetaboliteconditions <- function(dset, fun) {
    stopifnot(!missing(dset), !missing(fun))
    map_df(dset, ~as.data.frame(.x), .id = "setup") %>% 
        select(setup, starts_with("NMR_")) %>%
        rename_at(vars(-setup), ~bionames(.)) %>%
        group_by(setup) %>%
        mutate_at(., vars(-setup), ~ifelse(fun(.), 1, 0)) %>%
        summarize_all(sum) %>%
        gather(metabolite, outliers, -setup) %>%
        spread(setup, outliers) %>%
        filter(crosssectional + longitudinal > 0)
}

cleanlongitudinal <- function(df) {
    mutate(df, testee = gsub(".*_", "", sampleid)) %>%
        group_by(testee) %>%
        mutate(events = max(row_number())) %>%
        dplyr::ungroup() %>%
        dplyr::filter(events == 2) %>%
        dplyr::select(-events, -testee)
}

dropattributes <- function(x) {
    attributes(x) <- NULL
    x
}
