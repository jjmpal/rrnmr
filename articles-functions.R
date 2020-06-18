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

betacip <- function(df,
                    exclude = c("estimate", "conf.low", "conf.high", "qval", "pval", "statistic", "std.error"),
                    percent = FALSE) {
    dplyr::mutate(df, mean_ci = format.estimate_ci(estimate, conf.low, conf.high, percent)) %>%
        { if ("qval" %in% colnames(df)) mutate(., p.value = pub.p(qval)) else . } %>%
        dplyr::select(., colnames(.)[!colnames(.) %in% exclude])
}

format.estimate_ci <- function(estimate, conf.low, conf.high, percent = FALSE) {
    if(percent) sprintf("%.1f%% (%.1f–%.1f%%)",  estimate*100, conf.low*100, conf.high*100)
    else sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high)
}

standardnames <- function(x, prefix = "NMR_") {
    toupper(paste0(prefix, gsub("[/-]", "_", x)))
}

bioproperty <- function(metabolites, property = "name") {
    md <- as.data.frame(ggforestplot::df_NG_biomarker_metadata) 
    metabolites %>%
        purrr::map_chr(function(id) {
            id <- gsub("NMR_", "", id)
            names <- md %>% filter(machine_readable_name == id)
            if (property == "name")
                return(descriptivenames(names$abbreviation[1], names$description[1]))
            else
                return(names[[property]][1])
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

descriptivenames <- function(abbr, desc, length = 24) {
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
