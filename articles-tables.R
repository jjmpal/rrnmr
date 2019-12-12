characteristics.names <- function(onlyvars = FALSE) {
    if (onlyvars == TRUE)
        return(c("age", "female", "bmi", "sys", "dias", "htn", "smoker",
                 "diabetes", "leisure", "bptreat", "koltreat"))
    list("n" = "N",
         "age" = "Age, y (SD)",
         "female" = "Female, N (%)",
         "bmi" = "BMI, kg/m² (SD)",
         "sys" = "Systolic BP, mmHg (SD)",
         "dias" = "Diastolic BP, mmHg (SD)",
         "htn" = "Hypertension, N (%)",
         "smoker" = "Current smoker, N (%)",
         "diabetes" = "Diabetes mellitus, N (%)",
         "leisure" = "Exercise, N (%)",
         "bptreat" = "Antihypertensive medication, N (%)",
         "koltreat" = "Lipid medication, N (%)",
         "1"  =  "  Light",
         "2"  =  "  Moderate",
         "3"  =  "  Heavy",
         "4"  = "  Competitive")


}

characteristicsTableFull <- function(dset) {
    table.sub <- characteristicsTable(dset, "cohort")
    table.tot <- characteristicsTable(dset)
    cbind(table.tot, table.sub)
}

characteristicsTable <- function(dset, strata) {
    nostrata <- missing(strata)
    characteristics <- tableone::CreateTableOne(
                                     strata = strata,
                                     data = dset,
                                     vars = characteristics.names(TRUE),
                                     factorVars = getfactorvariables(dset,
                                                                     characteristics.names(TRUE)),
                                     test = !missing(strata))
    print(characteristics,
                      exact = "stage",
                      quote = FALSE,
                      noSpaces = TRUE,
                      printToggle = FALSE,
                      digits = 1,
                      pDigits = 3,
                      contDigits=1)  %>%
        as.data.frame %>%
        tibble::rownames_to_column(var = "rowname") %>%
        format(justify = "left", trim = TRUE) %>%
        mutate(rowname = characteristics.names()[gsub("^ *([A-Za-z_0-9]+).*", "\\1", rowname)]) %>%
        { if (!nostrata) select(., -test) else . }
}

tableone <- function(characteristics) {
    flextable(data = characteristics) %>%
#        set_header_labels(rowname = "Characteristics",
#                          Overall = paste0("Cases, n=", dim(dset)[1])) %>%
        flextable::width(j = 1, width = 3) %>%
        flextable::width(j = 2, width = 1.2) %>%
        flextable::border(border = fp_border(width=0), part="body") %>%
        flextable::border(border = fp_border(width=0), part="header") %>%
        flextable::border(part="header", border.bottom = fp_border(width=1)) %>%
    flextable::border(i = nrow(characteristics),
                      part="body", border.bottom = fp_border(width=1)) %>%
        flextable::bold(bold = FALSE, part = "header") %>%
        flextable::bold(bold = FALSE, part = "body") %>%
        flextable::fontsize(size = 12, part = "header") %>%
        flextable::fontsize(size = 12, part = "body") %>%
        flextable::align(align = "center", part = "all") %>%
        flextable::align(align = "left", part = "header", j = 1) %>%
        flextable::align(align = "left", part = "body", j = 1)
}


oddsratiotable <- function(df, modelres, regularisation, htn = "HT") {
    riskmodel.results  <- createoddstable(df, modelres, names(regularisation), htn = htn)

    regruns.hits <- purrr::map_df(regularisation, ~as.data.frame(.x), .id="id") %>%
        dplyr::filter(grepl("mzid", term)) %>%
        group_by(id) %>%
        summarize(n = n()) %>%
        spread(id, n)
    
    typology.tbl4 <- data.frame(
        col_keys = colnames(riskmodel.results),
        what = c( "", "", 
                 rep(paste("Risk score based on", regruns.hits$fwdbonf, "eicosanoids selected",
                           "using forward selection with Bonferroni-based",
                           "threshold"), 3),
                 rep(paste("Risk score based on", regruns.hits$lasso ,"eicosanoids selected using LASSO"), 3),
                 rep(paste("Risk score based on", regruns.hits$fwdaic ,"eicosanoids selected using forward",
                           "selection aiming at maximum model AIC"), 3)),
        measure = c("", "N of individuals", 
                    rep(c("N with HTN", 
                          "Unadjusted OR (95% CI)", 
                          "Adjusted OR (95% CI)"), 3)),
        stringsAsFactors = FALSE)

    typologyformatter(
        data = riskmodel.results, 
        typology = typology.tbl4,
        font = 10) %>%
        flextable::width(j = 1, width = 0.8) %>%
        flextable::width(j = 2, width = 0.8) %>%
        flextable::width(j = c(3, 6, 9), width = 0.5) %>%
        flextable::width(j = c(4:5, 7:8, 10:11), width = 1.3)
}

smalloddsratiotable <- function(df, modelresbp, modelreshtn, regularisation, htn = "HT") {
    riskmodel.results <- cbind(createoddstable(df, modelresbp, "fwdbonf", htn = htn) %>%
                               rename(fwdbonf.unadjusted.mean_ci_sbp = fwdbonf.unadjusted.mean_ci,
                                      fwdbonf.adjusted.mean_ci_sbp = fwdbonf.adjusted.mean_ci),
                               createoddstable(df, modelreshtn, "fwdbonf", htn = htn)[4:5])

    typology.tbl4 <- data.frame(
        col_keys = colnames(riskmodel.results),
        what = c( "", "", "", 
                 rep(paste("Odds for systolic blood pressure"), 2),
                 rep(paste("Odds for hypertension"), 2)),
        measure = c("", "N of individuals", 
                    "N with HTN",
                    rep(c("Unadjusted OR (95% CI)", 
                          "Adjusted OR (95% CI)"), 2)),
        stringsAsFactors = FALSE)

    typologyformatter(
        data = riskmodel.results, 
        typology = typology.tbl4,
        font = 10) %>%
        flextable::width(j = 1, width = 0.8) %>%
        flextable::width(j = 2, width = 0.8) 
#        flextable::width(j = c(3, 6), width = 0.5) %>%
#        flextable::width(j = c(4:5, 7:8), width = 1.3)
}

myspread <- function(ret, list = c2l("mean_ci", "p.value"), key = "response", by = "term") {
    lapply(list, function(x) ret %>%
                             dplyr::select(term, key, x) %>%
                             tidyr::spread(key, x) %>%
                             dplyr::rename_at(vars(-term), ~paste0(., "_", x))) %>%
        purrr::reduce(full_join, by = by)  %>%
        dplyr::select(term, noquote(order(colnames(.))))
}


mykable <- function(df) {
    knitr::kable(df, format = ifelse(exists("rmdgeneration"), "html", "markdown")) %>%
        { if (exists("rmdgeneration") == TRUE)
              kable_styling(., bootstrap_options = "striped", full_width = F)
          else
              .}
}

improveprob.tidy <- function(x) {
    data.frame(index = c("NRI", "IDI"),
               estimate = c(x$nri, x$idi),
               conf.low = c(x$nri - qnorm(1- 0.05/2)*x$se.nri,
                            x$idi - qnorm(1- 0.05/2)*x$se.idi),
               conf.high = c(x$nri + qnorm(1- 0.05/2)*x$se.nri,
                             x$idi + qnorm(1- 0.05/2)*x$se.idi))
}

roc.tidy <- function(df, models) {
    lapply(models, function(model) {
        roc.ret <- pROC::roc(myoutcomes(df), fitted(model))
        roc.ci <- pROC::ci(roc.ret)
        attributes(roc.ci) <- NULL
        data.frame(estimate = roc.ci[2],
                   conf.low = roc.ci[1],
                   conf.high = roc.ci[3])
    }) %>%
        purrr::map_df(~as.data.frame(.x), .id="model")


}

myresulttable <- function(list) {
    lapply(list, function(x) {
            broom::tidy(x) %>%
            mutate(conf.low = estimate - qnorm(0.975) * std.error,
                   conf.high = estimate + qnorm(0.975) * std.error,
                   mean_ci = sprintf("%.2f (%.2f–%.2f)%s",
                                     estimate,
                                     conf.low,
                                     conf.high,
                                     gtools::stars.pval(p.value))) %>%
            select(term, mean_ci)
    }) %>%
        purrr::map_df(., ~as.data.frame(.x), .id = "model") %>%
        spread(model, mean_ci)
} 
