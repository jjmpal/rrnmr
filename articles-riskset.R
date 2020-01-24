getriskset <- function(df, regularisation, by = "HAVTUN") {
    stopifnot(typeof(regularisation) == "list", typeof(regularisation[[1]]) == "list")
    rf <- lapply(c2l(names(regularisation)), function(x) 
        riskscore(df, regularisation[[x]], prefix = paste0(x, "_"), id = by))
    Reduce(function(x,y) full_join(x, y, by=by, all = TRUE), c(list(df), rf))
}

riskscore <- function(dset, betas, nriskclass=4, prefix = "", id = "HAVTUN", subset = "nmr_") {
    betas.nmr <- dplyr::filter(betas, grepl(subset, term))
    dset %>% mutate(risk = apply(.[betas.nmr$term], 1, weighted.sum, betas.nmr$estimate),
                    riskclass = relevel(factor(dplyr::ntile(risk, nriskclass)), ref='1'),
                    riskpersd = scale(risk)) %>%
        select(id, riskclass, riskpersd) %>%
        dplyr::rename(!!paste0(prefix, "riskclass") := riskclass,
                      !!paste0(prefix, "riskpersd") := riskpersd) 
}

weighted.sum <- function(x, w=rep(1, length(x))) {
    sum(x*w)
} 
