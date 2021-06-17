## Operators

`%difference%` <- function(a, b) {
    a[!a %in% b]
}

`%intersect%` <- function(a, b) {
    intersect(a, b)
}


`%union%` <- function(a, b) {
    c(a, b)
}

## General functions

tidywithresponse <- function(model) {
    dependant <- ifelse(isS4(model),
                        all.vars(model@call)[1],
                        all.vars(model$call)[1])
    tidy(model, exponentiate = FALSE) %>%
        tibble::add_column(response = dependant, .before = 1)
}

myformula <- function(response, term, covariates, interaction, mixed) {
    terms <- term %union% covariates %>%
        { if(response %in% c("htn", "htn.female", "htn.male", "htn.old", "htn.young")) . %difference% "bptreat" else . }
    if (!missing(interaction)) {
        terms <- terms %union% paste0(term, ":", interaction)
    }
    if (!missing(mixed)) {
        if (length(names(mixed)) > 0)
            terms <- terms %union% sprintf("(%s|%s)", mixed, names(mixed))
        else
            terms <- terms %union% sprintf("(1|%s)", mixed) 
    }
    sprintf("%s ~ %s", response, paste(terms, collapse = " + "))
}

loop.results <- function(..., filterstr = "NMR_", idtoterm = FALSE) {
    purrr::map_df(c(...), tidywithresponse, .id = "model") %>%
        { if (idtoterm) dplyr::mutate(., term = model) else . } %>%
        dplyr::select(-model) %>%
        dplyr::filter(grepl(filterstr, term)) %>%
        dplyr::mutate(conf.low = estimate - qnorm(0.975) * std.error,
                      conf.high = estimate + qnorm(0.975) * std.error) %>%
        group_by(., response) %>%
        mutate(qval = p.adjust(p.value, method="BH")) %>%
        ungroup 
}


## Mixed-effect models

loop.lmer <- function(dset,
                      response,
                      mixed = "sampleid",
                      loops,
                      covariates = c()) {
    models <- parallel::mclapply(c2l(loops), function(loop) {
        fo <- myformula(response, loop, covariates, mixed = mixed)
        ret <- lmerTest::lmer(formula = formula(fo), data = dset)
        ret@call <- str2lang(sprintf("lmerModLmerTest(%s)", fo))
        ret
    }, mc.cores = 1)
}

loop.binomialmixed <- function(dset,
                    response,
                    loops,
                    mixed = c("sampleid"),
                    covariates = c()) {
    parallel::mclapply(c2l(loops), function(loop) {
        fo <- myformula(response, loop, covariates, mixed = mixed)
        ret <- lme4::glmer(formula = as.formula(fo),
                           family = binomial(link = "logit"),
                           data = dset)
        ret@call <- str2lang(sprintf("lme4::glmer(%s)", fo))
        ret
    }, mc.cores = 1)
}

## COX

loop.cox <- function(dset,
                    response,
                    loops,
                    covariates = c()) {
    models <- parallel::mclapply(c2l(loops), function(loop) {
        fo <- myformula(response, loop, covariates)
        ret <- coxph(formula = as.formula(fo), ties = "breslow", data = dset)
        ret$call <- as.formula(fo)
    }, mc.cores = 1)
}

## Linear and logistic models

loop.lm <- function(dset,
                    response,
                    loops,
                    covariates = c()) {
    models <- parallel::mclapply(c2l(loops), function(loop) {
        fo <- myformula(response, loop, covariates)
        ret <- stats::lm(formula = as.formula(fo), data = dset, na.action = na.omit)
        ret$call <- as.formula(fo)
        ret
    }, mc.cores = 1)
}

loop.rma <- function(dset,
                    response,
                    loops,
                    covariates = c()) {
    lapply(c2l(loops), function(loop) {
        ret_by_cohorts <- dset %>%
            group_by(cohort) %>%
            do(single.lm(., response, loop, covariates) %>%
               tidy %>%
               filter(term == loop))
        rma(yi = ret_by_cohorts$estimate, sei = ret_by_cohorts$std.error, method="FE")
    })
}

single.lm <- function(dset, response, metabolite, covariates = c()) {
    fo <- myformula(response, metabolite, covariates)
    ret <- stats::lm(formula = as.formula(fo), data = dset, na.action = na.omit)
    ret$call <- as.formula(fo)
    ret
}

loop.gamma <- function(dset,
                    response,
                    loops,
                    covariates = c()) {
    stopifnot(!missing(dset), !missing(response), !missing(loops))
    mclapply(c2l(loops), function(loop) {
        fo <- myformula(response, loop, covariates)
        ret <- stats::glm(formula = as.formula(fo),
                          family = inverse.gaussian(link = "identity"),
                          data = dset)
        ret$call <- as.formula(fo)
        ret
    }, mc.cores = 1)
}

loop.binomial <- function(dset,
                    response,
                    loops,
                    covariates = c()) {
    stopifnot(!missing(dset), !missing(response), !missing(loops))
    lapply(c2l(loops), function(loop) {
        fo <- myformula(response, loop, covariates)
        ret <- stats::glm(formula = as.formula(fo),
                          family=binomial(link='logit'),
                          data = dset,
                          na.action = na.omit)
        ret$call <- as.formula(fo)
        ret
    })
}


check.proportionality <- function(ordered, multinomial, loops) {
    parallel::mclapply(c2l(loops), function(loop) {
        ord_dof <- attributes(logLik(ordered[[loop]]))$df
        mul_dof <- attributes(logLik(multinomial[[loop]]))$df
        ord_dev <- ordered[[loop]]$deviance
        mul_dev <- multinomial[[loop]]$deviance
        pchisq(abs(ord_dev - mul_dev), abs(ord_dof - mul_dof), lower.tail = FALSE)
    }, mc.cores = 1) %>%
        purrr::map_df(., ~data.frame("p.value" = .x), .id = "model")
}

loop.residuals <- function(..., ignore = c()) {
    models <- c(...)
    lapply(models[names(models) %difference% ignore], function(model) {
        term <- ifelse(isS4(model),
                       all.vars(model@call)[2],
                       all.vars(model$call)[2])
        ggplot(broom::augment(model), aes(.fitted,.resid)) +
            geom_point(shape = ".") +
            geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
            geom_hline(yintercept = 0, linetype = "dotted") +
            scale_x_continuous(name = "Fitted Values") +
            scale_y_continuous(name = "Residual") +
            ggtitle(sprintf("%s", bioproperty(term))) +
            theme_bw()
    })
}

loop.qq <- function(..., ignore = c()) {
    models <- c(...)
    lapply(models[names(models) %difference% ignore], function(model) {
        term <- ifelse(isS4(model),
                       all.vars(model@call)[2],
                       all.vars(model$call)[2])
        ggplot(data = broom::augment(model), mapping = aes(sample = .resid)) +
            geom_abline(intercept = 0, size = 1, slope=1, color = "blue") +
            stat_qq(geom = "line", size = 1) +
            scale_x_continuous("Fitted Values") +
            scale_y_continuous("Residual") +
            coord_fixed() +
            ggtitle(sprintf("%s", bioproperty(term))) +
            theme_bw()
    })
}


results.table <- function(df,
                          exponentiate = c("htn", "htn3", "htn.followup"),
                          percent = FALSE,
                          drop = FALSE,
                          nosubclass = TRUE) {
    dplyr::mutate(df,
                  group = bioproperty(term, "group"),
                  term = bioproperty(term),
                  estimate = ifelse(response %in% exponentiate, exp(estimate), estimate),
                  conf.low = ifelse(response %in% exponentiate, exp(conf.low), conf.low),
                  conf.high = ifelse(response %in% exponentiate, exp(conf.high), conf.high)) %>%
        { if (nosubclass) dplyr::filter(., group != "Lipoprotein subclasses") else . } %>%
        betacip(percent = percent) %>%
        myspread() %>%
        { if (drop) filter_at(., vars(contains("p.value")), any_vars(. < 0.05)) else .}
}

# Others

filternullcohorts <- function(dset, var, missinglevel = -0.5) {
    exclude <- dset %>%
        group_by(cohort) %>%
        summarize(value = max(!! rlang::sym(var)) < missinglevel) %>%
        filter(value == TRUE) %>%
        pull(cohort)
    if (length(exclude) > 0)
        message(sprintf("Excluding for %s cohorts %s.", var, paste(exclude, collapse = ",")))
    dset %>% dplyr::filter(!cohort %in% exclude)
}

myprcomp <- function(matrix) {
    df <- matdf(matrix)
    prcomp.ret <- prcomp(~.,
           data = dplyr::select(df, starts_with("NMR_")),
           center = TRUE,
           scale = TRUE,
           na.action = na.omit)
    prcomp.loading <- matdf(prcomp.ret$rotation)
    prcomp.contribution <- prcomp.ret$sdev^2/sum(prcomp.ret$sdev^2)
    prcomp.axes <- full_join(df, matdf(prcomp.ret$x), by = "rowname") %>%
        prcomprenamer %>%
        mutate_at(vars(starts_with("PC")), scale)
    list("axes" = prcomp.axes,
         "contribution " = prcomp.contribution,
         "model" = prcomp.ret,
         "loading" = prcomp.loading)
}

mypredict <- function(mod, dset) {
    full_join(dset %>% select(-starts_with("NMR_")) %>% matdf,
              predict(mod$model, dset %>% select(starts_with("NMR_"))) %>% matdf,
              by = "rowname") %>%
        select(-rowname)
}

prcomprenamer <- function(df) {
    as.data.frame(df) %>%
        rename_at(vars(starts_with("PC")),
                  ~sprintf("PC%03i", as.numeric(gsub("PC", "", .))))
}

matdf <- function(matrix) {
    as.data.frame(matrix) %>% tibble::rownames_to_column(var = "rowname")
}

myoutcomes <- function(df, response = "htn.followup") {
    df %>% pull(!!response) %>% as.numeric - 1
}

myglm <- function(vars, df, response = "htn.followup") {
    fo <- sprintf("%s ~ %s", response, paste0(vars, collapse = "+"))
    glm(as.formula(fo), family = "binomial", data = df, x=TRUE)
}

correlationmatrix  <- function(dset, vars, method = "spearman", adjust = "fdr") {
    dat <- dset %>% dplyr::select(one_of(vars[["id"]])) %>% na.omit()
    colnames(dat) <- vars[["group_name"]]
    psych::corr.test(as.matrix(dat), method = method, adjust = adjust)
}



mynri <- function(std, new, cut = 0, niter = 100, alpha = 0.05) {
    invisible(capture.output(ret <- nribin(mdl.std = std,
                                           mdl.new = new,
                                           cut = cut,
                                           niter = niter,
                                           alpha = alpha,
                                           msg = FALSE,
                                           updown = 'diff')))
    return(ret)
    ret$nri %>%
        as.data.frame %>%
        tibble::rownames_to_column(var = "term") %>%
        mutate_if(is.numeric, round, 3) %>%
        filter(grepl("NRI", term)) %>%
        select(-Std.Error) %>%
        rename(estimate = Estimate, conf.low = Lower, conf.high = Upper)
}

myreclassification <-  function (data, response = "htn.followup", predrisk1, predrisk2, z = qnorm(0.975))  {
    y <- Hmisc::improveProb(x1 = predrisk1, x2 = predrisk2, y = myoutcomes(data, response = response))
    return(data.frame(model = c("NRI", "NRI+", "NRI-", "IDI"),
                      estimate = c(y$nri, y$nri.ev, y$nri.ne, y$idi),
                      conf.low = c(y$nri - z*y$se.nri,
                                   y$nri.ev - z*y$se.nri.ev,
                                   y$nri.ne - z*y$se.nri.ne,
                                   y$idi - z*y$se.idi),
                      conf.high = c(y$nri + z*y$se.nri,
                                    y$nri.ev + z*y$se.nri.ev,
                                    y$nri.ne + z*y$se.nri.ne,
                                    y$idi + z*y$se.idi),
           p.value = 2*pnorm(-abs(c(y$z.nri,
                                    y$z.nri.ev,
                                    y$z.nri.ne,
                                    y$z.idi)))))
}

getincidental <- function(dset, keephypertensive = TRUE) {
    dset %>%
        filter(cohort %in% c("F2007", "T2000"), (htn == 0 | keephypertensive)) %>% 
        left_join(.,
                  dset %>%
                  dplyr::filter(cohort %in% c("F2014", "T2011")) %>%
                  dplyr::select(sampleid, htn, sys, dias, bptreat),
                  by = "sampleid",
                  suffix = c("", ".followup")) %>%
        mutate(dias.diff.pct = (dias.followup - dias)/dias,
               sys.diff.pct = (sys.followup - sys)/sys,
               agegroup = ifelse(age < median(age), "young", "old"))
}

getLongitudal <- function(dset, keephypertensive = FALSE) {
    ids <- dset %>%
        filter(cohort %in% c("F2007", "T2000"), htn == 0) %>%
        pull(sampleid)
    dset %>% dplyr::filter(sampleid %in% ids)
}

stratify.hypertension <- function(dset, value, rename, covariate) {
    dset %>%
        filter(htn == value) %>%
        rename(!!quo_name(rename) := sys)
}
