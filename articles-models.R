## Operators

`%difference%` <- function(a, b) {
    a[!a %in% b]
}

`%union%` <- function(a, b) {
    c(a, b)
}

## General functions

mytidy <- function(model) {
    dependant <- ifelse(isS4(model),
                        all.vars(model@call)[1],
                        all.vars(model$call)[1])
    tidy(model, exponentiate = FALSE) %>%
        tibble::add_column(response = dependant, .before = 1)
}

myformula <- function(response, term, covariates, interaction, mixed) {
    terms <- term %union% covariates
    if (!missing(interaction)) {
        terms <- terms %union% paste0(term, ":", interaction)
    }
    if (!missing(mixed)) {
        terms <- terms %union% sprintf("(1|%s)", mixed)
    }
    sprintf("%s ~ %s", response, paste(terms, collapse = " + "))
}

loop.results <- function(..., filterstr = "NMR_") {
    purrr::map_df(c(...), ~as.data.frame(mytidy(.x))) %>%
        dplyr::filter(grepl(filterstr, term)) %>%
        dplyr::mutate(pval = pnorm(abs(statistic), lower.tail = FALSE) * 2,
                      qval = p.adjust(pval, method="BH"),
                      conf.low = estimate - qnorm(0.975) * std.error,
                      conf.high = estimate + qnorm(0.975) * std.error)
}


## Mixed-effect models

loop.lmer <- function(dset,
                    response,
                    loops,
                    covariates = c()) {
    models <- parallel::mclapply(c2l(loops), function(loop) {
        fo <- myformula(response, loop, covariates, mixed = "sampleid")
        ret <- lmerTest::lmer(formula = formula(fo), data = dset)
        ret@call <- str2lang(sprintf("lmerModLmerTest(%s)", fo))
        ret
    }, mc.cores = min(length(loops), 8))
}

loop.binomialmixed <- function(dset,
                    response,
                    loops,
                    covariates = c()) {
    models <- parallel::mclapply(c2l(loops), function(loop) {
        fo <- myformula(response, loop, covariates, mixed = "sampleid")
        ret <- lme4::glmer(formula = as.formula(fo),
                           family = binomial(link = "logit"),
                           data = dset)
        ret@call <- str2lang(sprintf("lme4::glmer(%s)", fo))
        ret
    }, mc.cores = min(length(loops), 8))
}

## Linear and logistic models

loop.lm <- function(dset,
                    response,
                    loops,
                    covariates = c()) {
    models <- parallel::mclapply(c2l(loops), function(loop) {
        fo <- myformula(response, loop, covariates)
        ret <- stats::lm(formula = as.formula(fo), data = dset)
        ret$call <- as.formula(fo)
        ret
    }, mc.cores = min(length(loops), 8))
}

loop.binomial <- function(dset,
                    response,
                    loops,
                    covariates = c()) {
    stopifnot(!missing(dset), !missing(response), !missing(loops))
    lapply(c2l(loops), function(loop) {
        fo <- myformula(response, loop, covariates)
        ret <- stats::glm(formula = as.formula(fo), family=binomial(link='logit'), data = dset)
        ret$call <- as.formula(fo)
        ret
    })#, mc.cores = min(length(loops), 8))
}


check.proportionality <- function(ordered, multinomial, loops) {
    parallel::mclapply(c2l(loops), function(loop) {
        ord_dof <- attributes(logLik(ordered[[loop]]))$df
        mul_dof <- attributes(logLik(multinomial[[loop]]))$df
        ord_dev <- ordered[[loop]]$deviance
        mul_dev <- multinomial[[loop]]$deviance
        pchisq(abs(ord_dev - mul_dev), abs(ord_dof - mul_dof), lower.tail = FALSE)
    }, mc.cores = min(length(loops), 8)) %>%
        purrr::map_df(., ~data.frame("p.value" = .x), .id = "model")
}

loop.residuals <- function(...) {
    lapply(c(...), function(model) {
        term <- ifelse(isS4(model),
                       all.vars(model@call)[2],
                       all.vars(model$call)[2])
        ggplot(data.frame(fitted = fitted(model),
                          resid = resid(model)),
               aes(x = fitted, y = resid)) +
        geom_jitter(shape = ".") +
            geom_smooth(se = FALSE, method = "loess", formula = "y ~ x") +
            geom_hline(yintercept = 0, linetype = "dotted") +
            scale_x_continuous("Fitted Values") +
            scale_y_continuous("Residual") +
            ggtitle(sprintf("%s", bionames(term))) +
            theme_bw()
    })
}

results.table <- function(df, exponentiate = c("htn", "htn3")) {
    dplyr::mutate(df,
                  term = bionames(term),
                  estimate = ifelse(response %in% exponentiate, exp(estimate), estimate),
                  conf.low = ifelse(response %in% exponentiate, exp(conf.low), conf.low),
                  conf.high = ifelse(response %in% exponentiate, exp(conf.high), conf.high)) %>%
        betacip() %>%
        myspread()
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
    prcomp.ret <- prcomp(df %>% dplyr::select(starts_with("NMR_")))
    prcomp.loading <- matdf(prcomp.ret$rotation)
    prcomp.contribution <- prcomp.ret$sdev^2/sum(prcomp.ret$sdev^2)
    prcomp.axes <- full_join(df, matdf(prcomp.ret$x), by = "rowname") %>%
             prcomprenamer
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

spearmancorrelation  <- function(dset, vars) {
  dat <- dset %>% dplyr::select(vars)
  colnames(dat) <- colnames(dat) %>% bionames()
  cor <- cor(dat, method = 'spearman')
}

mynri <- function(std, new, cut = 0, niter = 100, alpha = 0.05) {
    invisible(capture.output(ret <- nribin(mdl.std = std,
                                           mdl.new = new,
                                           cut = cut,
                                           niter = niter,
                                           alpha = alpha,
                                           msg = FALSE,
                                           updown = 'diff')))
    ret$nri %>%
        as.data.frame %>%
        tibble::rownames_to_column(var = "term") %>%
        mutate_if(is.numeric, round, 3) %>%
        filter(grepl("NRI", term)) %>%
        select(-Std.Error) %>%
        rename(estimate = Estimate, conf.low = Lower, conf.high = Upper)
}