get_cv_dset <- function(dset, nfolds) {
    dset <- dset %>% tibble::rowid_to_column("id")
    dset_train <- dset %>% sample_frac((nfolds-1)/nfolds)
    dset_test <- dplyr::anti_join(dset, dset_train, by = 'id')
    list(train = dset_train %>% select(-id),
         test = dset_test %>% select(-id))
}

my_prepare_matrix <- function(dset) {
    my_label <- dset %>% pull(sys)
    my_matrix <- dset %>% select(-any_of(c("sys", "fold"))) %>% data.matrix
    xgb.DMatrix(my_matrix, label = my_label)
}

my_leaveoneout_crossvalidation <- function(x,
                                         params,
                                         dset,
                                         debug = FALSE) {
    matrix_train <- dset %>% dplyr::filter(fold != x["fold"]) %>% my_prepare_matrix
    matrix_test <- dset %>% dplyr::filter(fold == x["fold"]) %>% my_prepare_matrix

    if (debug) {
        sprintf("Single objective function running training N=%i and test N=%i",
                nrow(matrix_train), nrow(matrix_test)) %>%
            message
    }
    allparams <- modifyList(params,
                         list(eta = x["eta"],
                              max_depth = x["max_depth"],
                              min_child_weight = x["min_child_weight"],
                              gamma = x["gamma"],
                              subsample = x["subsample"],
                              colsample_bytree = x["colsample_bytree"]))
    
    mywatchlist <- list(train = matrix_train, test = matrix_test)
    
    xgboost_loss_function(params = allparams,
                          data =  matrix_train,
                          nrounds = x["nrounds"],
                          watchlist = mywatchlist)
}

my_n_fold_crossvalidation <- function(x,
                                         params,
                                         dset,
                                         nfolds = 5,
                                         debug = FALSE) {
    dsets <- get_cv_dset(dset, nfolds)
    matrix_train <- dsets[["train"]] %>% my_prepare_matrix
    matrix_test <- dsets[["test"]] %>% my_prepare_matrix
    
    if (debug) {
        sprintf("Single objective function running training N=%i and test N=%i",
                nrow(matrix_train), nrow(matrix_test)) %>%
            message
    }
    allparams <- modifyList(params,
                         list(eta = x["eta"],
                              max_depth = x["max_depth"],
                              min_child_weight = x["min_child_weight"],
                              gamma = x["gamma"],
                              subsample = x["subsample"],
                              colsample_bytree = x["colsample_bytree"]))
    
    mywatchlist <- list(train = matrix_train, test = matrix_test)
    
    xgboost_loss_function(params = allparams,
                          data =  matrix_train,
                          nrounds = x["nrounds"],
                          watchlist = mywatchlist)
}

xgboost_loss_function <- function(data, watchlist, nrounds, params) {
    ret <- xgb.train(data = data,
              watchlist = watchlist,
              params = params,
              nrounds = nrounds,
              verbose = FALSE)
    ret[["evaluation_log"]] %>% pull(test_rmse) %>% tail(n=1)
}


my_optimal_xboost <- function(dset,
                               params,
                               nrounds,
                               early_stepping_rounds = 25,
                              nthreads = 40) {
    xgboost(params = params,
            data = my_prepare_matrix(dset),
            nrounds = nrounds,
            nthread = nthreads,
            prediction = FALSE,
            showsd = TRUE,
            early_stopping_rounds = early_stepping_rounds,
            verbose = 0,
            print_every_n = 500,)
}

# If no design is given by the user, mlrMBO will generate a maximin Latin Hypercube Design of size 4 times the number of the black-box function’s parameters.
do_bayes <- function(n_design = NULL, opt_steps = NULL, of = obj.fun) {
    des <- ParamHelpers::generateDesign(n = n_design,
                                        par.set = getParamSet(of),
                                        fun = lhs::randomLHS)
    control <- mlrMBO::setMBOControlTermination(makeMBOControl(), iters = opt_steps)
    learner <- makeLearner("regr.km",
                           predict.type = "se",
                           covtype = "matern3_2",
                           control = list(trace = FALSE))
    mlrMBO::mbo(fun = of,
                design = des,
                learner = learner,
                control = control, 
                show.info = TRUE)
}

get_true_and_predicted_feature <- function(dset, model) {
    sys_true <- dset %>% pull("sys")
    sys_pred <- predict(model, my_prepare_matrix(dset))
    tibble(sys = sys_true, pred = sys_pred)
}

summarize_run <- function(dset, model, plot = FALSE) {
    df <- get_true_and_predicted_feature(dset, model)
    if (plot) {
        df %>%
            arrange(sys) %>%
            tibble::rowid_to_column("x") %>%
            ggplot(aes(x = x)) +
            geom_line(aes(y = sys), , colour="green") +
            geom_line(aes(y = pred) , colour="red")
    } else {
        df %>%
            mutate(diff = sys - pred) %>%
            summarise(diff_mean = mean(diff),
                      diff_sd = sd(diff),
                      diff_min = min(diff),
                      diff_max = max(diff),
                      mae = caret::MAE(sys, pred),
                      rmse = caret::RMSE(sys, pred),
                      r2 = 1 - var(sys - pred)/var(sys))
    }
}

get_nparameters <- function(params) {
    length(params[["pars"]])
}

my_pdplot <- function(dset, model, vars, palette) {
    dset_pdp <- dset %>% select(-any_of(c("sys", "fold"))) %>% data.matrix

    p <- ggplot()
    labels <- tibble(x = numeric(), y = numeric(), color = character(), desc = character())
    for (i in seq_along(vars)) {
        plot_data <- pdp::partial(model,
                                  pred.var = vars[[i]],
                                  train = dset_pdp,
                                  type = "regression",
                                  plot = FALSE) %>%
            as_tibble %>%
            rename(y = yhat, x := !!quo_name(vars[[i]])) %>%
            mutate(color = vars[[i]],
                   desc = new_colnames[[vars[[i]]]])
        labels <- add_row(labels, closest_to_point(plot_data, point = 2.0))
        p <-  p +
            geom_line(data = plot_data, aes(x = x, y = y, color = color))
    }

    p <- p +
        geom_text_repel(
            data = labels,
            aes(x = x, y = y, label = desc, color = color),
            size = 4,
            hjust = 0,
            direction = "y",
            nudge_x = 0.3,
            segment.color = NA,
            show.legend = FALSE)

    p +
        scale_x_continuous(breaks = seq(-2,2,1),
                           expand = expansion(add = c(0.5, 3.5)),
                           limits = c(-2, 2)) +
        scale_y_continuous(breaks = seq(110, 180, 2), limits = c(130, 140)) +
        scale_color_manual(guide = "none",
                           values = palette) +
        xlab("1-SD change in metabolic measure") +
        ylab("Systolic BP") +
        theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
}

my_multi_pdplot <- function(dset, model, vars, palette, extra = "female") {
    dset_pdp <- dset %>% select(-any_of(c("sys", "fold"))) %>% data.matrix

    p <- ggplot()
    labels <- tibble(x = numeric(), y = numeric(), color = character(), desc = character())
    for (i in seq_along(vars)) {
        plot_data <- pdp::partial(model,
                                  pred.var = c(vars[[i]], extra),
                                  train = dset_pdp,
                                  type = "regression",
                                  plot = FALSE) %>%
            as_tibble %>%
            rename(y = yhat, x := !!quo_name(vars[[i]])) %>%
            mutate(color = vars[[i]],
                   desc = new_colnames[[vars[[i]]]])
        labels <- add_row(labels, closest_to_point(plot_data, point = 2.0))
        p <-  p +
            geom_line(data = plot_data, aes(x = x, y = y, color = color))
    }

    p <- p +
        geom_text_repel(
            data = labels,
            aes(x = x, y = y, label = desc, color = color),
            size = 4,
            hjust = 0,
            direction = "y",
            nudge_x = 0.3,
            segment.color = NA,
            show.legend = FALSE)

    p +
        scale_x_continuous(breaks = seq(-2,2,1),
                           expand = expansion(add = c(0.5, 3.5)),
                           limits = c(-2, 2)) +
        scale_y_continuous(breaks = seq(110, 180, 2), limits = c(130, 140)) +
        scale_color_manual(guide = "none",
                           values = palette) +
        xlab("1-SD change in metabolic measure") +
        ylab("Systolic BP") +
        theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
}


get_longitudinal_dset <- function(dset) {
    dset %>%
        filter(cohort %in% c("F2007", "T2000")) %>%
        rename(sys_old = sys) %>%
        left_join(.,
                  dset %>%
                  dplyr::filter(cohort %in% c("F2014", "T2011")) %>%
                  dplyr::select(sampleid, sys),
                  by = "sampleid") %>%
        mutate(sys_diff = sys - sys_old)
}

fix_label_heights <- function(plot) {
#    g <- ggplotGrob(plot)
    g <- ggplot_gtable(ggplot_build(plot))
    index <- which(sapply(g$grobs, function(x) x$name == "strip"))
    g$grobs <- lapply(seq_along(g$grobs), function(.x) {
        if(.x %in% index) {
            g$grobs[[.x]]$heights <- unit(3,"mm")
        } 
        g$grobs[[.x]]
    } )
    return(g)
}

fix_ggplot_labels <- function(plot) {
    pg <- ggplotGrob(plot)
    for(i in which(grepl("strip-t", pg$layout$name))){
        pg$grobs[[i]]$layout$clip <- "off"
    }
    return(pg)
}

nmetabolites <- function(dset) {
    dset %>% colnames %>% grepl("^NMR", .) %>% sum
}

get_xgb_top_gain <- function(...) {
    list(...) %>%
        map_dbl(~ xgb.importance(model = .x) %>% pull(Gain) %>% max) %>%
        max
}

round_up_decimal <- function(x) {
    ceiling(x * 10)/10
}

