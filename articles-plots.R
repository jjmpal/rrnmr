myforestplot <- function(dset,
                         responses,
                         logodds = FALSE,
                         xlim = c(-1, 3),
                         legend.position = c(-1.15, -3.8),
                         panelone = c("Cholesterol",
                                      "Lipoprotein particle sizes",
                                      "Apolipoproteins",
                                      "Fatty acids",
                                      "Glycolysis related metabolites"),
                         paneltwo = c("Glycerides and phospholipids",
                                      "Fluid balance",
                                      "Ketone bodies",
                                      "Amino acids",
                                      "Inflammation"),
                         file,
                         dims = c(18, 22),
                         percentage = FALSE,
                         scale = bpscale,
                         preprocess = forestplotpreprocess) {
    stopifnot(!missing(responses), !missing(file))
    stopifnot(responses %in% unique(dset$response))
    
    df.metadata <- preprocess(dset, responses)
    groups <- df.metadata %>% pull(nmrgroup) %>% unique

    ggplot_multi <- myforest.loopgroups(df.metadata,
                                        tail(panelone %intersect% groups, n = 1),
                                        tail(paneltwo %intersect% groups, n = 1),
                                        logodds = logodds,
                                        legend.position = legend.position,
                                        percentage = percentage,
                                        xlim = xlim,
                                        scale = scale,
                                        omitlegend = length(responses) == 1)

    heightadj <- 2/length(responses)

    ag <- arrangeGrob(forest.order.and.scale(ggplot_multi, panelone %intersect% groups, heightadj),
                      forest.order.and.scale(ggplot_multi, paneltwo %intersect% groups, heightadj),
                      ncol = 2)

    ggsave(file = file, plot = ag, width = dims[1], height = dims[2], dpi = 300, unit = "cm")
}

forestplotpreprocess <- function(dset, responses) {
    dset %>%
        mutate(nmrname = bioproperty(term),
               nmrgroup = bioproperty(term, "group")) %>%
        filter(response %in% responses) %>%
        arrange(nmrgroup)
}


subclasspreprocess <- function(dset, responses) {    
    dset %>%
        mutate(nmrname = bioproperty(term, length = 0),
               nmrgroup = bioproperty(term, "group")) %>%
        filter(response %in% responses, nmrgroup == "Lipoprotein subclasses") %>%
        mutate(nmrgroup = case_when(grepl("IDL", nmrname) ~ "Intermediate-density lipoprotein",
                                 grepl("VLDL", nmrname) ~ "Very low-density lipoprotein",
                                 grepl("LDL", nmrname) ~ "Low-density lipoprotein",
                                 grepl("HDL", nmrname) ~ "High-density lipoprotein",
                                 TRUE ~ "General lipoprotein subclasses"))
}

forestorder <- function(list, sizelevels = c("XXL","XL", "L", "M", "IDL", "S", "XS", "XXS")) {
    data.frame(name = list) %>%
        mutate(level = gsub("-.*", "", name),
               found = level %in% sizelevels,
               level = factor(ifelse(found, level, "custom"),
                              levels = "custom" %union% sizelevels)) %>%
        arrange(level, name) %>%
        pull(name) %>%
        as.character
}

forest.customorder <- function(x) {
    ret <- c("Total cholesterol",
             "Total free cholesterol",
      "VLDL cholesterol",
      "LDL cholesterol",
      "HDL cholesterol",
      "HDL2 cholesterol",
      "HDL3 cholesterol",
      "Esterified-C",
      "Remnant-C",
      "VLDL size",
      "LDL size",
      "HDL size",
      "Total fatty acids",
      "Degree of unsaturation",
      "Omega-3 %",
      "DHA %",
      "Omega-6 %",
      "LA %",
      "PUFA %",
      "MUFA %",
      "SFA %",
      "TG/PG",
      "Total triglycerides",
      "Triglycerides in VLDL",
      "Triglycerides in LDL",
      "Triglycerides in HDL",
      "Total cholines",
      "Alanine",
      "Glutamine",
      "Glycine",
      "Histidine",
      "Isoleucine (branched)",
      "Leucine (branched)",
      "Valine (branched)",
      "Phenylalanine (aromatic)",
      "Tyrosine (aromatic)",
      "Glucose",
      "Pyruvate",
      "Lactate",
      "Citrate",
      "Glycerol",
      forestorder(x)) %>%
        unique %intersect% x %>%
      rev
    message(paste(ret, collapse = ","))
    ret
}

forest.order.and.scale <- function(df, groups, heightadj) {
    df %>%
        filter(nmrgroup %in% groups) %>%
        arrange(match(nmrgroup, groups)) %>%
        mutate(adjustheights = purrr::map2(gg_groups,
                                           rel_heights,
                                           ~myforestplot.equalheights(.x, .y, heightadj))) %>%
        pull(adjustheights) %>%
        append(., list(egg::.dummy_gtable)) %>%
        do.call(gtable_rbind, .)
}

bpscale <- function() {
 scale_colour_manual(name = "Blood pressure variables",
                     labels = c("sys" = "Systole (mmHg)",
                                "log(sys)" = "Systole",
                                       "dias" = "Diastole (mmHg)",
                                       "htn" = "Hypertension",
                                       "sys.diff.pct" = "Percentage change in systole",
                                       "dias.diff.pct" = "Percantage change in diastole",
                                       "sys.followup" = "Systole (mmHg)",
                                "htn.followup" = "Hypertension",
                                "followup" = "Hypertension"),
                     values = c("sys" = "red",
                                "log(sys)" = "red",
                                       "dias" = "blue",
                                       "htn" = "black",
                                       "sys.diff.pct" = "red",
                                       "sys.followup" = "black",
                                       "dias.diff.pct" = "blue",
                                "htn.followup" = "black",
                                "followup" = "black"))
}

agescale <- function() {
     list(scale_shape_manual(name = "Sex",
                       labels = c("sys.young" = "Younger than median age, systole",
                                  "dias.young" = "Younger than median age, diastole",
                                  "sys.old" = "Older than median age, systole",
                                  "dias.old" = "Older than median age, systole",
                                  "htn.young" = "Younger than median age",
                                  "htn.old" = "Older than median age"),
                       values = c("htn.young" = 25,
                                  "htn.old" = 24,
                                  "sys.young" = 25,
                                  "dias.young" = 25,
                                  "sys.old" = 24,
                                  "dias.old" = 24)),
          scale_color_manual(name = "Sex",
                             labels = c("sys.young" = "Younger than median age, systole",
                                  "dias.young" = "Younger than median age, diastole",
                                  "sys.old" = "Older than median age, systole",
                                  "dias.old" = "Older than median age, systole",
                                  "htn.young" = "Younger than median age",
                                  "htn.old" = "Older than median age"),
                             values = c("htn.young" = "black",
                                        "htn.old" = "black",
                                        "sys.young" = "red",
                                        "dias.young" = "blue",
                                        "sys.old" = "red",
                                        "dias.old" = "blue")))
}


sexscale <- function() {
     list(scale_shape_manual(name = "Sex",
                       labels = c("sys.female" = "Women, systole",
                                  "dias.female" = "Women, diastole",
                                  "sys.male" = "Men, systole",
                                  "dias.male" = "Men, systole",
                                  "htn.female" = "Women",
                                  "htn.male" = "Men"),
                       values = c("htn.female" = 24,
                                  "htn.male" = 25,
                                  "sys.female" = 24,
                                  "dias.female" = 24,
                                  "sys.male" = 25,
                                  "dias.male" = 25)),
          scale_color_manual(name = "Sex",
                             labels = c("sys.female" = "Women, systole",
                                  "dias.female" = "Women, diastole",
                                  "sys.male" = "Men, systole",
                                  "dias.male" = "Men, systole",
                                  "htn.female" = "Women",
                                  "htn.male" = "Men"),
                             values = c("htn.female" = "black",
                                        "htn.male" = "black",
                                        "sys.female" = "red",
                                        "dias.female" = "blue",
                                        "sys.male" = "red",
                                        "dias.male" = "blue")))
}

cohortscale <- function() {
    scale_colour_manual(name = "Cohort",
                        labels = c("DILGOM", "HEALTH", "1997", "2000",
                                   "2002", "2007", "2012"),
                        breaks = c("DILGOM", "HEALTH", "1997", "T2000",
                                   "2002", "2007", "2012"),
                        values = c("DILGOM" = "gray60",
                                   "HEALTH" = "black",
                                   "1997" = "gray0",
                                   "T2000" = "gray40",
                                   "2002" = "gray60",
                                   "2007" = "gray80",
                                   "2012" = "gray90"))
}

myforest.nestedplot <- function(data,
                                group,
                                logodds,
                                percentage,
                                scale,
                                xlim,
                                label = "Effect size (95%-CI) per 1-SD metabolite concentration",
                                legend.position = c(-1.15, -3.8),
                                omitlegend,
                                test,
                                panelonelast,
                                paneltwolast) {
    ggforestplot::forestplot(df = data,
                             name = nmrname,
                             estimate = estimate,
                             se = std.error,
                             logodds = logodds,
                             pvalue = qval,
                             title = sprintf("%-50s", group),
                             psignif = 0.05,
                             xlim = xlim,
                             shape = response,
                             colour = response) +
        scale_y_discrete(limits = forest.customorder(data$nmrname)) +
        xlab(ifelse(logodds,
                    myforest.labelformatter.logodds(label),
                    myforest.labelformatter(label, percentage))) +
        { if (logodds) 
              scale_x_log10(breaks = scales::pretty_breaks(n = 4))
          else
              if (percentage)
                  scale_x_continuous(breaks = scales::pretty_breaks(n = 5),
                                     labels = scales::percent_format(accuracy = 1))
          else scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) } +
        scale() +
        guides(colour = guide_legend(),
               shape = guide_legend()) +
        guides(colour = guide_legend(override.aes = list(size = 2))) +
        ggplot2::theme(plot.title = element_text(size = 11,
                                                 face = "italic",
                                                 margin = unit(c(0, 1, 1, 0), "mm")),
                       axis.title.x = element_text(size = 8)) +
        { if (group == paneltwolast && !omitlegend)
              ggplot2::theme(legend.position = legend.position,
                             legend.direction = "vertical",
                             legend.justification = "left",
                             legend.title = element_text(size = 11, face = "italic"))
          else
              ggplot2::theme(legend.position = "none") } +
        { if (!group %in% c(panelonelast, paneltwolast)) 
              ggplot2::theme(axis.text.x = element_blank(),
                             axis.title.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             plot.margin = unit(c(1, 2, 1, 2), "mm"))
        }
}

myforest.labelformatter.logodds <- function(label, width = 35) {
    label %>%
        gsub("Effect size", "Odds ratio", .) %>%
        strwrap(width = width) %>%
        paste0(collapse = "\n")
}

myforest.labelformatter <- function(label, percentage, width = 35) {
    label %>%
        ifelse(percentage, gsub("Effect size", "Percentage change", .), .) %>%
        strwrap(width = width) %>%
        paste0(collapse = "\n")
}

  
myforest.loopgroups <- function(dset,
                     panelonelast,
                     paneltwolast,
                     xlim = c(-1, 3),
                     legend.position = c(-1.15, -3.8),
                     logodds = FALSE,
                     percentage = percentage,
                     scale,
                     omitlegend = FALSE) {
    nest(dset, -nmrgroup, .key = "data") %>%
        mutate(gg_groups =
                   purrr::map2(data,
                               nmrgroup,
                               ~ suppressMessages(myforest.nestedplot(.x,
                                                     .y,
                                                     logodds = logodds,
                                                     percentage = percentage,
                                                     scale = scale,
                                                     xlim = xlim,
                                                     omitlegend = omitlegend,
                                                     legend.position = legend.position,
                                                     panelonelast = panelonelast,
                                                     paneltwolast = paneltwolast))),
               rel_heights =  purrr::map(data,  ~nrow(.)) %>% unlist())
}

myforestplot.equalheights <- function(plot, height, heightadj) {
    g <- gtable_frame(ggplotGrob(plot),
                      height = unit(heightadj*3*height, 'mm'),
                      width = unit(4, 'cm'))
    temp <- g$grobs[[1]]
    g$grobs[[1]] <- g$grobs[[4]]
    g$grobs[[4]] <- temp
    return(g)
}

myforestplot.heightadj <- function(ggplot_multi, height = 1) {
    purrr::map2(ggplot_multi$gg_groups,
                ggplot_multi$rel_heights, 
                function(p,h) {
                    g <- gtable_frame(ggplotGrob(p),
                                      height = unit(height*3*h, 'mm'),
                                      width = unit(4, 'cm'))
                    temp <- g$grobs[[1]]
                    g$grobs[[1]] <- g$grobs[[4]]
                    g$grobs[[4]] <- temp
                    return(g)
                })
}

myoutlier <- function(list) {
    list < quantile(list, 0.25, na.rm = TRUE) - IQR(list, na.rm = TRUE) * 1.5 |
        list > quantile(list, 0.75, na.rm = TRUE) + IQR(list, na.rm = TRUE) * 1.5 |
        (list > quantile(list, 0.40, na.rm = TRUE) &
         list < quantile(list, 0.60, na.rm = TRUE)) # buggy jitterdodge
}

myboxplot <- function(df, vars = c(), rename = FALSE, normalize = FALSE) {
    df %>%
        dplyr::select(vars, htn, cohort) %>%
        { if (normalize) mutate_at(., vars(vars), list(scale)) else .} %>%
        { if (rename) rename_at(., vars(vars), funs(bioproperty(.))) else . } %>%
        gather(covariate, value, -htn, -cohort) %>%
        group_by(covariate, cohort) %>%
        mutate(outlier = myoutlier(value)) %>%
        ungroup %>%
        ggplot(aes(x = covariate, y = value, color = cohort)) +
        ylim(-10, 10) +
        geom_point(data = function(x) dplyr::filter_(x, ~ outlier),
                   position = position_jitterdodge(jitter.width = 0.1),
                   alpha = 0.4, shape = ".") +
        geom_boxplot(lwd=0.1, outlier.shape = NA, notch = TRUE) +
        theme_classic() +
        theme(legend.position = "bottom",
              axis.title.x = element_blank(), 
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1))
}

pca.title <- function(pca, axis = 1, decimals = 1) {
    sprintf("PCA axis %i (%s%%)", axis, round(pca[axis]*100, decimals))
}

pca.scatterplots <- function(pca) {
    pca$axes %>%
        dplyr::select(starts_with("PC"), htn) %>%
        gather(key, value, -PC001, -htn) %>% 
        ggplot2::ggplot(ggplot2::aes_string(x = "PC001", y = "value", col = "htn")) +
        facet_wrap(~key, ncol = 10, scales = "free_y") +
     geom_point(shape = ".") +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 3), lim = c(-10, 10)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
        theme_classic() +
        theme(legend.position = "none")
}

pca.plot <- function(pca, correlation) {
    ggplot(pca$axes, aes(x = PC001, y = PC002)) +
        stat_summary_hex(aes(z = scale(sys, scale = FALSE)),
                         fun = median,
                         binwidth = c(0.33, 0.33)) +
        scale_fill_gradient2(name = "Deviation from\nmedian systolic\nBP in mmHG",
                             low = "darkblue",
                             mid = "white",
                             high = "darkred",
                             midpoint = 0,
                             breaks = c(-40, 0, 40),
                             limits = c(-40, 40)) +
        coord_equal() +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 5),
                           name = pca.title(pca$contribution, axis = 1)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                           name = pca.title(pca$contribution, axis = 2)) +
        theme_classic() +
        theme(legend.background = element_blank(),
              legend.title = element_text(size = 9),
              legend.text = element_text(size = 8),
              legend.position = c(1.15, 0.50),
              legend.key.height = unit(4, "mm"),
              legend.key.width = unit(3, "mm"),
              legend.direction = "vertical",
              plot.margin = unit(c(0, 26, 0, 2), "mm"))
}


pca.loading <- function(pca,
                        axis = "PC1",
                        ylim = c(-0.6, 0.4),
                        tag = list("PC1" = "A.", "PC2" = "B.")) {
    ggplot(pca %>% rename(pca_axis = !!axis),
           aes(x = reorder(nmrname, pca_axis),
               y = pca_axis)) +
        geom_bar(stat = 'identity', fill = "gray80") +
        geom_text(y = 0.005,
                  aes(label = nmrgroup),
                  size = 2,
                  angle = 0,
                  fontface="italic",
                  hjust = 0) +
        ylab(sprintf("PCA loading for %s", axis)) +
        xlab("") +
        coord_flip() +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0),
                           lim = ylim,
                           breaks = scales::pretty_breaks(n = 5)) +
        theme_classic() +
        theme(aspect.ratio = 4,
              legend.title = element_blank(),
              legend.background=element_blank(),
              legend.position=c(0.92, .940),
              legend.text = element_text(size = 6),
              legend.key.size = unit(3, "mm"),
              axis.title.y = element_text(size = 9),
              axis.text.x = element_text(size = 7)) +
         labs(tag = tag[[axis]])
}


diagnosticqqplot <- function(dset, vars, size = 16) {
    dplyr::select(dset, vars) %>%
        rename_all(~bioproperty(.)) %>%
        gather(key, value) %>%
        ggplot2::ggplot(ggplot2::aes(sample = value)) +
        facet_wrap(~key, ncol = 8) +
        stat_qq(size = 0.1) +
        stat_qq_line() +
        theme_classic() +
        theme(legend.position = "none",
              strip.background = element_blank(),
              strip.placement = "bottom",
              strip.text.x = element_text(size = size))
}

pca.kde2d <- function(df, htn = 0, n = 500) {
    x <- pca.df$axes %>% filter(htn == !!htn) %>% pull(PC001)
    y <- pca.df$axes %>% filter(htn == !!htn) %>% pull(PC002)
    xrange <- pca.df$axes %>% pull(PC001) %>% range
    yrange <- pca.df$axes %>% pull(PC002) %>% range
    
    kde2d(x, y, n = n, lims = c(xrange, yrange))
}

mycorrplot <- function(dset, metabolites) {
    corr <- spearmancorrelation(dset, metabolites)

    ggcorrplot(corr = corr,
               hc.order = TRUE,
               type = "full",
               insig = "blank",
               pch = 4,
               pch.col = "gray",
               show.legend = TRUE) +
        scale_fill_gradientn(name = "Spearman\ncorrelation",
                             colors = c("blue", "white", "red"),
                             breaks = c(-1, -0.5, 0, 0.5, 1),
                             limits = c(-1.1, 1.1)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.key.width = unit(4, "mm"),
          legend.key.height = unit(10, "mm"),
          axis.text.x = element_text(angle = 90, hjust = 1))
}

myrocplot <- function(dset, glm.min, glm.max) {
    ggroc(list("min" = pROC::roc(myoutcomes(dset), fitted(glm.min)),
               "max" = pROC::roc(myoutcomes(dset), fitted(glm.max))),
          size = 0.5) +
        geom_abline(intercept = 1) +
        coord_fixed() +
        scale_x_reverse(expand = c(0, 0), breaks = scales::pretty_breaks(n = 3)) +
        scale_y_continuous(expand = c(0, 0), breaks = scales::pretty_breaks(n = 3)) +
        scale_colour_manual(name = "Model variables",
                        labels = c("min" = "Clinical",
                                   "max" = "Clinical and\nmetabolites"),
                        values = c("min" = "blue",
                                   "max" = "red")) +
        xlab("Specifity") +
        ylab("Sensitivity") +
        theme_classic() +
        theme(legend.position = c(0.8, 0.2),
              legend.title = element_blank())
}

plot.splines <- function(...) {
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
            ggtitle(sprintf("%s", bioproperty(term))) +
            theme_bw()
    })
}

