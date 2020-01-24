forest.args <- function(list, range) {
    selected <- c(lapply(range, function(x) list[[x]]), list(egg::.dummy_gtable))
    do.call(gtable_rbind, selected)
}

merge.metadata <- function(dset, select = c("nmrid", "nmrname", "nmrgroup"), term = "term") {
    ggforestplot::df_NG_biomarker_metadata %>%
        mutate(nmrid = lapply(c(alternative_names), function(x) standardnames(x)),
               nmrname = ifelse(nchar(description) > 25, abbreviation, description),
               nmrgroup = group) %>%
        dplyr::select_(.dots = select) %>%
        unnest(nmrid) %>%
        distinct() %>%
        left_join(dset, ., by = setNames("nmrid", term)) %>%
        arrange(nmrgroup)
}

bpscale <- function() {
 scale_colour_manual(name = "Blood pressure variables",
                            labels = c("sys" = "Systole",
                                       "dias" = "Diastole",
                                       "htn" = "Hypertension",
                                       "htn.followup" = "Hypertension",
                                       "young" = "Young",
                                       "old" = "Old"),
                            values = c("sys" = "blue",
                                       "dias" = "red",
                                       "htn" = "black",
                                       "htn.followup" = "black",
                                       "young" = "red",
                                       "old" = "blue"))
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
 
foresttheme <- function(g, logodds = FALSE, colorscale = bpscale, omitlegend = FALSE) {
    g +
        { if (logodds) 
              scale_x_log10(breaks = scales::pretty_breaks(n = 4))
          else
              scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) } +
        colorscale() +
        guides(colour = guide_legend(override.aes = list(size=4, shape = 19))) +
        ggplot2::theme(legend.position = c(-0.5, -2.5),
                       legend.direction = "vertical",
                       legend.justification = "left",
                       plot.title = element_text(size = 11,
                                                    face = "italic",
                                                    margin = unit(c(0, 1, 1, 0), "mm")),
                       axis.title.x = element_text(size = 8)) +
    { if (omitlegend) ggplot2::theme(legend.position = "none") }
}

forestnolegend <- function() {
    ggplot2::theme(legend.position = "none")
}

forestremovelabels <- function() {
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(1, 2, 1, 2), "mm"))
}

nestplot <- function(dset,
                     first,
                     second,
                     xlim = c(-1, 3),
                     logodds = FALSE,
                     label = "Effect size (95%-CI) per 1-SD\nmetabolite concentration",
                     colorscale = bpscale,
                     omitlegend = FALSE) {
    nest(dset, -nmrgroup, .key = "data") %>%
        mutate(gg_groups = purrr::map2(data,
                                       nmrgroup,
                                       ~ (ggforestplot::forestplot(
                                                          df = .x,
                                                          name = nmrname,
                                                          estimate = estimate,
                                                          se = std.error,
                                                          logodds = logodds,
                                                          pvalue = qval,
                                                          title = sprintf("%-40s", .y),
                                                          psignif = 0.05,
                                                          colour = response,
                                                          xlim = xlim) +
                                          xlab(label)) %>%
                                           foresttheme(.,
                                                       logodds = logodds,
                                                       colorscale = colorscale,
                                                       omitlegend = omitlegend)),
               gg_groups = case_when(row_number() == second[length(second)] ~
                                         gg_groups,
                                     row_number() == first[length(first)] ~
                                         purrr::map(gg_groups, ~ . + forestnolegend()),
                                     TRUE ~
                                         purrr::map(gg_groups, ~ . + forestremovelabels())),
               rel_heights =  purrr::map(data,  ~ nrow(.)) %>% unlist())
}

myforestplot <- function(dset,
                         responses,
                         logodds = FALSE,
                         xlim = c(-1, 3),
                         first = c(3, 10, 2, 4, 8, 7),
                         second = c(6, 1, 9, 5),
                         file,
                         dims = c(18, 24),
                         colorscale = bpscale) {
    stopifnot(!missing(responses), !missing(file))
    df.metadata <- dset %>%
        mutate(nmrname = bionames(term),
               nmrgroup = biogroups(term)) %>%
        filter(response %in% responses, !grepl("%", nmrname)) %>%
        arrange(nmrgroup)

    ggplot_multi <- nestplot(df.metadata,
                             first,
                             second,
                             logodds = logodds,
                             xlim = xlim,
                             colorscale = colorscale,
                             omitlegend = length(responses) == 1)


    height <- 2/length(responses)
    
    lg <- purrr::map2(ggplot_multi$gg_groups,
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

    ag <- arrangeGrob(forest.args(lg, first),
                      forest.args(lg, second),
                      ncol = 2)

    ggsave(file = file, plot = ag, width = dims[1], height = dims[2], dpi = 300, unit = "cm")
}



myoutlier <- function(list) {
    list < quantile(list, 0.25) - IQR(list) * 1.5 |
        list > quantile(list, 0.75) + IQR(list) * 1.5 |
        (list > quantile(list, 0.40) & list < quantile(list, 0.60)) # buggy jitterdodge
}

myboxplot <- function(df, vars = c(), rename = FALSE, normalize = FALSE) {
    df %>%
        dplyr::select(vars, htn, cohort) %>%
        { if (normalize) mutate_at(., vars(vars), list(scale)) else .} %>%
        { if (rename) rename_at(., vars(vars), funs(bionames(.))) else . } %>%
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
        theme(legend.position = "none")
}

pca.plot <- function(pca, correlation) {
    ggplot(pca$axes, aes(x = PC001, y = PC002)) +
    stat_summary_hex(aes(z = scale(sys, scale = FALSE)), fun = median, binwidth = c(1.0, 1.0)) +
    scale_fill_gradient2(name = "Deviation from\nmedian systolic\nBP in mmHG",
                         low = "darkblue",
                         mid = "white",
                         high = "darkred",
                         midpoint = 0) +
    coord_equal() +
    xlab(pca.title(pca$contribution, axis = 1)) +
        ylab(pca.title(pca$contribution, axis = 2)) +
    theme_classic() +
    theme(legend.background=element_blank(),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 7),
          legend.position = c(0.97, 0.82),
          legend.key.height = unit(4, "mm"),
          legend.key.width = unit(3, "mm"),
          legend.direction = "vertical")
}

pca.loading <- function(pca,
                        axis = "PC1",
                        ylim = c(-0.4, 0.4),
                        tag = list("PC1" = "A.", "PC2" = "B.")) {
    ggplot(pca %>% rename(pca_axis = !!axis),
           aes(x = reorder(nmrname, pca_axis),
               y = pca_axis)) +
        geom_bar(stat = 'identity', fill = "pink") +
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
                           breaks = seq(min(ylim), max(ylim), 0.2)) +
        theme_classic() +
        theme(legend.title = element_blank(),
              legend.background=element_blank(),
              legend.position=c(0.92, .940),
              legend.text = element_text(size = 6),
              legend.key.size = unit(3, "mm"),
              axis.title.y = element_text(size = 9),
              axis.text.x = element_text(size = 7)) +
         labs(tag = tag[[axis]])
}


diagnosticqqplot <- function(dset, vars) {
    dplyr::select(dset, vars) %>%
        rename_all(~bionames(.)) %>%
        gather(key, value) %>%
        ggplot2::ggplot(ggplot2::aes(sample = value)) +
        facet_wrap(~key, ncol = 8) +
        stat_qq(size = 0.1) +
        stat_qq_line() +
        theme(legend.position = "none")
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
