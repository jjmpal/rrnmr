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

`%uniqunion%` <- function(a, b) {
    all <- c(a, b)
    all[!duplicated(all)]
}

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
    if(percent) sprintf("%.1f%% (%.1f to %.1f%%)",  estimate*100, conf.low*100, conf.high*100)
    else sprintf("%.2f (%.2f to %.2f)", estimate, conf.low, conf.high)
}

standardnames <- function(x, prefix = "NMR_") {
    toupper(paste0(prefix, gsub("[/-]", "_", x)))
}

bioproperty <- function(metabolites, property = "name", length = 26) {
    md <- as.data.frame(ggforestplot::df_NG_biomarker_metadata) %>%
        mutate(description = case_when(description == "Isoleucine" ~ "Isoleucine (branched)",
                                       description == "Leucine" ~ "Leucine (branched)",
                                       description == "Valine" ~ "Valine (branched)",
                                       description == "Phenylalanine" ~ "Phenylalanine (aromatic)",
                                       description == "Tyrosine" ~ "Tyrosine (aromatic)",
                                       TRUE ~ description))
    metabolites %>%
        purrr::map_chr(function(id) {
            id <- gsub("NMR_", "", id)
            names <- md %>% filter(machine_readable_name == id)
            if (property == "name")
                return(descriptivenames(names$abbreviation[1],
                                        names$description[1],
                                        length = length))
            else
                return(names[[property]][1])
        })
}

               

human_readable_names <- function(df, length = 40) {
    ggforestplot::df_NG_biomarker_metadata %>%
        mutate(machine_readable_name = paste0("NMR_", machine_readable_name),
               label = ifelse(nchar(description) > length,
                                    name,
                                    description)) %>%
        left_join(tibble(id = colnames(df)), ., by = c("id" = "machine_readable_name")) %>%
        mutate(label = case_when(id == "sys" ~ "Systolic BP",
                                id == "sys_old" ~ "Baseline systolic BP",
                                id == "dias" ~ "Diastolic BP",
                                id == "age" ~ "Age",
                                id == "bmi" ~ "BMI",
                                id == "leisure" ~ "Exercise",
                                id == "diabetes" ~ "Diabetes mellitus",
                                id == "female" ~ "Female",
                                id == "smoker" ~ "Current smoker",
                                id == "koltreat" ~ "Lipid medication",
                                id == "bptreat" ~ "Antihypertensive medication",
                                id == "Esterified_C" ~ "Total estrefied cholesterol",
                                id == "NMR_Free_C" ~ "Total free cholesterol",
                                TRUE ~ label)) %>%
        select(id, label)
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

getsubclasses <- function(list) {
    data.frame(term = list) %>%
        mutate(group = bioproperty(term, "group")) %>%
        filter(group == "Lipoprotein subclasses") %>%
        pull(term)
}

mkdir <- function(...) {
    vmkdir <- Vectorize(dir.create, vectorize.args = "path")
    vmkdir(path = c(...), showWarnings = FALSE) %>%
        invisible
}

rmdir <- function(...) {
    vrmdir <- Vectorize(unlink, vectorize.args = "x")
    vrmdir(x = c(...), recursive=TRUE) %>%
        invisible
}

submodel <- function(dset, filter, newname, oldname) {
    dset %>%
        dplyr::filter(!!enquo(filter)) %>%
        dplyr::rename(!!enquo(newname) := !!enquo(oldname))
}

mysavefactory <- function(arguments) {
    force(arguments)
    function(object, name, path = "rds/") {
        filename <- sprintf("%s%s_%s.rds", path, name, argstostring(arguments))
        saveRDS(object, file = filename)
  }
}

myloadfactory <- function(arguments) {
    force(arguments)
    function(name, path = "rds/") {
        filename <- sprintf("%s%s_%s.rds", path, name, argstostring(arguments))
        message(filename)
        if (file.exists(filename))
            return(readRDS(file = filename))
        else
            return(NULL)
  }
}

## myload <- function(name, path = "rds/") {
##     files <- list.files(path = "rds", pattern = paste0("^", name), full.names = TRUE)
##     i <- menu(files, title="Select file to load?")
##     message("Loading file ", files[[i]])
##     readRDS(files[[i]])
## }

closest_to_point <- function(df, point = 2.1) {
    if (point > 0)
        df %>% filter(x <= point + 0.001) %>% tail(n=1)
    else
        df %>% filter(x >= point - 0.001) %>% head(n=1)
}


normalize_plots <- function(...) {
    plotlist <- list(...)
    widths_list <- c()
    for (i in seq_along(plotlist)) {
        widths_list[[i]] <- ggplotGrob(plotlist[[i]])[["widths"]][2:5]
    }
    
    maxWidth <- do.call(grid::unit.pmax, widths_list)

    for (i in seq_along(plotlist)) {
        g <- ggplotGrob(plotlist[[i]])
        g[["widths"]][2:5] <- as.list(maxWidth)
        plotlist[[i]] <- g
    }
    return(plotlist)
}

list_get <- function(..., list = list) {
    lapply(c(...), function(x) list[[x]])
}

myggsavefactory <- function(datetime) {
    force(datetime)

    function(plot, file, draw = TRUE, ...) {
        filename_with_time <- sprintf("%s-%s.%s",
                                      tools::file_path_sans_ext(file),
                                      datetime,
                                      tools::file_ext(file))
        ggplot2::ggsave(filename = filename_with_time, plot = plot, ...)
        if (knitr::is_html_output())
        cat("![](", filename_with_time, ")")
    }
}

argstostring <- function(list) {
    list[order(names(list))] %>%
        purrr::list_modify(help = NULL) %>%
        paste(names(.), ., sep="=",collapse="_" )
}
