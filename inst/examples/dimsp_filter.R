#' ---
#' title: "Metabolomics data filtering"
#' author: "Wanchang Lin"
#' date: "`r Sys.Date()`"
#' output:
#'   BiocStyle::html_document:
#'     toc_depth: 3
#'     number_section: false
#'     toc_float: false
#'   BiocStyle::pdf_document:
#'     keep_tex: true
#'     toc: true
#'     toc_depth: 3
#'     number_section: false
#'     citation_package: natbib
#'     latex_engine: xelatex
#' always_allow_html: true
#' geometry: margin=1in
#' fontsize: 11pt
#' ---

#' <!--
#' # Metabolomics data filtering
#' > 23-11-2021, Tue: use mtExtra. do not split data. so no merge. <br/>
#' > 24-11-2021, Wed: test again and apply: <br/> 
#' >   Rscript -e 'knitr::spin("dimsp_filter.R")' <br/>
#' > 25-11-2021, Thu: Use a small data set. <br/>
#' > wl-01-08-2024, Thu: Use this
#' > Rscript -e "rmarkdown::render('dimsp_filter.R', BiocStyle::html_document())"
#' -->

#+ common, include=F
rm(list = ls(all = TRUE))
options(help_type = "html")
source("_common.R")
pkgs <- c("mt", "mtExtra", "tidyverse")
invisible(lapply(pkgs, library, character.only = TRUE))

## ---- Parameter Setting ----
#' ## Parameter Setting

#' parameter setting: put into a list
opt <- list(
  ## Input
  peak_file = here::here("extdata", "dims_peak.tsv"),

  ## QC filtering
  qc = TRUE,
  ## wl-01-04-2019, Mon: It is critical to choose rsd threshold. The
  ## distribution of rsd could be a guardian.
  qc_rsd_thres = 60.0,
  qc_mv_filter = TRUE,
  qc_mv_qc_sam = TRUE,
  qc_mv_thres = 0.30,

  ## blank filtering
  bl = TRUE,
  bl_method = "mean",
  bl_factor = 1,
  bl_mv_filter = FALSE,
  bl_mv_thres = 0.30,

  ## MV filtering on samples
  mv = FALSE,
  mv_thres = 0.30
)

## ---- Data preparation ----
#' ## Data preparation

#' Load peak table
peak <- read_tsv(opt$peak_file)

#' get only numerics
dat <- peak[, -1]

#' get group information for filtering
groups <- names(dat)
groups <- gsub("\\d+_", "", groups, perl = TRUE)
groups <- gsub("pool", "qc", groups, perl = TRUE)
groups

#' transpose to make rows as replicates and columns as features.
dat <- t(dat)
colnames(dat) <- peak[, 1, drop = TRUE]

#' change zero as NA
dat <- as_tibble(mv.zene(dat))
dat

## ---- Missing value and RSD checking ----
#' ## Missing value and RSD checking

#' combine data matrix with filtering group
data <- as_tibble(cbind(grp = groups, dat))
# data <- bind_cols(grp = groups, dat) 
data

#' function to select only numeric columns for computering rsd or mv.
df_summ <- function(dat, method = mean, ...) {
  dat %>% select(where(is.numeric)) %>% map_dfr(method) 
}

#' check RSD
vec_method <- rsd    
val <- data %>% group_by(grp) %>% group_modify(~ df_summ(., method = vec_method))
## val <- data %>% group_by(grp) %>% do(df_summ(., method = vec_method))
## val <- plyr::ddply(data, plyr::.(grp), function(x) {
##   x <- x[, -1, drop = FALSE]
##   rsd(x)
##   # mv_perc(x)
## }) %>% as_tibble()

val <- pivot_longer(val, cols = !grp) %>% filter(!is.na(value))

ggplot(val, aes(x = value)) +
    # geom_histogram(binwidth = .5, colour = "black", fill = "white") +
    geom_histogram(colour = "black", fill = "white") +
    facet_grid(grp ~ .) +
    ggtitle("RSD")
ggplot(val, aes(x = grp, y = value)) + geom_boxplot() + ggtitle("RSD")

#' check percentage of missing values
val <- data %>% 
  group_by(grp) %>% 
  group_modify(~ df_summ(., method = mv_perc)) %>%
  pivot_longer(cols = !grp) %>%
  filter(!is.na(value))

ggplot(val, aes(x = value)) +
    geom_histogram(colour = "black", fill = "white") +
    facet_grid(grp ~ .) +
    ggtitle("MV")
ggplot(val, aes(x = grp, y = value)) + geom_boxplot() + ggtitle("MV")

## ---- Feature filtering ----
#' ## Feature filtering

#' qc filtering
if (opt$qc) {
  dat <- qc_filter(dat, groups, thres_rsd = opt$qc_rsd_thres,
                   f_mv = opt$qc_mv_filter, f_mv_qc_sam = opt$qc_mv_qc_sam,
                   thres_mv = opt$qc_mv_thres)
}
dim(dat)

#' blank filtering
if (opt$bl) {
  dat <- blank_filter(dat, groups, method = opt$bl_method,
                      factor = opt$bl_factor, f_mv = opt$bl_mv_filter,
                      thres_mv = opt$bl_mv_thres)
}
dim(dat)

#' wl-28-11-2018, Wed: MV filtering can be done in qc_filter,
#'  blank_filter or mv_filter. Note that mv filtering on sample must be
#'  performed once. Otherwise even mv imputation does not work for large
#'  portion of missing values in some variables.
#'
#' wl-28-03-2019, Thu: make sure mv filtering must be performed once on
#'  **sample**
(mv_qc <- (opt$qc) && ("qc" %in% levels(groups)) && (opt$qc_mv_filter) && !(opt$qc_mv_qc_sam))
(mv_bl <- (opt$bl) && ("blank" %in% levels(groups)) && (opt$bl_mv_filter))
(mv <- opt$mv || !(mv_qc || mv_bl))

#' mv filtering
if (mv) {
  idx_sample <- groups %in% "sample"
  idx <- mv_filter(dat[idx_sample, , drop = FALSE], thres = opt$mv_thres)$idx
  dat <- dat[, idx, drop = FALSE]
}
dim(dat)
dat
