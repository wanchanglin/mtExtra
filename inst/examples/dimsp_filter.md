


# Metabolomics data filtering




## Setting
load R packages 


```r
pkgs <- c("mt", "mtExtra", "tidyverse")
invisible(lapply(pkgs, library, character.only = TRUE))
#> -- Attaching packages ----------------------- tidyverse 1.3.1 --
#> v ggplot2 3.3.5     v purrr   0.3.4
#> v tibble  3.1.5     v dplyr   1.0.7
#> v tidyr   1.1.4     v stringr 1.4.0
#> v readr   2.0.2     v forcats 0.5.1
#> -- Conflicts -------------------------- tidyverse_conflicts() --
#> x dplyr::filter()     masks stats::filter()
#> x dplyr::lag()        masks stats::lag()
#> x stringr::str_trim() masks mtExtra::str_trim()
```

parameter setting


```r
opt <- list(
  #' Input
  peak_file = here::here("extdata", "dims_peak.tsv"),

  #' QC filtering
  qc = TRUE,
  #' wl-01-04-2019, Mon: It is critical to choose rsd threshold. The
  #' distribution of rsd could be a guardian.
  qc_rsd_thres = 60.0,
  qc_mv_filter = TRUE,
  qc_mv_qc_sam = TRUE,
  qc_mv_thres = 0.30,

  #' blank filtering
  bl = TRUE,
  bl_method = "mean",
  bl_factor = 1,
  bl_mv_filter = FALSE,
  bl_mv_thres = 0.30,

  #' MV filtering on samples
  mv = FALSE,
  mv_thres = 0.30
)
```


## Data preparation
Load peak table


```r
peak <- read_tsv(opt$peak_file)
#> Rows: 100 Columns: 15
#> -- Column specification ----------------------------------------
#> Delimiter: "\t"
#> chr  (1): name
#> dbl (14): 01_sample, 02_sample, 03_sample, 04_sample, 05_sample, 06_sample, ...
#> 
#> i Use `spec()` to retrieve the full column specification for this data.
#> i Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

get only numerics


```r
dat <- peak[, -1]
```

get group information for filtering


```r
groups <- names(dat)
groups <- gsub("\\d+_", "", groups, perl = TRUE)
groups <- gsub("pool", "qc", groups, perl = TRUE)
groups
#>  [1] "sample" "sample" "sample" "sample" "sample" "sample" "sample" "sample"
#>  [9] "sample" "sample" "qc"     "qc"     "blank"  "blank"
```

transpose to make rows as replicates and columns as features.


```r
dat <- t(dat)
colnames(dat) <- peak[, 1, drop = TRUE]
```

change zero as NA


```r
dat <- as_tibble(mv.zene(dat))
dat
#> # A tibble: 14 x 100
#>    `TG_6:0_[M+NH4]1+` `TG_6:0_[M+Na]1+` `LPC_2:0_[M+H]1+` `LPC-2O_3:0_[M+H]1+`
#>                 <dbl>             <dbl>             <dbl>                <dbl>
#>  1             23680.              84.9              NA                  456. 
#>  2             18618.             193.               NA                   74.6
#>  3             22791.              NA                NA                   71.7
#>  4             26643.             151.               97.3                161. 
#>  5             22989.              53.2              NA                   NA  
#>  6             23964.             104.               NA                  141. 
#>  7             22829.              NA                NA                   NA  
#>  8             19742.              NA                NA                  107. 
#>  9             22408.              NA                NA                   NA  
#> 10             19530.              NA                NA                   NA  
#> 11             11295.             149.               NA                  301. 
#> 12             11773.              NA                98.5                825. 
#> 13              4716.             232.               NA                   NA  
#> 14             17952.              99.9              NA                  376. 
#> # ... with 96 more variables: LPC_3:1_[M+H]1+ <dbl>, PC-O_3:0_[M+H]1+ <dbl>,
#> #   LPC-2O_4:0_[M+H]1+ <dbl>, Cer_20:1_[M+H-H2O]1+ <dbl>,
#> #   PC-O_4:0_[M+H]1+ <dbl>, LPC_O-5:0_[M+H]1+ <dbl>, PC_4:0_[M+H]1+ <dbl>,
#> #   LPC_5:0_[M+H]1+ <dbl>, Cer_20:1_[M+H]1+ <dbl>, MG_16:0_[M+NH4]1+ <dbl>,
#> #   Cholesterol_[M+H-2H2O]1+ <dbl>, LPC_6:0_[M+H]1+ <dbl>,
#> #   MG_17:0_[M+NH4]1+ <dbl>, Oxocholesterol_[M+H-2H2O]1+ <dbl>,
#> #   Hydroxycholesterol_[M+H-2H2O]1+ <dbl>, Cholesterol_[M+H-H2O]1+ <dbl>, ...
```


## Missing value and RSD checking
combine data matrix with filtering group


```r
data <- as_tibble(cbind(grp = groups, dat))
# data <- bind_cols(grp = groups, dat) 
data
#> # A tibble: 14 x 101
#>    grp    `TG_6:0_[M+NH4]1+` `TG_6:0_[M+Na]1+` `LPC_2:0_[M+H]1~ `LPC-2O_3:0_[M+~
#>    <chr>               <dbl>             <dbl>            <dbl>            <dbl>
#>  1 sample             23680.              84.9             NA              456. 
#>  2 sample             18618.             193.              NA               74.6
#>  3 sample             22791.              NA               NA               71.7
#>  4 sample             26643.             151.              97.3            161. 
#>  5 sample             22989.              53.2             NA               NA  
#>  6 sample             23964.             104.              NA              141. 
#>  7 sample             22829.              NA               NA               NA  
#>  8 sample             19742.              NA               NA              107. 
#>  9 sample             22408.              NA               NA               NA  
#> 10 sample             19530.              NA               NA               NA  
#> 11 qc                 11295.             149.              NA              301. 
#> 12 qc                 11773.              NA               98.5            825. 
#> 13 blank               4716.             232.              NA               NA  
#> 14 blank              17952.              99.9             NA              376. 
#> # ... with 96 more variables: LPC_3:1_[M+H]1+ <dbl>, PC-O_3:0_[M+H]1+ <dbl>,
#> #   LPC-2O_4:0_[M+H]1+ <dbl>, Cer_20:1_[M+H-H2O]1+ <dbl>,
#> #   PC-O_4:0_[M+H]1+ <dbl>, LPC_O-5:0_[M+H]1+ <dbl>, PC_4:0_[M+H]1+ <dbl>,
#> #   LPC_5:0_[M+H]1+ <dbl>, Cer_20:1_[M+H]1+ <dbl>, MG_16:0_[M+NH4]1+ <dbl>,
#> #   Cholesterol_[M+H-2H2O]1+ <dbl>, LPC_6:0_[M+H]1+ <dbl>,
#> #   MG_17:0_[M+NH4]1+ <dbl>, Oxocholesterol_[M+H-2H2O]1+ <dbl>,
#> #   Hydroxycholesterol_[M+H-2H2O]1+ <dbl>, Cholesterol_[M+H-H2O]1+ <dbl>, ...
```

function to select only numeric columns for computering rsd or mv.


```r
df_summ <- function(dat, method = mean, ...) {
  dat %>% select(where(is.numeric)) %>% map_dfr(method) 
}
```

check RSD


```r
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
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
ggplot(val, aes(x = grp, y = value)) + geom_boxplot() + ggtitle("RSD")
```

<img src="figure/unnamed-chunk-10-1.png" title="plot of chunk unnamed-chunk-10" alt="plot of chunk unnamed-chunk-10" style="display: block; margin: auto;" /><img src="figure/unnamed-chunk-10-2.png" title="plot of chunk unnamed-chunk-10" alt="plot of chunk unnamed-chunk-10" style="display: block; margin: auto;" />

check percentage of missing values


```r
val <- data %>% 
  group_by(grp) %>% 
  group_modify(~ df_summ(., method = mv_perc)) %>%
  pivot_longer(cols = !grp) %>%
  filter(!is.na(value))

ggplot(val, aes(x = value)) +
    geom_histogram(colour = "black", fill = "white") +
    facet_grid(grp ~ .) +
    ggtitle("MV")
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
ggplot(val, aes(x = grp, y = value)) + geom_boxplot() + ggtitle("MV")
```

<img src="figure/unnamed-chunk-11-1.png" title="plot of chunk unnamed-chunk-11" alt="plot of chunk unnamed-chunk-11" style="display: block; margin: auto;" /><img src="figure/unnamed-chunk-11-2.png" title="plot of chunk unnamed-chunk-11" alt="plot of chunk unnamed-chunk-11" style="display: block; margin: auto;" />


## Feature filtering
qc filtering


```r
if (opt$qc) {
  dat <- qc_filter(dat, groups, thres_rsd = opt$qc_rsd_thres,
                    f_mv = opt$qc_mv_filter, f_mv_qc_sam = opt$qc_mv_qc_sam,
                    thres_mv = opt$qc_mv_thres)
}
dim(dat)
#> [1] 14 26
```

blank filtering


```r
if (opt$bl) {
  dat <- blank_filter(dat, groups, method = opt$bl_method,
                      factor = opt$bl_factor, f_mv = opt$bl_mv_filter,
                      thres_mv = opt$bl_mv_thres)
}
dim(dat)
#> [1] 14 17
```

wl-28-11-2018, Wed: MV filtering can be done in qc_filter,
 blank_filter or mv_filter. Note that mv filtering on sample must be
 performed once. Otherwise even mv imputation does not work for large
 portion of missing values in some variables.

wl-28-03-2019, Thu: make sure mv filtering must be performed once on
 **sample**


```r
(mv_qc <- (opt$qc) && ("qc" %in% levels(groups)) && (opt$qc_mv_filter) && !(opt$qc_mv_qc_sam))
#> [1] FALSE
(mv_bl <- (opt$bl) && ("blank" %in% levels(groups)) && (opt$bl_mv_filter))
#> [1] FALSE
(mv <- opt$mv || !(mv_qc || mv_bl))
#> [1] TRUE
```

mv filtering


```r
if (mv) {
  idx_sample <- groups %in% "sample"
  idx <- mv_filter(dat[idx_sample, , drop = FALSE], thres = opt$mv_thres)$idx
  dat <- dat[, idx, drop = FALSE]
}
dim(dat)
#> [1] 14 14
dat
#> # A tibble: 14 x 14
#>    `TG_6:0_[M+NH4]1+` `Cer_20:1_[M+H]1+` `Cholesterol_[M+H-H2~ `MG_18:1_[M+NH4]~
#>                 <dbl>              <dbl>                 <dbl>             <dbl>
#>  1             23680.              1431.               762819.            16913.
#>  2             18618.               718.               977877.            12147.
#>  3             22791.               254.              1187192.             1550.
#>  4             26643.               467.              1462344.             2217.
#>  5             22989.               211.              1435597.             2827.
#>  6             23964.               463.              1136028.             3157.
#>  7             22829.               585.              1249338.             8240.
#>  8             19742.              2148.              1767738.             8960.
#>  9             22408.              1381.              1683483.             6319.
#> 10             19530.              1279.              1476509.             2893.
#> 11             11295.              1506.              1927678.            37331.
#> 12             11773.              1203.              1039088.            16054.
#> 13              4716.                NA                 91983.             2153.
#> 14             17952.               512.                 2819.             5025.
#> # ... with 10 more variables: MG_22:6_[M+NH4]1+ <dbl>, V45 <dbl>, V46 <dbl>,
#> #   V47 <dbl>, LPG_12:0_[M+NH4]1+ <dbl>,
#> #   PC-O_13:1_[M+H]1+ / LPE_16:1_[M+H]1+ <dbl>, LPG_13:0_[M+NH4]1+ <dbl>,
#> #   LPG-O_16:0_[M+H]1+ <dbl>, LPE_18:4_[M+H]1+ <dbl>, LPG_14:0_[M+NH4]1+ <dbl>
```

