## wl-09-11-2021, Tue: gather all general functions from 2015

## ------------------------------------------------------------------------
#' Filtering variable based on local false discovery rate
#' 
#' Filter data based on local false discovery rate. This function uses `z.2`
#' of `locfdr`.
#' 
#' @param x a data metirx
#' @param plot an interger for plotting. 0 gives no plots. 1 gives single 
#'   plot showing the histogram of zz and fitted densities f and p0*f0.
#' @param thres a user defined threshold for filtering. The default is NULL, 
#' which use local FDR as threshold for filtering.
#' @param ... other parameters to be passed to `locfdr`. 
#' @return a list of with contents: \itemize{
#'  \item dat the filtered data matrix.
#'  \item idx a vector of filtering index.
#'  \item thres threshold used for filtering.
#' }
#' 
#' @details 
#' From R package `locfdr` vignette: 
#'   `z.2` is the interval along the zz-axis outside of which `fdr(z) < 0.2`,
#'   the locations of the yellow triangles in the histogram plot. If no
#'   elements of `zz` on the left or right satisfy the criterion, the
#'   corresponding element of z.2 is NA, and the corresponding triangle does
#'   not appear.
#' @seealso [locfdr()] 
#' @importFrom locfdr locfdr
#' @export 
## wl-09-12-2020, Wed: Filter data based on z.2 of 'locfdr'
## wl-16-12-2020, Wed: add fixed threshold in case failurer of 'locfdr'
## wl-26-11-2021, Fri: Should check potential multiple arguments for `plot`
##  and `nulltype`.
locfdr_filter <- function(x, plot = 1, thres = NULL, ...) {
  if (!is.null(thres)) {
    if (length(thres) == 1 && thres <= 0) {
      stop("Single threshold value must be larger than zero")
    }
    if ((length(thres) > 1) && (thres[1] >= thres[2])) {
      stop("The threshold range should be [lower, upper]")
    }
    if (length(thres) == 1) {
      thres <- c(-thres, thres)
    } else {
      thres <- thres[1:2]
    }
  } else {
    vec <- as.vector(as.matrix(x))
    vec <- vec[!is.na(vec)]

    fdr <- locfdr::locfdr(vec, nulltype = 1, plot = plot, ...)
    thres <- fdr$z.2
  }

  ## Filtering function: Keep vector having at least one significant elements.
  sig_cut <- function(x, lower, upper) {
    if (sum(x <= lower) > 0 | sum(x >= upper) > 0) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }

  ## filter data
  idx <- apply(x, 2, sig_cut, thres[1], thres[2])
  x <- x[, idx, drop = F]

  res <- list(dat = x, idx = idx, thres = thres)
  return(res)
}

## ------------------------------------------------------------------------
#' Symbilise data set
#' 
#' Symbilise data values as (-1, 0, 1) based on a range
#' 
#' @param x an vector or data matrix.
#' @param thres an vector with lower and higher values.
#' @return a symbolic data set.
#' @export 
## wl-10-12-2020, Thu: convert data with (-1, 0, 1) based on threshold.  
dat_symb <- function(x, thres) {
  
  if (length(thres) == 1 && thres <= 0) {
    stop("Single threshold value must be larger than zero")
  }
  if ((length(thres) > 1) && (thres[1] >= thres[2])) {
    stop("The threshold range should be [lower, upper]")
  }
  if (length(thres) == 1) {
    thres <- c(-thres, thres)
  } else {
    thres <- thres[1:2]
  }

  ## wl-26-11-2021, Fri: same resukts but the second one is simple and clear
  if (F) {
    data_l <- -1 * (x <= thres[1]) # left tail
    data_u <- 1 * (x >= thres[2])  # right tail
    x <- data_l + data_u
  } else {
    x[(x > thres[1]) & (x < thres[2])] <- 0
    x[x >= thres[2]] <- 1
    x[x <= thres[1]] <- -1
  }

  return(x)
}

## ------------------------------------------------------------------------
#' Transform data
#' 
#' Perform data set transformation
#'
#' @param x a matrix, data frame or vector.
#' @param method transformation method, including: "center", "auto", "range",
#'   "pareto", "vast", "level", "log", "log10", "sqrt" and "asinh". 
#' @param add a shift value for log transformation.
#' @return transformed data
#' @references 
#'  Berg, R., Hoefsloot, H., Westerhuis, J., Smilde, A. and Werf, M. (2006),
#'  Centering, scaling, and transformations: improving the biological
#'  information content of metabolomics data, \emph{BMC Genomics}, 7:142
#' @name trans
#' @examples
#'  data(iris)
#'  ## transform an vector
#'  vec <- iris[, 1]
#'  dat_trans(vec, method = "auto")
#'  vec_trans(vec, method = "auto")
#'  ## transform a data frame
#'  mat <- iris[, 1:4]
#'  dat_trans(mat, method = "log")
#'  ## transform data frame under condition
#'  plyr::ddply(iris, ("Species"), function(x, method) {
#'    dat_trans(x[, 1:4], method = method)
#'  }, method = "range")
NULL

## ------------------------------------------------------------------------
#' @export 
#' @rdname trans 
#' @order 1
## wll-29-07-2015: Transform a data matrix
## Note: The programming structure is from ststs's 'sd' function.
dat_trans <- function(x, method = "auto", add = 1) {
  if (is.matrix(x)) {
    res <- apply(x, 2, vec_trans, method = method, add = add)
  } else if (is.vector(x)) {
    res <- vec_trans(x, method = method, add = add)
  } else if (is.data.frame(x)) {
    res <- sapply(x, vec_trans, method = method, add = add)
  } else {
    res <- vec_trans(as.vector(x), method = method, add = add)
  }

  return(res)
}

## ------------------------------------------------------------------------
#' @export 
#' @rdname trans
#' @order 2
## wll-29-07-2015: transform an vector
vec_trans <- function(x, method = "auto", add = 1) {
  method <- match.arg(method, c(
    "center", "auto", "range", "pareto", "vast",
    "level", "log", "log10", "sqrt", "asinh"
  ))

  ## basic functions
  me <- function(x) mean(x, na.rm = T)
  se <- function(x) sd(x, na.rm = T)
  mx <- function(x) max(x, na.rm = T)
  mn <- function(x) min(x, na.rm = T)
  ## sm  <- function(x) sum(x,na.rm=T)

  res <- switch(method,
    "center"  = (x - me(x)),
    "auto"    = (x - me(x)) / se(x),
    "range"   = (x - me(x)) / (mx(x) - mn(x)),
    "pareto"  = (x - me(x)) / sqrt(se(x)),
    "vast"    = (x - me(x)) * me(x) / se(x)^2,
    "level"   = (x - me(x)) / me(x),
    "log"     = log(x + add),
    "log10"   = log10(x + add),
    "sqrt"    = sqrt(x),
    "asinh"   = asinh(x)
  )
  return(res)
}

## ----------------------------------------------------------------------- 
#' MS/NMR data filtering based on RSD of "qc" data
#' 
#' MS/NMR data filtering based on RSD of "qc" data
#' 
#' @param x a data matrix.
#' @param y a character string with contents of "sample", "qc" and "blank".
#' @param thres_rsd threshold of RSD on QC. Features less than this
#'   threshold will be kept.
#' @param f_mv  a flag indicating whether or not to performance missing value
#'   filtering on either "sample" or "qc" data.
#' @param f_mv_qc_sam a flag for filtering using percentage of missing values
#'  on "qc" or "sample". TRUE is for "qc".
#' @param thres_mv threshold of missing values. Features less than this
#'  threshold will be kept.
#' @return a filtered data.
#' @details This filter process takes two steps. First, the missing values 
#'  filterting is performed on either "qc" or "sample". Then RSD-based 
#'  filtering is applied to "qc" data.
#' @family variable filters
#' @export 
## wl-06-11-2018, Tue: Feature filtering based on QC's RSD
## wl-14-11-2018, Wed: add flag to missing value filtering
qc_filter <- function(x, y, thres_rsd = 20, f_mv = TRUE,
                      f_mv_qc_sam = FALSE, thres_mv = 0.30) {
  ## 1) filtering based on missing values: sample or qc.
  if (f_mv) {
    if (f_mv_qc_sam) {
      tmp <- y %in% "qc" 
    } else {
      tmp <- y %in% "sample" 
    }
    idx <- mv_filter(x[tmp, , drop = FALSE], thres = thres_mv)$idx
    x <- x[, idx, drop = FALSE]
  }

  ## 2) filtering based rsd of "qc"
  tmp <- y %in% "qc" 
  idx <- rsd_filter(x[tmp, , drop = FALSE], thres = thres_rsd)$idx
  x <- x[, idx, drop = FALSE]

  return(x)
}

## ----------------------------------------------------------------------- 
#' MS/NMR data filtering based on "blank" data
#' 
#' MS/NMR data filtering based on "blank" data
#' 
#' @param x a data matrix.
#' @param y a character string with contents of "sample", "qc" and "blank".
#' @param method method for stats. Support "mean", "median" and "max".
#' @param factor multiplier for blank stats
#' @param f_mv a flag indicating whether perform missing value filtering on
#'   "sample" data.
#' @param thres_mv threshold of missing values on QC. Features less than this
#'    threshold will be kept.
#' @return  a filtered data. 
#' @details This function provides an option to perform missing value filtering
#'   on "sample" data. 
#' @family variable filters
#' @export 
## wl-06-11-2018, Tue: Feature filtering based on blank
## wl-14-11-2018, Wed: change order of missing value filtering
blank_filter <- function(x, y, method = c("mean", "median", "max"),
                         factor = 1, f_mv = TRUE, thres_mv = 0.30) {
  method <- match.arg(method)
  idx_sample <- y %in% "sample" 
  idx_blank  <- y %in% "blank"

  ## 1) filtering based on missing values of "sample".
  if (f_mv) {
    idx <- mv_filter(x[idx_sample, , drop = FALSE], thres = thres_mv)$idx
    x <- x[, idx, drop = FALSE]
  }

  ## 2) filtering based on characteristics of blank intensities: mean, median
  ##    or max
  stats.blank <- apply(x[idx_blank, , drop = FALSE], 2, method, na.rm = TRUE)
  stats.blank <- factor * stats.blank
  stats.sample <- apply(x[idx_sample, , drop = FALSE], 2, method, na.rm = TRUE)

  ## keep features with sample stats are larger than blank
  idx <- stats.sample >= stats.blank
  idx[is.na(idx)] <- FALSE
  ## Also keep features whose values are NA in blank
  idx.1 <- is.na(stats.blank)
  ## take union
  idx <- idx | idx.1
  
  ## update data set
  x <- x[, idx, drop = FALSE]

  return(x)
}

## ------------------------------------------------------------------------
#' Filtering variable based on variability
#'
#' Perform variable filtering based on feature variability measurement such as
#' standard deviation(SD) and interquantile range(IQR).
#'
#' @param x a matrix or data frame.
#' @param method variability measurement method, such as `IQR` and `sd`.
#' @param na.rm	a logical value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param thres a numeric value between 0 and 1 for the threshold of
#'  `quantile`. Features whose variability value is large than this threshold
#'  will be kept.
#' @return a list of with contents: \itemize{
#'  \item dat the filtered data matrix
#'  \item idx a vector of filtering index.
#' }
#' @examples
#' set.seed(100)
#' x <- matrix(rnorm(20 * 10), ncol = 10)
#' res <- var_filter(x, method = "sd", thres = 0.25)
#' sum(res$idx)
#' @family variable filters
#' @export
## wl-04-10-2016, Tue: Variable filtering based on IQR.
## wl-19-01-2017, Thu: fix a bug. Use drop=F
## wl-15-10-2020, Thu: return index as well
## wl-22-11-2021, Mon: change function name and modify
var_filter <- function(x, method = "IQR", na.rm = FALSE, thres = 0.25) {

  method <-
    if (is.function(method)) {
      method
    } else if (is.character(method)) {
      get(method)
    } else {
      eval(method)
    }

  vars <- apply(x, 2, method, na.rm = na.rm)

  if (0 < thres && thres < 1) {
    quant <- quantile(vars, probs = thres)
    selected <- !is.na(vars) & vars > quant
  } else {
    stop("Quantile threshold has to be between 0 and 1.")
  }

  x <- x[, selected, drop = F]
  return(list(dat = x, idx = selected))
}

## -----------------------------------------------------------------------
#' Filtering variable based on the percentage of missing values
#'
#' This function calculates the percentage of missing values and keeps those
#' features with missing values percentage less than the designed threshold.
#'
#' @param x a data matrix. The columns are features.
#' @param thres threshold of missing values. Features less than this
#'   threshold will be kept. Value has to be between 0 and 1.
#' @return a list of with contents: \itemize{
#'  \item dat the filtered data matrix
#'  \item idx a logical vector of index for keeping features.
#' }
#' @family variable filters
#' @export
## wl-14-06-2011: Filter features based on the percentage of missing values
## wl-17-06-2021, Thu: several version but this one is simple. Need to test
## wl-06-11-2018, Tue: feature filter index based on missing values
mv_filter <- function(x, thres = 0.3) {

  if (!(is.matrix(x) || is.data.frame(x))) {
    stop("\n'data' is not matrix or data frame.\n")
  }

  thres <- thres * nrow(x)

  ## get number of Na, NaN and zero in each of feature/variable
  count <- apply(x, 2, function(y) {
    tmp <- sum(is.na(y) | is.nan(y) | (y == 0))
  })
  ## na.mat <- is.na(data)
  ## count  <- apply(na.mat,2,sum)

  ## index of features whose number of MVs are less than threshold
  idx <- count <= thres

  x <- x[, idx, drop = F]
  return(list(dat = x, idx = idx))
}

## -----------------------------------------------------------------------
#' Filtering variable based on relative standard deviation (RSD)
#' 
#' Filtering variable based on relative standard deviation (RSD).
#' 
#' @param x a data frame where columns are features.
#' @param thres threshold of RSD. Features less than this threshold will be
#'    kept.
#' @return a list of with contents: \itemize{
#'  \item dat the filtered data matrix
#'  \item idx a logical vector of index for keeping features.
#' }
#' @family variable filters
#' @export 
## wl-06-11-2018, Tue: feature filter index based on RSD
rsd_filter <- function(x, thres = 20) {
  res <- rsd(x)
  idx <- res < thres
  idx[is.na(idx)] <- FALSE

  ## summary(res)
  ## tmp <- hist(res, plot = F)$counts
  ## hist(res,
  ##   xlab = "rsd", ylab = "Counts", col = "lightblue",
  ##   ylim = c(0, max(tmp) + 10)
  ## )

  x <- x[, idx, drop = F]
  return(list(dat = x, idx = idx))
}

## ------------------------------------------------------------------------
#' Select random sample with stratification
#' 
#' Select random sample with stratification from a binary group. 
#' 
#' @param grp a character string indicating two group information.
#' @param len  length of random selection.
#' @param strat a logical value indicated the sampling should be stratified.
#' @return returns a selected index.
#' @details multiple group is not supported at the current stage.
#' @examples 
#' cls <- iris[, 5, drop = TRUE]
#' cls <- cls[cls == "setosa" | cls == "versicolor"]
#' cls <- droplevels(cls)
#' ind <- samp_strat(cls, 6, strat = TRUE)
#' cls[ind]
#' ind <- samp_strat(cls, 5, strat = TRUE)
#' cls[ind]
#' @export 
## lwc-05-10-2010: randomly sample with stratification.
samp_strat <- function(grp, len, strat = TRUE) {
  ## to-do: should check the validity of arguments.
  grp <- droplevels(factor(grp))
  n <- length(grp)
  if (strat) {
    ## idx <- sample(1:n,n,replace=F)  #' shuffle
    ## grp <- grp[idx]

    idx.g <- lapply(levels(grp), function(x) which(grp == x))
    ## idx.g <- lapply(levels(grp), function(x) grep(x,grp, fixed=T))

    ratio <- table(grp)[1]
    ratio <- ratio / n
    ## only consider two class. should be more general.
    len.1 <- trunc(len * ratio) 
    len.2 <- len - len.1

    ind <- mapply(function(x, y) {
      sample(x, y, replace = F)
    }, idx.g, c(len.1, len.2))
    
    ## wl-19-11-2021, Fri: fix a bug
    if (is.matrix(ind)) {
      ind <- as.vector(ind)
    } else {
     ind <- do.call(c, ind)
    }
    ind <- sample(ind, length(ind), replace = F) # shuffle
  } else {
    ind <- sample(1:n, len, replace = F)
  }
  return(ind)
}

## ------------------------------------------------------------------------
#' Select random samples from each group
#'  
#' Select random samples from each group.
#' 
#' @param x  a character string or factor indicating class info.
#' @param k  number of samples selected from each class.
#' @param n  number of replicates.
#' @return returns a selected index.
#' @export 
#' @examples  
#'  cls   <- iris[, 5, drop = TRUE]
#'  (tmp <- samp_sub(cls,k = 6, n = 2))
#'  cls[tmp[[1]]]
#'  table(cls[tmp[[1]]])
## lwc-02-02-2007: randomly select k samples from each factor/class
## NOTE: It should extend to chose different number of samples for different
##       classes.
samp_sub <- function(x, k, n = 10) {
  ## if (!is.factor(y)) stop("y is not of class factor")
  x <- droplevels(as.factor(x))
  idx <- 1:length(x)
  g <- levels(x)
  ng <- length(g)

  nidx <- list()
  for (j in 1:ng) {
    nidx <- c(nidx, list(idx[which(x == g[j])]))
  }

  out.idx <- list()
  for (i in 1:n) {
    kidx <- c()
    for (j in 1:ng) {
      kidx <- c(kidx, sample(nidx[[j]], k))
    }
    kidx <- sample(kidx) # shuffling
    out.idx <- c(out.idx, list(kidx))
  }
  return(out.idx)
}

## ------------------------------------------------------------------------
#' Univariate outlier detection
#' 
#' Perform outlier detection using univariate method.
#' 
#' @param x a numeric vector.
#' @param method method for univariate outlier detection. Only `boxplot`, 
#'   `mean` and `median` are supported.
#' @return returns a logical vector.
#' @references 
#'   Wilcox R R, Fundamentals of Modern Statistical Methods: Substantially 
#'   Improving Power and Accuracy, Springer 2010 (2nd edition), pages 31-35.
#' @details
#'  - `mean`: the absolute difference between that observation and the sample 
#'    mean is more than 2 times of SD.
#'  - `median`: the absolute difference between the observation and the sample
#'    median is larger than 2 times of the Median Absolute Deviation divided
#'    by 0.6745.
#'  - `boxplot`: either smaller than the 1st quartile minus 1.5 times of IQR, 
#'     or larger than the 3rd quartile plus 1.5 times of IQR.
#' @family  outlier detectors
#' @examples 
#' x <- c(2, 3, 4, 5, 6, 7, 8, 9, 50, 50)
#' outl_det_u(x, "boxplot")
#' @export 
## wl-19-09-2020, Sat: Univariate outlier detection.
##   Modified from R package GmAMisc.
outl_det_u <- function(x, method = c("boxplot", "median", "mean")) {
  method <- match.arg(method)
  if (method == "mean") {
    outlier <- abs(x - mean(x)) > 2 * sd(x)
  }
  if (method == "median") {
    med <- median(x)
    mad <- median(abs(med - x))
    outlier <- abs(x - med) > 2 * (mad / 0.6745)
  }
  if (method == "boxplot") {
    q1 <- quantile(x, 0.25)
    q3 <- quantile(x, 0.75)
    iqr <- q3 - q1
    outlier <- x < q1 - 1.5 * iqr | x > q3 + 1.5 * iqr
  }
  return(outlier)
}

## ------------------------------------------------------------------------
#' Multivariate outlier detection
#' 
#' Perform multivariate outlier detection.
#' 
#' @param x a data matrix.
#' @param method methods for resistant estimation of multivariate location 
#'   and scatter. Only `mve`, `mcd` and `classical` are supported.
#' @param conf.level a confidential level.
#' @return retuns a logical vector.
#' @importFrom MASS cov.rob
#' @seealso [cov.rob()] for "Resistant Estimation of Multivariate Location
#'  and Scatter"
#' @family  outlier detectors
#' @export 
#' @examples  
#' set.seed(134)
#' x <- cbind(rnorm(80), rnorm(80), rnorm(80))
#' y <- cbind(rnorm(10, 5, 1), rnorm(10, 5, 1), rnorm(10, 5, 1))
#' x <- rbind(x, y)
#' outl <- outl_det_m(x, method = "mcd", conf.level = 0.95)
## lwc-30-01-2013: Multivariate outlier detection
## lwc-04-02-2013: Major changes.
## wll-01-12-2015: the results bewtween pca.outlier and outl_det_m are
##   different. Need to be careful.
## Note: "mve" and "mcd" are based on approximate search. User need to set up
##       random seed by set.seed. For details, see ?cov.rob
outl_det_m <- function(x, method = "mcd", conf.level = 0.95) {
  method <- match.arg(method, c("mve", "mcd", "classical"))
  x <- as.matrix(x)
  covr <- MASS::cov.rob(x, method = method) # MASS. NAs are not allowed.
  dis <- sqrt(mahalanobis(x, center = covr$center, cov = covr$cov))
  cutoff <- sqrt(qchisq(conf.level, ncol(x)))

  outlier <- dis > cutoff
  # outlier <- which(outlier)
  # if (!is.null(names(outlier))) outlier <- names(outlier)

  return(outlier)
}

## ------------------------------------------------------------------------
#' Batch shifting
#' 
#' Remove batch effect withing each block.
#' 
#' @param x a data matrix.
#' @param y  a categorical data for batch/block information.
#' @param method method for shifting.
#' @return returns a shifted data matrix.
#' @references 
#'   Silvia Wagner, et.al, Tools in Metabonomics: An Integrated Validation
#'   Approach for LC-MS Metabolic Profiling of Mercapturic Acids in Human
#'   Urine Anal. Chem., 2007, 79 (7), pp 2918-2926, DOI: 10.1021/ac062153w
#' @export 
## lwc-07-07-2011: batch shifting: remove mean withing each batch/block
batch_shift <- function(x, y, method = "mean") {
  x <- as.data.frame(x)

  g.mean <- sapply(x, function(x) tapply(x, y, method, na.rm = T))
  g.mean <- sapply(1:ncol(x), function(i) g.mean[, i][y])
  x <- x - g.mean

  return(x)
}

## ------------------------------------------------------------------------
#' Get the count number of features
#' 
#' Calculate the features counts based on the frequency of multiple selectors.
#' 
#' @param fs.ord data matrix for multiple feature order lists.
#' @param top.k top feature number to be processed.
#' @return a data matrix of feature counts.
#' @export 
## lwc-10-11-2021, Wed: get the count number of feature selectors in top-k
##  feature orders.
feat_count <- function(fs.ord, top.k = 30) {

  fs_top <- fs.ord[1:top.k, ]

  fs_top <- fs_top %>%
    as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "variable")

  fs_top <- fs_top %>%
    pivot_wider(names_from = variable, values_from = variable,
                id_cols = value, values_fn = length, values_fill = 0) %>%
    rowwise() %>%
    mutate(count = rowSums(across(where(is.numeric)))) %>%
    ## mutate(count = sum(c_across(where(is.numeric)))) %>%
    ungroup() %>%
    arrange(desc(count))

  return(fs_top)
}

## ------------------------------------------------------------------------
#' Get the count number of features
#' 
#' Calculate the features counts based on the frequency of multiple selectors.
#' 
#' @param fs.ord data matrix for multiple feature order lists.
#' @param top.k top feature number to be processed.
#' @return a data matrix of feature counts.
#' @noRd 
## wl-19-01-2016, Sat: get the count number of feature selectors in top-k
##  feature orders.
## wl-16-11-2021, Tue: 'reshape' works fine but 'reshape2::dcast()'
##  margins = "grand_col"  has problem.
feat_count_1 <- function(fs.ord, top.k = 30) {
  fs.ord <- fs.ord[1:top.k, ]

  ## Use melt and cast to get frequency table
  # tmp <- reshape::melt(fs.ord, id = NULL)
  # cons <- reshape::cast(tmp, value ~ variable,
  #   fun.aggregate = length,
  #   margins = "grand_col"
  # )
  ## FIX-ME: 'grand_col' does not work here. 
  tmp <- melt(fs.ord, id = NULL)
  cons <- dcast(tmp, value ~ variable,
    fun.aggregate = length,
    margins = "grand_col"
  )
  head(cons)
  names(cons)[which(names(cons) == "(all)")] <- "count"
  cons <- arrange(cons, desc(count))
}

## ------------------------------------------------------------------------
#' Relative standard deviation
#'
#' Caculate Relative Standard Deviation(RSD). RSD is also known as the 
#' coefficient of variation (CV)
#'
#' @param x an vector, matrix or data frame.
#' @param na.rm	a logical value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @return  RSD value multiplied by 100.
#' @details 
#' Some notes: 
#'   - The higher the CV, the greater the dispersion in the variable. 
#'   - The CV is often presented as the given ratio multiplied by 100
#'   - Basically CV<10 is very good, 10-20 is good, 20-30 is acceptable, 
#'     and CV>30 is not acceptable. 
#' @examples
#' data(iris)
#' rsd(iris[, 1:4])
#' 
#' ## group rsd
#' \dontrun{
#' df_summ <- function(dat, method = mean, ...) {
#'   dat %>% select(where(is.numeric)) %>% map_dfr(method) 
#' }
#' 
#' val <- iris %>% 
#'   group_by(Species) %>% 
#' 	group_modify(~ df_summ(., method = rsd)) %>%
#'   pivot_longer(cols = !Species) %>% filter(!is.na(value))
#' 
#' ggplot(val, aes(x = value)) +
#'     geom_histogram(colour = "black", fill = "white") +
#'     facet_grid(Species ~ .)
#' ggplot(val, aes(x = Species, y = value)) + geom_boxplot()
#' 
#' ## The plotting of missing value is similiar like this.
#' }
#' @export
## lwc-02-06-2011: Relative standard deviation.
## wl-22-11-2021, Mon: 
rsd <- function(x, na.rm = TRUE) {
  if (is.matrix(x)) {
    apply(x, 2, rsd, na.rm = na.rm)
  } else if (is.vector(x)) {
    100 * (sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm))
  } else if (is.data.frame(x)) {
    sapply(x, rsd, na.rm = na.rm)
  } else {
    100 * sd(as.vector(x), na.rm = na.rm) / mean(as.vector(x), na.rm = na.rm)
  }
}

## ------------------------------------------------------------------------
#' Missing value percentage
#' 
#' Calculate missing value percentage.
#'
#' @param x an vector, matrix or data frame.
#' @return missing value percentage.
#' @export 
## wl-24-11-2021, Wed: extract from 'mv.stats' in 'mt'. 
##   For plot of rsd and mv using tidyverse, see 'debug_filter_1.R' in 
##   'r_data/cam/dismp'
mv_perc <- function(x) {
  if (is.matrix(x)) {
    apply(x, 2, mv_perc)
  } else if (is.vector(x)) {
    round(sum(is.na(x) | is.nan(x)) / length(x), digits = 3)
  } else if (is.data.frame(x)) {
    sapply(x, mv_perc)
  } else {
    round(sum(is.na(as.vector(x)) | is.nan(as.vector(x))) / length(as.vector(x)), digits = 3)
  }
}

## ------------------------------------------------------------------------
#' Wrapper function for plotting classification restults
#' 
#' This function plots accuracy, AUC and margin (aam) of classification 
#' results from package \code{mt}.
#' 
#' @param aam a data matrix of classification results.
#' @param fig_title a string of figure title
#' @return an object of class `ggplot2`
#' @export 
## lwc-06-05-2011: Wrapper function for plotting the aam results
##    (Accuracy, AUC and Margin).
## wl-09-11-2021, Tue: implement with ggplot2
plot_aam <- function(aam, fig_title = "Accuracy, AUC and Margin") {
  z <- melt(aam)
  z <- z[complete.cases(z), ] ## in case NAs
  names(z) <- c("classifier", "assessment", "value", "data")
  z$classifier <- factor(z$classifier)

  aam.p <-
    ggplot(z, aes(y = value, x = data, color = classifier,
                  group = classifier)) +
    geom_line(aes(linetype = classifier)) +
    geom_point(aes(shape = classifier)) +
    ggtitle(fig_title) +
    facet_wrap(~ assessment) +
    coord_flip()
  # aam.p <-
  #     dotplot(factor(data, levels = rev(unique.default(data))) ~ value | assessment,
  #       data = z, groups = classifier, as.table = T,
  #       layout = c(length(unique(z$assessment)), 1),
  #       par.settings = list(
  #         superpose.line = list(lty = c(1:7)),
  #         superpose.symbol = list(pch = rep(1:25))
  #       ),
  #       type = "o", ## comment this line to get original dot plot
  #       auto.key = list(lines = TRUE, space = "bottom", columns = nlevels(z$classifier)),
  #       xlab = "", main = fig_title, ...
  #     )
  return(aam.p)
}

## ------------------------------------------------------------------------
#' Plot adjusted p-values
#' 
#' Plot the adjusted p-values using `ggplot2`
#' 
#' @param pval_list a data matix or a list of data matrix of p-value correction
#' @return an object of class `ggplot2`
#' @export  
## wl-09-11-2021, Tue: plot adjusted p-values
## wl-11-11-2021, Thu: deal with single matrix/data.frame
plot_pval <- function(pval_list) {

  if (!is.list(pval_list)) pval_list <- list(pval_list)

  tmp <- lapply(pval_list, function(x) {
    rownames(x) <- 1:nrow(x)
    ## wl-03-10-2021, Sun: the output of melt will be different if x is
    ##  either matrix or data frame.
    # x <- as.data.frame(x)
    x
  })

  ## wl-09-11-2021, Tue: actually melt is easy to use
  ## tmp <- reshape2::melt(tmp)
  ## p.pval <- xyplot(value ~ Var1 | L1,
  ##   data = tmp, groups = Var2,
  ##   par.settings = list(
  ##     superpose.line = list(lty = c(1:7)),
  ##     superpose.symbol = list(pch = rep(1:25))
  ##   ),
  ##   as.table = T, type = "l",
  ##   par.strip.text = list(cex = 1.0), ylim = c(-0.005, 1.0),
  ##   ylab = "P-values", xlab = "Index of feature ranks",
  ##   main = paste0(DF, "P-value rejection numbers"),
  ##   auto.key = list(lines = TRUE, points = F, space = "bottom", columns = 4),
  ##   panel = function(x, y, ...) {
  ##     panel.xyplot(x, y, ...)
  ##     panel.abline(h = alpha, col = "black", lty = 2)
  ##   }
  ## )

  ## wl-09-11-2021, Tue: equvalent to 'melt' by tidyverse. May be a temp
  ## function.
  ## wl-16-11-2021, Tue: use '.data' to prevent warning in 'rcmd check'
  tmp <- lapply(tmp, function(x) {
    x <- as_tibble(x) %>%
    rownames_to_column(var = "rn") %>%
    mutate(rn = as.numeric(.data$rn))
  })
  tmp <- tmp %>%
    bind_rows(.id = "data") %>%
    pivot_longer(-c(data, rn), names_to = "variable") %>%
    filter(complete.cases(.))

  ## wl-24-05-2021, Mon: combine two legends with the same legend title
  p <- ggplot(tmp, aes(x = rn, y = value)) +
    geom_line(aes(linetype = variable, color = variable)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    labs(x = "Index of variable", y = "p-values",
         color = "p-values type", linetype = "p-values type")
  ## p + facet_wrap(~ data, scales = "fixed")
  p.pval <- p + facet_grid(data ~ ., scales = "fixed")
  return(p.pval)
}

## ------------------------------------------------------------------------
#' Correlation hierarchical cluster analysis
#' 
#' Perform hierarchical cluster analysis based on correlation.
#' 
#' @param mat a data matrix
#' @param cutoff a threshold for correlation anaalysis.
#' @param use a string giving a method for computing covariances. For details, 
#'   see `cor`.
#' @param method method for correlation `cor`.
#' @param fig_title title for plotting
#' @param size text font size of `ggplot2`.
#' @param rotate a logical indicating wheter to rotate teh plot.
#' @return retuns a list of a `ggplot2` object for clustering and a list of 
#'   cluster centres.
#' @importFrom ggdendro ggdendrogram dendro_data
#' @export 
## wl-05-10-2021, Tue: correlation hierarchical cluster analysis
## Note:
##  1.) use 'ggdendro' but more eddforts are needed.
##  2.) modify from 'cor.hcl' in R package 'mt'
cor_hcl <- function(mat, cutoff = 0.75,
                    use = "pairwise.complete.obs",
                    method = "pearson",
                    fig_title = "Cluster Dendrogram",
                    size = 3,
                    rotate = FALSE) {
  # require(ggdendro)
  co <- cor(mat, use = use, method = method)
  hc <- hclust(as.dist(1 - co))

  # get dendrogram
  if (T) {  # Using ggdendro directly from object class hclust
    p <- ggdendrogram(hc, rotate = rotate, theme_dendro = F, size = size)
    p <- p +
      geom_hline(yintercept = 1 - cutoff, color = "red", linetype = "dashed") +
      labs(title = fig_title, x = "", y = "1 - Correlation")
  } else {  # get dendrogram manually
    dd <- dendro_data(hc, type = "rectangle") # "rectangle" or "triangle"
    # dd <- dendro_data(as.dendrogram(hc), type = "rectangle")
    p <- ggplot(dd$segments) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_text(data = dd$labels, aes(x, y, label = label),
                hjust = 1, angle = 90, size = size)
    p <- p +
      geom_hline(yintercept = 1 - cutoff, color = "red", linetype = "dashed") +
      labs(title = fig_title, x = "", y = "1 - Correlation")
    if (rotate) {
      p <- p + coord_flip()
    }
  }

  ## get clusters which correlation coefficients larger than cutoff
  id <- cutree(hc, h = 1 - cutoff)
  res <- lapply(unique(id), function(x) {
    cent <- mat[, id == x, drop = FALSE]
    res <- if (ncol(cent) < 2) NA else cor(cent, use = use, method = method)
  })

  ## shrink the list
  id <- sapply(res, function(x) {
    if (!any(is.na(x))) TRUE else FALSE
  })
  res <- res[id]
  names(res) <- paste("Cluster", 1:length(res), sep = "_")

  # res$all <- co
  return(list(dendro_plot = p, clusters = res))
}

## ------------------------------------------------------------------------
#' Update data set by a subset of features
#' 
#' Update data set after feature selection. 
#' 
#' @param dat.all a list of metobolomics data including indensity data, peak 
#'   data and meta data
#' @param ord a subset of selected freature orders
#' @export 
## lwc-03-03-2010: Update data set by a subset of features or variables.
upd_data <- function(dat.all, ord) {
  dn <- names(ord)

  data <- lapply(dn, function(i) {
    dat.all[[i]]$dat <- dat.all[[i]]$dat[, ord[[i]], drop = F]
    dat.all[[i]]
  })
  names(data) <- dn
  ## debug:print out the dimension
  print(sapply(data, function(x) dim(x$dat)))
  data
}

## -----------------------------------------------------------------------
#' Correlation network analysis
#' 
#' Perform network analysis for correlation coefficient.
#' 
#' @param mat a data matrix for correlation analysis
#' @param use,method parameters of `cor`.
#' @param thres correlation coeeficient threshold for network analysis. Only 
#'   keep those with coeeficient larger than `thres`.
#' @param fig_title a chracter string for gifure title.
#' @return a list of `ggplot2` plots. 
#' @importFrom ggraph create_layout ggraph geom_edge_link geom_node_point
#'    scale_edge_colour_discrete geom_node_text
#'    facet_edges facet_nodes
#' @importFrom igraph graph_from_data_frame delete_vertices cluster_louvain
#'   E<- V<-
#' @export  
## wl-07-10-2021, Thu: Correlation network analysis
## wl-08-10-2021, Fri: Call function 'net_graph'
##  Note: This function is different from 'corr_net' which is for
##        Bipartite-mode correlation network.
cor_net <- function(mat,
                    use = "pairwise.complete.obs",
                    method = "pearson",
                    thres = 0.6,
                    fig_title = "Correlation network") {

  co <- cor(mat, use = use, method = method)

  ## wl-05-11-2020, Thu: create igraph object
  edge <- sym2long(co, tri = "upper")
  idx <- abs(edge$var) >= thres
  edge <- edge[idx, ]
  edge <- mutate(edge, mod = ifelse(var > 0, "pos", "neg"))

  ## only positive or negative?
  if (F) {
    if (T) {
      idx <- edge$var > 0 ## positive
    } else {
      idx <- edge$var < 0 ## negative
    }
    edge <- edge[idx, ]
  }

  ## get graph object
  if (F) {   ## direct method
    ## get net graph
    g <- graph_from_data_frame(d = edge, directed = F, vertices = NULL)
    ## delete isolated vertices
    g <- delete_vertices(g, which(igraph::degree(g) < 1))
    ## Community detection
    com <- cluster_louvain(g, weights = NULL)
    mem <- as.factor(membership(com))
    V(g)$clus <- mem
  } else {  ##  call function
    g <- net_graph(edge_df = edge, node_df = NULL)
  }

  lay <- create_layout(g, layout = "fr")
  p <- ggraph(lay) +
    geom_edge_link(alpha = .25, aes(color = as.factor(mod))) +
    geom_node_point(size = 2) +
    scale_edge_colour_discrete(name = "Edge") +
    ggtitle(fig_title)
  ## theme_graph()

  ## clusters based on network community detection
  p1 <- p + geom_node_point(size = 2, aes(color = as.factor(clus))) +
    scale_colour_discrete(name = "Network Community")

  ## with text
  p2 <- p1 + geom_node_text(aes(label = name), repel = TRUE)

  ## facet with positive and negative
  p3 <- p2 + facet_edges(~mod)

  ## facet with clusters
  p4 <- p2 + facet_nodes(~clus)

  res <- list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)
  return(res)
}

## -----------------------------------------------------------------------
#' Bipartite/two-mode correlation network
#' 
#' Perform bipartite/two-mode correlation network analysis.
#' 
#' @param co_mat a correlation coffeficient matrix.
#' @param thres correlation coeeficient threshold for network analysis.
#' @param dn a character string for correlation name.
#' @return a list of `ggplot2` plots.
#' @importFrom igraph graph_from_data_frame bipartite_mapping V E
#'   cluster_edge_betweenness cluster_fast_greedy degree membership V<- E<-
#' @importFrom ggraph create_layout ggraph facet_nodes geom_node_point 
#'   geom_edge_link geom_node_text scale_edge_color_manual geom_edge_arc
#' @export 
## wl-07-05-2021, Fri: Bipartite/two-mode correlation network
## wl-19-06-2021, Thu: change family font from 'Arial' to 'sans'. So no
##  'extrafont' is needed for PDF output.
bi_cor_net <- function(co_mat, thres = 0.6, dn = NULL) {
  g_dat <- melt(co_mat)
  if (!is.null(dn)) {
    names(g_dat) <- c(dn, "corr")
  } else {
    names(g_dat)[3] <- "corr"
  }

  ## wl-29-04-2021, Thu: differenct data sets need different threshold
  ## ind <- g_dat$corr > thres
  ind <- abs(g_dat$corr) > thres
  mat <- g_dat[ind, ]
  dim(mat)

  ## ---- Network analysis with igraph/ggraph ----
  g         <- graph_from_data_frame(mat, directed = F)
  V(g)$type <- bipartite_mapping(g)$type    ## convert into bipartite
  V(g)$data <- ifelse(V(g)$type, names(g_dat)[2], names(g_dat[1]))
  V(g)$size <- degree(g)
  E(g)$edge <- ifelse(E(g)$corr > 0, "pos", "neg")
  ## edge_attr_names(g)
  ## vertex_attr_names(g)

  ## Community detection
  ## com <- cluster_edge_betweenness(g)         ## wl-09-06-2021: not good
  com <- cluster_fast_greedy(g)                 ## based on modularity score
  ## com <- cluster_louvain(g, weights = NULL)  ## wl-09-06-2021: problem
  clus <- as.factor(membership(com))

  ## set edge colour platte
  ## edge_col <- c("neg" = "#EEB422", "pos" = "#3D3D3D")
  edge_col <- c("pos" = "red", "neg" = "black")

  ## plot network. Control vertex: platform or clus?
  lay <- create_layout(g, layout = "fr")    ## kk, graphopt
  (np <-
    ggraph(lay) +
    theme(text = element_text(family = "sans")) +
    scale_edge_color_manual(values = edge_col) +
    geom_edge_arc(strength = 0.2, width = 0.5, alpha = 0.15, aes(color = edge)) +
    ## geom_edge_link(alpha = .25, aes(color = edge)) +
    geom_node_point(aes(color = clus, size = size), shape = 20)
    ## geom_node_point(aes(color = data, size = size), shape = 20)
  )

  ## (np1 <- np + geom_node_text(aes(label = name, size = size), repel = T))
  ## (np1 <- np + geom_node_text(aes(label = name, size = 2), repel = T))
  (np1 <- np + geom_node_text(aes(label = name, color = clus, size = 2),
                              repel = T))
  (np2 <- np1 + facet_nodes(clus))

  ## ---- Community detection plot with dendrogram ----
  ## convert communities to 'hclust' or 'dendrogram'
  hc <- as.hclust(com)

  ## get number of groups for cutree/rect.hclust
  k <- length(com)

  ## plot dendrogram
  hcdata <- dendro_data_k(hc, k = k)
  (dp1 <- plot_ggdendro(hcdata, fan = TRUE, label.size = 3,
                        nudge.label = 0.02, expand.y = 0.4))

  (dp2 <- plot_ggdendro(hcdata, direction = "lr",  ## 'tb'
                        label.size = 3, branch.size = 0.5, expand.y = 0.2))
  return(list(net1 = np1, net2 = np2, dend1 = dp1, dend2 = dp2))
}

## ------------------------------------------------------------------------
#' Correlation between two data sets
#' 
#' Performs correlation between two data sets
#' 
#' @param dat_1,dat_2 two data matrix for correlation analysis. Note that they
#'   must have the same row numbers.
#' @param partial a logical value indicating whether to perform partial 
#'   correlation analysis.
#' @param method correlation method.
#' @param ... further parameters to correlation method.
#' @return a mstrix of correlation.
#' @export 
#' @details This function is very specific. Will change for more general 
#'   purpose.
## wl-04-05-2021, Tue: Correlation between two data sets. Note
##   that samples must be the same. This function is not general.
## wl-17-05-2021, Mon: more general
## wl-03-06-2021, Thu: move out filtering
## wl-04-06-2021, Fri: add partial correlation
dat_cor <- function(dat_1, dat_2, partial = FALSE,
                     method = c("pearson", "kendall", "spearman"), ...) {

  ## make sure the rows are the same for the two data sets
  idtp_1 <- dat_1$meta$id_tp
  idtp_2 <- dat_2$meta$id_tp

  (idtp <- intersect(idtp_1, idtp_2))
  setdiff(idtp_1, idtp)
  setdiff(idtp_2, idtp)

  ## update both data sets
  idx_1 <- unlist(sapply(idtp, function(x) which(idtp_1 == x)))
  idx_2 <- unlist(sapply(idtp, function(x) which(idtp_2 == x)))
  ## or use grep with exact match
  ## idx_1 <- unlist(sapply(idtp, function(x) grep(paste0("^",x,"$"), idtp_1, perl = T)))
  ## idx_2 <- unlist(sapply(idtp, function(x) grep(paste0("^",x,"$"), idtp_2, perl = T)))

  dat_1$data <- dat_1$data[idx_1,]
  dat_2$data <- dat_2$data[idx_2,]

  ## get the correlation matrix.
  if (partial) {
    co <- pcor_dat(dat_1$data, dat_2$data, method = method)
  } else {
    co <- cor(dat_1$data, dat_2$data, method = method, ...)
  }
  return(co)
}

## ------------------------------------------------------------------------
#' Partial correlation analysis
#' 
#' Perform partial correlation analysis of two data matrix
#' 
#' @param x,y two data matrix for correlation anlysys.
#' @param method correlation method.
#' @return returns a correlation matrix.
#' @importFrom ppcor pcor
#' @export 
## wl-04-06-2021, Fri: partial correlation between two data frame
pcor_dat <- function(x, y, method = c("pearson", "kendall", "spearman")) {
  if (!(is.matrix(x) || is.data.frame(x))) {
    stop("supply a matrix-like 'x'")
  }
  if (!(is.matrix(y) || is.data.frame(y))) {
    stop("supply a matrix-like 'y'")
  }
  if (nrow(x) != nrow(y)) {
    stop(" x and y must have the same number of rows")
  }

  dat <- cbind(x, y)
  res <- ppcor::pcor(dat, method = method)

  ## only get correlation coefficients
  res <- res$estimate[-c(ncol(x) + 1:ncol(dat)), -c(1:ncol(x))]
  rownames(res) <- colnames(x)
  colnames(res) <- colnames(y)

  return(res)
}

## -----------------------------------------------------------------------
#' Heatmap with dendrograms with ggplot2
#' 
#' Plot heatmap of a data matrix using `ggplot2`. This funcion is modofied 
#' from https://bit.ly/2UUnY2L.
#' 
#' @param mat a data mstrix to be plotted.
#' @param row.dend plot row dendrogram or not.
#' @param col.dend plot column dendrogram ot not.
#' @param row.dend.right a ligical value to indivate the position of row 
#'   dendrogram.
#' @param colors a vector of colours for heatmap.
#' @param font.size label font size.
#' @param x.rot plot rotate degree.
#' @param legend.title lengend title.
#' @param dist.method distance method.
#' @param clust.method cluster method.
#' @param dend.line.size dendrogram line size.
#' @return retuns an object of class `ggplot2`.
#' @importFrom cowplot axis_canvas insert_yaxis_grob ggdraw
#' @importFrom ggdendro dendro_data segment
#' @export 
#' @examples 
#' library(ggplot2)
#' gg_heat_dend(mtcars)
## wl-24-11-2020, Tue: Heatmap with dendrograms with ggplot2
## wl-15-06-2021, Tue: add 'row_den_left'
## wl-15-10-2021, Fri: Review
gg_heat_dend <- function(mat,
                         row.dend = T,
                         col.dend = T,
                         row.dend.right = TRUE,
                         colors = c("red", "white", "blue"),
                         ## colors = c("yellow", "red"),
                         font.size = 10,
                         x.rot = 60,
                         legend.title = "value",
                         dist.method = "euclidean",
                         clust.method = "complete",
                         dend.line.size = 0.5) {

  ## data_m <- tibble::rownames_to_column(mat) %>% reshape2::melt()
  mat <- as.data.frame(mat)
  data_m <- cbind(rn = rownames(mat), mat)
  data_m <- reshape2::melt(data_m, id.vars = "rn")

  # Cluster rows
  if (row.dend) {
    dd.row <- as.dendrogram(hclust(dist(mat, method = dist.method),
                                   method = clust.method))
    row.ord <- order.dendrogram(dd.row)
    ordered_row_names <- row.names(mat[row.ord, ])
    data_m$rn <- factor(data_m$rn, levels = ordered_row_names)
  }

  # Cluster columns
  if (col.dend) {
    dd.col <- as.dendrogram(hclust(dist(t(mat), method = dist.method),
                                   method = clust.method))
    col.ord <- order.dendrogram(dd.col)
    ordered_col_names <- colnames(mat[, col.ord])
    data_m$variable <- factor(data_m$variable, levels = ordered_col_names)
  }

  ## wl-17-06-2021, Thu: remove vjust for good
  p_heat <- ggplot(data_m, aes(x = variable, y = rn, fill = value)) +
    theme_minimal() + geom_tile() + xlab("") + ylab("") +
    theme(axis.line = element_line(size = 0),
          text = element_text(size = font.size),
          axis.text.x = element_text(angle = x.rot, hjust = 1))

  ## adding colors
  p_heat <- p_heat +
    scale_fill_gradient2(low = colors[1], mid = colors[2], high = colors[3],
                         midpoint = mean(data_m$value),
                         space = "Lab", name = legend.title)

  if (row.dend) {
    dendro_data_row <- ggdendro::dendro_data(dd.row, type = "rectangle")

    if (row.dend.right) {
      dendro_row <-
        cowplot::axis_canvas(p_heat, axis = "y", coord_flip = TRUE) +
        geom_segment(data = ggdendro::segment(dendro_data_row),
                     aes(y = y, x = x, xend = xend, yend = yend),
                     size = dend.line.size) +
        coord_flip()

      p <- cowplot::insert_yaxis_grob(p_heat, dendro_row,
                                      grid::unit(0.2, "null"),
                                      position = "right")
    } else {
      dendro_row <-
        cowplot::axis_canvas(p_heat, axis = "y", coord_flip = TRUE) +
        geom_segment(data = ggdendro::segment(dendro_data_row),
                     aes(y = -y, x = x, xend = xend, yend = -yend),
                     size = dend.line.size) +
        coord_flip()

      p_heat <- p_heat + scale_y_discrete(position = "right")
      p <- cowplot::insert_yaxis_grob(p_heat, dendro_row,
                                      grid::unit(0.2, "null"),
                                      position = "left")
    }
  } else {
    p <- p_heat
  }

  if (col.dend) {
    dendro_data_col <- ggdendro::dendro_data(dd.col, type = "rectangle")
    dendro_col <- cowplot::axis_canvas(p_heat, axis = "x") +
      geom_segment(data = ggdendro::segment(dendro_data_col),
                   aes(x = x, y = y, xend = xend, yend = yend),
                   size = dend.line.size)

    p <- cowplot::insert_xaxis_grob(p, dendro_col, grid::unit(0.2, "null"),
                                    position = "top")
  }

  p <- cowplot::ggdraw(p)
  return(p)
}

## ------------------------------------------------------------------------
#' Heatmap with dendrogram on both sides
#' 
#' Plot heatmap with dendrogram on both sides using lattice
#' 
#' @param mat a data matrix for plotting.
#' @param x.rot the rotate degree.
#' @param col colours for heatmap.
#' @param cex font size.
#' @param xlab,ylab,main character strings for `xlab`, `ylab` and figure title.
#' @param ... further parameters for `latticw`.
#' @return an object of class `lattice`.
#' @importFrom lattice levelplot panel.fill panel.levelplot
#' @importFrom latticeExtra dendrogramGrob
#' @export 
#' @examples 
#' heat_dend(mtcars)
## wl-04-09-2015, Fri: Heatmap with dendrogram on both sides using lattice
## wl-15-06-2021, Tue: Makes heatmap consistent with data matrix and remove
## non-lattice parts
heat_dend <- function(mat, x.rot = 60,
                      col = c("red", "white", "blue"),
                      cex = 0.5, xlab = "", ylab = "", main = "", ...) {
  mat <- as.matrix(mat)
  row.dd <- as.dendrogram(hclust(dist(mat)))
  row.ord <- order.dendrogram(row.dd)

  col.dd <- as.dendrogram(hclust(dist(t(mat))))
  col.ord <- order.dendrogram(col.dd)

  col <- colorRampPalette(col)

  ph <-
    levelplot(t(mat[row.ord, col.ord]),
      aspect = "fill",
      col.regions = col, ## colorRampPalette(c("yellow", "red")),
      xlab = xlab, ylab = ylab, main = main,
      scales = list(x = list(rot = x.rot), cex = cex),
      colorkey = list(space = "bottom"),
      ## wll-10-09-2015: Use panel.fill() to fill the background with
      ## your 'NA' color.From Deepayan Sarkar
      panel = function(...) {
        panel.fill(col = "gray") # "black")
        panel.levelplot(...)
      },

      ## colorkey = list(space = "bottom", labels=list(cex=cex)),
      legend =
        list(
          right =
            list(
              fun = dendrogramGrob,
              args =
                list(
                  x = row.dd, ord = row.ord,
                  side = "right",
                  size = 10
                )
            ),
          top =
            list(
              fun = dendrogramGrob,
              args =
                list(
                  x = col.dd, ord = col.ord,
                  side = "top",
                  type = "rectangle", size = 5
                )
            )
        ),
      ...
    )

  return(ph)
}

## -----------------------------------------------------------------------
#' Vector summary
#' 
#' Calculate the statistical summary of a vector.
#' 
#' @param x a numveric vector.
#' @param na.rm remove NA or not.
#' @param conf.interval a numeric value for confidenice interval.
#' @return 
#'   retuns a vector of summary consisting:
#'   - number of vector length
#'   - vector mean
#'   - vector standard derivation   
#'   - standard error of mean
#'   - confidence interval
#' @details Can be used for error bar plotting. Modify from
#'   https://bit.ly/3onsqot
#' @family vector stats functions
#' @export 
## wl-18-05-2021, Tue: stats of a vector.
vec_stats <- function(x, na.rm = FALSE, conf.interval = .95) {

  ## Handle NA's: if na.rm==T, don't count them
  length2 <- function(x, na.rm = FALSE) {
    if (na.rm) sum(!is.na(x)) else length(x)
  }
  n    <- length2(x, na.rm = na.rm)
  mean <- mean(x, na.rm = na.rm)
  sd   <- sd(x, na.rm = na.rm)
  se   <- sd / sqrt(n)                 ## standard error of mean

  ## Confidence interval multiplier for standard error
  ## Calculate t-statistic for confidence interval:
  ## e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ci_mult <- qt(conf.interval / 2 + .5, n - 1)
  ci <- se * ci_mult

  c(n = n, mean = mean, sd = sd, se = se, ci = ci)
}

## ------------------------------------------------------------------------
#' Vector statistics for error bar plotting
#' 
#' Calculate vector's standard derivation, standard error of mean and 
#' confidencd interval.
#' 
#' @param x an vector.
#' @param bar a chracter string, supporting "SD", "SE" and "CI".
#' @return retuns an vector including lower, center and upper values.
#' @examples
#' vec_segment(iris[,1])
#' mat <- reshape2::melt(iris)
#' plyr::ddply(mat, plyr::.(Species,variable), function(x,bar) {
#'   vec_segment(x$value, bar = bar)
#' }, bar = "SD")
#' @family vector stats functions
#' @importFrom stats t.test
#' @export 
## lwc-13-05-2013: Calculate error bar statistics
## lwc-08-01-2014: add some examples.
vec_segment <- function(x, bar = c("SD", "SE", "CI")) {
  bar <- match.arg(bar)

  centre <- mean(x, na.rm = T)

  if (bar == "SD") {
    stderr <- sd(x, na.rm = T) # Standard derivation (SD)
    lower <- centre - stderr
    upper <- centre + stderr
  } else if (bar == "SE") {   # Standard error(SE) of mean
    stderr <- sd(x, na.rm = T) / sqrt(sum(!is.na(x)))
    ## stderr <- sqrt(var(x, na.rm = T)/length(x[complete.cases(x)]))
    lower <- centre - stderr
    upper <- centre + stderr
  } else if (bar == "CI") {   # Confidence interval (CI), here 95%.
    conf <- t.test(x)$conf.int
    lower <- conf[1]
    upper <- conf[2]
  } else {
    stop("'method' invalid")
  }

  res <- c(lower = lower, centre = centre, upper = upper)
  return(res)
}

## --------------------------------------------------------------------------
#' Create an igraph object
#' 
#' Create an `igraph` object
#' 
#' @param edge_df graph eddge matrix.
#' @param node_df graph vertex matrix.
#' @return an `igraph` object with some vertex statistics and community
#'  detection center.
#' @importFrom igraph graph_from_data_frame delete_vertices set_vertex_attr
#' @importFrom igraph edge_attr edge_attr_names vertex_attr vertex_attr_names
#'   gorder gsize V E 
#' @importFrom igraph cluster_fast_greedy cluster_edge_betweenness
#'   cluster_louvain membership components
#' @export 
## wl-05-11-2020, Thu: create an igraph object
net_graph <- function(edge_df, node_df = NULL) {

  ## Create network
  g <- graph_from_data_frame(d = edge_df, directed = F, vertices = node_df)
  ## delete isolated vertices
  g <- delete_vertices(g, which(igraph::degree(g) < 1))

  ## get node metrics
  stats <- vertex_stats(g)

  ## set vertex attributes. Note that the apply does not work in this occasion.
  dn <- names(stats)
  for (i in dn) {
    g <- set_vertex_attr(g, i, index = V(g), stats[[i]])
  }
  ## edge_attr_names(g)
  ## vertex_attr_names(g)
  ## vertex_attr(g)
  ## edge_attr(g)
  ## gorder(g)                 ## vcount(g)
  ## gsize(g)                  ## ecount(g)

  ## Community detection
  # com <- cluster_edge_betweenness(g)
  # com <- cluster_louvain(g, weights = NULL)
  com <- cluster_fast_greedy(g)        ## based on modularity score
  mem <- membership(com)
  # mem <- as.factor(membership(com))
  V(g)$clus <- as.integer(mem)

  ## wl-24-10-2020, Sat: Find connected components
  ## cl <- components(g)
  ## plot(g, layout = layout, vertex.color = cl$membership + 1L)

  return(g)
}

## ------------------------------------------------------------------------
#' Graph level metrics
#' 
#' Return graph level metrics.
#' 
#' @param g an `igraph` object.
#' @return a vector of graph level metrics.
#' @importFrom igraph transitivity diameter edge_density reciprocity 
#'   assortativity_degree cluster_fast_greedy modularity mean_distance radius
#'   centr_degree centr_clo centr_betw centr_eigen
#' @export 
## wl-22-03-2017: Graph level metrics
graph_stats <- function(g) {
  ## num.nod <- gorder(g)                              ## nodes
  ## num.edg <- gsize(g)                               ## edges
  ## num.neg <- (sum(E(g)$corr < 0)/num.edg) * 100     ## positive edges
  ## num.pos <- (sum(E(g)$corr >= 0)/num.edg) * 100    ## negative edges

  ## Transitivity: measures the probability that the adjacent
  ## vertices of a vertex are connected. This is sometimes also called the
  ## **clustering coefficient**.  Graph is treated as an undirected network.
  tran <- transitivity(g, type = "global")

  ## Diameter: The diameter of a graph is the length of the longest geodesic.
  diam <- diameter(g, directed = F, weights = NA)

  ## Density: The proportion of present edges from all possible ties.
  dens <- edge_density(g, loops = F)

  ## Reciprocity: The proportion of reciprocated ties (for a directed
  ## network).  The measure of reciprocity defines the proportion of mutual
  ## connections, in a directed graph. It is most commonly defined as the
  ## probability that the opposite counterpart of a directed edge is also
  ## included in the graph.
  reci <- reciprocity(g)

  ## Assortativity (homophily): The tendency of nodes to connect to others
  ## who are similar on some variable.
  asso <- assortativity_degree(g, directed = F)

  ## Modularity: The modularity of a graph with respect to some division (or
  ## vertex types) measures how good the division is, or how separated are
  ## the different vertex types from each other.  Modularity calculates the
  ## modularity of a graph with respect to the given membership vector.
  cfg <- cluster_fast_greedy(g)
  modu <- modularity(cfg)

  ## Average path length: Calculating the shortest paths between all pairs
  ## of vertices (both ways for directed graphs). This function does not
  ## consider edge weights currently and uses a breadth-first search.
  ave.path <- mean_distance(g) ## "mean_distance"

  ## The eccentricity of a vertex is its shortest path distance from the
  ## farthest other node in the graph. The smallest eccentricity in a graph
  ## is called its *radius*.
  radi <- radius(g)

  ## centralization of four centrality
  c_deg <- centr_degree(g)$centralization
  c_clo <- centr_clo(g, mode = "all")$centralization
  c_btw <- centr_betw(g, directed = FALSE)$centralization
  c_eig <- centr_eigen(g, directed = FALSE)$centralization

  ## graph descriptive
  graph.stats <- c(
    Transitivity = tran,
    Diameter = diam,
    Density = dens,
    Reciprocity = reci,
    Assortativity = asso,
    Modularity = modu,
    AveragePathLength = ave.path,
    Radius = radi,
    CentralizationDegree = c_deg,
    CentralizationCloseness = c_clo,
    CentralizationBetweenness = c_btw,
    CentralizationEigen = c_eig
  )

  return(graph.stats)
}

## ------------------------------------------------------------------------
#' Node descriptive metrics
#' 
#' Return node descriptive metrics.
#' 
#' @param g an `igraph` object.
#' @importFrom igraph distances gorder diameter transitivity betweenness
#'   eccentricity knn cluster_edge_betweenness cluster_fast_greedy
#'   membership
#' @export 
## wl-08-02-2017: Node descriptive
## wl-09-02-2017: add more descriptive
## wl-01-03-2017: add the membership of community detection for grouping.
## wl-22-03-2017: move graph-level part into 'graph.stats' and change as
## 'vertex.stats'.
vertex_stats <- function(g) {

  ## wl-09-02-2017: radiality. Taken from package `CePa`
  radiality <- function(graph, mode = c("all", "in", "out")) {
    mode <- mode[1] ## wl: why only "all"?

    sp <- distances(graph, mode = mode, weights = NA)
    n <- gorder(graph)
    diam <- diameter(graph, directed = ifelse(mode == "all", FALSE, TRUE))
    s <- apply(sp, 1, function(x) {
      if (all(x == Inf)) {
        return(0)
      }
      else {
        return(sum(diam + 1 - x[x != Inf]) / (n - 1))
      }
    })
    return(s)
  }

  ## Distance matrix
  dis <- distances(g) ## or use `distance`
  ## Need to convert inf to NA for calculating mean
  dis[mapply(is.infinite, dis)] <- NA

  ## Transitivity in vertex level, also called clustering coefficients
  transi <- transitivity(g, type = "local")
  transi[is.nan(transi)] <- 0

  ## Vertex desc
  ## wl-10-03-2017: It is very strange that the values of closeness are too
  ##  small. Have tried out, in, all and total, but all these values are
  ##  small.
  node.stats <- data.frame(
    ## Abandance = V(g)$abd,             ## specific. not in node.
    AverageShortestPathLength = rowMeans(dis, na.rm = T),
    BetweennessCentrality = betweenness(g, directed = FALSE),

    ## wl-05-11-2020, Thu: igraph complain the following line
    ## ClosenessCentrality = closeness(g, mode = "all", normalized = F),

    ClusteringCoefficient = transi,
    DegreeAll = igraph::degree(g),
    Eccentricity = eccentricity(g),
    NeighborhoodConnectivity = knn(g)$knn,
    ## K-NN neighbour degree. DegreeOfNeighbours
    ## wl-16-02-2017: For the information for treating knn as neighbourhood
    ## Connectivity see https://goo.gl/tRV2tg. Beware that it is different from
    ## degree substantially.

    Radiality = radiality(g)
    ## PageRank                  = page_rank(g)$vector,
    ## WeightedDegree            = strength(g),
    ## EigenvectorCentrality     = eigen_centrality(g,directed=FALSE)$vector,
    ## PowerCentrality           = power_centrality(g),
    ## HubCentrality             = hub_score(g)$vector,
    ## AuthorityCentrality       = authority_score(g)$vector,
    ## AlphaCentrality           = alpha_centrality(g),
    ## SubgraphCentrality        = subgraph_centrality(g)
  )

  ## wl-01-03-2017: Community detection
  ## Note: need to debug for potential bugs such as failing of computing of
  ## community detection for extremely small networks.
  ceb <- cluster_edge_betweenness(g)
  cfg <- cluster_fast_greedy(g)
  ceb_mem <- as.vector(membership(ceb))
  cfg_mem <- as.vector(membership(cfg))

  node.stats <- cbind(node.stats,
    EdgeBetweennessClass = ceb_mem,
    GreedyModularityClass = cfg_mem
  )

  ## sort based on degree and abundance
  ## node.stats <- arrange_row(node.stats, desc(Degree), desc(Abundance))

  return(node.stats)
}

## ------------------------------------------------------------------------
#' Get correlation coefficient and p-values
#' 
#' Return correlation coefficient and p-values.
#' 
#' @param x a data frame or matrix for correlation analysis column-wise.
#' @param cor.method method for correlation
#' @param adj.method p-value correction method
#' @param ... other parameter for correlation.
#' @return a list of with contents: \itemize{
#'  \item r correlation coefficient
#'  \item p statistics matrix, in which the lower triangular is p-values 
#'    and the upper triangular is adjusted p-values
#' }
#' @details 
#'   This file is modified from 'cor.table' of package 'picante'
#'   and 'corr.test' of package 'psych'.
#'   The original implementation is from Bill Venables, the author of R great
#'   book MASS. For details, see
#'   https://stat.ethz.ch/pipermail/r-help/2001-November/016201.html
#' @export 
#' @examples 
#' co <- cor_tab(mtcars, cor.method = "spearman", adj.method = "BH")
#' names(co)
#' 
#' \dontrun{
#' library(pysch)
#' co <- corr.test(mtcars, method="spearman",adjust="BH")
#' # From pysch: For symmetric matrices, p values adjusted for multiple tests
#' # are reported above the diagonal.
#' }
## wl-23-06-2015: Get correlation coefficient and p-values
cor_tab <- function(x,
                    cor.method = c("pearson", "kendall", "spearman"),
                    adj.method = c("holm", "hochberg", "hommel",
                                   "bonferroni", "BH", "BY", "fdr", "none"),
                    ...) {
  R <- cor(x, method = cor.method, ...)
  df <- dim(x)[1] - 2

  if (T) {
    t <- R * sqrt(df / (1 - R^2))
    P <- 2 * (1 - pt(abs(t), df)) ## from corr.test: two-tailed
    ## P <- 2*pt(t, df)           ## from cor.table: right-tailed (greater)
    ## P[R>0] <- 2*pt(t[R>0], df,lower.tail=FALSE)
  } else { ## from Bill Venables
    F <- R^2 * df / (1 - R^2)
    P <- 1 - pf(F, 1, df)
  }
  diag(P) <- NA


  ## get adjusted p-values
  idx <- upper.tri(P, diag = FALSE)
  pval <- P[idx]
  padj <- p.adjust(pval, method = adj.method)
  P[upper.tri(P, diag = FALSE)] <- padj

  list(r = R, p = P)
}

## -----------------------------------------------------------------------
#' Convert a symmetric table(short format) to long format
#' 
#' Convert a symmetric table(short format) to long format
#' 
#' @param x a symmetric matrix-like data set.
#' @param tri triangular being used.
#' @return returns a data frame of pair-wise comparison.
#' @examples 
#' co <- cor_tab(mtcars, cor.method = "spearman", adj.method = "BH")
#' names(co)
#' corr <- sym2long(co$r, tri = "upper")
#' pval <- sym2long(t(co$p), tri = "upper")
#' padj <- sym2long(co$p, tri = "upper")
#' tmp  <- data.frame(corr, pval, padj)
#' @export 
## wl-24-06-2015: Convert a symmetric table(short format) to long format
sym2long <- function(x, tri = c("upper", "lower")) {
  tri <- match.arg(tri)

  if (tri == "upper") {
    ind <- lower.tri(x)
    x[ind] <- NA
    diag(x) <- NA
    x <- x[-nrow(x), -1, drop = F]
  } else if (tri == "lower") {
    ind <- upper.tri(x)
    x[ind] <- NA
    diag(x) <- NA
    x <- x[-1, -ncol(x), drop = F]
  } else { ## It never reaches here because of match.arg.
    stop("Invalid method")
  }

  if (F) { 
    res <- reshape2::melt(x)
    res <- res[complete.cases(res), ]
    colnames(res) <- c("com1", "com2", "var")
  } else {
    idx <- which(!is.na(x), arr.ind = T)
    fs1 <- rownames(x)[idx[, 1]]
    fs2 <- colnames(x)[idx[, 2]]
    res <- data.frame(
      com1 = fs1, com2 = fs2, var = x[idx],
      stringsAsFactors = FALSE
    )
  }
  return(res)
}

## ------------------------------------------------------------------------
#' Scale vector to a defined range
#' 
#' Scale vector to a defined range.
#'
#' @param x a numeric vector
#' @param range a vector with two values: lower and higher.
#' @return a scaled vector.
#' @export 
#' @examples 
#' set.seed(100)
#' x <- rnorm(10)
#' range_scale(x, range = c(10, 20))
## wl-08-07-2015: Map values to a defined range. 
## wl-05-05-2016: Rescale vector between min and max.
##   similar function `rescale` in package `network`
range_scale <- function(x, range = c(0, 1)) {

  if (range[2] <= range[1]) {
    stop("range should be low and high")
  }

  x_range <- range(x, na.rm = TRUE)
  ## check if all values are the same
  if (!diff(x_range)) {
    stop("Numeric vector has the same values.")
  }
  
  min <- x_range[1]
  max <- x_range[2]
  low  <- range[1] 
  high <- range[2]

  res <- ((high - low) * (x - min)) / (max - min) + low
  res
}

## --------------------------------------------------------------------------
#' Trim white spaces
#' 
#' Trim white spaces.
#' 
#' @param string a character string to be processed.
#' @return returns a trimmed string.
#' @export 
## lwc-29-04-2013: trim white spaces
str_trim <- function(string) {
  string <- gsub("^[ \t]+|[ \t]+$", "", string) ## trim white spaces
  return(string)
}

## ------------------------------------------------------------------------
#' Calculate the percentage of non digits
#' 
#' Calculate the percentage of non digits.
#' 
#' @param mat a data matrix.
#' @return returns an vector.
#' @export  
## wll-15-09-2015: Calculate the percentage of non digits
non_digit <- function(mat) {
  mat <- as.data.frame(mat)
  func <- function(x) round(sum(is.na(x) | is.nan(x)) / length(x), digits = 3) * 100
  res <- apply(mat, 2, func)
  return(res)
}

## ------------------------------------------------------------------------
#' Calculate the non-zero percentage
#' 
#' Calculate the non-zero percentage.
#' 
#' @param mat a data matrix.
#' @return returns an vector.
#' @export  
## wll-09-09-2015: Calculate the non-zero percentage
## wll-22-10-2015: minor changes
## ------------------------------------------------------------------------
non_zero <- function(mat) {
  mat <- as.data.frame(mat)
  res <- sapply(mat, function(x) {
    round(sum(x > .Machine$double.eps) / length(x), digits = 3) * 100
  })
  return(res)
}

## ------------------------------------------------------------------------
#' Calculate the zero percentage
#' 
#' Calculate the zero percentage
#' 
#' @param mat a data matrix.
#' @return returns an vector.
#' @export  
## wll-09-09-2015: Calculate the zero percentage
## wll-22-10-2015: minor changes
is_zero <- function(mat) {
  mat <- as.data.frame(mat)
  res <- sapply(mat, function(x) {
    round(sum(x < .Machine$double.eps) / length(x), digits = 3) * 100
  })
  return(res)
}

## --------------------------------------------------------------------------
#' Get row index of missing values in data frame
#' 
#' Get row index of missing values in data frame.
#' 
#' @param df  data frame being used to check index of missing values.
#' @param vars vector of variable names in `df` for missing values checking.
#' @return returns an vector of missing value index.
#' @export 
## lwc-25-08-2011: get row index of missing values in data frame.
df_na_idx <- function(df, vars) {
  tmp <- df[, vars]
  idx <- complete.cases(tmp)
  return(idx)
}

## ------------------------------------------------------------------------
#' Row binding of a list of data matrix
#' 
#' Row bind of a list of matrix or data frame with the same dimension.
#' 
#' @param list a list of data matrix.
#' @return returns a data mstrix.
#' @examples  
#' df1 <- data.frame(x = c(1, 3), y = c(2, 4))
#' df2 <- data.frame(x = c(5, 7), y = c(6, 8))
#' df  <- list(df1 = df1, df2 = df2)
#' rbind_df(df)
#' 
#' tmp <- vector(mode = "list", length = 3)
#' for (i in 1:3) {
#'   tmp[[i]] <- data.frame(
#'     a = sample(letters, 5, rep = TRUE),
#'     b = rnorm(5), c = rnorm(5)
#'   )
#' }
#' names(tmp) <- c("abc", "def", "ghi")
#' tmp
#' 
#' do.call("rbind", tmp)
#' dplyr::bind_rows(tmp, .id = "column_label")
#' plyr::ldply(tmp, data.frame, .id = "column_label")
#' rbind_df(tmp)
#' @export 
## wll-26-01-2016: 'rbind' matrix or data frame with the same dimension
## from a list.
rbind_df <- function(list) {
  ## convert to data frame
  list <- lapply(list, mat2df)

  ## need to use 'stopifnot' for handling exception
  df <- do.call(rbind, list)
  nam <- rownames(df)

  ## replace the first '.' with '-' to indicate two variables
  nam <- sapply(nam, function(x) {
    x <- gsub("^([^.]+)\\.", "\\1-", x) ## regular expression
  })

  ## use function in 'reshape'
  tmp <- colsplit(nam, "-", c("Var1", "Var2"))

  df <- cbind(tmp, df)
  rownames(df) <- NULL
  df
}

## ------------------------------------------------------------------------
#' Convert matrix to df
#' 
#' Convert matrix to df
#' 
#' @param x a matrix
#' @return returns a data frame.
#' @details from `.matrix_to_df` of package `plyr`.
#' @export 
## lwc-28-01-2016: Convert matrix to df.
mat2df <- function(x) {
  cnames <- colnames(x)
  if (is.null(cnames)) cnames <- rep("", ncol(x))
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  colnames(x) <- cnames
  x
}

## ------------------------------------------------------------------------
#' Transpose a numeric data frame (with/without NAs)
#' 
#' Transpose a numeric data frame (with/without NAs)
#' 
#' @param x a data frame
#' @return returns a transposed data frame.
#' @details It is used mostly for numeric data frame with NAs and return a
##   numeric data frame.
#' @export 
#' @examples  
#' x <- data.frame(group = c(1,1,2,NA,2), val = c(6,4,6,3.1,NA))
#' x <- df_t(x)
## wl-01-06-2011: Transpose a numeric data frame (with/without NAs)
df_t <- function(x) {
  x <- data.frame((t(x)), stringsAsFactors = F)

  dn <- dimnames(x)
  res <- sapply(x, function(i) as.numeric(i))
  dimnames(res) <- dn

  res <- as.data.frame(res)

  ## lwc-04-11-2011: try type.convert
  ## res <- colwise(as.numeric)(x)
  return(res)
}

## ------------------------------------------------------------------------
#' Get dimension matrix of a list of data frame
#' 
#' Retuns a matrix for dimension information of a list of data frame.
#' 
#' @param mat_list  a list of data frame.
#' @return a data frame.
#' @export 
## wll-25-01-2017: Wrapper function for dim matrix of a list of data frame
dim_mat <- function(mat_list) {
  res <- as.data.frame(t(sapply(mat_list, dim)))
  names(res) <- c("row", "col")
  res
}

## ------------------------------------------------------------------------
#' Sort a data frame
#' 
#' sort a date fram and keep the row names.
#' 
#' @param df a data frame.
#' @param ... further parameters to `order`.
#' @return returns a sorted data frame.
#' @export 
#' @details This function is modified from `arrange` of R package `plyr`.
#'   `arrange` deliberately remove rownames, but sometimes it is worth to 
#'   keep them.
## wll-08-02-2017: 'arrange' deliberately remove rownames, but sometimes it
## is worth to keep them.
arrange_row <- function(df, ...) {
  stopifnot(is.data.frame(df))
  ord <- eval(substitute(order(...)), df, parent.frame())
  if (length(ord) != nrow(df)) {
    stop("Length of ordering vectors don't match data frame size",
      call. = FALSE
    )
  }
  ## unrowname(df[ord, , drop = FALSE])  wll: comment this line
  df[ord, , drop = FALSE]
}

## ------------------------------------------------------------------------
#' Convert a list of unequal vectors to a data frame
#' 
#' Convert a list of unequal vectors to a data frame.
#' 
#' @param x a list of vector.
#' @return returns a data matrix.
#' @importFrom plyr rbind.fill
#' @examples 
#' lst <- list(data.frame(a = 1, b = 2), data.frame(a = 2, c = 3, d = 5))
#' do.call(plyr::rbind.fill, lst)
#' vec <- list(var1 = c(a = 1, b = 2), var2 = c(a = 2, c = 3, d = 5))
#' vec2dat(vec)
#' @export 
## wll-09-11-2015: convert a list of unequal vectors to a data frame
vec2dat <- function(x) {

  ## to-do: check if x is a list
  res <- lapply(x, function(y) as.data.frame(t(as.data.frame(y))))
  res <- do.call(plyr::rbind.fill, res)
  ## Beware that no such function 'cbind.fil' in 'plyr'.
	res
}

## -------------------------------------------------------------------------
#' Save PS and TIFF file
#' 
#' Save PS and TIFF file
#' 
#' @param handle figure handle from lattice
#' @param PRE directory and prefix of file name
#' @param fn file name
#' @return no return.
#' @importFrom lattice trellis.device
#' @noRd 
#' @keywords internal
## lwc-23-11-2013: save PS and TIFF file
ps_tiff <- function(handle, PRE = NULL, fn = "temp") {
  trellis.device(
    device = "postscript", color = TRUE,
    file = paste(PRE, fn, ".ps", sep = "_"),
    colormodel = "rgb", width = 10, height = 10, paper = "special",
    ## family = "ComputerModern",
    horizontal = FALSE, onefile = T
  )
  plot(handle)
  dev.off()

  tiff(
    filename = paste(PRE, fn, ".tiff", sep = "_"), width = 10, height = 10,
    units = "in", res = 300, compression = "lzw"
  )
  plot(handle)
  dev.off()
}

## --------------------------------------------------------------------------
#' Extension of `ggdendro` with cluster centre.
#' 
#' Retrun a `ggdendro` extension with cluster centres.
#' 
#' @param hc an object of `hclust`.
#' @param k the number of clusters.
#' @return returns `ggplot2` plot.
#' @details See `ggdendro` extension from https://bit.ly/2QOPSyJ for details.
#' @noRd 
## wl-06-05-2021, Thu: ggdendro extension. From https://bit.ly/2QOPSyJ
## dendro_data_k(): takes a k argument, an integer, specifying the
## number of desired clusters.
dendro_data_k <- function(hc, k) {

  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)

  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }

  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust

  hcdata
}

## --------------------------------------------------------------------------
#' Plot `ggdendro` with clusters
#' 
#' Plot `ggdendro` with clusters.
#' 
#' @param hcdata a data mstrix return from `dendro_data_k`.
#' @param direction where to plot dendrogram.
#' @param fan the orientation.
#' @param branch.size,label.size,nudge.lable,expand.y a set of parameters to 
#'   control `ggdendro`.
#' @return returns a `ggplot2` plot of dendrogram.
#' @importFrom ggdendro label
#' @noRd  
plot_ggdendro <- function(hcdata,
                          direction   = c("lr", "rl", "tb", "bt"),
                          fan         = FALSE,
                          scale.color = NULL,
                          branch.size = 1,
                          label.size  = 3,
                          nudge.label = 0.01,
                          expand.y    = 0.1) {

  direction <- match.arg(direction) # if fan = FALSE
  ybreaks   <- pretty(segment(hcdata)$y, n = 5)
  ymax      <- max(segment(hcdata)$y)

  ## branches
  p <- ggplot() +
    geom_segment(data         =  segment(hcdata),
                 aes(x        =  x,
                     y        =  y,
                     xend     =  xend,
                     yend     =  yend,
                     linetype =  factor(line),
                     colour   =  factor(clust)),
                 lineend      =  "round",
                 show.legend  =  FALSE,
                 size         =  branch.size)

  ## orientation
  if (fan) {
    p <- p +
      coord_polar(direction = -1) +
      scale_x_continuous(breaks = NULL,
                         limits = c(0, nrow(ggdendro::label(hcdata)))) +
      scale_y_reverse(breaks = ybreaks)
  } else {
    p <- p + scale_x_continuous(breaks = NULL)
    if (direction %in% c("rl", "lr")) {
      p <- p + coord_flip()
    }
    if (direction %in% c("bt", "lr")) {
      p <- p + scale_y_reverse(breaks = ybreaks)
    } else {
      p <- p + scale_y_continuous(breaks = ybreaks)
      nudge.label <- - (nudge.label)
    }
  }

  ## labels
  set_labels_params <- function(nbLabels,
                                direction = c("tb", "bt", "lr", "rl"),
                                fan       = FALSE) {
    if (fan) {
      angle       <-  360 / nbLabels * 1:nbLabels + 90
      idx         <-  angle >= 90 & angle <= 270
      angle[idx]  <-  angle[idx] + 180
      hjust       <-  rep(0, nbLabels)
      hjust[idx]  <-  1
    } else {
      angle       <-  rep(0, nbLabels)
      hjust       <-  0
      if (direction %in% c("tb", "bt"))  angle <- angle + 45
      if (direction %in% c("tb", "rl"))  hjust <- 1
    }
    list(angle = angle, hjust = hjust, vjust = 0.5)
  }

  labelParams <- set_labels_params(nrow(hcdata$labels), direction, fan)
  hcdata$labels$angle <- labelParams$angle

  p <- p +
    geom_text(data         =  ggdendro::label(hcdata),
               aes(x       =  x,
                   y       =  y,
                   label   =  label,
                   colour  =  factor(clust),
                   angle   =  angle),
               vjust       =  labelParams$vjust,
               hjust       =  labelParams$hjust,
               nudge_y     =  ymax * nudge.label,
               size        =  label.size,
               show.legend =  FALSE)

  ## colors and limits
  if (!is.null(scale.color)) {
    p <- p + scale_color_manual(values = scale.color)
  }

  ylim <- -round(ymax * expand.y, 1)
  p    <- p + expand_limits(y = ylim)

  p
}

## ------------------------------------------------------------------------
#' Split a vector into chunks
#' 
#' Split a vector into chunks.
#' 
#' @param x an vector.
#' @param n length of chunk.
#' @return return the split chunk.
#' @export  
#' @details For details, see https://bit.ly/2SM4m2G
#' @examples 
#' x <- 1:10
#' n <- 3
#' chunk(x, n)
## wl-01-06-2021, Tue: Split a vector into chunks.
chunk <- function(x,n) {
  res <- split(x, ceiling(seq_along(x) / n))
}

chunk_1 <- function(x,n) {
  split(x, ceiling(seq_along(x)/n))
}

chunk_2 <- function(x,n) {
  split(x, cut(seq_along(x), n, labels = FALSE))
}

chunk_3 <- function(x, n, force.number.of.groups = TRUE, len = length(x),
                    groups = trunc(len/n), overflow = len%%n) {
  if(force.number.of.groups) {
    f1 <- as.character(sort(rep(1:n, groups)))
    f <- as.character(c(f1, rep(n, overflow)))
  } else {
    f1 <- as.character(sort(rep(1:groups, n)))
    f <- as.character(c(f1, rep("overflow", overflow)))
  }

  g <- split(x, f)

  if(force.number.of.groups) {
    g.names <- names(g)
    g.names.ordered <- as.character(sort(as.numeric(g.names)))
  } else {
    g.names <- names(g[-length(g)])
    g.names.ordered <- as.character(sort(as.numeric(g.names)))
    g.names.ordered <- c(g.names.ordered, "overflow")
  }

  return(g[g.names.ordered])
}

## ------------------------------------------------------------------------
#' PCA plot
#'
#' Plot PCA of a matrix or data frame with base R function.
#'
#' @param x a matrix for plotting. Should have row names
#' @param y a factor or character vector specifying the group information of
#'   the row.
#' @param scale a logical flag to indicate whether or not scale \code{x}.
#' @param ep.plot plot ellipse or not.
#' @param ... further parameters passed to `plot`
#'
#' @return a base R graphics object.
#' @details
#'   The `rownames` of `x` will be shown in the plot. `y` is used to indicate
#'   the group of row of `x`.
#' @author Wanchang Lin
#' @examples
#' data(iris)
#' pca_plot(iris[, 1:4], iris[, 5], ep.plot = TRUE)
#' @keywords multivariate
#' @export
## wl-04-11-2015, Wed: PCA plot with rownames of matrix.
## wl-09-11-2021, Tue: Slightly modified from mt::pca.plot.
pca_plot <- function(x, y = NULL, scale = TRUE, ep.plot = FALSE, ...) {

  ## wl-09-11-2021, Tue: fix a bug
  if (is.null(rownames(x))) {
    rownames(x) <- 1:nrow(x)
  }

  ## wll-12-12-2008: Plot ellipse
  elli.plot <- function(x, y, ...) {
    Var <- var(cbind(x, y))
    Mean <- cbind(mean(x), mean(y))
    Elli <- ellipse(Var, centre = Mean, level = 0.975)
    lines(Elli[, 1], Elli[, 2], ...)
  }

  ## wll-29-03-2008: Compute the PCA scores and proportion of variance
  pca.comp <- function(x, scale = FALSE, pcs = 1:2, ...) {
    pca <- prcomp(x, scale. = scale, ...)
    vars <- pca$sdev^2       # i.e. eigenvalues/variance
    vars <- vars / sum(vars) # Proportion of Variance
    names(vars) <- colnames(pca$rotation)
    vars <- round(vars * 100, 2)
    ## dfn  <- paste(names(vars)," (",vars[names(vars)],"%)",sep="")
    dfn <- paste(names(vars), ": ", vars[names(vars)], "%", sep = "")
    x <- data.frame(pca$x)
    x <- x[, pcs]
    vars <- vars[pcs]
    dfn <- dfn[pcs]
    ## names(x) <- dfn

    return(list(scores = x, vars = vars, varsn = dfn))
  }

  ## x <- as.matrix(x)
  pca <- pca.comp(x, scale = scale, pcs = 1:2, ...)

  val <- pca$scores
  val <- val[c("PC2", "PC1")] # Swap position

  plot(val,
    type = "n", cex = 1.0, cex.lab = 1.0, cex.axis = 1.0, cex.main = 1.0,
    ylab = paste("PC1", " (", pca$vars[1], "%)", sep = ""),
    xlab = paste("PC2", " (", pca$vars[2], "%)", sep = ""),
    ...
  )

  if (is.null(y)) {
    text(val[, 1], val[, 2], rownames( x), cex = 0.7, ...)
  } else {
    y <- factor(y)
    text(val[, 1], val[, 2], rownames(x), cex = 0.7, col = unclass(y), ...)

    if (ep.plot) {
      tmp <- as.factor(as.numeric(unclass(y)))
      for (i in levels(tmp)) {
        idx <- tmp == i
        elli.plot(val[idx, 1], val[idx, 2], col = i)
      }
    }
  }
  invisible(NULL)
}

## ------------------------------------------------------------------------
#' @description `mtExtra` provides some general functions and more specific 
#'  functions to R package `mt` for metabolomics data analysis.
#'
#' @section Main functions:
#' The main mtExtra provides more functions for metabolomics data analysis.
#' These functions are statistical analysis and plot methods with `ggplot2`
#' and `lattice`. It uses `tidyverse` as well as `reshape2` and `plyr` 
#' packages.
#'
#' @section Package context:
#' This package follows the principles of the "tidyverse" as mush as possible.
#' It also uses `melt` in `reshape2` if `tidyr` is complicated in some
#' circumstances.
#'
#' @importFrom ellipse ellipse
#' @importFrom graphics lines text
#' @importFrom stats prcomp quantile sd var mahalanobis median qchisq
#' @import tidyr dplyr purrr tibble ggplot2
#' @importFrom reshape2 melt dcast colsplit
#' @importFrom grDevices colorRampPalette dev.off tiff
#' @importFrom stats as.dendrogram as.dist as.hclust complete.cases cor
#'   cutree dist hclust line order.dendrogram p.adjust pf pt qt
#' @importFrom utils data head
#' @importFrom rlang .data
#' @keywords internal
"_PACKAGE"
