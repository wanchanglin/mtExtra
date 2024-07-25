
## -----------------------------------------------------------------------
#' Wrapper function for QC-RLSC
#'
#' Wrapper function for QC-RLSC
#'
#' @param dat A data frame with samples (row) and variables (column).
#' @param cls.qc A vector with string of "qc" and "sample".
#' @param cls.bl A vector with string of batch indicators.
#' @param method Data scaling method. Support "subtract" and "divide"
#' @param intra A logical value indicating whether signal correction is
#'   performed inside each batch.
#' @param opti A logical value indicating whether or not 'span' parameters
#'   are optimised.
#' @param log10 A logical value indicating whether log10 transformation for
#'   the data set or not.
#' @param outl A logical value indicating whether or not QC outlier
#'   detection is employed.
#' @param shift A logical value indicating whether or not batch shift is
#' applied after signal correction.
#' @param ... Other parameter for 'loess'.
#' @return  A corrected data frame.
#' @export
## wl-19-07-2024, Fri: wrapper function for QC-RLSC
qc_rlsc_wrap <- function(dat, cls.qc, cls.bl,
                         method = c("subtract", "divide"),
                         intra = FALSE, opti = TRUE, log10 = TRUE,
                         outl = TRUE, shift = TRUE, ...) {

  ## log transformation
  if (log10) {
    dat[dat ==0] <- NA
    dat <- log10(dat)
  }

  ## QC outlier detection
  if  (outl) {
    dat <- sapply(dat, function(x){ #' x <- dat[, 6, drop = T]
      qc_ind <- grepl("qc", cls.qc, ignore.case =  TRUE, perl = TRUE)
      ## get median of qc data
      qc_dat <- x[qc_ind]
      qc_median <- median(qc_dat, na.rm = TRUE)
      ## assign other data as NA for QC outlier detection
      tmp <- x
      tmp[!qc_ind] <- NA
      ## QC outlier detection
      out_ind <- outl_det_u(tmp)
      ## asisgn outlier as qc median
      x[out_ind] <- qc_median
      return(x)
    }) %>% as_tibble()
  }

  ## QC-RLSC
  if (!intra) {
    res <- qc_rlsc(dat, cls.qc, method = method, opti = opti)
  } else { # wl-16-08-2023, Wed: do signal correction inside each batch
    res <- lapply(levels(cls.bl), function(x){
      idx <- cls.bl %in% x
      tmp <- qc_rlsc(dat[idx,], cls.qc[idx], method = method, opti = opti)
    })
    res <- bind_rows(res)
  }

  ## batch shift. sensitive to missing values
  if (shift) {
    res <- batch_shift(res, cls.bl, overall_average = T) %>% as_tibble()
  }

  ## inverse log10 transformation
  if (log10) {
    res <- 10^res %>% as_tibble()
  }

  return(res)
}

## -----------------------------------------------------------------------
#' QC based robust LOESS signal correction (QC-RLSC)
#'
#' QC based robust LOESS (locally estimated scatterplot smoothing) signal
#' correction (QC-RLSC)
#'
#' @param x A data frame with samples (row) and variables (column).
#' @param y A vector with string of "qc" and "sample".
#' @param method Data scaling method.
#' @param opti A logical value indicating whether or not optimise 'span'
#' @param ... Other parameter for 'loess'.
#' @return  A corrected data frame.
#' @references
#'  Dunn et al. Procedures for large-scale metabolic profiling of serum and
#'  plasma using gas chromatography and liquid chromatography coupled to
#'  mass spectrometry. Nature Protocols 6, 1060–1083 (2011)
#' @export
## wl-14-08-2023, Mon: QC-RLSC
##   Note that the variables are divided by predicted values using qc-based
##   'loess'. Have nothing to do with batch information. Also see:
##     1.) "statTarget"(https://bit.ly/454KIiT)
##     2.) "Rcpm" (https://github.com/ricoderks/Rcpm/)
##   see 'staTarget_shiftcor_v2.R'  and 'Rcpm_qc_rlsc.R'
## wl-20-06-2024, Thu: Add 'method'. The default method is 'subtract'
##   which gets the data back on the original scale.
## wl-02-07-2024, Tue: 'losses' fails in vector which includes only missing
##   values. Should filter based on each batch or lower down threshold
## wl-08-07-2024, Mon: call 'loess_gcv' for optimisation span
qc_rlsc <- function(x, y, method = c("subtract", "divide"), opti = TRUE,
                    ...) {
  method <- match.arg(method)

  ## order for qc
  ind <- grep("qc", y, ignore.case =  TRUE, perl = TRUE)
  if (length(ind) == 0) stop("No QC samples")

  ## order for interpolation (all samples and qcs)
  ord <- 1:length(y)

  ## wl-02-07-2024, Tue: problem if x[ind, i] are all missing values.
  nc <- ncol(x)
  smo <- sapply(1:nc, function(i) { #' i = 60 #' cat(i, "\n")
    ## apply loess to QCs
    if (opti) { # optimise span
      loe <- loess_gcv(ind, x[ind, i, drop = TRUE],
                       span.range = c(.05, .95), ...)
    } else {    # default span: 0.75
      loe <- loess(x[ind, i, drop = TRUE] ~ ind, ...)
    }
    if (T) {              # predict all (sample and qc)
      yf <- predict(loe, ord)
    } else {              # approximate all the samples
      yf <- approx(x = ind, y = loe$fitted, xout = ord)$y
    }
  })

  ## get average values of QC
  mn <- colMeans(x[ind, , drop = F], na.rm = TRUE)
  mn <- matrix(mn, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)

  ## normalise data set
  res <-
    switch(method,
      "subtract" = (x - smo) + mn,
      "divide"   = (x / smo) * mn
    )

  return(as_tibble(res))
}

## -----------------------------------------------------------------------
#' Optimise LOESS's span parameter using generalized cross-validation (GCV)
#'
#' Optimise LOESS's span parameter using generalized cross-validation (GCV)
#'
#' @param x, y vectors giving the coordinates of the points in the scatter
#'   plot.
#' @param span.range a vector for span optimisation range.
#' @param ... Other parameter for 'loess'.
#' @return  A optimised loess model.
#' @importFrom stats optimize resid update
#' @noRd
#' @keywords internal
## wl-08-07-2024, Mon: Optimal loess model with GCV
## Modified from https://bit.ly/3zBo3Qn
loess_gcv <- function(x, y, span.range = c(.05, .95), ...){

  ## ---------------------------------------------------------------------
  ## GCV for loess.
  ## Modified from code by Michael Friendly
  ## http://tolstoy.newcastle.edu.au/R/help/05/11/15899.html
  loessGCV <- function(x) {
    if (!(inherits(x, "loess")))
      stop("Error: argument must be a loess object")
    span <- x$pars$span
    n <- x$n
    traceL <- x$trace.hat
    sigma2 <- sum(resid(x)^2) / (n - 1)
    gcv <- n * sigma2 / (n - traceL)^2
    result <- list(span = span, gcv = gcv)
    result
  }

  ## ---------------------------------------------------------------------
  ## Wrapper function to link optimize() with loessGCV()
  bestLoess <- function(model, spans = c(.05, .95)) {
    f <- function(span) {
      mod <- update(model, span = span)
      loessGCV(mod)[["gcv"]]
    }
    result <- optimize(f, spans)
    result
  }

  ## ---------------------------------------------------------------------
  ## default span = 0.75
  loe <- loess(y ~ x, ...)

  ## optimise model span
  opti <- bestLoess(loe, spans = span.range)

  ## optimal model
  loe_gcv <- update(loe, span = opti$minimum)

  return(loe_gcv)
}

## ------------------------------------------------------------------------
#' Quality control–based robust LOESS signal correction (QC-RLSC)
#'
#' Signal correction based on quality control.
#'
#' @param x a data matrix.
#' @param y a categorical data for batch/block information.
#' @param ... other parameter for 'loess' such as 'span' (the parameter
#' which controls the degree of smoothing. Default: 0.75) and
#' 'degree' (the degree of the polynomials to be used. Default: 2).
#' @return  a corrected data matrix.
#' @importFrom stats loess predict approx
#' @noRd
#' @keywords internal
#' @references 
#'  Warwick B Dunn, et.al, Procedures for large-scale metabolic profiling of
#'  serum and plasma using gas chromatography and liquid chromatography
#'  coupled to mass spectrometry, Nature Protocols, 6:1060–1083 (2011)
## wl-14-08-2023, Mon: QC-RLSC 
## Note than the variables are divided by predicted values using qc-based
## 'loess'. Do not use any block/batch information. 
## wl-21-08-2023, Mon: use 'grep'
qc_rlsc_1 <- function(x, y, ...) {

  ## order for QCs
  ind <- grep("qc", y, ignore.case =  TRUE, perl = TRUE)
  if (length(ind) == 0) stop("No QC samples")

  ## order for interpolation (all samples and QCs)
  ord <- 1:length(y)

  ## number of variables
  nc <- ncol(x)

  res <- x
  for (i in 1:nc) { #' i = 3
    ## apply loess to QCs
    loe <- loess(x[ind, i, drop = TRUE] ~ ind, ...)
    ## yf <- predict(loe, ind)     # only predict qc
    if (F) {                       # predict all (sample and qc)
      yf <- predict(loe, ord)
    } else {                       # approximate all the samples
      yf <- approx(x = ind, y = loe$fitted, xout = ord)$y
    }
    res[, i] <- x[, i, drop = TRUE] / yf
  }
  return(res)
}

## -----------------------------------------------------------------------
#' Data filtering based on "qc" missing values
#'
#' @param x a data matrix.
#' @param y a character string with contents of "sample", "qc" and "blank".
#' @param thres threshold of missing values. Features less than this
#'  threshold will be kept.
#' @return a list of with contents: \itemize{
#'  \item dat the filtered data matrix
#'  \item idx a logical vector of index for keeping features.
#' }
#' @export
## wl-14-06-2011, Tue: Filter features based on missing values in QC
## wl-01-07-2024, Mon: Review and minor change. different from 'qc_filter'
##   in 'mtExtra'.
mv_filter_qc <- function(x, y, thres = 0.3) {
  tmp <- grep("qc", y, ignore.case =  TRUE, perl = TRUE)
  idx <- mv_filter(x[tmp, , drop = FALSE], thres = thres)$idx
  x <- x[, idx, drop = FALSE]
  return(list(dat = x, idx = idx))
}

## -----------------------------------------------------------------------
#' Data filtering based on quality assurance in QC samples
#'
#' Data filtering based on quality assurance in QC samples
#'
#' @param x A data matrix.
#' @param y A character string with contents of "sample", "qc" and "blank".
#' @param tol A torrance for quality assurance
#' @param thres A threshold for quality assurance
#' @return a list of with contents: \itemize{
#'  \item dat the filtered data matrix
#'  \item idx a logical vector of index for keeping features.
#' }
#' @export
## wl-31-05-2011, Tue: Filter based on quality assurance in QC samples.
## wl-21-06-2024, Fri: Minor changes
qc_pass <- function(x, y, tol = 0.3, thres = 0.4) {
  ind <- grep("qc", y, ignore.case =  TRUE, perl = TRUE)
  if (length(ind) == 0) stop("No QC samples")

  x.qc <- x[ind, , drop = F]
  nc <- ncol(x.qc)
  nr <- nrow(x.qc)

  mean.qc <- colMeans(x.qc, na.rm = TRUE)
  upper_tol <- mean.qc + mean.qc * tol
  lower_tol <- mean.qc - mean.qc * tol

  ## get the sum of count
  count <- sapply(1:nc, function(i) {
    tmp <- x.qc[, i]
    ## tmp[is.na(x)] <- 0.0
    n.upp <- tmp > upper_tol[i]
    n.low <- tmp < lower_tol[i]
    res <- sum((n.upp | n.low), na.rm = TRUE)
    ## res <- sum(n.upp, na.rm=TRUE) + sum(n.low, na.rm=TRUE)
    res <- res / nr
  })

  idx <- (count <= thres)
  x <- x[, idx, drop = FALSE]

  return(list(dat = x, idx = idx))
}
