## wl-09-11-2021, Tue: gather all general functions from 2015

#' MS/NMR data filtering
## wl-03-06-2021, Thu: wrapper for filtering MS/NMR data
#'
#' Wrapper function for filtering MS/NMR data.
#'
#' @param x a list of metabolomics data, including intensity data and peak
#'   information.
#' @param method the method for filtering.
#' @param ... further parameters to be passed to `method`.
#' @return a filtered list of metabolomics data.
#' @export
dat_filter <- function(x, method = "var_filter", ...) {

  method <-
    if (is.function(method)) {
      method
    } else if (is.character(method)) {
      get(method)
    } else {
      eval(method)
    }

  #' filtering features
  idx <- method(x$data, ...)$idx
  x$data <- x$data[, idx]
  x$peak <- x$peak[idx, ]

  return(x)
}

#' Filtering variable based on IQR
## wl-04-10-2016, Tue: Variable filtering(column-wise) based on IQR.
#'
#' Perform variable filtering(column-wise) based on IQR. This function is
#' slightly modified from `varFilter` in `genefilter`.
#'
#' @param x a matrix or data frame.
#' @param var.func filtering method.
#' @param var.cutoff a numeric value for the threshold.
#' @param filterByQuantile	a logical indicating whether `var.cutoff` is to be
#'   a quantile of all `var.func` values or as an absolute value.
#'
#' @return a list with contents of filtered data and the filtered index.
#' @export
iqr_filter <- function(x, var.func = "IQR", var.cutoff = 0.25,
                       y = TRUE) {
  vars <- apply(x, 2, var.func)

  if (filterByQuantile) {
    if (0 < var.cutoff && var.cutoff < 1) {
      quant <- quantile(vars, probs = var.cutoff)
      selected <- !is.na(vars) & vars > quant
    } else {
      stop("Cutoff Quantile has to be between 0 and 1.")
    }
  } else {
    selected <- !is.na(vars) & vars > var.cutoff
  }
  x <- x[, selected, drop = F]

  return(list(dat = x, idx = selected))
}

#' ------------------------------------------------------------------------
#' wl-04-05-2021, Tue: filtering variables with sd values close to zero
#' wl-03-06-2021, Thu: return filtering data as well
sd_zero_filter <- function(x, na.rm = FALSE) {
  #' take off the columns with the same values.
  idx <- which(apply(x, 2, sd, na.rm = na.rm) > .Machine$double.eps)
  x <- x[, idx, drop = F]
  return(list(dat = x, idx = idx))
}

#' ------------------------------------------------------------------------
#' wl-26-02-2013: Filtering by SD
#' - Calculate the standard deviation for every feature across all samples.
#'   - That gives one SD value per feature
#' - Select, say, 5-50% of the features with the highest SD
#' - That leaves the features that vary a lot
#'   - Between groups
#'   - Randomply fluctuating expression
#' wl-17-06-2021, Thu: need to test and debug
sd_filter <- function(x, sig.cutoff = 0.95, na.rm = FALSE) {
  sds <- apply(x, 2, sd, na.rm = na.rm)
  cutsd <- quantile(sds, sig.cutoff)
  idx <- sds > cutsd
  x <- x[, idx, drop = F]
  return(list(dat = x, idx = idx))
}

#' -----------------------------------------------------------------------
#' wl-14-06-2011: Filter features based on the percentage of missing values
#' wl-17-06-2021, Thu: several version but this one is simple. Need to test
mv_filter <- function(x, thres = 0.3) {

  if (!(is.matrix(x) || is.data.frame(x))) {
    stop("\n'data' is not matrix or data frame.\n")
  }

  thres <- thres * nrow(x)

  #' get number of Na, NaN and zero in each of feature/variable
  count <- apply(x, 2, function(y) {
    tmp <- sum(is.na(y) | is.nan(y) | (y == 0))
  })
  #' na.mat <- is.na(data)
  #' count  <- apply(na.mat,2,sum)

  #' index of features whose number of MVs are less than threshold
  idx <- count <= thres

  x <- x[, idx, drop = F]
  return(list(dat = x, idx = idx))
}

#' Relative standard deviation.
## lwc-02-06-2011: Relative standard deviation.
#'
#' Perform Relative Standard Deviation (RSD).
#'
#' @param x an vector, matrix or data frame.
#' @param na.rm	a logical value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @return  RSD values
#' @examples
#' data(iris)
#' rsd(iris[, 1:4])
#' @export
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

#' PCA plot
## wl-04-11-2015, Wed: PCA plot with rownames of matrix.
## wl-09-11-2021, Tue: Slightly modified from mt::pca.plot.
#'
#' Plot PCA of a matrix or data frame with base R function.
#'
#' @param x a matrix for plotting. Should have row names
#' @param y a factor or character vector specifiing the group information of
#'   the row.
#' @param scale a logical flag to indicate whether or not scale \code{x}.
#' @param ep.plot plot ellipse or not.
#' @param ... further parameters passed to `plot`
#'
#' @return a base R graphics object.
#' @author Wanchang Lin
#' @examples
#' data(iris)
#' pca_plot(iris[, 1:4], iris[, 5], ep.plot = TRUE)
#' @keywords multivariate
#' @export
#' @importFrom ellipse ellipse
#' @importFrom graphics lines text
#' @importFrom stats prcomp quantile sd var
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
    vars <- pca$sdev^2 #' i.e. eigenvalues/variance
    vars <- vars / sum(vars) #' Proportion of Variance
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
  val <- val[c("PC2", "PC1")] #' Swap position

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

#' IQR-based Variable filtering
## wl-04-10-2016, Tue: Variable filtering(row-wise) based on IQR.
##  Modified from `varFilter` in `genefilter`.
## wl-19-01-2017, Thu: fix a bug. Use drop=F
## wl-15-10-2020, Thu: return index as well
#'
#' Perform variable filtering(row-wise) based on IQR.
#'
#' @param x a matrix for processing.
#' @param var.func a method for filtering
#' @param var.cutoff a threshold for cut-off
#' @param filterByQuantile a binary value to control filter by quantitle.
#' @return a list of with contents: \itemize{
#'  \item x the filtered data matrix
#'  \item idx a vector of filtring index.
#' }
#'
#' @export
#' @examples
#' x <- matrix(rnorm(20*3), ncol = 3)
#' var_filter(x)
var_filter <- function(x, var.func = "IQR", var.cutoff = 0.25,
                       filterByQuantile = TRUE) {
  vars <- apply(x, 1, var.func)

  if (filterByQuantile) {
    if (0 < var.cutoff && var.cutoff < 1) {
      quant <- quantile(vars, probs = var.cutoff)
      selected <- !is.na(vars) & vars > quant
    } else {
      stop("Cutoff Quantile has to be between 0 and 1.")
    }
  } else {
    selected <- !is.na(vars) & vars > var.cutoff
  }
  x <- x[selected, , drop = F]

  return(list(dat = x, idx = selected))
}

#' vector normalisation
## wl-22-09-2020, Tue: vector normalisation
#'
#' Normalise a numeric vector
#'
#' @param x a vector
#' @param method the normalisation method to be used.
#' @param scale a logical to scale `x` or not
#' @return returns a normalised vector
#' @export
#' @examples
#' x <- c(2, 3, 4, 5, 6, 7, 8, 9, 50, 50)
#' vec_norm(x, method = "median", scale = TRUE)
vec_norm <- function(x, method = "median", scale = TRUE) {
  method <- match.arg(method, c("median", "mean"))
  method <- get(method)
  center <- method(x, na.rm = T)
  x <- x - center
  if (scale) {
    x <- x / sd(x, na.rm = T)
  }
  return(x)
}

