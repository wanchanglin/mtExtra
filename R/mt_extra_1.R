## wl-09-11-2021, Tue: gather all general functions from 2015

#' PCA plot
## wl-04-11-2015, Wed: PCA plot with rownames of matrix.
## wl-09-11-2021, Tue: Slightly modified from mt::pca.plot.
#'
#' Plot PCA of a matrix or data frame with base R function
#'
#' @param x a matrix for plotting. Should have row names
#' @param y a factor or character vector specifiing the group information of
#' the row.
#' @param scale a logical flag to indicate whether or not scale \code{x}.
#' @param ep.plot plot ellipse or not.
#'
#' @return a base R graphics object.
#' @keywords multivariate
#' @export
#' @examples
#' data(iris)
#' pca_plot(iris[, 1:4], iris[, 5], ep.plot = TRUE)
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
#' Variable filtering(row-wise) based on IQR.
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
#' @param cale a logical to scale `x` or not
#' @return returns a normalised vector
#' @export
#' @examples
#' x <- c(2, 3, 4, 5, 6, 7, 8, 9, 50, 50)
#' vec_norm(x, method = "median", scale = T)
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

