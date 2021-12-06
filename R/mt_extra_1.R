
## ------------------------------------------------------------------------
#' Data matrix summary
#' 
#' Summarise a data matrix. This is a wrapper function for any summary 
#' function on an vector.
#' 
#' @param x a matrix-like object.
#' @param method summary method for an vector.
#' @param ... further parameters for `method`.
#' @return a summarised table.
#' @export 
#' @examples 
#' library(dplyr)
#' library(tidyr)
#' library(purrr)
#' iris %>% group_by(Species) %>% group_modify(~ dat_summ(., method = mean))
#' iris %>% group_by(Species) %>% group_modify(~ dat_summ(., method = vec_stats))
#' iris %>% group_by(Species) %>% do(dat_summ(., method = vec_segment))
#' iris %>% dat_summ(method = sd)
## wl-20-10-2021, Wed: Data matrix summary using tidyverse
## wl-24-11-2021, Wed: use '.id' for rownames
## wl-29-11-2021, Mon: 'method'  a vector or data matrix. 
dat_summ <- function(x, method = mean, ...) {
  x %>% 
    as_tibble() %>% 
    select(where(is.numeric)) %>% 
    map_dfr(method, .id = "var") 
}

## -----------------------------------------------------------------------
#' Vector summary
#' 
#' Calculate the statistical summary of a vector.
#' 
#' @param x a numeric vector.
#' @param na.rm remove NA or not.
#' @param conf.interval a numeric value for confidence interval.
#' @return 
#'   a vector of summary consisting:
#'   - number of vector length
#'   - vector mean
#'   - vector standard derivation   
#'   - standard error of mean
#'   - confidence interval
#' @details Can be used for error bar plotting. Modify from
#'   https://bit.ly/3onsqot
#' @family vector stats functions
#' @examples  
#' library(dplyr)
#' library(tidyr)
#' library(purrr)
#' 
#' iris %>% dat_summ(method = vec_stats)
#' iris %>% group_by(Species) %>% group_modify(~ dat_summ(., method = vec_stats))
#' iris %>%
#'   pivot_longer(cols = !Species, names_to = "var") %>% 
#'   group_nest(Species, var) %>%
#'   mutate(map_dfr(data, ~ vec_stats(.x$value))) %>%
#'   select(!data)
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
#' confidence interval.
#' 
#' @param x an vector.
#' @param bar a character string, supporting "SD", "SE" and "CI".
#' @return an vector including lower, center and upper values.
#' @examples
#' library(plyr)
#' library(reshape2)
#' library(dplyr)
#' library(tidyr)
#' library(purrr)
#' 
#' vec_segment(iris[,1])
#' 
#' ## Use 'plyr' and 'reshape2' for group 
#' mat <- melt(iris)
#' ddply(mat, .(Species,variable), function(x,bar) {
#'   vec_segment(x$value, bar = bar)
#' }, bar = "SD")
#' 
#' ## Use 'tidyverse' for group 
#' iris %>%
#'   pivot_longer(cols = !Species, names_to = "variable") %>% 
#'   group_by(Species, variable) %>%
#'   nest() %>%
#'   mutate(map_dfr(.x = data, .f = ~ vec_segment(.x$value))) %>%
#'   select(!data)
#'	
#' iris %>%
#'   pivot_longer(cols = !Species, names_to = "variable") %>% 
#'   group_nest(Species, variable) %>%
#'   mutate(map_dfr(data, ~ vec_segment(.x$value))) %>%
#'   select(!data)
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

## ------------------------------------------------------------------------
#' Transform data
#' 
#' Perform data set transformation
#'
#' @param x a matrix, data frame or vector.
#' @param method transformation method, including: "center", "auto", "range",
#'   "pareto", "vast", "level", "log", "log10", "sqrt" and "asinh". 
#' @param na.rm	a logical value indicating whether NA values should be 
#'   stripped before the computation proceeds.
#' @param add a shift value for log transformation.
#' @return transformed data
#' @references 
#'  Berg, R., Hoefsloot, H., Westerhuis, J., Smilde, A. and Werf, M. (2006),
#'  Centering, scaling, and transformations: improving the biological
#'  information content of metabolomics data, \emph{BMC Genomics}, 7:142
#' @name trans
#' @examples
#'  data(iris)
#' 
#'  ## transform an vector
#'  vec <- iris[, 1]
#'  dat_trans(vec, method = "auto")
#'  vec_trans(vec, method = "auto")
#' 
#'  ## transform a data frame
#'  mat <- iris[, 1:4]
#'  dat_trans(mat, method = "log")
#' 
#'  ## transform data frame under different conditions
#'  plyr::ddply(iris, ("Species"), function(x, method) {
#'    dat_trans(x[, 1:4], method = method)
#'  }, method = "range")
#' 
#' ## use 'tidyverse'
#' library(dplyr)
#' library(tidyr)
#' library(purrr)
#' 
#' ## transform whole data set
#' iris %>% mutate(across(where(is.numeric), ~ vec_trans(., method = "range")))
#' 
#' ## transform data set within groups
#' iris %>% 
#'   group_by(Species) %>%
#'   mutate(across(where(is.numeric), ~ vec_trans(., method = "range")))
NULL

## ------------------------------------------------------------------------
#' @export 
#' @rdname trans 
#' @order 1
## wll-29-07-2015: Transform a data matrix
## Note: The programming structure is from ststs's 'sd' function.
dat_trans <- function(x, method = "auto", na.rm = TRUE, add = 1) {
  if (is.matrix(x)) {
    res <- apply(x, 2, vec_trans, method = method, na.rm = na.rm, add = add)
  } else if (is.vector(x)) {
    res <- vec_trans(x, method = method, na.rm = na.rm, add = add)
  } else if (is.data.frame(x)) {
    res <- sapply(x, vec_trans, method = method, na.rm = na.rm, add = add)
  } else {
    res <- vec_trans(as.vector(x), method = method, na.rm = na.rm, add = add)
  }

  return(res)
}

## ------------------------------------------------------------------------
#' @export 
#' @rdname trans
#' @order 2
## wll-29-07-2015: transform an vector
vec_trans <- function(x, method = "auto", na.rm = TRUE, add = 1) {
  method <- match.arg(method, c(
    "center", "auto", "range", "pareto", "vast",
    "level", "log", "log10", "sqrt", "asinh"
  ))

  ## basic functions
  me <- function(x) mean(x, na.rm = na.rm)
  se <- function(x) sd(x, na.rm = na.rm)
  mx <- function(x) max(x, na.rm = na.rm)
  mn <- function(x) min(x, na.rm = na.rm)
  ## sm  <- function(x) sum(x, na.rm = na.rm)

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
#' Melt a numeric data matrix to long format
#'
#' Reshape a matrix or data frame to long format with row names and column
#' names in two columns.
#' 
#' @param x a matrix or data frame
#' @return  a `tibble` object
#' @details `reshape2::melt` keeps the rownames when melting a matrix, but
#'   not when melting a data frame. This function keeps rownames for both 
#'   matrix and data frame.
#' @export
#' @examples
#' class(mtcars)     # data.frame
#' reshape2::melt(mtcars)
#' reshape2::melt(as.matrix(mtcars))
#' dat2long(mtcars)
## wl-04-12-2021, Sat: write this function with tidyr
dat2long <- function(x) {
  x %>% 
    as.data.frame() %>%
    rownames_to_column(var = "row_name") %>%
    pivot_longer(-c("row_name"), names_to = "col_name")
    # pivot_longer(!row_name, names_to = "col_name")
    # pivot_longer(!.data$row_name, names_to = "col_name")
}

## -----------------------------------------------------------------------
#' Convert a symmetric table(short format) to long format
#' 
#' Convert a symmetric table(short format) to long format
#' 
#' @param x a symmetric matrix-like data set.
#' @param tri triangular being used.
#' @return  a data frame of pair-wise comparison.
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

  ## res <- reshape2::melt(x)
  ## res <- res[complete.cases(res), ]
  ## colnames(res) <- c("com1", "com2", "var")

  idx <- which(!is.na(x), arr.ind = T)
  fs1 <- rownames(x)[idx[, 1]]
  fs2 <- colnames(x)[idx[, 2]]
  res <- data.frame(
    com1 = fs1, com2 = fs2, var = x[idx],
    stringsAsFactors = FALSE
  )
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

## ------------------------------------------------------------------------
#' Wrapper function for plotting classification results
#' 
#' This function plots accuracy, AUC and margin (aam) of classification 
#' results from package \code{mt}.
#' 
#' @param aam_list a data matrix of classification results.
#' @param fig_title a string of figure title
#' @return an object of class `ggplot2`
#' @examples 
#' aam <- mtExtra:::aam
#' plot_aam(aam)
#' @seealso `aam.mcl` in R package `mt` for how to get accuracy, AUC and 
#'   margin. 
#' @export 
## lwc-06-05-2011: Wrapper function for plotting the aam results
##    (Accuracy, AUC and Margin).
## wl-09-11-2021, Tue: implement with ggplot2
## wl-03-12-2021, Fri: replace melt with functions in tidyverse
plot_aam <- function(aam_list, fig_title = "Accuracy, AUC and Margin") {

  if (!is.list(aam_list)) aam_list <- list(aam_list)

  ## lst <- reshape2::melt(aam)
  ## lst <- lst[complete.cases(lst), ]

  lst <- lapply(aam_list, function(x) {
    x <- as.data.frame(x)
    x <- rownames_to_column(x, var = "rn")
  }) %>%
    bind_rows(.id = "data") %>%
    pivot_longer(-c("data", "rn"), names_to = "var") %>%
    drop_na() %>%
    rename("classifier" = rn, "assessment" = var)

  aam.p <-
    ggplot(lst, aes(y = value, x = data, color = classifier,
                    group = classifier)) +
    geom_line(aes(linetype = classifier)) +
    geom_point(aes(shape = classifier)) +
    ggtitle(fig_title) +
    facet_wrap(~ assessment) +
    coord_flip()
  return(aam.p)
}

## ------------------------------------------------------------------------
#' Plot adjusted p-values
#' 
#' Plot the adjusted p-values using `ggplot2`
#' 
#' @param pval_list a data matrix or a list of data matrix of p-value correction
#' @return an object of class `ggplot2`
#' @export  
#' @examples 
#' pval <- mtExtra:::pval
#' plot_pval(pval)
## wl-03-10-2021, Sun: the output of melt will be different if x is
##  either matrix or data frame.
## wl-09-11-2021, Tue: plot adjusted p-values
## wl-11-11-2021, Thu: deal with single matrix/data.frame
## wl-02-12-2021, Thu: use 'usethis::use_data(pval, y, internal = TRUE)' to 
##  get 'sysdata.rda'
plot_pval <- function(pval_list) {

  if (!is.list(pval_list)) pval_list <- list(pval_list)

  tmp <- lapply(pval_list, function(x) {
    rownames(x) <- 1:nrow(x)
    x
  })
  ## tmp <- reshape2::melt(tmp)

  ## wl-16-11-2021, Tue: use '.data' to prevent warning in 'rcmd check'
  tmp <- lapply(tmp, function(x) {
    x <- as_tibble(x) %>%
    rownames_to_column(var = "rn") %>%
    mutate(rn = as.numeric(.data$rn))
  })
  tmp <- tmp %>%
    bind_rows(.id = "data") %>%
    pivot_longer(-c("data", "rn"), names_to = "variable") %>%
    drop_na()

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

## -----------------------------------------------------------------------
#' Heatmap with dendrograms with ggplot2
#' 
#' Plot heatmap of a data matrix using `ggplot2`. This function is modified 
#' from https://bit.ly/2UUnY2L.
#' 
#' @param mat a data mstrix to be plotted.
#' @param row.dend plot row dendrogram or not.
#' @param col.dend plot column dendrogram ot not.
#' @param row.dend.right a logical value to indicate the position of row 
#'   dendrogram.
#' @param colors a vector of colours for heatmap.
#' @param font.size label font size.
#' @param x.rot plot rotate degree.
#' @param legend.title legend title.
#' @param dist.method distance method.
#' @param clust.method cluster method.
#' @param dend.line.size dendrogram line size.
#' @return an object of class `ggplot2`.
#' @importFrom cowplot axis_canvas insert_yaxis_grob ggdraw
#' @importFrom ggdendro dendro_data segment
#' @export 
#' @examples 
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

  mat <- as.data.frame(mat)
  # data_m <- cbind(rn = rownames(mat), mat)
  # data_m <- reshape2::melt(data_m, id.vars = "rn")
  data_m <- rownames_to_column(mat, var = "rn")
  data_m <- data_m %>% 
   pivot_longer(cols = !rn, names_to = "variable", values_to = "value")

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
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt dcast colsplit
#' @importFrom grDevices colorRampPalette dev.off tiff
#' @importFrom stats as.dendrogram as.dist as.hclust complete.cases cor
#'   cutree dist hclust line order.dendrogram p.adjust pf pt qt
#' @importFrom utils data head
#' @keywords internal
"_PACKAGE"

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr place-holder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
## wl-02-12-2021, Thus: get from running 'usethis::use_pipe()'
NULL

utils::globalVariables(c(
  'where',
  'angle',
  'classifier',
  'clus',
  'clust',
  'edge',
  'mod',
  'name',
  'rn',
  'size',
  'value',
  'variable',
  'x',
  'xend',
  'y',
  'yend'
))

##  1) dat_summ
##  2) vec_stats
##  3) vec_segment
##  4) dat_trans
##  5) vec_trans
##  6) dat2long
##  7) sym2long
##  8) range_scale
##  9) plot_aam
## 10) plot_pval
## 11) gg_heat_dend
## 12) heat_dend
## 13) pca_plot
