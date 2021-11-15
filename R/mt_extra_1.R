## wl-09-11-2021, Tue: gather all general functions from 2015

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
## lwc-19-01-2016: get the count number of feature selectors in top-k
##  feature orders.
## Use 'reshape' (not 'reshape2'). 'reshape2::dcast()' margins = "grand_col"
##  has problem.
feat_count_1 <- function(fs.ord, top.k = 30) {
  fs.ord <- fs.ord[1:top.k, ]

  ## Use melt and cast to get frequency table
  tmp <- melt(fs.ord, id = NULL)
  # tmp <- reshape::melt(fs.ord, id = NULL)
  # cons <- reshape::cast(tmp, value ~ variable,
  cons <- dcast(tmp, value ~ variable,
    fun.aggregate = length,
    margins = "grand_col"
  )
  head(cons)
  names(cons)[which(names(cons) == "(all)")] <- "count"
  cons <- arrange(cons, desc(count))
}

## ------------------------------------------------------------------------
#' MS/NMR data filtering
#'
#' Wrapper function for filtering MS/NMR data.
#'
#' @param x a list of metabolomics data, including intensity data and peak
#'   information.
#' @param method the method for filtering.
#' @param ... further parameters to be passed to `method`.
#' @return  a filtered list of metabolomics data.
#' @export
## wl-03-06-2021, Thu: wrapper for filtering MS/NMR data
dat_filter <- function(x, method = "var_filter", ...) {

  method <-
    if (is.function(method)) {
      method
    } else if (is.character(method)) {
      get(method)
    } else {
      eval(method)
    }

  ## filtering features
  idx <- method(x$data, ...)$idx
  x$data <- x$data[, idx]
  x$peak <- x$peak[idx, ]

  return(x)
}

## ------------------------------------------------------------------------
#' Filtering variable based on IQR
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
#' @return a list of with contents: \itemize{
#'  \item x the filtered data matrix
#'  \item idx a vector of filtering index.
#' }
#' @examples
#' x <- matrix(rnorm(20*3), ncol = 3)
#' iqr_filter(x)
#' @export
## wl-04-10-2016, Tue: Variable filtering(column-wise) based on IQR.
## wl-19-01-2017, Thu: fix a bug. Use drop=F
## wl-15-10-2020, Thu: return index as well
iqr_filter <- function(x, var.func = "IQR", var.cutoff = 0.25, 
                       filterByQuantile = TRUE) {
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

## ------------------------------------------------------------------------
#' Filtering variables with sd values close to zero
#'
#' Filter variables with sd values close to zero
#'
#' @param x an vector, matrix or data frame.
#' @param na.rm	a logical value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @return a list of with contents: \itemize{
#'  \item x the filtered data matrix
#'  \item idx a vector of filtering index.
#' }
#' @export
## wl-04-05-2021, Tue: filtering variables with sd values close to zero
## wl-03-06-2021, Thu: return filtering data as well
sd_zero_filter <- function(x, na.rm = FALSE) {
  ## take off the columns with the same values.
  idx <- which(apply(x, 2, sd, na.rm = na.rm) > .Machine$double.eps)
  x <- x[, idx, drop = F]
  return(list(dat = x, idx = idx))
}

## ------------------------------------------------------------------------
#' Filtering variable based SD
#'
#' This function calculates the standard deviation for every feature across all
#' samples and keep the features with large quantile controlled by `sig.cutoff`.
#'
#' @param x an vector, matrix or data frame.
#' @param sig.cutoff a numeric value (default 0.95) for `quantile`.
#' @param na.rm	a logical value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @return a list of with contents: \itemize{
#'  \item x the filtered data matrix
#'  \item idx a vector of filtering index.
#' }
#' @export
## wl-26-02-2013: Filtering by SD
## - Calculate the standard deviation for every feature across all samples.
##   - That gives one SD value per feature
## - Select, say, 5-50% of the features with the highest SD
## - That leaves the features that vary a lot
##   - Between groups
##   - Randomply fluctuating expression
## wl-17-06-2021, Thu: need to test and debug
sd_filter <- function(x, sig.cutoff = 0.95, na.rm = FALSE) {
  sds <- apply(x, 2, sd, na.rm = na.rm)
  cutsd <- quantile(sds, sig.cutoff)
  idx <- sds > cutsd
  x <- x[, idx, drop = F]
  return(list(dat = x, idx = idx))
}

## -----------------------------------------------------------------------
#' Filtering features based on the percentage of missing values
#'
#' This function calculates the percentage of missing values and keeps those
#' features with missing values percentage less than the designed threshold.
#'
#' @param x an vector, matrix or data frame.
#' @param thres a threshold to drop off the variables.
#' @return a list of with contents: \itemize{
#'  \item x the filtered data matrix
#'  \item idx a vector of filtering index.
#' }
#' @export
## wl-14-06-2011: Filter features based on the percentage of missing values
## wl-17-06-2021, Thu: several version but this one is simple. Need to test
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

## ------------------------------------------------------------------------
#' Relative standard deviation
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
## lwc-02-06-2011: Relative standard deviation.
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
#' Vector normalisation
#'
#' Perform normalisation on a numeric vector.
#'
#' @param x a vector
#' @param method the normalisation method to be used.
#' @param scale a logical to scale `x` or not.
#' @return returns a normalised vector.
#' @export
#' @examples
#' x <- c(2, 3, 4, 5, 6, 7, 8, 9, 50, 50)
#' vec_norm(x, method = "median", scale = TRUE)
## wl-22-09-2020, Tue: vector normalisation
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
##  pval.ad - a list of matrix of p-values and their correction
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
  tmp <- lapply(tmp, function(x) {
    x <- as_tibble(x) %>%
    rownames_to_column() %>%
    mutate(rowname = as.numeric(rowname))
  })
  tmp <- tmp %>%
    bind_rows(.id = "data") %>%
    pivot_longer(-c(data, rowname), names_to = "variable") %>%
    filter(complete.cases(.))

  ## wl-24-05-2021, Mon: combine two legends with the same legend title
  p <- ggplot(tmp, aes(x = rowname, y = value)) +
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
## wl-07-10-2021, Thu: Correlation network analysis
## wl-08-10-2021, Fri: Call function 'net_graph'
##  Note: This function is different from 'corr_net' which is for
##        Bipartite-mode correlation network.
##  mat - data matrix for correlation analysis
##  thres - correlation coeeficient threshold for network analysis
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
## wl-07-05-2021, Fri: Bipartite/two-mode correlation network
## wl-19-06-2021, Thu: change family font from 'Arial' to 'sans'. So no
##  'extrafont' is needed for PDF output.
##  co_mat - correlation coffeficient matrix
##  thres - correlation coeeficient threshold for network analysis
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
#' @param colors a vector of colours for heatmap.
#' @param font.size label font size.
#' @param x.rot plot rotate degree.
#' @param legend.title lengend title.
#' @param dist.metode distance method.
#' @param clust.method cluster method.
#' @param dend.line.size dendrogram line size.
#' @return retuns an object of class `ggplot2`.
#' @importFrom cowplot axis_canvas insert_yaxis_grob ggdraw
#' @export 
#' @examples 
#' library(ggplot2)
#' gg_heat_dend(mtcars)
## wl-24-11-2020, Tue: Heatmap with dendrograms with ggplot2
## Modified from https://bit.ly/2UUnY2L <br/>
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
  data_m <- cbind(rowname = rownames(mat), mat)
  data_m <- reshape2::melt(data_m, id.vars = "rowname")

  # Cluster rows
  if (row.dend) {
    dd.row <- as.dendrogram(hclust(dist(mat, method = dist.method),
                                   method = clust.method))
    row.ord <- order.dendrogram(dd.row)
    ordered_row_names <- row.names(mat[row.ord, ])
    data_m$rowname <- factor(data_m$rowname, levels = ordered_row_names)
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
  p_heat <- ggplot(data_m, aes(x = variable, y = rowname, fill = value)) +
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
## wl-04-09-2015, Fri: Heatmap with dendrogram on both sides using lattice
## wl-15-06-2021, Tue: Makes heatmap consistent with data matrix and remove
## non-lattice parts
## heat_dend(mtcars)
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
## wl-18-05-2021, Tue: stats of a vector.
##  Used for error bar plotting. Modify from https://bit.ly/3onsqot
vec_ci <- function(x, na.rm = FALSE, conf.interval = .95) {

  ## Handle NA's: if na.rm==T, don't count them
  length2 <- function(x, na.rm = FALSE) {
    if (na.rm) sum(!is.na(x)) else length(x)
  }
  n    <- length2(x, na.rm = na.rm)
  mean <- mean(x, na.rm = na.rm)
  sd   <- sd(x, na.rm = na.rm)
  se   <- sd / sqrt(n)                 ## standard error of the mean

  ## Confidence interval multiplier for standard error
  ## Calculate t-statistic for confidence interval:
  ## e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ci_mult <- qt(conf.interval / 2 + .5, n - 1)
  ci <- se * ci_mult

  c(n = n, mean = mean, sd = sd, se = se, ci = ci)
}

## --------------------------------------------------------------------------
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
## wl-08-02-2017: Node descriptive
## wl-09-02-2017: add more descriptive
## wl-01-03-2017: add the membership of community detection for grouping.
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
  ## node.stats <- arrange.row(node.stats, desc(Degree), desc(Abundance))

  return(node.stats)
}

## ------------------------------------------------------------------------
## wl-23-06-2015: Get correlation coefficient and p-values
## Note:
##   This file is modified from 'cor.table' of package 'picante'
## and 'corr.test' of package 'psych'.
##   The original implementation is from Bill Venables, the author of R great
## book MASS. For details, see
## https://stat.ethz.ch/pipermail/r-help/2001-November/016201.html
## Arguments:
##  x - a data frame or matrix for correlation analysis column-wise
##  cor.method - method for correlation
##  adj.method - p-value correction method
##  ... - other parameter for correlation.
## Values:
##  r - correlation coefficient
##  p - statistics matrix, in which the lower triangular is p-values and the
##      upper triangular is adjusted p-values
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
## wl-24-06-2015: Convert a symmetric table(short format) to long format
## Arguments:
##   x     - A symmetric matrix-like data set
##   tri   - Triangular being used
## Returns:
##   A data frame of pair-wise comparison
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

  if (F) { ## require(reshape2)
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
## wl-09-10-2021, Sat: scale vector with standarization
std_scale <- function(x, na.rm = TRUE) {
  res <- (x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm)
  return(res)
}

## ------------------------------------------------------------------------
## wll-05-05-2016: Rescale vector between min and max
## From 'rescale' package 'network'
mm_scale <- function(nchar, low, high) {
  min_d <- min(nchar)
  max_d <- max(nchar)
  rscl <- ((high - low) * (nchar - min_d)) / (max_d - min_d) + low
  rscl
}

## ------------------------------------------------------------------------
## wll-12-01-2016: calculate the ellipse values
elli <- function(x, y, conf.level = 0.95) {
  Var <- var(cbind(x, y))
  Mean <- cbind(mean(x), mean(y))
  Elli <- ellipse::ellipse(Var, centre = Mean, level = conf.level)
}

## --------------------------------------------------------------------------
## lwc-29-04-2013: trim white spaces
str_trim <- function(string) {
  string <- gsub("^[ \t]+|[ \t]+$", "", string) ## trim white spaces
  return(string)
}

## ------------------------------------------------------------------------
## wll-15-09-2015: Calculate the percentage of non digits
non_digit <- function(mat) {
  mat <- as.data.frame(mat)
  func <- function(x) round(sum(is.na(x) | is.nan(x)) / length(x), digits = 3) * 100
  res <- apply(mat, 2, func)
  return(res)
}

## ------------------------------------------------------------------------
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
## wll-09-09-2015: Calculate the non-zero percentage
## wll-22-10-2015: minor changes
is_zero <- function(mat) {
  mat <- as.data.frame(mat)
  res <- sapply(mat, function(x) {
    round(sum(x < .Machine$double.eps) / length(x), digits = 3) * 100
  })
  return(res)
}

## --------------------------------------------------------------------------
## lwc-25-08-2011: get row index of missing values in data frame.
## Arguments:
##   df   -  data frame being used to check index of missing values.
##   vars -  vector of variable names in df for missing values checking.
df_na_idx <- function(df, vars) {
  tmp <- df[, vars]
  idx <- complete.cases(tmp)
  return(idx)
}

## ------------------------------------------------------------------------
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
## lwc-28-01-2016: Convert matrix to df.
##  Note: from '.matrix_to_df' of package 'plyr'
mat2df <- function(.data) {
  cnames <- colnames(.data)
  if (is.null(cnames)) cnames <- rep("", ncol(.data))
  .data <- as.data.frame(.data, stringsAsFactors = FALSE)
  colnames(.data) <- cnames
  .data
}

## ------------------------------------------------------------------------
## wl-01-06-2011: Transpose a numeric data frame (with/without NAs)
## Note: It is used mostly for numeric data frame with NAs and return a
## numeric data frame
## Usages:
## x <- data.frame(group = c(1,1,2,NA,2), val = c(6,4,6,3.1,NA))
## x <- df.t(x)
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
## wll-25-01-2017: Wrapper function for dim matrix of a list of data frame
dim_mat <- function(mat) {
  res <- as.data.frame(t(sapply(mat, dim)))
  names(res) <- c("row", "col")
  res
}

## ------------------------------------------------------------------------
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
## wll-09-11-2015: convert a list of unequal vectors to a data frame
vec2dat <- function(x) {

  ## to-do: check if x is a list
  ## require(plyr)
  res <- lapply(x, function(y) as.data.frame(t(as.data.frame(y))))
  res <- do.call(rbind.fill, res)
}

## -------------------------------------------------------------------------
## lwc-23-11-2013: save PS and TIFF file
## Arguments:
##  handle -- figure handle from lattice
##  PRE    -- directory and prefix of file name
##  fn     -- file name
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
    file = paste(PRE, fn, ".tiff", sep = "_"), width = 10, height = 10,
    units = "in", res = 300, compression = "lzw"
  )
  plot(handle)
  dev.off()
}

## --------------------------------------------------------------------------
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
                         limits = c(0, nrow(label(hcdata)))) +
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
    geom_text(data         =  label(hcdata),
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
## wl-01-06-2021, Tue: Split a vector into chunks.
##  see https://bit.ly/2SM4m2G
## x <- 1:10
## n <- 3
## chunk_1(x, n)
## chunk_2(x, n)
## chunk_3(x, n)
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
#' @title Extra functions for metabolomics data analysis
#' @name mtExtra
#'
#' @description `mtExtra` has more functions for metabolomics data analysis.
#'   It should use ' with R package \code{mt}.
#'
#' @section Main functions:
#' The main mtExtra provides more functions for metabolomics data analysis.
#' These functions are statistical analysis and plot methods with `ggplot2` and
#' `lattice`. It uses `tidyverse` as well as `reshape2` and `plyr` packages.
#'
#' @section Package context:
#' This package follows the principles of the "tidyverse" as mush as possible.
#' It also uses `melt` in `reshape2` if `tidyr` is complicated in some
#' circumstances.
#'
#' @importFrom ellipse ellipse
#' @importFrom graphics lines text
#' @importFrom stats prcomp quantile sd var
#' @import tidyr dplyr purrr
#' @importFrom reshape2 melt
#' @importFrom grDevices colorRampPalette dev.off tiff
#' @importFrom stats as.dendrogram as.dist as.hclust complete.cases cor
#'   cutree dist hclust line order.dendrogram p.adjust pf pt qt
#' @importFrom utils data head
#' @docType package
#' @aliases mtExtra mtExtra-package
NULL