## wl-09-11-2021, Tue: gather all general functions from 2015

## ------------------------------------------------------------------------
#' Caculate correlation coefficient and p-values
#' 
#' Calculate correlation coefficient and p-values.
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
#' @examples 
#' cor_hcl(mtcars)
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
#' @examples 
#' res <- cor_net(mtcars)
#' names(res)
#' res
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
  ## g_dat <- reshape2::melt(co_mat)
  ## g_dat <- reshape2:::melt.matrix(co_mat)
  g_dat <- dat2long(co_mat)
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
#' - Keep the variables which have at least one significant element. The 
#'   significant element is defined as larger than the lower of threshold or
#'   less than the upper of theshold. 
#' - Threshold can be given by user or be estimated by `locfdr`, i.e. the 
#'   returned `z.2` as a threshold. It is not gurranted that `locfdr` returns
#'   `z.2`. If not, user must provide this value.
#' - From R package `locfdr` vignette: 
#'   `z.2` is the interval along the zz-axis outside of which `fdr(z) < 0.2`,
#'   the locations of the yellow triangles in the histogram plot. If no
#'   elements of `zz` on the left or right satisfy the criterion, the
#'   corresponding element of z.2 is NA, and the corresponding triangle does
#'   not appear.
#'  
#' @seealso [locfdr()] 
#' @family variable filters
#' @importFrom locfdr locfdr
#' @examples 
#' library(dplyr)
#' library(tidyr)
#' library(purrr)
#' library(readr)
#' 
#' ## get ionomics data
#' dat <- read_csv("https://github.com/wanchanglin/ionflow/raw/master/extra/paper_ko.csv")
#' dim(dat)
#' 
#' ## missing valuee filling with mean
#' dat <- dat %>% 
#'   mutate(across(where(is.numeric), function(x) {
#'     m <- mean(x, na.rm = TRUE)
#'     x[is.na(x)] <- m
#'     x
#'   }))
#' dat
#' 
#' res <- locfdr_filter(t(dat[, -1]), plot = 1)
#' res$thres
#' 
#' ## filter data
#' dat <- dat[res$idx, , drop = FALSE]
#' 
#' ## symbolise data
#' dat_sym <- dat %>% 
#'   mutate(across(where(is.numeric), ~ dat_symb(., thres = res$thres)))
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
#' @details 
#'   `thres` can be estimated by [locfdr_filter()]. 
#' @seealso [locfdr_filter()] and [locfdr()] for `z.2`. 
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
#' Caculate the count number of features
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
    dplyr::mutate(count = rowSums(across(where(is.numeric)))) %>%
    ## mutate(count = sum(c_across(where(is.numeric)))) %>%
    ungroup() %>%
    arrange(desc(count))

  return(fs_top)
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
#' library(ggplot2)
#' library(dplyr)
#' library(tidyr)
#' library(purrr)
#' 
#' data(iris)
#' rsd(iris[, 1:4])
#' 
#' ## group rsd
#' val <- iris %>% 
#'   group_by(Species) %>% 
#' 	 group_modify(~ dat_summ(., method = rsd)) %>%
#'   pivot_longer(cols = !Species) %>% filter(!is.na(value))
#' 
#' ggplot(val, aes(x = value)) +
#'   geom_histogram(colour = "black", fill = "white") +
#'   facet_grid(Species ~ .)
#' ggplot(val, aes(x = Species, y = value)) + geom_boxplot()
#' 
#' ## The plotting of missing value is similiar like this.
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

##  1) cor_tab
##  2) cor_hcl
##  3) cor_net
##  4) bi_cor_net
##  6) pcor_dat
##  7) qc_filter
##  8) blank_filter
##  9) locfdr_filter
## 10) dat_symb
## 11) var_filter
## 12) mv_filter
## 13) rsd_filter
## 14) outl_det_u
## 15) outl_det_m
## 16) batch_shift
## 17) feat_count
## 18) rsd
## 19) mv_perc
## 20) net_graph
## 21) graph_stats
## 22) vertex_stats
## 23) samp_strat
## 24) samp_sub
