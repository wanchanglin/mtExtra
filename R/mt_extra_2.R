
## --------------------------------------------------------------------------
#' Trim white spaces
#' 
#' Trim head and hail white spaces of a string.
#' 
#' @param string a character string to be processed.
#' @return a trimmed string.
#' @export 
## lwc-29-04-2013: trim white spaces
trim_str <- function(string) {
  string <- gsub("^[ \t]+|[ \t]+$", "", string) ## trim white spaces
  return(string)
}

## ------------------------------------------------------------------------
#' Calculate the percentage of non digits
#' 
#' Calculate the percentage of non digits.
#' 
#' @param mat a data matrix.
#' @return  an vector.
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
#' @return  an vector.
#' @export  
## wll-09-09-2015: Calculate the non-zero percentage
## wll-22-10-2015: minor changes
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
#' @return  an vector.
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
#' @return  an vector of missing value index.
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
#' @return  a data mstrix.
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
## Note: compare with 'save.tab' in package 'mt': `save.tab' do no require
##  the same dimension but 'rbind_df' does.
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
#' @return  a data frame.
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
#' @return  a transposed data frame.
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
#' Get a matrix for dimension information of a list of data frame.
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
#' Sort a date fram and keep the row names.
#' 
#' @param df a data frame.
#' @param ... further parameters to `order`.
#' @return  a sorted data frame.
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
#' @return  a data matrix.
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

  if (force.number.of.groups) {
    g.names <- names(g)
    g.names.ordered <- as.character(sort(as.numeric(g.names)))
  } else {
    g.names <- names(g[-length(g)])
    g.names.ordered <- as.character(sort(as.numeric(g.names)))
    g.names.ordered <- c(g.names.ordered, "overflow")
  }

  return(g[g.names.ordered])
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
#' Return a `ggdendro` extension with cluster centres.
#' 
#' @param hc an object of `hclust`.
#' @param k the number of clusters.
#' @return  `ggplot2` plot.
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
#' @return  a `ggplot2` plot of dendrogram.
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
#' Update data set by a subset of features
#' 
#' Update data set after feature selection. 
#' 
#' @param dat.all a list of metabolomics data including intensity data, peak 
#'   data and meta data
#' @param ord a subset of selected feature orders
#' @export 
## lwc-03-03-2010: Update data set by a subset of features or variables.
## lwc-10-08-2010: Fix a bug
## lwc-17-01-2011: Only consider dat.all corresponding to feature order
##   provided by ord
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

##  1) str_trim
##  2) non_digit
##  3) non_zero
##  4) is_zero
##  5) df_na_idx
##  6) rbind_df
##  7) mat2df
##  8) df_t
##  9) dim_mat
## 10) arrange_row
## 11) vec2dat
## 12) chunk
## 13) chunk_1
## 14) chunk_2
## 15) chunk_3
## 16) ps_tiff
## 17) dendro_data_k
## 18) plot_ggdendro
## 19) upd_data
