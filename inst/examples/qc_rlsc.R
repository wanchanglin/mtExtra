#+ common, include=F
## > 15-08-2023, Tue: test code for qc_rlsc. <br/>
## > 17-08-2023, Thu: finalise. <br/>
rm(list = ls(all = TRUE))
source("_common.R")

#' # Quality control–based robust LOESS signal correction (QC-RLSC)
#'

## ---- Setting ----
#' ## Setting

#' Load R packages 
pkgs <- c("mt", "mtExtra", "tidyverse", "readxl")
invisible(lapply(pkgs, library, character.only = TRUE))

## ----- Read data -----
#' ## Read data

PATH <- here::here("extdata", "data_qcrlsc.xlsx")
xls  <- PATH %>%
  excel_sheets() %>% 
  set_names() %>% 
  map(read_excel, path = PATH)

#' Check the data
names(xls)
t(sapply(xls, dim))

meta <- xls$meta
meta
data <- xls$data %>%
  mutate_if(is.character, as.numeric)
data

#' Group information for data screening
(cls.qc <- factor(meta$type))
(cls.bl <- factor(meta$block))

## ----- Missing value filtering and filling -----
#' ## Missing value filtering and filling

#' Check percentage of missing value
tail(sort(mv_perc(data)), 20)

dat <- mv_filter(data, thres = 0.3)$dat %>%
  mv.fill(method = "median", ze_ne = T) %>% as_tibble()
dat

#' Data screening before signal correction
pcaplot(dat, cls.qc, pcs = c(2, 1), ep = 1)
pcaplot(dat, cls.bl, pcs = c(2, 1), ep = 1)
plot(pcalda(dat, cls.bl), dimen = c(1:2), ep = 2)
plot(plslda(dat, cls.bl), dimen = c(1:2), ep = 2)

## ----- Batch correction by QC-RLSC  -----
#' ## Batch correction by QC-RLSC

if (T) {
  res <- qc_rlsc(dat, cls.qc)
} else { # wl-16-08-2023, Wed: do signal correction inside each batch
  res <- lapply(levels(cls.bl), function(x){
    #' x <- levels(cls.bl)[[1]]
    idx <- cls.bl %in% x
    tmp <- qc_rlsc(dat[idx,], cls.qc[idx])
  })
  res <- bind_rows(res)
}

#' Data screening after singal correction
pcaplot(res, cls.bl, pcs = c(2, 1), ep = 1)
plot(pcalda(res, cls.bl), dimen = c(1:2), ep = 2)
plot(plslda(res, cls.bl), dimen = c(1:2), ep = 2)

## ----- OPTION: batch shift -----
#' ## OPTION: batch shift

#' Remove batch mean
res <- batch_shift(res, cls.bl)

#' Data screening after batch shift
pcaplot(res, cls.bl, pcs = c(2, 1), ep = 1)
plot(pcalda(res, cls.bl), dimen = c(1:2), ep = 2)
plot(plslda(res, cls.bl), dimen = c(1:2), ep = 2)
