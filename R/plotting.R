# Functions to plot geoBAMr objects


#' Plot a geoBAMr-created object
#'
#' @param bamobj Currently supported objects are bamdata objects (created by
#'   \code{bam_data()}) and bamval objects (created by \code{bam_validate()}).
#'
#' @export
bam_plot <- function(bamobj) {
  UseMethod("bam_plot")
}

#' Plot a bamdata object
#'
#' @importFrom reshape2 melt
#' @importFrom stats setNames
#' @import ggplot2
#' @importFrom rlang .data
#' @param bamobj A bamdata object, as produced by \code{bam_data()}
#'
#' @export
bam_plot.bamdata <- function(bamobj) {
  piece <- c("w", "s", "dA")

  if (sum(bamobj$Sobs != 1, na.rm = TRUE) == 0) {
    piece <- "w"
  }
  nx <- bamobj$nx
  w_df <- as.data.frame(t(bamobj$Wobs)) %>%
    setNames(1:nx)
  s_df <- as.data.frame(t(bamobj$Sobs)) %>%
    setNames(1:nx)
  dA_df <- as.data.frame(t(bamobj$dAobs)) %>%
    setNames(1:nx)
  # browser()
  w_df$time <- s_df$time <- dA_df$time <- 1:bamobj$nt
  # browser()
  sw <- suppressWarnings
  data_long <- sw(dplyr::bind_rows(w = melt(w_df, id.vars = "time",
                                            variable.name = "xs"),
                                   s = melt(s_df, id.vars = "time",
                                           variable.name = "xs"),
                                  dA = melt(dA_df, id.vars = "time",
                                            variable.name = "xs"),
                                 .id = "variable"))
  data_long$xs <- as.numeric(as.character(data_long$xs))
  plotdata <- data_long[data_long[["variable"]] %in% piece, ]

  out <- ggplot(plotdata, aes(x = .data$time, y = .data$value)) +
    geom_line(aes(color = .data$xs, group = .data$xs)) +
    scale_color_gradient() +
    facet_wrap(~variable, scales = "free_y", ncol = 1) +
    xlab("time") + ylab("value") +
    guides(color = guide_colorbar("location"))

  out
}

#' Plot a bamval object to show predictive performance
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @param bamobj a \code{bamval} object, as produced by \code{bam_validate}
#' @export
bam_plot.bamval <- function(bamobj) {
  valdata <- bamobj$valdata
  ggplot(valdata, aes(x = .data$qobs, y = .data$qpred)) +
    geom_point() +
    geom_abline(aes(intercept = 0, slope = 1))
}

#' Plot flow time series from BAM inference
#'
#' @param fit A stanfit object, as returned from \code{bam_estimate()}
#' @param qobs An optional vector giving observed flow for comparison
#' @importFrom dplyr "%>%"
#' @importFrom rlang .data
#' @export
bam_hydrograph <- function(fit, qobs = NULL) {

  nchains <- length(fit@stan_args)
  qpred <- lapply(1:nchains, function(x) bam_qpred(fit, x)) %>%
    setNames(paste0("chain", 1:nchains)) %>%
    dplyr::bind_rows(.id = "series") %>%
    reshape2::melt(id.vars = c("series", "time"),
                   measure.vars = c("mean", "conf.low", "conf.high"),
                   variable.name = "stat", value.name = "flow") %>%
    dplyr::mutate(stat = as.character(stat))

  out <- ggplot(qpred, aes(x = .data$time, y = .data$flow, color = .data$stat)) +
    geom_line(aes(linetype = .data$series))

  if (!is.null(qobs)) {
    obsdf <- data.frame(time = 1:max(qpred$time),
                        flow = qobs, series = "observed", stat = NA)
    out <- out +
      geom_line(aes(x = .data$time, y = .data$flow, linetype = .data$series),
                data = obsdf) +
      scale_linetype_manual(values = c(2:(nchains + 1), 1))
  }


  out
}

