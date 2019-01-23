#' Overlay the traceplots of MCMC output of different models/algorithms
#'
#' @param df a data frame of the MCMC iterations.
#' @param fit_str a vector of character string, possibly defined in get_par_df().
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot geom_line aes theme labs
#' @export
plot_trace <- function(df, fit_str = c("0", "1")) {
    df0 <- dplyr::filter(df, .data$fit %in% fit_str)
    gg0 <- ggplot2::ggplot(df0) +
        ggplot2::geom_line(ggplot2::aes(.data$iteration, .data$value, col = .data$fit)) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(x = "iteration", y = "value")
    gg0
}

#' Overlay the posterior densities of MCMC output of different models/algorithms
#'
#' @param df a data frame of the MCMC iterations.
#' @param fit_str a vector of character string, possibly defined in get_par_df().
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot geom_density aes theme labs
#' @export
plot_dpost <- function(df, fit_str = c("0", "1")) {
    df0 <- dplyr::filter(df, .data$fit %in% fit_str)
    gg0 <- ggplot2::ggplot(df0) +
        ggplot2::geom_density(ggplot2::aes(.data$value, col = .data$fit, lty = .data$fit), lwd = 1.2) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(x = "value")
    gg0
}

#' Plot retweet count
#'
#' @param df a data frame with at least the following variables: id_str and t_ij.
#' @param t_0 start of data collection; POSIXct.
#' @param t_inf end of data collection relative to t_0; a positive scalar.
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom lubridate seconds
#' @importFrom ggplot2 ggplot geom_line aes scale_x_datetime scale_y_log10 theme element_text labs
#' @export
plot_retweet_count <- function(df, t_0, t_inf) {
    df0 <- dplyr::mutate(dplyr::ungroup(dplyr::mutate(dplyr::group_by(df, .data$id_str), count = rank(.data$t_ij))), t = t_0 + lubridate::seconds(.data$t_ij))
    gg0 <- ggplot2::ggplot(df0) +
        ggplot2::geom_line(ggplot2::aes(.data$t, .data$count, col = .data$id_str)) +
        ggplot2::scale_x_datetime(limits = c(t_0, t_0 + seconds(t_inf)), expand = c(0.1, 0)) +
        ggplot2::scale_y_log10(breaks = 10L^(-1L:4L), minor_breaks = NULL) +
        ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5)) +
        ggplot2::labs(x = NULL, y = "retweet count")
    gg0
}
