#' Simulate retweet times according to hybrid process
#'
#' @param T duration of simulation period; a positive scalar.
#' @param x ``mean-centred'' follower count; a real-valued scalar.
#' @param beta coefficient for the linear term x; a real-valued scalar.
#' @param kappa coefficient for the quadratic term x squared; a real-valued scalar.
#' @param lambda parameter for the power-law effect in the hybrid process; a scalar smaller than 1.
#' @param phi scale parameter for the hybrid process; a positive scalar.
#' @param psi shift in the power-law term; a non-negative scalar.
#' @param tau precision of the Gaussian distribution for the random error; a positive scalar.
#' @param theta parameter for the exponential decay in the hybrid process; a non-negative scalar.
#' @importFrom stats pgamma qgamma rnorm rpois runif setNames time
sim_hybrid <- function(T, x, beta, kappa, lambda, phi, psi, tau, theta) {
    omega <- 1.0 - lambda
    g1ml <- gamma(omega)
    gamma <- exp(beta * x + kappa * x * x + rnorm(1L, 0.0, tau ^ (-0.5))) # not the gamma function
    if (theta < 0.0) {
        stop("sim_hybrid: theta has to be non-negative.")
    }
    else if (theta == 0.0) {
        H <- phi * gamma * ((T + psi) ^ omega - psi ^ omega) / omega
    }
    else { # theta > 0.0
        c <- phi * gamma * g1ml / theta ^ omega * exp(theta * psi)
        H <- c * (pgamma(theta * (T + psi), omega) - pgamma(theta * psi, omega))
    }
    n <- rpois(1L, H)
    if (n == 0L) { # 0 RT
        t <- as.numeric(NA)
    }
    else { # >= 1 RT's
        u <- runif(n, 0.0, H)
        if (theta == 0.0) {
            t <- ((u / gamma * omega / phi) + (psi ^ omega)) ^ (1.0 / omega) - psi
        }
        else { # theta > 0.0
            t <- qgamma(u / c + pgamma(theta * psi, omega), omega) / theta - psi
        }
    }
    sort(t, na.last = TRUE)
}

#' Simulate hybrid process for multiple original tweets with respective follower counts and times of creation
#'
#' @param df.data a data frame with the following variables: id_str (character, unique identifier of original tweet), log_followers_count (logarithm of 1 + user follower count), t_i (time of creation of original tweet, relative to start of data collection).
#' @param df.summary a one-row data frame with the true values of the following parameters defined in sim_hybrid(): beta, kappa, lambda, phi, psi, tau, theta.
#' @param t_inf total duration of data collection / simulation period.
#' @importFrom dplyr mutate transmute
#' @importFrom purrr map2 map_lgl map_int
#' @importFrom magrittr %>% add
sim_and_mutate <- function(df.data, df.summary, t_inf) {
    df.data %>%
        dplyr::mutate(x_star = .data$log_followers_count - mean(.data$log_followers_count),
               t = map2(t_inf - .data$t_i, .data$x_star, sim_hybrid, beta = df.summary$beta, kappa = df.summary$kappa, lambda = df.summary$lambda, phi = df.summary$phi, psi = df.summary$psi, tau = df.summary$tau, theta = df.summary$theta),
               retweeted = !map_lgl(.data$t, identical, as.numeric(NA)),
               length = map_int(.data$t, length)) %>%
        transmute(id_str = .data$id_str,
###                  user_followers_count, # remove if not needed
                  retweet_count = ifelse(.data$retweeted, .data$length, 0L),
                  t_i = .data$t_i,
                  t_ij = map2(.data$t, .data$t_i, magrittr::add),
                  t_relative = .data$t,
                  retweeted = .data$retweeted,
                  log_followers_count = .data$log_followers_count,
                  log_retweet_count = log(.data$retweet_count + 1L))
}

#' Wrapper of sim_and_mutate(), by creating a list that aligns with l.original
#'
#' @param l.original a list containing two data frames df1.original and df0.range. The former can be used in the argument df.data in sim_and_mutate(), while the latter is a one-row data frame with (at least) the variable t_inf.
#' @param df.summary a one-row data frame with the true values of the following parameters defined in sim_hybrid(): beta, kappa, lambda, phi, psi, tau, theta.
#' @param seed a single value, interpreted as an integer, to be used in set.seed().
#' @importFrom tidyr unnest
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @export
sim_wrapper <- function(l.original, df.summary, seed) {
    l0 <- l.original # clone & modify
    set.seed(seed)
    df1.original <- l0$df1.original %>%
        sim_and_mutate(df.summary, t_inf = l0$df0.range$t_inf)
    df0.original.a <- df1.original %>%
        tidyr::unnest() %>%
        dplyr::filter(is.na(.data$t_ij))
    df0.original.b <- df1.original %>%
        tidyr::unnest() %>%
        dplyr::filter(!is.na(.data$t_ij))
    invisible(
        list(
            df0.range = l0$df0.range,
            df0.original.a = df0.original.a,
            df0.original.b = df0.original.b,
            df1.original = df1.original
        )
    )
}
