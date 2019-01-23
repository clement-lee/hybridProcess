#' Get and augment the parameter data frame of an MCMC output
#'
#' @param df a data frame of the MCMC iterations of the parameters.
#' @param fit_chr a scalar character string for identifying the model/algorithm when combining with other outputs.
#' @importFrom dplyr mutate
#' @importFrom tibble as_tibble
#' @export
get_par_df <- function(df, fit_chr) {
    tibble::as_tibble(dplyr::mutate(df, iteration = seq(nrow(df)), fit = fit_chr))
}

#' Computation time (in hours) from MCMC output
#'
#' @param l an output from run_mcmc().
#' @export
computation_time <- function(l) {
    round(attributes(l)$duration_mcmc[["elapsed"]]/3600.0,1L)
}

#' Wrapper of the MCMC algorithms coded in C++
#'
#' @param date_str a character string representing the date.
#' @param seed a single value, interpreted as an integer, to be used in set.seed().
#' @param f_mcmc name of an MCMC algorithm coded in C++. Current choices are "mh_etas", "rj_etas" and "me_etas" (remove the quotes when using).
#' @param N desired chain length AFTER burn-in and thinning.
#' @param thin thinning in MCMC.
#' @param burnin burn-in in MCMC.
#' @param print_freq how frequent should the current values in the MCMC be printed.
#' @param a_theta_prop shape parameter of pseudoprior for theta.
#' @param b_theta_prop rate parameter of pseudoprior for theta.
#' @param p0 prior for model 0 (power-law process), as opposed to model 1 (hybrid process).
#' @param p01 jump probability from model 0 to 1 in RJMCMC.
#' @param p10 jump probability from model 1 to 0 in RJMCMC.
#' @param write boolean; should csv files of estimated and predicted retweet counts be written?
#' @param l0 a list containing two data frames df1.original and df0.range.
#' @param beta initial value of beta, coefficient for the linear term.
#' @param kappa initial value of kappa, coefficient for the quadratic term.
#' @param lambda initial value of lambda, parameter for the power-law effect in the hybrid process.
#' @param phi initial value of phi, scale parameter for the hybrid process.
#' @param psi initial value of psi, shift in the power-law term.
#' @param tau initial value of tau, precision of the Gaussian distribution for the random error.
#' @param theta initial value of theta, parameter for the exponential decay in the hybrid process.
#' @param s_beta initial proposal standard deviation for the Metropolis step for beta.
#' @param s_kappa initial proposal standard deviation for the Metropolis step for kappa.
#' @param s_lambda initial proposal standard deviation for the Metropolis step for lambda.
#' @param s_psi initial proposal standard deviation for the Metropolis step for psi.
#' @param s_theta initial proposal standard deviation for the Metropolis step for theta.
#' @param s_e_init initial proposal standard deviation for the Metropolis step for each random error.
#' @param rebuild boolean; should the MCMC be re-run if an output (an rds files) already exists?
#' @importFrom glue glue
#' @useDynLib hybridProcess
#' @importFrom Rcpp sourceCpp
#' @importFrom readr write_rds
#' @export
run_mcmc <- function(date_str, seed, f_mcmc, N, thin, burnin, print_freq, a_theta_prop, b_theta_prop, p0, p01, p10, write, l0, beta, kappa, lambda, phi, psi, tau, theta, s_beta, s_kappa, s_lambda, s_psi, s_theta, s_e_init, rebuild = FALSE) {
    ## wrapper for running an MCMC algo. to data
    if (identical(f_mcmc, mh_etas)) {
        if (p0 != 1.0 && p0 != 0.0) {
            stop("run_mcmc: p0 has to be 0.0 or 1.0 for mh_etas.")
        }
        fit_str <- ifelse(p0 == 0.0, "mcmc_1", "mcmc_0")
    }
    else if (identical(f_mcmc, rj_etas)) {
        fit_str <- "rjmcmc"
    }
    else if (identical(f_mcmc, ms_etas)) {
        fit_str <- "gvs_cc"
    }
    else if (identical(f_mcmc, me_etas)) {
        fit_str <- "bridge"
    }
    else {
        stop("run_mcmc: f_mcmc has to be one of mh_etas, rj_etas, ms_etas or me_etas.")
    }
    prefix <- glue::glue("results/{date_str}_{fit_str}")
    name_rds <- glue::glue("{prefix}_single.rds")
    if (!file.exists(name_rds) || rebuild) {
        set.seed(seed)
        filename_par <- glue::glue("{prefix}_par_single.csv")
        filename_m <- glue::glue("{prefix}_m_single.csv")
        filename_p <- glue::glue("{prefix}_p_single.csv")
        t0.fit <- system.time({
            l0.fit <- f_mcmc(beta, kappa, lambda, phi, psi, tau, theta, l0$df1.original$log_followers_count, l0$df1.original$t_i, l0$df1.original$retweet_count, l0$df0.original.b$t_relative, l0$df0.range$t_inf, s_beta, s_kappa, s_lambda, s_psi, s_theta, s_e_init, N, thin, burnin, print_freq, a_theta_prop, b_theta_prop, p0, p01, p10, write, filename_par, filename_m, filename_p)
        })
        attr(l0.fit, "datetime_data") <- l0$df0.range$t_0
        attr(l0.fit, "duration_mcmc") <- t0.fit
        readr::write_rds(l0.fit, name_rds)
    }
    invisible(0L)
}

