// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace Rcpp;
using namespace std;
using namespace arma;

// 00) Prelim
void nan_to_minus_infinity(double & x) {
  // turn NaN x to -Inf
  if (isnan(x)) {
    x = -INFINITY;
  }
}

void update(double & par_curr, const double par_prop, double & lpost_curr, const double lpost_prop, double & s, const int i, const int burnin, const double factor = 10.0) {
  // M-H update
  const bool accept_reject = log(runif(1)[0]) < lpost_prop - lpost_curr;
  par_curr = accept_reject ? par_prop : par_curr;
  lpost_curr = accept_reject ? lpost_prop : lpost_curr;
  if (i < burnin) {
    s = sqrt(s * s + (accept_reject ? 3.0 : (-1.0)) * s * s / factor / sqrt(i + 1.0));
  }
}

void update_vec(NumericVector & par_curr, const NumericVector par_prop, NumericVector & lpost_curr, const NumericVector lpost_prop, NumericVector & s, const int i, const int burnin, const double factor = 10.0) {
  // vectorised version of update - template in the future?
  const int n = par_curr.size();
  const LogicalVector accept_reject = log(runif(n)) < lpost_prop - lpost_curr;
  par_curr = ifelse(accept_reject, par_prop, par_curr);
  lpost_curr = ifelse(accept_reject, lpost_prop, lpost_curr);
  if (i < burnin) {
    s = sqrt(s * s + ifelse(accept_reject, 3.0, -1.0) * s * s / factor / sqrt(i + 1.0));
  }
}

void stat_pars(const int n,
               running_stat<double> & M_stat, const int M_curr,
               running_stat<double> & beta_stat, const double beta_curr,
               running_stat<double> & kappa_stat, const double kappa_curr,
               running_stat<double> & lambda_stat, const double lambda_curr,
               running_stat<double> & phi_stat, const double phi_curr,
               running_stat<double> & psi_stat, const double psi_curr,
               running_stat<double> & tau_stat, const double tau_curr,
               running_stat<double> & theta_stat, const double theta_curr,
               running_stat<double> & e_first_stat, running_stat<double> & e_last_stat, const NumericVector e_curr) {
  // record the pars in the running_stat objects
  M_stat((double) M_curr);
  beta_stat(beta_curr);
  kappa_stat(kappa_curr);
  lambda_stat(lambda_curr);
  phi_stat(phi_curr);
  psi_stat(psi_curr);
  tau_stat(tau_curr);
  theta_stat(theta_curr);
  e_first_stat(e_curr[0]);
  e_last_stat(e_curr[n-1]);
}

void print_all(const int i, const int burnin, const int n,
               const int M_curr, const running_stat<double> M_stat,
               const double beta_curr, const running_stat<double> beta_stat, const double s_beta,
               const double kappa_curr, const running_stat<double> kappa_stat, const double s_kappa,
               const double lambda_curr, const running_stat<double> lambda_stat, const double s_lambda,
               const double phi_curr, const running_stat<double> phi_stat,
               const double psi_curr, const running_stat<double> psi_stat, const double s_psi,
               const double tau_curr, const running_stat<double> tau_stat,
               const double theta_curr, const running_stat<double> theta_stat, const double s_theta,
               const NumericVector e_curr, const running_stat<double> e_first_stat, const running_stat<double> e_last_stat) {
  // print chain values at pre-determined iteration
  Rcout << "Model (curr): " << M_curr << endl;
  Rcout << "Model (mean): " << M_stat.mean() << endl;
  Rcout << "beta (curr): " << beta_curr << endl;
  Rcout << "beta (mean): " << beta_stat.mean() << endl;
  if (i < burnin) {
    Rcout << "s_beta: " << s_beta << endl;
  }
  Rcout << "kappa (curr): " << kappa_curr << endl;
  Rcout << "kappa (mean): " << kappa_stat.mean() << endl;
  if (i < burnin) {
    Rcout << "s_kappa: " << s_kappa << endl;
  }
  Rcout << "lambda (curr): " << lambda_curr << endl;
  Rcout << "lambda (mean): " << lambda_stat.mean() << endl;
  if (i < burnin) {
    Rcout << "s_lambda: " << s_lambda << endl;
  }
  Rcout << "phi (curr): " << phi_curr << endl;
  Rcout << "phi (mean): " << phi_stat.mean() << endl;
  Rcout << "psi (curr): " << psi_curr << endl;
  Rcout << "psi (mean): " << psi_stat.mean() << endl;
  if (i < burnin) {
    Rcout << "s_psi: " << s_psi << endl;
  }
  Rcout << "tau (curr): " << tau_curr << endl;
  Rcout << "tau (mean): " << tau_stat.mean() << endl;
  Rcout << "theta (curr): " << theta_curr << endl;
  Rcout << "theta (mean): " << theta_stat.mean() << endl;
  if (i < burnin && M_curr == 1) {
    Rcout << "s_theta: " << s_theta << endl;
  }
  Rcout << "e_first (curr): " << e_curr[0] << endl;
  Rcout << "e_first (mean): " << e_first_stat.mean() << endl;
  Rcout << "e_last (curr): " << e_curr[n-1] << endl;
  Rcout << "e_last (mean): " << e_last_stat.mean() << endl;
  Rcout << endl;
}

DataFrame df_par(const IntegerVector M_vec,
                 const NumericVector beta_vec,
                 const NumericVector kappa_vec,
                 const NumericVector lambda_vec,
                 const NumericVector phi_vec,
                 const NumericVector psi_vec,
                 const NumericVector tau_vec,
                 const NumericVector theta_vec,
                 const NumericVector e_first_vec,
                 const NumericVector e_last_vec,
                 const NumericVector chisq_actual_vec,
                 const NumericVector chisq_estimated_vec,
                 const NumericVector chisq_predicted_vec) {
  // data frame of parameters
  return DataFrame::create(Named("M") = M_vec,
                           Named("beta") = beta_vec,
                           Named("kappa") = kappa_vec,
                           Named("lambda") = lambda_vec,
                           Named("phi") = phi_vec,
                           Named("psi") = psi_vec,
                           Named("tau") = tau_vec,
                           Named("theta") = theta_vec,
                           Named("e_first") = e_first_vec,
                           Named("e_last") = e_last_vec,
                           Named("chisq_actual") = chisq_actual_vec,
                           Named("chisq_estimated") = chisq_estimated_vec,
                           Named("chisq_predicted") = chisq_predicted_vec);
}

DataFrame df_fitted(const NumericVector x,
                    const NumericVector t_original,
                    const IntegerVector m_retweets,
                    const running_stat_vec<vec> m_estimated_stat,
                    const running_stat_vec<vec> m_predicted_stat,
                    const running_stat_vec<vec> e_stat,
                    const NumericVector s_e) {
  // data frame of original data & fitted values
  NumericVector
    m_estimated_mean   = wrap(m_estimated_stat.mean()),
    m_estimated_stddev = wrap(m_estimated_stat.stddev()),
    m_estimated_min    = wrap(m_estimated_stat.min()),
    m_estimated_max    = wrap(m_estimated_stat.max()),
    m_predicted_mean   = wrap(m_predicted_stat.mean()),
    m_predicted_stddev = wrap(m_predicted_stat.stddev()),
    m_predicted_min    = wrap(m_predicted_stat.min()),
    m_predicted_max    = wrap(m_predicted_stat.max()),
    e_mean   = wrap(e_stat.mean()),
    e_stddev = wrap(e_stat.stddev());
  return DataFrame::create(Named("log_followers_count") = x,
                           Named("t_i") = t_original,
                           Named("retweet_count") = m_retweets,
                           // above: to align with pre-mcmc labels
                           // below: post-mcmc, so adhere to labels here
                           Named("m_estimated_mean")   = m_estimated_mean,
                           Named("m_estimated_stddev") = m_estimated_stddev,
                           Named("m_estimated_min")    = m_estimated_min,
                           Named("m_estimated_max")    = m_estimated_max,
                           Named("m_predicted_mean")   = m_predicted_mean,
                           Named("m_predicted_stddev") = m_predicted_stddev,
                           Named("m_predicted_min")    = m_predicted_min,
                           Named("m_predicted_max")    = m_predicted_max,
                           Named("e_mean") = e_mean,
                           Named("e_stddev") = e_stddev,
                           Named("s_e") = s_e);
}  

DataFrame df_init_sds(const double beta,
                      const double kappa,
                      const double lambda,
                      const double phi,
                      const double psi,
                      const double tau,
                      const double theta,
                      const double s_beta,
                      const double s_kappa,
                      const double s_lambda,
                      const double s_psi,
                      const double s_theta,
                      const double s_e_init) {
  // put all initial values and proposal standard deviations into a data frame
  return DataFrame::create(Named("beta") = NumericVector::create(beta),
                           Named("kappa") = NumericVector::create(kappa),
                           Named("lambda") = NumericVector::create(lambda),
                           Named("phi") = NumericVector::create(phi),
                           Named("psi") = NumericVector::create(psi),
                           Named("tau") = NumericVector::create(tau),
                           Named("theta") = NumericVector::create(theta),
                           Named("s_beta") = NumericVector::create(s_beta),
                           Named("s_kappa") = NumericVector::create(s_kappa),
                           Named("s_lambda") = NumericVector::create(s_lambda),
                           Named("s_psi") = NumericVector::create(s_psi),
                           Named("s_theta") = NumericVector::create(s_theta),
                           Named("s_e_init") = NumericVector::create(s_e_init));
}

DataFrame df_scalars(const int N,
                     const int thin,
                     const int burnin,
                     const int print_freq,
                     const double a_theta_prop,
                     const double b_theta_prop,
                     const double p0,
                     const double p01,
                     const double p10,
                     const bool write,
                     const double t_infinity,
                     const running_stat_vec<vec> m_estimated_stat,
                     const running_stat_vec<vec> m_predicted_stat,
                     const running_stat_vec<vec> e_stat) {
  // put all MCMC running & misc scalars into a data frame
  return DataFrame::create(Named("N") = IntegerVector::create(N),
                           Named("thin") = IntegerVector::create(thin),
                           Named("burnin") = IntegerVector::create(burnin),
                           Named("print_freq") = IntegerVector::create(print_freq),
                           Named("a_theta_prop") = NumericVector::create(a_theta_prop),
                           Named("b_theta_prop") = NumericVector::create(b_theta_prop),
                           Named("p0") = NumericVector::create(p0),
                           Named("p01") = NumericVector::create(p01),
                           Named("p10") = NumericVector::create(p10),
                           Named("write") = LogicalVector::create(write),
                           Named("t_infinity") = NumericVector::create(t_infinity),
                           Named("m_estimated_count") = m_estimated_stat.count(),
                           Named("m_predicted_count") = m_predicted_stat.count(),
                           Named("e_count") = e_stat.count());
}

DataFrame df_hyper(const double mu_beta,
                   const double tau_beta,
                   const double mu_kappa,
                   const double tau_kappa,
                   const double a_lambda,
                   const double b_lambda,
                   const double a_phi,
                   const double b_phi,
                   const double a_psi,
                   const double b_psi,
                   const double a_tau,
                   const double b_tau,
                   const double a_theta,
                   const double b_theta) {
  // put all hyperparameters into a data frame
  return DataFrame::create(Named("mu_beta") = NumericVector::create(mu_beta),
                           Named("tau_beta") = NumericVector::create(tau_beta),
                           Named("mu_kappa") = NumericVector::create(mu_kappa),
                           Named("tau_kappa") = NumericVector::create(tau_kappa),
                           Named("a_lambda") = NumericVector::create(a_lambda),
                           Named("b_lambda") = NumericVector::create(b_lambda),
                           Named("a_phi") = NumericVector::create(a_phi),
                           Named("b_phi") = NumericVector::create(b_phi),
                           Named("a_psi") = NumericVector::create(a_psi),
                           Named("b_psi") = NumericVector::create(b_psi),
                           Named("a_tau") = NumericVector::create(a_tau),
                           Named("b_tau") = NumericVector::create(b_tau),
                           Named("a_theta") = NumericVector::create(a_theta),
                           Named("b_theta") = NumericVector::create(b_theta));
}

List list_all(const IntegerVector M_vec,
              const NumericVector beta_vec,
              const NumericVector kappa_vec,
              const NumericVector lambda_vec,
              const NumericVector phi_vec,
              const NumericVector psi_vec,
              const NumericVector tau_vec,
              const NumericVector theta_vec,
              const NumericVector e_first_vec,
              const NumericVector e_last_vec,
              const NumericVector chisq_actual_vec,
              const NumericVector chisq_estimated_vec,
              const NumericVector chisq_predicted_vec,
              const NumericVector x,
              const NumericVector t_original,
              const IntegerVector m_retweets,
              const running_stat_vec<vec> m_estimated_stat,
              const running_stat_vec<vec> m_predicted_stat,
              const running_stat_vec<vec> e_stat,
              const NumericVector s_e,
              const double beta,
              const double kappa,
              const double lambda,
              const double phi,
              const double psi,
              const double tau,
              const double theta,
              const double s_beta,
              const double s_kappa,
              const double s_lambda,
              const double s_psi,
              const double s_theta,
              const double s_e_init,
              const int N,
              const int thin,
              const int burnin,
              const int print_freq,
              const double a_theta_prop,
              const double b_theta_prop,
              const double p0,
              const double p01,
              const double p10,
              const bool write,
              const double t_infinity,
              const double mu_beta,
              const double tau_beta,
              const double mu_kappa,
              const double tau_kappa,
              const double a_lambda,
              const double b_lambda,
              const double a_phi,
              const double b_phi,
              const double a_psi,
              const double b_psi,
              const double a_tau,
              const double b_tau,
              const double a_theta,
              const double b_theta) {
  // put all parameters & required output of an MCMC algo in a list of dfs
  DataFrame par = df_par(M_vec, beta_vec, kappa_vec, lambda_vec, phi_vec, psi_vec, tau_vec, theta_vec, e_first_vec, e_last_vec, chisq_actual_vec, chisq_estimated_vec, chisq_predicted_vec),
    fitted = df_fitted(x, t_original, m_retweets, m_estimated_stat, m_predicted_stat, e_stat, s_e),
    init_sds = df_init_sds(beta, kappa, lambda, phi, psi, tau, theta, s_beta, s_kappa, s_lambda, s_psi, s_theta, s_e_init),
    scalars = df_scalars(N, thin, burnin, print_freq, a_theta_prop, b_theta_prop, p0, p01, p10, write, t_infinity, m_estimated_stat, m_predicted_stat, e_stat),
    hyper = df_hyper(mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_tau, b_tau, a_theta, b_theta);
  List output = List::create(Named("par") = par, Named("fitted") = fitted, Named("init_sds") = init_sds, Named("scalars") = scalars, Named("hyper") = hyper);
  return output;
}  



// 01) Modelling originals or retweets of 1 single original w/o latent gammas
//' Evaluate the log-likelihood under the hybrid process
//'
//' @param par vector of length 3, containing (in that order) the power parameter lambda, the exponent parameter theta, and the scale parameter.
//' @param x vector of cumulative/absolute event times.
//' @export
// [[Rcpp::export]]
const double llik_nhpp_hybrid(const NumericVector par, const NumericVector x) {
  // log-likelihood of NHPP w/ hybrid intensity
  // x is cum./abs. event times (!= interarrivals)
  if (par.size() != 3) {
    stop("llik_nhpp_hybrid: length of par has to be 3.");
  }
  const double lambda = par[0], theta = par[1], gamma = par[2];
  double llik = 0.0;
  if (lambda > 1.0 || theta < 0.0 || gamma <= 0.0) {
    llik = -INFINITY;
  }
  else {
    const int n = x.size();
    const double T = max(x), omega = 1.0 - lambda;
    llik = (double) n * log(gamma) - theta * sum(x) - lambda * sum(log(x));
    if (theta == 0.0) {
      llik -= gamma / omega * pow(T, omega);
    }
    else { // theta > 0.0
      llik -= gamma * pgamma(NumericVector::create(theta * T), omega, 1.0)[0] * exp(lgamma(omega) - omega * log(theta));
    }
  }
  nan_to_minus_infinity(llik);
  return llik;
}

//' Evaluate the log-likelihood under a special case of hybrid process when the power parameter is 0
//'
//' @param par vector of length 2, containing (in that order) the exponent parameter theta, and the scale parameter.
//' @param x vector of cumulative/absolute event times.
//' @export
// [[Rcpp::export]]
const double llik_nhpp_exp(const NumericVector par, const NumericVector x) {
  // log-likelihood of NHPP w/ exp. decay intensity i.e. llik_nhpp_hybrid() w/ lambda = 0.0
  if (par.size() != 2) {
    stop("llik_nhpp_exp: length of par has to be 2.");
  }
  const NumericVector par0 = NumericVector::create(0.0, par[0], par[1]);
  return llik_nhpp_hybrid(par0, x);
}

//' Evaluate the log-likelihood under the power law process, equivalently a special case of hybrid process when the exponent parameter is 0
//'
//' @param par vector of length 2, containing (in that order) the power parameter lambda, and the scale parameter.
//' @param x vector of cumulative/absolute event times.
//' @export
// [[Rcpp::export]]
const double llik_nhpp_power(const NumericVector par, const NumericVector x) {
  // log-likelihood of NHPP w/ power intensity i.e. llik_nhpp_hybrid() w/ theta = 0.0
  if (par.size() != 2) {
    stop("llik_nhpp_power: length of par has to be 2.");
  }
  const NumericVector par0 = NumericVector::create(par[0], 0.0, par[1]);
  return llik_nhpp_hybrid(par0, x);
}



// 02) Modelling retweets with latent gammas
const NumericVector lcum_intensity_etas(const double beta, const double kappa, const double lambda, const double phi, const double psi, const double theta, const NumericVector x_star, const NumericVector e, const NumericVector t_diff) {
  // log cumulative intensity
  const int n = t_diff.size();
  if (x_star.size() != n || e.size() != n) {
    stop("lcum_intensity_etas: x_star, e & t_diff have to be of same length.");
  }
  const double omega = 1.0 - lambda;
  NumericVector lvec(n);
  if (theta == 0.0) {
    lvec = log(pow(t_diff + psi, omega) - pow(psi, omega)) - log(omega) + log(phi); // vector + scalar + scalar
  }
  else {
    lvec = log(pgamma(theta * (t_diff + psi), omega, 1.0) - pgamma(NumericVector::create(theta * psi), omega, 1.0)[0]) + lgamma(omega) - omega * log(theta) + log(phi) + theta * psi; // vector + scalar + scalar + scalar + scalar
  }
  lvec += beta * x_star + kappa * x_star * x_star + e;
  return ifelse(lvec != lvec, -INFINITY, lvec);
}

void save_all(const int j, const int n,
              IntegerVector & M_vec, const int M_curr,
              NumericVector & beta_vec, const double beta_curr,
              NumericVector & kappa_vec, const double kappa_curr,
              NumericVector & lambda_vec, const double lambda_curr,
              NumericVector & phi_vec, const double phi_curr,
              NumericVector & psi_vec, const double psi_curr,
              NumericVector & tau_vec, const double tau_curr,
              NumericVector & theta_vec, const double theta_curr,
              running_stat_vec<vec> & e_stat, const NumericVector e_curr, const NumericVector x_star, const NumericVector t_diff,
              NumericVector & e_first_vec, NumericVector & e_last_vec,
              running_stat_vec<vec> & m_estimated_stat, running_stat_vec<vec> & m_predicted_stat,
              NumericVector & chisq_actual_vec, NumericVector & chisq_estimated_vec, NumericVector & chisq_predicted_vec,
              const IntegerVector m_retweets,
              NumericVector & m_estimated,
              NumericVector & m_predicted) {
  // save parameters & output m_estimated & m_predicted
  // model pars
  M_vec[j] = M_curr;
  beta_vec[j] = beta_curr;
  kappa_vec[j] = kappa_curr;
  lambda_vec[j] = lambda_curr;
  phi_vec[j] = phi_curr;
  psi_vec[j] = psi_curr;
  tau_vec[j] = tau_curr;
  theta_vec[j] = theta_curr;
  e_stat(as<vec>(e_curr));
  e_first_vec[j] = e_curr[0];
  e_last_vec[j] = e_curr[n-1];
  // 2) actual, estimated & predicted #RT
  const NumericVector m_actual = (NumericVector) m_retweets;
  m_estimated = exp(lcum_intensity_etas(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, x_star, e_curr, t_diff));
  for (int i = 0; i < n; i++) {
    m_predicted[i] = (double) rpois(1, m_estimated[i])[0];
  }
  m_estimated_stat(as<vec>(m_estimated));
  m_predicted_stat(as<vec>(m_predicted));
  // 3) goodness-of-fit measures
  const NumericVector e_0(n, 0.0),
    m_0 = exp(lcum_intensity_etas(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, x_star, e_0, t_diff)); // #RT w/o errors
  const double E_ln = exp(0.5 / tau_curr),
    denom = exp(2.0 / tau_curr) - exp(1.0 / tau_curr),
    chisq_actual = sum(pow(m_actual / m_0 - E_ln, 2.0)) / denom,
    chisq_estimated = sum(pow(m_estimated / m_0 - E_ln, 2.0)) / denom, // or sum(pow(exp(e_curr) - E_ln, 2.0)) / denom
    chisq_predicted = sum(pow(m_predicted / m_0 - E_ln, 2.0)) / denom;
  chisq_actual_vec[j] = chisq_actual;
  chisq_estimated_vec[j] = chisq_estimated;
  chisq_predicted_vec[j] = chisq_predicted;
}

const NumericVector lden_e_etas(const double beta,
                                const double kappa,
                                const double lambda,
                                const double phi,
                                const double psi,
                                const double theta,
                                const NumericVector x,
                                const NumericVector e,
                                const NumericVector t_original,
                                const IntegerVector m_retweets,
                                const double t_infinity) {
  // log density that involves e
  if (is_true(any(t_original >= t_infinity))) {
    stop("lden_e_etas: all t_original have to be smaller than t_infinity.");
  }
  const int n = t_original.size();
  if (x.size() != n || e.size() != n || m_retweets.size() != n) {
    stop("lden_e_etas: x, e, t_original & m_retweets have to all be of same length.");
  }
  if (is_true(any(m_retweets < 0))) {
    stop("lden_e_etas: m_retweets (#RT) should be all non-negative");
  }
  if (is_true(any(x < 0.0))) {
    stop("lden_e_etas: x = log(#followers+1L) is non-negative by definition, even though a covariate in Gaussian regression is usually unrestricted.");
  }
  if (theta < 0.0) {
    stop("lden_e_etas: theta has to be non-negative.");
  }
  const NumericVector x_star = x - mean(x),
    t_diff = t_infinity - t_original,
    lvec_intensity = lcum_intensity_etas(beta, kappa, lambda, phi, psi, theta, x_star, e, t_diff),
    m_log_gamma = (NumericVector) m_retweets * (beta * x_star + kappa * x_star * x_star + e),
    lvec = m_log_gamma - exp(lvec_intensity);
  return ifelse(lvec != lvec, -INFINITY, lvec);
}

// [[Rcpp::export]]
const double llik_etas(const double beta,
                       const double kappa,
                       const double lambda,
                       const double phi,
                       const double psi,
                       const double theta,
                       const NumericVector x,
                       const NumericVector e,
                       const NumericVector t_original,
                       const IntegerVector m_retweets,
                       const NumericVector t_relative,
                       const double t_infinity) {
  // log-likelihood for NHPP modelling of retweets w/ regression by follower count
  // a) beta, kappa: coeff.s of lin. pred. for #RT upper limit of an original
  // b) lambda, theta: power & exp. decay rates in NHPP intensity
  // c) phi, psi: additional scale & location parameters
  // d) x & e: covariate i.e. log(followers_count+1L) & latent errors
  // e) t_original: times of originals rel. to when data collection begins
  // f) m_retweets: #RT for each original, sum to len. of t_relative
  // g) t_relative: times of RTs relative to those of their originals
  // h) t_infinity: time data collection ends rel. to when data collection begins; equivalently len. of obs. period
  if (is_true(any(t_relative <= 0.0))) {
    stop("llik_etas: all t_relative have to be positive.");
  }
  const int m = t_relative.size(); // # RTs
  if (sum(m_retweets) != m) {
    stop("llik_etas: sum of m_retweets has to equal length of t_relative.");
  }
  double llik;
  const NumericVector x_star = x - mean(x);
  if (lambda > 1.0 || phi < 0.0 || psi < 0.0 || theta < 0.0) {
    llik = -INFINITY;
  }
  else {
    const NumericVector lden = lden_e_etas(beta, kappa, lambda, phi, psi, theta, x, e, t_original, m_retweets, t_infinity);
    llik = sum(lden) - sum(t_relative) * theta - sum(log(t_relative + psi)) * lambda + (double) m * log(phi);
  }
  nan_to_minus_infinity(llik);
  return llik;
}

const double lpost_etas(const double beta,
                        const double kappa,
                        const double lambda,
                        const double phi,
                        const double psi,
                        const double theta,
                        const NumericVector x,
                        const NumericVector e,
                        const NumericVector t_original,
                        const IntegerVector m_retweets,
                        const NumericVector t_relative,
                        const double t_infinity,
                        const double mu_beta,
                        const double tau_beta,
                        const double mu_kappa,
                        const double tau_kappa,
                        const double a_lambda,
                        const double b_lambda,
                        const double a_phi,
                        const double b_phi,
                        const double a_psi,
                        const double b_psi,
                        const double a_theta,
                        const double b_theta) {
  double lpost = llik_etas(beta, kappa, lambda, phi, psi, theta, x, e, t_original, m_retweets, t_relative, t_infinity);
  if (lpost != -INFINITY) {
    // no need to do anything else if already -Inf
    const NumericVector hyper = NumericVector::create(tau_beta, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_theta, b_theta);
    if (is_true(any(hyper <= 0.0))) {
      lpost = -INFINITY;
    }
    else {
      lpost += dnorm(NumericVector::create(beta), mu_beta, 1.0 / sqrt(tau_beta), true)[0] +
        dnorm(NumericVector::create(kappa), mu_kappa, 1.0 / sqrt(tau_kappa), true)[0] +
        dgamma(NumericVector::create(1.0 - lambda), a_lambda, 1.0 / b_lambda, true)[0] +
        dgamma(NumericVector::create(phi), a_phi, 1.0 / b_phi, true)[0] +
        dgamma(NumericVector::create(psi), a_psi, 1.0 / b_psi, true)[0];
      if (theta > 0.0) {
        lpost += dgamma(NumericVector::create(theta), a_theta, 1.0 / b_theta, true)[0]; // theta_prop is never 0.0 in model 1; however, this matters in jump/selection algos!
      }
    }
    nan_to_minus_infinity(lpost);
  }
  return lpost;
}

//' Run the Metropolis-within-Gibbs algorithm for the hierarchical model of hybrid processes
//'
//' @param beta,kappa,lambda,phi,psi,tau,theta initial values of the parameters.
//' @param x,t_original,m_retweets,t_relative numeric vectors of follower counts, creation times of original tweets, number of retweets, and times of retweets relative to their corresponding original tweets.
//' @param t_infinity numeric scalar, duration of data collection.
//' @param s_beta,s_kappa,s_lambda,s_psi,s_theta,s_e_init initial values of proposal standard deviations.
//' @param N desired chain length AFTER burn-in and thinning.
//' @param thin thinning in MCMC.
//' @param burnin burn-in in MCMC.
//' @param print_freq how frequent should the current values in the MCMC be printed.
//' @param a_theta_prop shape parameter of pseudoprior for theta.
//' @param b_theta_prop rate parameter of pseudoprior for theta.
//' @param p0 prior for power-law process as opposed to hybrid process.
//' @param p01 jump probability from model 0 to 1 in RJMCMC.
//' @param p10 jump probability from model 1 to 0 in RJMCMC.
//' @param write boolean; should csv files of estimated and predicted retweet counts be written?
//' @param filename_par,filename_m,filename_p names of the csv files for the chains of the parameters, the estimated retweet counts and the predicted retweet counts, respectively; ignored if write is FALSE.
//' @param mu_beta,tau_beta,mu_kappa,tau_kappa,a_lambda,b_lambda,a_phi,b_phi,a_psi,b_psi,a_tau,b_tau,a_theta,b_theta hyperparameters for the priors of beta, kappa, lambda, phi, psi, tau, and theta, respectively.
//' @export
// [[Rcpp::export]]
List mh_etas(const double beta,
             const double kappa,
             const double lambda,
             const double phi,
             const double psi,
             const double tau,
             const double theta,
             const NumericVector x,
             const NumericVector t_original,
             const IntegerVector m_retweets,
             const NumericVector t_relative,
             const double t_infinity,
             double s_beta,
             double s_kappa,
             double s_lambda,
             double s_psi,
             double s_theta,
             const double s_e_init,
             const int N,
             const int thin,
             const int burnin,
             const int print_freq,
             const double a_theta_prop,
             const double b_theta_prop,
             const double p0,
             const double p01,
             const double p10,
             const bool write,
             std::string filename_par = "default_par.csv",
             std::string filename_m = "default_m.csv",
             std::string filename_p = "default_p.csv",
             const double mu_beta = 0.0,
             const double tau_beta = 1.0e-4,
             const double mu_kappa = 0.0,
             const double tau_kappa = 1.0e-4,
             const double a_lambda = 1.0,
             const double b_lambda = 1.0e-3,
             const double a_phi = 1.0,
             const double b_phi = 1.0e-3,
             const double a_psi = 1.0,
             const double b_psi = 1.0e-3,
             const double a_tau = 1.0,
             const double b_tau = 1.0e-3,
             const double a_theta = 1.0,
             const double b_theta = 1.0e-3) {
  // Metropolis-within-Gibbs algo. for individal models
  // 01) initialise & check
  if (phi <= 0.0 || psi < 0.0 || tau <= 0.0 || theta <= 0.0) {
    stop("mh_etas: Initial values of phi, psi, tau & theta can't be non-positive.");
  }
  if (lambda > 1.0) {
    stop("mh_etas: Initial value of lambda can't be >= 1.0.");
  }
  const NumericVector prop_sds = NumericVector::create(s_beta, s_kappa, s_lambda, s_psi, s_theta, s_e_init);
  if (is_true(any(prop_sds <= 0.0))) {
    stop("mh_etas: Initial proposal standard deviations can't be non-positive.");
  }
  const NumericVector hyper = NumericVector::create(tau_beta, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_tau, b_tau, a_theta, b_theta, a_theta_prop, b_theta_prop);
  if (is_true(any(hyper <= 0.0))) {
    stop("mh_etas: All hyperparameters except mu_beta & mu_kappa have to be positive.");
  }
  if (burnin <= 0) {
    stop("mh_etas: burnin has to be positive. (to be fixed)");
  }
  if (p0 > 1.0 || p0 < 0.0) {
    stop("mh_etas: p0 has to be b/w 0 & 1 inclusive.");
  }
  const int n = t_original.size(), m = t_relative.size();
  // 02) quantities for updating & saving
  int M_curr = (runif(1)[0] < p0) ? 0 : 1;
  double beta_curr = beta, beta_prop, kappa_curr = kappa, kappa_prop, lambda_curr = lambda, lambda_prop, phi_curr = phi, psi_curr = psi, psi_prop, tau_curr = tau, theta_curr = (M_curr == 1) ? theta : 0.0, theta_prop, lpost_curr, lpost_prop; // the usuals
  double a, b; // for tau
  NumericVector x_star = x - mean(x), e_curr(n, 0.0), e_prop(n), s_e(n, s_e_init), lvec_curr(n), lvec_prop(n); // for e
  NumericVector m_estimated(n), m_predicted(n);
  const NumericVector t_diff = t_infinity - t_original; // != t_relative!
  running_stat<double> M_stat, beta_stat, kappa_stat, lambda_stat, phi_stat, psi_stat, tau_stat, theta_stat, e_first_stat, e_last_stat;
  running_stat_vec<vec> e_stat, m_estimated_stat, m_predicted_stat;
  NumericVector beta_vec(N), kappa_vec(N), lambda_vec(N), phi_vec(N), psi_vec(N), tau_vec(N), theta_vec(N), e_first_vec(N), e_last_vec(N), chisq_actual_vec(N), chisq_estimated_vec(N), chisq_predicted_vec(N);
  IntegerVector M_vec(N);
  ofstream output_par, output_p, output_m;
  if (write) {
    output_par.open(filename_par);
    output_par << "M, beta, kappa, lambda, phi, psi, tau, theta, e_first, e_last, chisq_actual, chisq_estimated, chisq_predicted" << endl;
    output_m.open(filename_m);
    output_p.open(filename_p);
    for (int k = 0; k < n-1; k++) {
      output_m << "m_" << k+1 << ",";
      output_p << "p_" << k+1 << ",";
    }
    output_m << "m_" << n << endl;
    output_p << "p_" << n << endl;
  }
  // 03) run
  auto lpost = [x, t_original, m_retweets, t_relative, t_infinity, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_theta, b_theta](const double beta, const double kappa, const double lambda, const double phi, const double psi, const double theta, const NumericVector e) {
    return lpost_etas(beta, kappa, lambda, phi, psi, theta, x, e, t_original, m_retweets, t_relative, t_infinity, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_theta, b_theta);
  };
  int i, j, k;
  for (i = 0; i < N * thin + burnin; i++) {
    // a) update tau
    a = a_tau + 0.5 * (double) n;
    b = b_tau + 0.5 * sum(e_curr * e_curr);
    tau_curr = rgamma(1, a, 1.0 / b)[0];
    lpost_curr = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
    // c) update beta
    beta_prop = rnorm(1, beta_curr, s_beta)[0];
    lpost_prop = lpost(beta_prop, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
    update(beta_curr, beta_prop, lpost_curr, lpost_prop, s_beta, i, burnin);
    // d) update kappa
    kappa_prop = rnorm(1, kappa_curr, s_kappa)[0];
    lpost_prop = lpost(beta_curr, kappa_prop, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
    update(kappa_curr, kappa_prop, lpost_curr, lpost_prop, s_kappa, i, burnin);
    // e) update lambda
    lambda_prop = rnorm(1, lambda_curr, s_lambda)[0];
    lpost_prop = lpost(beta_curr, kappa_curr, lambda_prop, phi_curr, psi_curr, theta_curr, e_curr);
    update(lambda_curr, lambda_prop, lpost_curr, lpost_prop, s_lambda, i, burnin);
    // f) update psi
    psi_prop = rnorm(1, psi_curr, s_psi)[0];
    lpost_prop = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_prop, theta_curr, e_curr);
    update(psi_curr, psi_prop, lpost_curr, lpost_prop, s_psi, i, burnin);
    // g) update theta
    if (M_curr == 1) {
      theta_prop = rnorm(1, theta_curr, s_theta)[0];
      lpost_prop = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_prop, e_curr);
      update(theta_curr, theta_prop, lpost_curr, lpost_prop, s_theta, i, burnin);
    }
    // h) update e
    e_prop = e_curr + rnorm(n) * s_e;
    lvec_prop = lden_e_etas(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, x, e_prop, t_original, m_retweets, t_infinity) + dnorm(sqrt(tau_curr) * e_prop, 0.0, 1.0, true);
    lvec_curr = lden_e_etas(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, x, e_curr, t_original, m_retweets, t_infinity) + dnorm(sqrt(tau_curr) * e_curr, 0.0, 1.0, true);
    update_vec(e_curr, e_prop, lvec_curr, lvec_prop, s_e, i, burnin);
    // i) update phi
    lvec_curr = lcum_intensity_etas(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, x_star, e_curr, t_diff) - log(phi_curr);
    a = a_phi + (double) m;
    b = b_phi + sum(exp(lvec_curr));
    phi_curr = rgamma(1, a, 1.0 / b)[0];
    // j) print & save
    stat_pars(n, M_stat, M_curr, beta_stat, beta_curr, kappa_stat, kappa_curr, lambda_stat, lambda_curr, phi_stat, phi_curr, psi_stat, psi_curr, tau_stat, tau_curr, theta_stat, theta_curr, e_first_stat, e_last_stat, e_curr);
    if ((i + 1) % print_freq == 0) {
      Rcout << "Iteration " << i + 1 << endl;
      print_all(i, burnin, n, M_curr, M_stat, beta_curr, beta_stat, s_beta, kappa_curr, kappa_stat, s_kappa, lambda_curr, lambda_stat, s_lambda, phi_curr, phi_stat, psi_curr, psi_stat, s_psi, tau_curr, tau_stat, theta_curr, theta_stat, s_theta, e_curr, e_first_stat, e_last_stat);
    }
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      save_all(j, n, M_vec, M_curr, beta_vec, beta_curr, kappa_vec, kappa_curr, lambda_vec, lambda_curr, phi_vec, phi_curr, psi_vec, psi_curr, tau_vec, tau_curr, theta_vec, theta_curr, e_stat, e_curr, x_star, t_diff, e_first_vec, e_last_vec, m_estimated_stat, m_predicted_stat, chisq_actual_vec, chisq_estimated_vec, chisq_predicted_vec, m_retweets, m_estimated, m_predicted);
      // write to files
      if (write) {
        output_par << M_curr << "," << beta_curr << "," << kappa_curr << "," << lambda_curr << "," << phi_curr << "," << psi_curr << "," << tau_curr << "," << theta_curr << "," << e_curr[0] << "," << e_curr[n-1] << "," << chisq_actual_vec[j] << "," << chisq_estimated_vec[j] << "," << chisq_predicted_vec[j] << endl;
        for (k = 0; k < n-1; k++) {
          output_m << m_estimated[k] << ",";
          output_p << m_predicted[k] << ",";
        }
        output_m << m_estimated[n-1] << endl;
        output_p << m_predicted[n-1] << endl;
      }
    }
  }
  // 04) save
  if (write) {
    output_par.close();
    output_m.close();
    output_p.close();
  }
  List output = list_all(M_vec, beta_vec, kappa_vec, lambda_vec, phi_vec, psi_vec, tau_vec, theta_vec, e_first_vec, e_last_vec, chisq_actual_vec, chisq_estimated_vec, chisq_predicted_vec, x, t_original, m_retweets, m_estimated_stat, m_predicted_stat, e_stat, s_e, beta, kappa, lambda, phi, psi, tau, theta, s_beta, s_kappa, s_lambda, s_psi, s_theta, s_e_init, N, thin, burnin, print_freq, a_theta_prop, b_theta_prop, p0, p01, p10, write, t_infinity, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_tau, b_tau, a_theta, b_theta);
  return output;
}



// 03) model selection algorithms
//' Run the Reversible-jump MCMC algorithm for the hierarchical model of hybrid processes
//'
//' @param beta,kappa,lambda,phi,psi,tau,theta initial values of the parameters.
//' @param x,t_original,m_retweets,t_relative numeric vectors of follower counts, creation times of original tweets, number of retweets, and times of retweets relative to their corresponding original tweets.
//' @param t_infinity numeric scalar, duration of data collection.
//' @param s_beta,s_kappa,s_lambda,s_psi,s_theta,s_e_init initial values of proposal standard deviations.
//' @param N desired chain length AFTER burn-in and thinning.
//' @param thin thinning in MCMC.
//' @param burnin burn-in in MCMC.
//' @param print_freq how frequent should the current values in the MCMC be printed.
//' @param a_theta_prop shape parameter of pseudoprior for theta.
//' @param b_theta_prop rate parameter of pseudoprior for theta.
//' @param p0 prior for power-law process as opposed to hybrid process.
//' @param p01 jump probability from model 0 to 1 in RJMCMC.
//' @param p10 jump probability from model 1 to 0 in RJMCMC.
//' @param write boolean; should csv files of estimated and predicted retweet counts be written?
//' @param filename_par,filename_m,filename_p names of the csv files for the chains of the parameters, the estimated retweet counts and the predicted retweet counts, respectively; ignored if write is FALSE.
//' @param mu_beta,tau_beta,mu_kappa,tau_kappa,a_lambda,b_lambda,a_phi,b_phi,a_psi,b_psi,a_tau,b_tau,a_theta,b_theta hyperparameters for the priors of beta, kappa, lambda, phi, psi, tau, and theta, respectively.
//' @export
// [[Rcpp::export]]
List rj_etas(const double beta,
             const double kappa,
             const double lambda,
             const double phi,
             const double psi,
             const double tau,
             const double theta,
             const NumericVector x,
             const NumericVector t_original,
             const IntegerVector m_retweets,
             const NumericVector t_relative,
             const double t_infinity,
             double s_beta,
             double s_kappa,
             double s_lambda,
             double s_psi,
             double s_theta,
             const double s_e_init,
             const int N,
             const int thin,
             const int burnin,
             const int print_freq,
             const double a_theta_prop,
             const double b_theta_prop,
             const double p0,
             const double p01,
             const double p10,
             const bool write,
             std::string filename_par = "default_par.csv",
             std::string filename_m = "default_m.csv",
             std::string filename_p = "default_p.csv",
             const double mu_beta = 0.0,
             const double tau_beta = 1.0e-4,
             const double mu_kappa = 0.0,
             const double tau_kappa = 1.0e-4,
             const double a_lambda = 1.0,
             const double b_lambda = 1.0e-3,
             const double a_phi = 1.0,
             const double b_phi = 1.0e-3,
             const double a_psi = 1.0,
             const double b_psi = 1.0e-3,
             const double a_tau = 1.0,
             const double b_tau = 1.0e-3,
             const double a_theta = 1.0,
             const double b_theta = 1.0e-3) {
  // rjmcmc for model selection
  // 01) initialise & check
  if (phi <= 0.0 || psi < 0.0 || tau <= 0.0 || theta <= 0.0) {
    stop("rj_etas: Initial values of phi, psi, tau & theta can't be non-positive.");
  }
  if (lambda > 1.0) {
    stop("rj_etas: Initial value of lambda can't be >= 1.0.");
  }
  const NumericVector prop_sds = NumericVector::create(s_beta, s_kappa, s_lambda, s_psi, s_theta, s_e_init);
  if (is_true(any(prop_sds <= 0.0))) {
    stop("rj_etas: Initial proposal standard deviations can't be non-positive.");
  }
  const NumericVector hyper = NumericVector::create(tau_beta, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_tau, b_tau, a_theta, b_theta, a_theta_prop, b_theta_prop);
  if (is_true(any(hyper <= 0.0))) {
    stop("rj_etas: All hyperparameters except mu_beta & mu_kappa have to be positive.");
  }
  if (p01 <= 0.0 || p01 >= 1.0 || p10 <= 0.0 || p10 >= 1.0) {
    stop("rj_etas: p01 & p10 have to be b/w 0 & 1 exclusive.");
  }
  if (p0 > 1.0 || p0 < 0.0) {
    stop("rj_etas: p0 has to be b/w 0 & 1 inclusive.");
  }
  const int n = t_original.size(), m = t_relative.size();
  // 02) quantities for updating & saving
  int M_curr = (runif(1)[0] < p0) ? 0 : 1;
  double beta_curr = beta, beta_prop, kappa_curr = kappa, kappa_prop, lambda_curr = lambda, lambda_prop, phi_curr = phi, psi_curr = psi, psi_prop, tau_curr = tau, theta_curr = (M_curr == 0) ? 0.0 : theta, theta_prop, lpost_curr, lpost_prop; // the usuals
  double a, b, // for tau in both models
    lpost_0, lpost_1, lsudo_theta, log_A_0, log_A_1, la; // for model choice
  NumericVector x_star = x - mean(x), e_curr(n, 0.0), e_prop(n), s_e(n, s_e_init), lvec_curr(n), lvec_prop(n);
  NumericVector m_estimated(n), m_predicted(n);
  const NumericVector t_diff = t_infinity - t_original; // != t_relative!
  running_stat<double> M_stat, beta_stat, kappa_stat, lambda_stat, phi_stat, psi_stat, tau_stat, theta_stat, e_first_stat, e_last_stat;
  running_stat_vec<vec> e_stat, m_estimated_stat, m_predicted_stat;
  NumericVector beta_vec(N), kappa_vec(N), lambda_vec(N), phi_vec(N), psi_vec(N), tau_vec(N), theta_vec(N), e_first_vec(N), e_last_vec(N), chisq_actual_vec(N), chisq_estimated_vec(N), chisq_predicted_vec(N);
  IntegerVector M_vec(N);
  ofstream output_par, output_p, output_m;
  if (write) {
    output_par.open(filename_par);
    output_par << "M, beta, kappa, lambda, phi, psi, tau, theta, e_first, e_last, chisq_actual, chisq_estimated, chisq_predicted" << endl;
    output_m.open(filename_m);
    output_p.open(filename_p);
    for (int k = 0; k < n-1; k++) {
      output_m << "m_" << k+1 << ",";
      output_p << "p_" << k+1 << ",";
    }
    output_m << "m_" << n << endl;
    output_p << "p_" << n << endl;
  }
  // 03) run
  auto lpost = [x, t_original, m_retweets, t_relative, t_infinity, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_theta, b_theta](const double beta, const double kappa, const double lambda, const double phi, const double psi, const double theta, const NumericVector e) {
    return lpost_etas(beta, kappa, lambda, phi, psi, theta, x, e, t_original, m_retweets, t_relative, t_infinity, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_theta, b_theta);
  };
  int i, j, k;
  double p; // jump prob.
  for (i = 0; i < N * thin + burnin; i++) {
    p = runif(1)[0];
    // A) std MH if M_prop == M_curr
    if ((M_curr == 0 && p > p01) || (M_curr == 1 && p > p10)) { 
      // a) update tau
      a = a_tau + 0.5 * (double) n;
      b = b_tau + 0.5 * sum(e_curr * e_curr);
      tau_curr = rgamma(1, a, 1.0 / b)[0];
      lpost_curr = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
      // c) update beta
      beta_prop = rnorm(1, beta_curr, s_beta)[0];
      lpost_prop = lpost(beta_prop, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
      update(beta_curr, beta_prop, lpost_curr, lpost_prop, s_beta, i, burnin);
      // d) update kappa
      kappa_prop = rnorm(1, kappa_curr, s_kappa)[0];
      lpost_prop = lpost(beta_curr, kappa_prop, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
      update(kappa_curr, kappa_prop, lpost_curr, lpost_prop, s_kappa, i, burnin);
      // e) update lambda
      lambda_prop = rnorm(1, lambda_curr, s_lambda)[0];
      lpost_prop = lpost(beta_curr, kappa_curr, lambda_prop, phi_curr, psi_curr, theta_curr, e_curr);
      update(lambda_curr, lambda_prop, lpost_curr, lpost_prop, s_lambda, i, burnin);
      // f) update psi
      psi_prop = rnorm(1, psi_curr, s_psi)[0];
      lpost_prop = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_prop, theta_curr, e_curr);
      update(psi_curr, psi_prop, lpost_curr, lpost_prop, s_psi, i, burnin);
      // g) update theta
      if (M_curr == 1) {
        theta_prop = rnorm(1, theta_curr, s_theta)[0];
        lpost_prop = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_prop, e_curr);
        update(theta_curr, theta_prop, lpost_curr, lpost_prop, s_theta, i, burnin);
      }
      // h) update e
      e_prop = e_curr + rnorm(n) * s_e;
      lvec_prop = lden_e_etas(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, x, e_prop, t_original, m_retweets, t_infinity) + dnorm(sqrt(tau_curr) * e_prop, 0.0, 1.0, true);
      lvec_curr = lden_e_etas(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, x, e_curr, t_original, m_retweets, t_infinity) + dnorm(sqrt(tau_curr) * e_curr, 0.0, 1.0, true);
      update_vec(e_curr, e_prop, lvec_curr, lvec_prop, s_e, i, burnin);
      // i) update phi
      lvec_curr = lcum_intensity_etas(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, x_star, e_curr, t_diff) - log(phi_curr);
      a = a_phi + (double) m;
      b = b_phi + sum(exp(lvec_curr));
      phi_curr = rgamma(1, a, 1.0 / b)[0];
    }
    // B) proposed jump
    else {
      if (M_curr == 0) {
        theta_prop = rgamma(1, a_theta_prop, 1.0 / b_theta_prop)[0];
        lpost_0 = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
        lpost_1 = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_prop, e_curr);
        lsudo_theta = dgamma(NumericVector::create(theta_prop), a_theta_prop, 1.0 / b_theta_prop, true)[0];
      }
      else {
        theta_prop = 0.0;
        lpost_0 = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_prop, e_curr);
        lpost_1 = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
        lsudo_theta = dgamma(NumericVector::create(theta_curr), a_theta_prop, 1.0 / b_theta_prop, true)[0];
      }
      log_A_0 = lpost_0 + log(0.0 + p0) + lsudo_theta;
      log_A_1 = lpost_1 + log(1.0 - p0);
      la = (2.0 * (double) M_curr - 1.0) * (log_A_0 + log(p01) - log_A_1 - log(p10));
      if (log(runif(1)[0]) < la) {
        theta_curr = theta_prop;
        M_curr = 1 - M_curr; // jump
      }
    }
    // C) print & save
    stat_pars(n, M_stat, M_curr, beta_stat, beta_curr, kappa_stat, kappa_curr, lambda_stat, lambda_curr, phi_stat, phi_curr, psi_stat, psi_curr, tau_stat, tau_curr, theta_stat, theta_curr, e_first_stat, e_last_stat, e_curr);
    if ((i + 1) % print_freq == 0) {
      Rcout << "Iteration " << i + 1 << endl;
      print_all(i, burnin, n, M_curr, M_stat, beta_curr, beta_stat, s_beta, kappa_curr, kappa_stat, s_kappa, lambda_curr, lambda_stat, s_lambda, phi_curr, phi_stat, psi_curr, psi_stat, s_psi, tau_curr, tau_stat, theta_curr, theta_stat, s_theta, e_curr, e_first_stat, e_last_stat);
    }
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      save_all(j, n, M_vec, M_curr, beta_vec, beta_curr, kappa_vec, kappa_curr, lambda_vec, lambda_curr, phi_vec, phi_curr, psi_vec, psi_curr, tau_vec, tau_curr, theta_vec, theta_curr, e_stat, e_curr, x_star, t_diff, e_first_vec, e_last_vec, m_estimated_stat, m_predicted_stat, chisq_actual_vec, chisq_estimated_vec, chisq_predicted_vec, m_retweets, m_estimated, m_predicted);
      // write to files
      if (write) {
        output_par << M_curr << "," << beta_curr << "," << kappa_curr << "," << lambda_curr << "," << phi_curr << "," << psi_curr << "," << tau_curr << "," << theta_curr << "," << e_curr[0] << "," << e_curr[n-1] << "," << chisq_actual_vec[j] << "," << chisq_estimated_vec[j] << "," << chisq_predicted_vec[j] << endl;
        for (k = 0; k < n-1; k++) {
          output_m << m_estimated[k] << ",";
          output_p << m_predicted[k] << ",";
        }
        output_m << m_estimated[n-1] << endl;
        output_p << m_predicted[n-1] << endl;
      }
    }
  }
  // 04) save
  if (write) {
    output_par.close();
    output_m.close();
    output_p.close();
  }
  List output = list_all(M_vec, beta_vec, kappa_vec, lambda_vec, phi_vec, psi_vec, tau_vec, theta_vec, e_first_vec, e_last_vec, chisq_actual_vec, chisq_estimated_vec, chisq_predicted_vec, x, t_original, m_retweets, m_estimated_stat, m_predicted_stat, e_stat, s_e, beta, kappa, lambda, phi, psi, tau, theta, s_beta, s_kappa, s_lambda, s_psi, s_theta, s_e_init, N, thin, burnin, print_freq, a_theta_prop, b_theta_prop, p0, p01, p10, write, t_infinity, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_tau, b_tau, a_theta, b_theta);
  return output;
}

//' Run the Gibbs variable selection algorithm for the hierarchical model of hybrid processes
//'
//' @param beta,kappa,lambda,phi,psi,tau,theta initial values of the parameters.
//' @param x,t_original,m_retweets,t_relative numeric vectors of follower counts, creation times of original tweets, number of retweets, and times of retweets relative to their corresponding original tweets.
//' @param t_infinity numeric scalar, duration of data collection.
//' @param s_beta,s_kappa,s_lambda,s_psi,s_theta,s_e_init initial values of proposal standard deviations.
//' @param N desired chain length AFTER burn-in and thinning.
//' @param thin thinning in MCMC.
//' @param burnin burn-in in MCMC.
//' @param print_freq how frequent should the current values in the MCMC be printed.
//' @param a_theta_prop shape parameter of pseudoprior for theta.
//' @param b_theta_prop rate parameter of pseudoprior for theta.
//' @param p0 prior for power-law process as opposed to hybrid process.
//' @param p01 jump probability from model 0 to 1 in RJMCMC.
//' @param p10 jump probability from model 1 to 0 in RJMCMC.
//' @param write boolean; should csv files of estimated and predicted retweet counts be written?
//' @param filename_par,filename_m,filename_p names of the csv files for the chains of the parameters, the estimated retweet counts and the predicted retweet counts, respectively; ignored if write is FALSE.
//' @param mu_beta,tau_beta,mu_kappa,tau_kappa,a_lambda,b_lambda,a_phi,b_phi,a_psi,b_psi,a_tau,b_tau,a_theta,b_theta hyperparameters for the priors of beta, kappa, lambda, phi, psi, tau, and theta, respectively.
//' @export
// [[Rcpp::export]]
List ms_etas(const double beta,
             const double kappa,
             const double lambda,
             const double phi,
             const double psi,
             const double tau,
             const double theta,
             const NumericVector x,
             const NumericVector t_original,
             const IntegerVector m_retweets,
             const NumericVector t_relative,
             const double t_infinity,
             double s_beta,
             double s_kappa,
             double s_lambda,
             double s_psi,
             double s_theta,
             const double s_e_init,
             const int N,
             const int thin,
             const int burnin,
             const int print_freq,
             const double a_theta_prop,
             const double b_theta_prop,
             const double p0,
             const double p01,
             const double p10,
             const bool write,
             std::string filename_par = "default_par.csv",
             std::string filename_m = "default_m.csv",
             std::string filename_p = "default_p.csv",
             const double mu_beta = 0.0,
             const double tau_beta = 1.0e-4,
             const double mu_kappa = 0.0,
             const double tau_kappa = 1.0e-4,
             const double a_lambda = 1.0,
             const double b_lambda = 1.0e-3,
             const double a_phi = 1.0,
             const double b_phi = 1.0e-3,
             const double a_psi = 1.0,
             const double b_psi = 1.0e-3,
             const double a_tau = 1.0,
             const double b_tau = 1.0e-3,
             const double a_theta = 1.0,
             const double b_theta = 1.0e-3) {
  // model selection algo. for hierarchical model of NHPP w/ reg.
  // 01) initialise & check
  if (phi <= 0.0 || psi < 0.0 || tau <= 0.0 || theta <= 0.0) {
    stop("ms_etas: Initial values of phi, psi, tau & theta can't be non-positive.");
  }
  if (lambda > 1.0) {
    stop("ms_etas: Initial value of lambda can't be >= 1.0.");
  }
  const NumericVector prop_sds = NumericVector::create(s_beta, s_kappa, s_lambda, s_psi, s_theta, s_e_init);
  if (is_true(any(prop_sds <= 0.0))) {
    stop("ms_etas: Initial proposal standard deviations can't be non-positive.");
  }
  const NumericVector hyper = NumericVector::create(tau_beta, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_tau, b_tau, a_theta, b_theta, a_theta_prop, b_theta_prop);
  if (is_true(any(hyper <= 0.0))) {
    stop("ms_etas: All hyperparameters except mu_beta & mu_kappa have to be positive.");
  }
  if (burnin <= 0) {
    stop("ms_etas: burnin has to be positive. (to be fixed)");
  }
  if (p0 > 1.0 || p0 < 0.0) {
    stop("ms_etas: p0 has to be b/w 0 & 1 inclusive.");
  }
  const int n = t_original.size(), m = t_relative.size();
  // 02) quantities for updating & saving
  int M_curr = (runif(1)[0] < p0) ? 0 : 1;
  double beta_curr = beta, beta_prop, kappa_curr = kappa, kappa_prop, lambda_curr = lambda, lambda_prop, phi_curr = phi, psi_curr = psi, psi_prop, tau_curr = tau, theta_curr = (M_curr == 0) ? 0.0 : theta, theta_prop, lpost_curr, lpost_prop; // the usuals
  double a, b, // for tau in both models
    lpost_0, lpost_1, lsudo_theta, log_A_0, log_A_1, ratio; // for model choice
  NumericVector x_star = x - mean(x), e_curr(n, 0.0), e_prop(n), s_e(n, s_e_init), lvec_curr(n), lvec_prop(n);
  NumericVector m_estimated(n), m_predicted(n);
  const NumericVector t_diff = t_infinity - t_original; // != t_relative!
  running_stat<double> M_stat, beta_stat, kappa_stat, lambda_stat, phi_stat, psi_stat, tau_stat, theta_stat, e_first_stat, e_last_stat;
  running_stat_vec<vec> e_stat, m_estimated_stat, m_predicted_stat;
  NumericVector beta_vec(N), kappa_vec(N), lambda_vec(N), phi_vec(N), psi_vec(N), tau_vec(N), theta_vec(N), e_first_vec(N), e_last_vec(N), chisq_actual_vec(N), chisq_estimated_vec(N), chisq_predicted_vec(N);
  IntegerVector M_vec(N);
  ofstream output_par, output_p, output_m;
  if (write) {
    output_par.open(filename_par);
    output_par << "M, beta, kappa, lambda, phi, psi, tau, theta, e_first, e_last, chisq_actual, chisq_estimated, chisq_predicted" << endl;
    output_m.open(filename_m);
    output_p.open(filename_p);
    for (int k = 0; k < n-1; k++) {
      output_m << "m_" << k+1 << ",";
      output_p << "p_" << k+1 << ",";
    }
    output_m << "m_" << n << endl;
    output_p << "p_" << n << endl;
  }
  // 03) run
  auto lpost = [x, t_original, m_retweets, t_relative, t_infinity, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_theta, b_theta](const double beta, const double kappa, const double lambda, const double phi, const double psi, const double theta, const NumericVector e) {
    return lpost_etas(beta, kappa, lambda, phi, psi, theta, x, e, t_original, m_retweets, t_relative, t_infinity, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_theta, b_theta);
  };
  int i, j, k;
  for (i = 0; i < N * thin + burnin; i++) {
    // a) update tau
    a = a_tau + 0.5 * (double) n;
    b = b_tau + 0.5 * sum(e_curr * e_curr);
    tau_curr = rgamma(1, a, 1.0 / b)[0];
    lpost_curr = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
    // c) update beta
    beta_prop = rnorm(1, beta_curr, s_beta)[0];
    lpost_prop = lpost(beta_prop, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
    update(beta_curr, beta_prop, lpost_curr, lpost_prop, s_beta, i, burnin);
    // d) update kappa
    kappa_prop = rnorm(1, kappa_curr, s_kappa)[0];
    lpost_prop = lpost(beta_curr, kappa_prop, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
    update(kappa_curr, kappa_prop, lpost_curr, lpost_prop, s_kappa, i, burnin);
    // e) update lambda
    lambda_prop = rnorm(1, lambda_curr, s_lambda)[0];
    lpost_prop = lpost(beta_curr, kappa_curr, lambda_prop, phi_curr, psi_curr, theta_curr, e_curr);
    update(lambda_curr, lambda_prop, lpost_curr, lpost_prop, s_lambda, i, burnin);
    // f) update psi
    psi_prop = rnorm(1, psi_curr, s_psi)[0];
    lpost_prop = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_prop, theta_curr, e_curr);
    update(psi_curr, psi_prop, lpost_curr, lpost_prop, s_psi, i, burnin);
    // g) update theta
    if (M_curr == 1) {
      theta_prop = rnorm(1, theta_curr, s_theta)[0];
      lpost_prop = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_prop, e_curr);
      update(theta_curr, theta_prop, lpost_curr, lpost_prop, s_theta, i, burnin);
    }
    // h) update e
    e_prop = e_curr + rnorm(n) * s_e;
    lvec_prop = lden_e_etas(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, x, e_prop, t_original, m_retweets, t_infinity) + dnorm(sqrt(tau_curr) * e_prop, 0.0, 1.0, true);
    lvec_curr = lden_e_etas(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, x, e_curr, t_original, m_retweets, t_infinity) + dnorm(sqrt(tau_curr) * e_curr, 0.0, 1.0, true);
    update_vec(e_curr, e_prop, lvec_curr, lvec_prop, s_e, i, burnin);
    // i) update phi
    lvec_curr = lcum_intensity_etas(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, x_star, e_curr, t_diff) - log(phi_curr);
    a = a_phi + (double) m;
    b = b_phi + sum(exp(lvec_curr));
    phi_curr = rgamma(1, a, 1.0 / b)[0];
    // j) update M
    if (M_curr == 0) {
      theta_prop = rgamma(1, a_theta_prop, 1.0 / b_theta_prop)[0];
      lpost_0 = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
      lpost_1 = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_prop, e_curr);
      lsudo_theta = dgamma(NumericVector::create(theta_prop), a_theta_prop, 1.0 / b_theta_prop, true)[0];
    }
    else {
      theta_prop = 0.0;
      lpost_0 = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_prop, e_curr);
      lpost_1 = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
      lsudo_theta = dgamma(NumericVector::create(theta_curr), a_theta_prop, 1.0 / b_theta_prop, true)[0];
    }
    log_A_0 = lpost_0 + log(0.0 + p0) + lsudo_theta;
    log_A_1 = lpost_1 + log(1.0 - p0);
    ratio = 1.0 / (1.0 + exp(log_A_1 - log_A_0));
    if (runif(1)[0] < ratio) {
      // M becomes 0; if M STAYS 0, do nothing
      if (M_curr == 1) {
        M_curr = 0;
        theta_curr = 0.0;
      }
    }
    else {
      // M becomes 1; if M STAYS 1, do nothing
      if (M_curr == 0) {
        M_curr = 1;
        theta_curr = theta_prop;
      }
    }
    // k) print & save
    stat_pars(n, M_stat, M_curr, beta_stat, beta_curr, kappa_stat, kappa_curr, lambda_stat, lambda_curr, phi_stat, phi_curr, psi_stat, psi_curr, tau_stat, tau_curr, theta_stat, theta_curr, e_first_stat, e_last_stat, e_curr);
    if ((i + 1) % print_freq == 0) {
      Rcout << "Iteration " << i + 1 << endl;
      Rcout << "log_A_0: " << log_A_0 << endl;
      Rcout << "log_A_1: " << log_A_1 << endl;
      print_all(i, burnin, n, M_curr, M_stat, beta_curr, beta_stat, s_beta, kappa_curr, kappa_stat, s_kappa, lambda_curr, lambda_stat, s_lambda, phi_curr, phi_stat, psi_curr, psi_stat, s_psi, tau_curr, tau_stat, theta_curr, theta_stat, s_theta, e_curr, e_first_stat, e_last_stat);
    }
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      save_all(j, n, M_vec, M_curr, beta_vec, beta_curr, kappa_vec, kappa_curr, lambda_vec, lambda_curr, phi_vec, phi_curr, psi_vec, psi_curr, tau_vec, tau_curr, theta_vec, theta_curr, e_stat, e_curr, x_star, t_diff, e_first_vec, e_last_vec, m_estimated_stat, m_predicted_stat, chisq_actual_vec, chisq_estimated_vec, chisq_predicted_vec, m_retweets, m_estimated, m_predicted);
      // write to files
      if (write) {
        output_par << M_curr << "," << beta_curr << "," << kappa_curr << "," << lambda_curr << "," << phi_curr << "," << psi_curr << "," << tau_curr << "," << theta_curr << "," << e_curr[0] << "," << e_curr[n-1] << "," << chisq_actual_vec[j] << "," << chisq_estimated_vec[j] << "," << chisq_predicted_vec[j] << endl;
        for (k = 0; k < n-1; k++) {
          output_m << m_estimated[k] << ",";
          output_p << m_predicted[k] << ",";
        }
        output_m << m_estimated[n-1] << endl;
        output_p << m_predicted[n-1] << endl;
      }
    }
  }
  // 04) save
  if (write) {
    output_par.close();
    output_m.close();
    output_p.close();
  }
  List output = list_all(M_vec, beta_vec, kappa_vec, lambda_vec, phi_vec, psi_vec, tau_vec, theta_vec, e_first_vec, e_last_vec, chisq_actual_vec, chisq_estimated_vec, chisq_predicted_vec, x, t_original, m_retweets, m_estimated_stat, m_predicted_stat, e_stat, s_e, beta, kappa, lambda, phi, psi, tau, theta, s_beta, s_kappa, s_lambda, s_psi, s_theta, s_e_init, N, thin, burnin, print_freq, a_theta_prop, b_theta_prop, p0, p01, p10, write, t_infinity, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_tau, b_tau, a_theta, b_theta);
  return output;
}

//' Run the Metropolised Gibbs variable selection algorithm for the hierarchical model of hybrid processes
//'
//' @param beta,kappa,lambda,phi,psi,tau,theta initial values of the parameters.
//' @param x,t_original,m_retweets,t_relative numeric vectors of follower counts, creation times of original tweets, number of retweets, and times of retweets relative to their corresponding original tweets.
//' @param t_infinity numeric scalar, duration of data collection.
//' @param s_beta,s_kappa,s_lambda,s_psi,s_theta,s_e_init initial values of proposal standard deviations.
//' @param N desired chain length AFTER burn-in and thinning.
//' @param thin thinning in MCMC.
//' @param burnin burn-in in MCMC.
//' @param print_freq how frequent should the current values in the MCMC be printed.
//' @param a_theta_prop shape parameter of pseudoprior for theta.
//' @param b_theta_prop rate parameter of pseudoprior for theta.
//' @param p0 prior for power-law process as opposed to hybrid process.
//' @param p01 jump probability from model 0 to 1 in RJMCMC.
//' @param p10 jump probability from model 1 to 0 in RJMCMC.
//' @param write boolean; should csv files of estimated and predicted retweet counts be written?
//' @param filename_par,filename_m,filename_p names of the csv files for the chains of the parameters, the estimated retweet counts and the predicted retweet counts, respectively; ignored if write is FALSE.
//' @param mu_beta,tau_beta,mu_kappa,tau_kappa,a_lambda,b_lambda,a_phi,b_phi,a_psi,b_psi,a_tau,b_tau,a_theta,b_theta hyperparameters for the priors of beta, kappa, lambda, phi, psi, tau, and theta, respectively.
//' @export
// [[Rcpp::export]]
List me_etas(const double beta,
             const double kappa,
             const double lambda,
             const double phi,
             const double psi,
             const double tau,
             const double theta,
             const NumericVector x,
             const NumericVector t_original,
             const IntegerVector m_retweets,
             const NumericVector t_relative,
             const double t_infinity,
             double s_beta,
             double s_kappa,
             double s_lambda,
             double s_psi,
             double s_theta,
             const double s_e_init,
             const int N,
             const int thin,
             const int burnin,
             const int print_freq,
             const double a_theta_prop,
             const double b_theta_prop,
             const double p0,
             const double p01,
             const double p10,
             const bool write,
             std::string filename_par = "default_par.csv",
             std::string filename_m = "default_m.csv",
             std::string filename_p = "default_p.csv",
             const double mu_beta = 0.0,
             const double tau_beta = 1.0e-4,
             const double mu_kappa = 0.0,
             const double tau_kappa = 1.0e-4,
             const double a_lambda = 1.0,
             const double b_lambda = 1.0e-3,
             const double a_phi = 1.0,
             const double b_phi = 1.0e-3,
             const double a_psi = 1.0,
             const double b_psi = 1.0e-3,
             const double a_tau = 1.0,
             const double b_tau = 1.0e-3,
             const double a_theta = 1.0,
             const double b_theta = 1.0e-3) {
  // metropolised cc for model selection
  // 01) initialise & check
  if (phi <= 0.0 || psi < 0.0 || tau <= 0.0 || theta <= 0.0) {
    stop("me_etas: Initial values of phi, psi, tau & theta can't be non-positive.");
  }
  if (lambda > 1.0) {
    stop("me_etas: Initial value of lambda can't be >= 1.0.");
  }
  const NumericVector prop_sds = NumericVector::create(s_beta, s_kappa, s_lambda, s_psi, s_theta, s_e_init);
  if (is_true(any(prop_sds <= 0.0))) {
    stop("me_etas: Initial proposal standard deviations can't be non-positive.");
  }
  const NumericVector hyper = NumericVector::create(tau_beta, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_tau, b_tau, a_theta, b_theta, a_theta_prop, b_theta_prop);
  if (is_true(any(hyper <= 0.0))) {
    stop("me_etas: All hyperparameters except mu_beta & mu_kappa have to be positive.");
  }
  if (p01 <= 0.0 || p01 >= 1.0 || p10 <= 0.0 || p10 >= 1.0) {
    stop("me_etas: p01 & p10 have to be b/w 0 & 1 exclusive.");
  }
  if (p0 > 1.0 || p0 < 0.0) {
    stop("me_etas: p0 has to be b/w 0 & 1 inclusive.");
  }
  const int n = t_original.size(), m = t_relative.size();
  // 02) quantities for updating & saving
  int M_curr = (runif(1)[0] < p0) ? 0 : 1;
  double beta_curr = beta, beta_prop, kappa_curr = kappa, kappa_prop, lambda_curr = lambda, lambda_prop, phi_curr = phi, psi_curr = psi, psi_prop, tau_curr = tau, theta_curr = (M_curr == 0) ? 0.0 : theta, theta_prop, lpost_curr, lpost_prop;
  double a, b,
    lpost_0, lpost_1, lsudo_theta, log_A_0, log_A_1, la;
  NumericVector x_star = x - mean(x), e_curr(n, 0.0), e_prop(n), s_e(n, s_e_init), lvec_curr(n), lvec_prop(n);
  NumericVector m_estimated(n), m_predicted(n);
  const NumericVector t_diff = t_infinity - t_original; // != t_relative!
  running_stat<double> M_stat, beta_stat, kappa_stat, lambda_stat, phi_stat, psi_stat, tau_stat, theta_stat, e_first_stat, e_last_stat;
  running_stat_vec<vec> e_stat, m_estimated_stat, m_predicted_stat;
  NumericVector beta_vec(N), kappa_vec(N), lambda_vec(N), phi_vec(N), psi_vec(N), tau_vec(N), theta_vec(N), e_first_vec(N), e_last_vec(N), chisq_actual_vec(N), chisq_estimated_vec(N), chisq_predicted_vec(N);
  IntegerVector M_vec(N);
  ofstream output_par, output_p, output_m;
  if (write) {
    output_par.open(filename_par);
    output_par << "M, beta, kappa, lambda, phi, psi, tau, theta, e_first, e_last, chisq_actual, chisq_estimated, chisq_predicted" << endl;
    output_m.open(filename_m);
    output_p.open(filename_p);
    for (int k = 0; k < n-1; k++) {
      output_m << "m_" << k+1 << ",";
      output_p << "p_" << k+1 << ",";
    }
    output_m << "m_" << n << endl;
    output_p << "p_" << n << endl;
  }
  // 03) run
  auto lpost = [x, t_original, m_retweets, t_relative, t_infinity, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_theta, b_theta](const double beta, const double kappa, const double lambda, const double phi, const double psi, const double theta, const NumericVector e) {
    return lpost_etas(beta, kappa, lambda, phi, psi, theta, x, e, t_original, m_retweets, t_relative, t_infinity, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_theta, b_theta);
  };
  int i, j, k;
  double p;
  for (i = 0; i < N * thin + burnin; i++) {
    if (true) { // update eta no matter jump or not - diff. w/ rj
      // a) update tau
      a = a_tau + 0.5 * (double) n;
      b = b_tau + 0.5 * sum(e_curr * e_curr);
      tau_curr = rgamma(1, a, 1.0 / b)[0];
      lpost_curr = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
      // c) update beta
      beta_prop = rnorm(1, beta_curr, s_beta)[0];
      lpost_prop = lpost(beta_prop, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
      update(beta_curr, beta_prop, lpost_curr, lpost_prop, s_beta, i, burnin);
      // d) update kappa
      kappa_prop = rnorm(1, kappa_curr, s_kappa)[0];
      lpost_prop = lpost(beta_curr, kappa_prop, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
      update(kappa_curr, kappa_prop, lpost_curr, lpost_prop, s_kappa, i, burnin);
      // e) update lambda
      lambda_prop = rnorm(1, lambda_curr, s_lambda)[0];
      lpost_prop = lpost(beta_curr, kappa_curr, lambda_prop, phi_curr, psi_curr, theta_curr, e_curr);
      update(lambda_curr, lambda_prop, lpost_curr, lpost_prop, s_lambda, i, burnin);
      // f) update psi
      psi_prop = rnorm(1, psi_curr, s_psi)[0];
      lpost_prop = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_prop, theta_curr, e_curr);
      update(psi_curr, psi_prop, lpost_curr, lpost_prop, s_psi, i, burnin);
      // g) update theta
      if (M_curr == 1) {
        theta_prop = rnorm(1, theta_curr, s_theta)[0];
        lpost_prop = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_prop, e_curr);
        update(theta_curr, theta_prop, lpost_curr, lpost_prop, s_theta, i, burnin);
      }
      // h) update e
      e_prop = e_curr + rnorm(n) * s_e;
      lvec_prop = lden_e_etas(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, x, e_prop, t_original, m_retweets, t_infinity) + dnorm(sqrt(tau_curr) * e_prop, 0.0, 1.0, true);
      lvec_curr = lden_e_etas(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, x, e_curr, t_original, m_retweets, t_infinity) + dnorm(sqrt(tau_curr) * e_curr, 0.0, 1.0, true);
      update_vec(e_curr, e_prop, lvec_curr, lvec_prop, s_e, i, burnin);
      // i) update phi
      lvec_curr = lcum_intensity_etas(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, x_star, e_curr, t_diff) - log(phi_curr);
      a = a_phi + (double) m;
      b = b_phi + sum(exp(lvec_curr));
      phi_curr = rgamma(1, a, 1.0 / b)[0];
    }
    // j) proposed jump
    p = runif(1)[0];
    if ((M_curr == 0 && p > p01) || (M_curr == 1 && p > p10)) {
      if (M_curr == 0) {
        theta_prop = rgamma(1, a_theta_prop, 1.0 / b_theta_prop)[0];
        lpost_0 = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
        lpost_1 = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_prop, e_curr);
        lsudo_theta = dgamma(NumericVector::create(theta_prop), a_theta_prop, 1.0 / b_theta_prop, true)[0];
      }
      else {
        theta_prop = 0.0;
        lpost_0 = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_prop, e_curr);
        lpost_1 = lpost(beta_curr, kappa_curr, lambda_curr, phi_curr, psi_curr, theta_curr, e_curr);
        lsudo_theta = dgamma(NumericVector::create(theta_curr), a_theta_prop, 1.0 / b_theta_prop, true)[0];
      }
      log_A_0 = lpost_0 + log(0.0 + p0) + lsudo_theta;
      log_A_1 = lpost_1 + log(1.0 - p0);
      la = (2.0 * (double) M_curr - 1.0) * (log_A_0 + log(p01) - log_A_1 - log(p10));
      if (log(runif(1)[0]) < la) {
        theta_curr = theta_prop;
        M_curr = 1 - M_curr; // jump
      }
    }
    // k) print & save
    stat_pars(n, M_stat, M_curr, beta_stat, beta_curr, kappa_stat, kappa_curr, lambda_stat, lambda_curr, phi_stat, phi_curr, psi_stat, psi_curr, tau_stat, tau_curr, theta_stat, theta_curr, e_first_stat, e_last_stat, e_curr);
    if ((i + 1) % print_freq == 0) {
      Rcout << "Iteration " << i + 1 << endl;
      print_all(i, burnin, n, M_curr, M_stat, beta_curr, beta_stat, s_beta, kappa_curr, kappa_stat, s_kappa, lambda_curr, lambda_stat, s_lambda, phi_curr, phi_stat, psi_curr, psi_stat, s_psi, tau_curr, tau_stat, theta_curr, theta_stat, s_theta, e_curr, e_first_stat, e_last_stat);
    }
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      save_all(j, n, M_vec, M_curr, beta_vec, beta_curr, kappa_vec, kappa_curr, lambda_vec, lambda_curr, phi_vec, phi_curr, psi_vec, psi_curr, tau_vec, tau_curr, theta_vec, theta_curr, e_stat, e_curr, x_star, t_diff, e_first_vec, e_last_vec, m_estimated_stat, m_predicted_stat, chisq_actual_vec, chisq_estimated_vec, chisq_predicted_vec, m_retweets, m_estimated, m_predicted);
      // write to files
      if (write) {
        output_par << M_curr << "," << beta_curr << "," << kappa_curr << "," << lambda_curr << "," << phi_curr << "," << psi_curr << "," << tau_curr << "," << theta_curr << "," << e_curr[0] << "," << e_curr[n-1] << "," << chisq_actual_vec[j] << "," << chisq_estimated_vec[j] << "," << chisq_predicted_vec[j] << endl;
        for (k = 0; k < n-1; k++) {
          output_m << m_estimated[k] << ",";
          output_p << m_predicted[k] << ",";
        }
        output_m << m_estimated[n-1] << endl;
        output_p << m_predicted[n-1] << endl;
      }
    }
  }
  // 04) save
  if (write) {
    output_par.close();
    output_m.close();
    output_p.close();
  }
  List output = list_all(M_vec, beta_vec, kappa_vec, lambda_vec, phi_vec, psi_vec, tau_vec, theta_vec, e_first_vec, e_last_vec, chisq_actual_vec, chisq_estimated_vec, chisq_predicted_vec, x, t_original, m_retweets, m_estimated_stat, m_predicted_stat, e_stat, s_e, beta, kappa, lambda, phi, psi, tau, theta, s_beta, s_kappa, s_lambda, s_psi, s_theta, s_e_init, N, thin, burnin, print_freq, a_theta_prop, b_theta_prop, p0, p01, p10, write, t_infinity, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_tau, b_tau, a_theta, b_theta);
  return output;
}
