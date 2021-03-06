// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// llik_nhpp_hybrid
const double llik_nhpp_hybrid(const NumericVector par, const NumericVector x);
RcppExport SEXP _hybridProcess_llik_nhpp_hybrid(SEXP parSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(llik_nhpp_hybrid(par, x));
    return rcpp_result_gen;
END_RCPP
}
// llik_nhpp_exp
const double llik_nhpp_exp(const NumericVector par, const NumericVector x);
RcppExport SEXP _hybridProcess_llik_nhpp_exp(SEXP parSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(llik_nhpp_exp(par, x));
    return rcpp_result_gen;
END_RCPP
}
// llik_nhpp_power
const double llik_nhpp_power(const NumericVector par, const NumericVector x);
RcppExport SEXP _hybridProcess_llik_nhpp_power(SEXP parSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(llik_nhpp_power(par, x));
    return rcpp_result_gen;
END_RCPP
}
// llik_etas
const double llik_etas(const double beta, const double kappa, const double lambda, const double phi, const double psi, const double theta, const NumericVector x, const NumericVector e, const NumericVector t_original, const IntegerVector m_retweets, const NumericVector t_relative, const double t_infinity);
RcppExport SEXP _hybridProcess_llik_etas(SEXP betaSEXP, SEXP kappaSEXP, SEXP lambdaSEXP, SEXP phiSEXP, SEXP psiSEXP, SEXP thetaSEXP, SEXP xSEXP, SEXP eSEXP, SEXP t_originalSEXP, SEXP m_retweetsSEXP, SEXP t_relativeSEXP, SEXP t_infinitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type e(eSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type t_original(t_originalSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type m_retweets(m_retweetsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type t_relative(t_relativeSEXP);
    Rcpp::traits::input_parameter< const double >::type t_infinity(t_infinitySEXP);
    rcpp_result_gen = Rcpp::wrap(llik_etas(beta, kappa, lambda, phi, psi, theta, x, e, t_original, m_retweets, t_relative, t_infinity));
    return rcpp_result_gen;
END_RCPP
}
// mh_etas
List mh_etas(const double beta, const double kappa, const double lambda, const double phi, const double psi, const double tau, const double theta, const NumericVector x, const NumericVector t_original, const IntegerVector m_retweets, const NumericVector t_relative, const double t_infinity, double s_beta, double s_kappa, double s_lambda, double s_psi, double s_theta, const double s_e_init, const int N, const int thin, const int burnin, const int print_freq, const double a_theta_prop, const double b_theta_prop, const double p0, const double p01, const double p10, const bool write, std::string filename_par, std::string filename_m, std::string filename_p, const double mu_beta, const double tau_beta, const double mu_kappa, const double tau_kappa, const double a_lambda, const double b_lambda, const double a_phi, const double b_phi, const double a_psi, const double b_psi, const double a_tau, const double b_tau, const double a_theta, const double b_theta);
RcppExport SEXP _hybridProcess_mh_etas(SEXP betaSEXP, SEXP kappaSEXP, SEXP lambdaSEXP, SEXP phiSEXP, SEXP psiSEXP, SEXP tauSEXP, SEXP thetaSEXP, SEXP xSEXP, SEXP t_originalSEXP, SEXP m_retweetsSEXP, SEXP t_relativeSEXP, SEXP t_infinitySEXP, SEXP s_betaSEXP, SEXP s_kappaSEXP, SEXP s_lambdaSEXP, SEXP s_psiSEXP, SEXP s_thetaSEXP, SEXP s_e_initSEXP, SEXP NSEXP, SEXP thinSEXP, SEXP burninSEXP, SEXP print_freqSEXP, SEXP a_theta_propSEXP, SEXP b_theta_propSEXP, SEXP p0SEXP, SEXP p01SEXP, SEXP p10SEXP, SEXP writeSEXP, SEXP filename_parSEXP, SEXP filename_mSEXP, SEXP filename_pSEXP, SEXP mu_betaSEXP, SEXP tau_betaSEXP, SEXP mu_kappaSEXP, SEXP tau_kappaSEXP, SEXP a_lambdaSEXP, SEXP b_lambdaSEXP, SEXP a_phiSEXP, SEXP b_phiSEXP, SEXP a_psiSEXP, SEXP b_psiSEXP, SEXP a_tauSEXP, SEXP b_tauSEXP, SEXP a_thetaSEXP, SEXP b_thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type t_original(t_originalSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type m_retweets(m_retweetsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type t_relative(t_relativeSEXP);
    Rcpp::traits::input_parameter< const double >::type t_infinity(t_infinitySEXP);
    Rcpp::traits::input_parameter< double >::type s_beta(s_betaSEXP);
    Rcpp::traits::input_parameter< double >::type s_kappa(s_kappaSEXP);
    Rcpp::traits::input_parameter< double >::type s_lambda(s_lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type s_psi(s_psiSEXP);
    Rcpp::traits::input_parameter< double >::type s_theta(s_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type s_e_init(s_e_initSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const int >::type print_freq(print_freqSEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta_prop(a_theta_propSEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta_prop(b_theta_propSEXP);
    Rcpp::traits::input_parameter< const double >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< const double >::type p01(p01SEXP);
    Rcpp::traits::input_parameter< const double >::type p10(p10SEXP);
    Rcpp::traits::input_parameter< const bool >::type write(writeSEXP);
    Rcpp::traits::input_parameter< std::string >::type filename_par(filename_parSEXP);
    Rcpp::traits::input_parameter< std::string >::type filename_m(filename_mSEXP);
    Rcpp::traits::input_parameter< std::string >::type filename_p(filename_pSEXP);
    Rcpp::traits::input_parameter< const double >::type mu_beta(mu_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type tau_beta(tau_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type mu_kappa(mu_kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type tau_kappa(tau_kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_lambda(a_lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_lambda(b_lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_phi(a_phiSEXP);
    Rcpp::traits::input_parameter< const double >::type b_phi(b_phiSEXP);
    Rcpp::traits::input_parameter< const double >::type a_psi(a_psiSEXP);
    Rcpp::traits::input_parameter< const double >::type b_psi(b_psiSEXP);
    Rcpp::traits::input_parameter< const double >::type a_tau(a_tauSEXP);
    Rcpp::traits::input_parameter< const double >::type b_tau(b_tauSEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta(a_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta(b_thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(mh_etas(beta, kappa, lambda, phi, psi, tau, theta, x, t_original, m_retweets, t_relative, t_infinity, s_beta, s_kappa, s_lambda, s_psi, s_theta, s_e_init, N, thin, burnin, print_freq, a_theta_prop, b_theta_prop, p0, p01, p10, write, filename_par, filename_m, filename_p, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_tau, b_tau, a_theta, b_theta));
    return rcpp_result_gen;
END_RCPP
}
// rj_etas
List rj_etas(const double beta, const double kappa, const double lambda, const double phi, const double psi, const double tau, const double theta, const NumericVector x, const NumericVector t_original, const IntegerVector m_retweets, const NumericVector t_relative, const double t_infinity, double s_beta, double s_kappa, double s_lambda, double s_psi, double s_theta, const double s_e_init, const int N, const int thin, const int burnin, const int print_freq, const double a_theta_prop, const double b_theta_prop, const double p0, const double p01, const double p10, const bool write, std::string filename_par, std::string filename_m, std::string filename_p, const double mu_beta, const double tau_beta, const double mu_kappa, const double tau_kappa, const double a_lambda, const double b_lambda, const double a_phi, const double b_phi, const double a_psi, const double b_psi, const double a_tau, const double b_tau, const double a_theta, const double b_theta);
RcppExport SEXP _hybridProcess_rj_etas(SEXP betaSEXP, SEXP kappaSEXP, SEXP lambdaSEXP, SEXP phiSEXP, SEXP psiSEXP, SEXP tauSEXP, SEXP thetaSEXP, SEXP xSEXP, SEXP t_originalSEXP, SEXP m_retweetsSEXP, SEXP t_relativeSEXP, SEXP t_infinitySEXP, SEXP s_betaSEXP, SEXP s_kappaSEXP, SEXP s_lambdaSEXP, SEXP s_psiSEXP, SEXP s_thetaSEXP, SEXP s_e_initSEXP, SEXP NSEXP, SEXP thinSEXP, SEXP burninSEXP, SEXP print_freqSEXP, SEXP a_theta_propSEXP, SEXP b_theta_propSEXP, SEXP p0SEXP, SEXP p01SEXP, SEXP p10SEXP, SEXP writeSEXP, SEXP filename_parSEXP, SEXP filename_mSEXP, SEXP filename_pSEXP, SEXP mu_betaSEXP, SEXP tau_betaSEXP, SEXP mu_kappaSEXP, SEXP tau_kappaSEXP, SEXP a_lambdaSEXP, SEXP b_lambdaSEXP, SEXP a_phiSEXP, SEXP b_phiSEXP, SEXP a_psiSEXP, SEXP b_psiSEXP, SEXP a_tauSEXP, SEXP b_tauSEXP, SEXP a_thetaSEXP, SEXP b_thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type t_original(t_originalSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type m_retweets(m_retweetsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type t_relative(t_relativeSEXP);
    Rcpp::traits::input_parameter< const double >::type t_infinity(t_infinitySEXP);
    Rcpp::traits::input_parameter< double >::type s_beta(s_betaSEXP);
    Rcpp::traits::input_parameter< double >::type s_kappa(s_kappaSEXP);
    Rcpp::traits::input_parameter< double >::type s_lambda(s_lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type s_psi(s_psiSEXP);
    Rcpp::traits::input_parameter< double >::type s_theta(s_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type s_e_init(s_e_initSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const int >::type print_freq(print_freqSEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta_prop(a_theta_propSEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta_prop(b_theta_propSEXP);
    Rcpp::traits::input_parameter< const double >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< const double >::type p01(p01SEXP);
    Rcpp::traits::input_parameter< const double >::type p10(p10SEXP);
    Rcpp::traits::input_parameter< const bool >::type write(writeSEXP);
    Rcpp::traits::input_parameter< std::string >::type filename_par(filename_parSEXP);
    Rcpp::traits::input_parameter< std::string >::type filename_m(filename_mSEXP);
    Rcpp::traits::input_parameter< std::string >::type filename_p(filename_pSEXP);
    Rcpp::traits::input_parameter< const double >::type mu_beta(mu_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type tau_beta(tau_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type mu_kappa(mu_kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type tau_kappa(tau_kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_lambda(a_lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_lambda(b_lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_phi(a_phiSEXP);
    Rcpp::traits::input_parameter< const double >::type b_phi(b_phiSEXP);
    Rcpp::traits::input_parameter< const double >::type a_psi(a_psiSEXP);
    Rcpp::traits::input_parameter< const double >::type b_psi(b_psiSEXP);
    Rcpp::traits::input_parameter< const double >::type a_tau(a_tauSEXP);
    Rcpp::traits::input_parameter< const double >::type b_tau(b_tauSEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta(a_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta(b_thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(rj_etas(beta, kappa, lambda, phi, psi, tau, theta, x, t_original, m_retweets, t_relative, t_infinity, s_beta, s_kappa, s_lambda, s_psi, s_theta, s_e_init, N, thin, burnin, print_freq, a_theta_prop, b_theta_prop, p0, p01, p10, write, filename_par, filename_m, filename_p, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_tau, b_tau, a_theta, b_theta));
    return rcpp_result_gen;
END_RCPP
}
// ms_etas
List ms_etas(const double beta, const double kappa, const double lambda, const double phi, const double psi, const double tau, const double theta, const NumericVector x, const NumericVector t_original, const IntegerVector m_retweets, const NumericVector t_relative, const double t_infinity, double s_beta, double s_kappa, double s_lambda, double s_psi, double s_theta, const double s_e_init, const int N, const int thin, const int burnin, const int print_freq, const double a_theta_prop, const double b_theta_prop, const double p0, const double p01, const double p10, const bool write, std::string filename_par, std::string filename_m, std::string filename_p, const double mu_beta, const double tau_beta, const double mu_kappa, const double tau_kappa, const double a_lambda, const double b_lambda, const double a_phi, const double b_phi, const double a_psi, const double b_psi, const double a_tau, const double b_tau, const double a_theta, const double b_theta);
RcppExport SEXP _hybridProcess_ms_etas(SEXP betaSEXP, SEXP kappaSEXP, SEXP lambdaSEXP, SEXP phiSEXP, SEXP psiSEXP, SEXP tauSEXP, SEXP thetaSEXP, SEXP xSEXP, SEXP t_originalSEXP, SEXP m_retweetsSEXP, SEXP t_relativeSEXP, SEXP t_infinitySEXP, SEXP s_betaSEXP, SEXP s_kappaSEXP, SEXP s_lambdaSEXP, SEXP s_psiSEXP, SEXP s_thetaSEXP, SEXP s_e_initSEXP, SEXP NSEXP, SEXP thinSEXP, SEXP burninSEXP, SEXP print_freqSEXP, SEXP a_theta_propSEXP, SEXP b_theta_propSEXP, SEXP p0SEXP, SEXP p01SEXP, SEXP p10SEXP, SEXP writeSEXP, SEXP filename_parSEXP, SEXP filename_mSEXP, SEXP filename_pSEXP, SEXP mu_betaSEXP, SEXP tau_betaSEXP, SEXP mu_kappaSEXP, SEXP tau_kappaSEXP, SEXP a_lambdaSEXP, SEXP b_lambdaSEXP, SEXP a_phiSEXP, SEXP b_phiSEXP, SEXP a_psiSEXP, SEXP b_psiSEXP, SEXP a_tauSEXP, SEXP b_tauSEXP, SEXP a_thetaSEXP, SEXP b_thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type t_original(t_originalSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type m_retweets(m_retweetsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type t_relative(t_relativeSEXP);
    Rcpp::traits::input_parameter< const double >::type t_infinity(t_infinitySEXP);
    Rcpp::traits::input_parameter< double >::type s_beta(s_betaSEXP);
    Rcpp::traits::input_parameter< double >::type s_kappa(s_kappaSEXP);
    Rcpp::traits::input_parameter< double >::type s_lambda(s_lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type s_psi(s_psiSEXP);
    Rcpp::traits::input_parameter< double >::type s_theta(s_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type s_e_init(s_e_initSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const int >::type print_freq(print_freqSEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta_prop(a_theta_propSEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta_prop(b_theta_propSEXP);
    Rcpp::traits::input_parameter< const double >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< const double >::type p01(p01SEXP);
    Rcpp::traits::input_parameter< const double >::type p10(p10SEXP);
    Rcpp::traits::input_parameter< const bool >::type write(writeSEXP);
    Rcpp::traits::input_parameter< std::string >::type filename_par(filename_parSEXP);
    Rcpp::traits::input_parameter< std::string >::type filename_m(filename_mSEXP);
    Rcpp::traits::input_parameter< std::string >::type filename_p(filename_pSEXP);
    Rcpp::traits::input_parameter< const double >::type mu_beta(mu_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type tau_beta(tau_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type mu_kappa(mu_kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type tau_kappa(tau_kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_lambda(a_lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_lambda(b_lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_phi(a_phiSEXP);
    Rcpp::traits::input_parameter< const double >::type b_phi(b_phiSEXP);
    Rcpp::traits::input_parameter< const double >::type a_psi(a_psiSEXP);
    Rcpp::traits::input_parameter< const double >::type b_psi(b_psiSEXP);
    Rcpp::traits::input_parameter< const double >::type a_tau(a_tauSEXP);
    Rcpp::traits::input_parameter< const double >::type b_tau(b_tauSEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta(a_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta(b_thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(ms_etas(beta, kappa, lambda, phi, psi, tau, theta, x, t_original, m_retweets, t_relative, t_infinity, s_beta, s_kappa, s_lambda, s_psi, s_theta, s_e_init, N, thin, burnin, print_freq, a_theta_prop, b_theta_prop, p0, p01, p10, write, filename_par, filename_m, filename_p, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_tau, b_tau, a_theta, b_theta));
    return rcpp_result_gen;
END_RCPP
}
// me_etas
List me_etas(const double beta, const double kappa, const double lambda, const double phi, const double psi, const double tau, const double theta, const NumericVector x, const NumericVector t_original, const IntegerVector m_retweets, const NumericVector t_relative, const double t_infinity, double s_beta, double s_kappa, double s_lambda, double s_psi, double s_theta, const double s_e_init, const int N, const int thin, const int burnin, const int print_freq, const double a_theta_prop, const double b_theta_prop, const double p0, const double p01, const double p10, const bool write, std::string filename_par, std::string filename_m, std::string filename_p, const double mu_beta, const double tau_beta, const double mu_kappa, const double tau_kappa, const double a_lambda, const double b_lambda, const double a_phi, const double b_phi, const double a_psi, const double b_psi, const double a_tau, const double b_tau, const double a_theta, const double b_theta);
RcppExport SEXP _hybridProcess_me_etas(SEXP betaSEXP, SEXP kappaSEXP, SEXP lambdaSEXP, SEXP phiSEXP, SEXP psiSEXP, SEXP tauSEXP, SEXP thetaSEXP, SEXP xSEXP, SEXP t_originalSEXP, SEXP m_retweetsSEXP, SEXP t_relativeSEXP, SEXP t_infinitySEXP, SEXP s_betaSEXP, SEXP s_kappaSEXP, SEXP s_lambdaSEXP, SEXP s_psiSEXP, SEXP s_thetaSEXP, SEXP s_e_initSEXP, SEXP NSEXP, SEXP thinSEXP, SEXP burninSEXP, SEXP print_freqSEXP, SEXP a_theta_propSEXP, SEXP b_theta_propSEXP, SEXP p0SEXP, SEXP p01SEXP, SEXP p10SEXP, SEXP writeSEXP, SEXP filename_parSEXP, SEXP filename_mSEXP, SEXP filename_pSEXP, SEXP mu_betaSEXP, SEXP tau_betaSEXP, SEXP mu_kappaSEXP, SEXP tau_kappaSEXP, SEXP a_lambdaSEXP, SEXP b_lambdaSEXP, SEXP a_phiSEXP, SEXP b_phiSEXP, SEXP a_psiSEXP, SEXP b_psiSEXP, SEXP a_tauSEXP, SEXP b_tauSEXP, SEXP a_thetaSEXP, SEXP b_thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type t_original(t_originalSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type m_retweets(m_retweetsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type t_relative(t_relativeSEXP);
    Rcpp::traits::input_parameter< const double >::type t_infinity(t_infinitySEXP);
    Rcpp::traits::input_parameter< double >::type s_beta(s_betaSEXP);
    Rcpp::traits::input_parameter< double >::type s_kappa(s_kappaSEXP);
    Rcpp::traits::input_parameter< double >::type s_lambda(s_lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type s_psi(s_psiSEXP);
    Rcpp::traits::input_parameter< double >::type s_theta(s_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type s_e_init(s_e_initSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const int >::type print_freq(print_freqSEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta_prop(a_theta_propSEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta_prop(b_theta_propSEXP);
    Rcpp::traits::input_parameter< const double >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< const double >::type p01(p01SEXP);
    Rcpp::traits::input_parameter< const double >::type p10(p10SEXP);
    Rcpp::traits::input_parameter< const bool >::type write(writeSEXP);
    Rcpp::traits::input_parameter< std::string >::type filename_par(filename_parSEXP);
    Rcpp::traits::input_parameter< std::string >::type filename_m(filename_mSEXP);
    Rcpp::traits::input_parameter< std::string >::type filename_p(filename_pSEXP);
    Rcpp::traits::input_parameter< const double >::type mu_beta(mu_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type tau_beta(tau_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type mu_kappa(mu_kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type tau_kappa(tau_kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_lambda(a_lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_lambda(b_lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_phi(a_phiSEXP);
    Rcpp::traits::input_parameter< const double >::type b_phi(b_phiSEXP);
    Rcpp::traits::input_parameter< const double >::type a_psi(a_psiSEXP);
    Rcpp::traits::input_parameter< const double >::type b_psi(b_psiSEXP);
    Rcpp::traits::input_parameter< const double >::type a_tau(a_tauSEXP);
    Rcpp::traits::input_parameter< const double >::type b_tau(b_tauSEXP);
    Rcpp::traits::input_parameter< const double >::type a_theta(a_thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type b_theta(b_thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(me_etas(beta, kappa, lambda, phi, psi, tau, theta, x, t_original, m_retweets, t_relative, t_infinity, s_beta, s_kappa, s_lambda, s_psi, s_theta, s_e_init, N, thin, burnin, print_freq, a_theta_prop, b_theta_prop, p0, p01, p10, write, filename_par, filename_m, filename_p, mu_beta, tau_beta, mu_kappa, tau_kappa, a_lambda, b_lambda, a_phi, b_phi, a_psi, b_psi, a_tau, b_tau, a_theta, b_theta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hybridProcess_llik_nhpp_hybrid", (DL_FUNC) &_hybridProcess_llik_nhpp_hybrid, 2},
    {"_hybridProcess_llik_nhpp_exp", (DL_FUNC) &_hybridProcess_llik_nhpp_exp, 2},
    {"_hybridProcess_llik_nhpp_power", (DL_FUNC) &_hybridProcess_llik_nhpp_power, 2},
    {"_hybridProcess_llik_etas", (DL_FUNC) &_hybridProcess_llik_etas, 12},
    {"_hybridProcess_mh_etas", (DL_FUNC) &_hybridProcess_mh_etas, 45},
    {"_hybridProcess_rj_etas", (DL_FUNC) &_hybridProcess_rj_etas, 45},
    {"_hybridProcess_ms_etas", (DL_FUNC) &_hybridProcess_ms_etas, 45},
    {"_hybridProcess_me_etas", (DL_FUNC) &_hybridProcess_me_etas, 45},
    {NULL, NULL, 0}
};

RcppExport void R_init_hybridProcess(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
