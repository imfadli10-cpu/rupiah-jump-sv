// =====================================================================
// svjx_sampler.cpp
// Fast draw of the latent log-variance path h (length n) from a Gaussian
// with tridiagonal precision matrix, used in the Gibbs sampler of the
// SV / SVJ / SVJ-X models (all-without-a-loop sampler; Rue 2001,
// McCausland et al. 2011, Kastner and Fruhwirth-Schnatter 2014).
//
// Posterior: Q h = b,  Q = Q_prior(phi, sigma2) + diag(1 / v_obs)
//            b = Q_prior * mu * 1 + (ytilde - m_obs) / v_obs
// Returns one draw h ~ N(Q^{-1} b, Q^{-1}).
// =====================================================================
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector draw_h_tridiag(NumericVector ytilde, NumericVector m_obs, NumericVector v_obs,
                             double mu, double phi, double sigma2, NumericVector z) {
  const int n = ytilde.size();
  if (n < 2) stop("need at least two observations");
  if (m_obs.size() != n || v_obs.size() != n || z.size() != n) stop("length mismatch");
  if (!(sigma2 > 0.0) || !(std::fabs(phi) < 1.0)) stop("invalid phi or sigma2");

  const double a = 1.0 / sigma2;
  std::vector<double> diag(n), off(n - 1, -a * phi), b(n);

  for (int t = 0; t < n; ++t) {
    double qd = (t == 0 || t == n - 1) ? a : a * (1.0 + phi * phi);
    double qb = (t == 0 || t == n - 1) ? a * (1.0 - phi) * mu
                                       : a * (1.0 - phi) * (1.0 - phi) * mu;
    diag[t] = qd + 1.0 / v_obs[t];
    b[t]    = qb + (ytilde[t] - m_obs[t]) / v_obs[t];
  }

  // Cholesky Q = L L', L lower bidiagonal with diagonal d and sub-diagonal l
  std::vector<double> d(n), l(n - 1);
  d[0] = std::sqrt(diag[0]);
  for (int t = 1; t < n; ++t) {
    l[t - 1] = off[t - 1] / d[t - 1];
    double s = diag[t] - l[t - 1] * l[t - 1];
    if (!(s > 0.0)) stop("precision matrix not positive definite");
    d[t] = std::sqrt(s);
  }

  // Solve L w = b (forward), then L' x = w + z (backward):
  // x = Q^{-1} b + L'^{-1} z  ~  N(Q^{-1} b, Q^{-1})
  std::vector<double> w(n);
  w[0] = b[0] / d[0];
  for (int t = 1; t < n; ++t) w[t] = (b[t] - l[t - 1] * w[t - 1]) / d[t];

  NumericVector h(n);
  h[n - 1] = (w[n - 1] + z[n - 1]) / d[n - 1];
  for (int t = n - 2; t >= 0; --t) h[t] = (w[t] + z[t] - l[t] * h[t + 1]) / d[t];
  return h;
}
