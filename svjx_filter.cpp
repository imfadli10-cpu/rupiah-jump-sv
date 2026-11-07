// =====================================================================
// svjx_filter.cpp
// Bootstrap particle filter for the SV / SVJ / SVJ-X models with the
// jump indicator and jump size integrated out analytically:
//
//   p(y_t | h_t) = p_t N(y_t - m; mu_xi, sigma_xi^2 + e^{h_t})
//                + (1 - p_t) N(y_t - m; 0, e^{h_t})      (m: daily drift)
//
// For each day t it returns, given information up to t-1:
//   logscore  log p(y_t | y_1:t-1)          (one-step-ahead log predictive density)
//   pit       P(Y_t <= y_t | y_1:t-1)       (probability integral transform)
//   q         one-step-ahead predictive quantiles at the levels in 'probs'
//   h_filt    E(h_t | y_1:t)                (filtered log-variance)
// The jump probability p_t is supplied by the caller (0 for the SV model).
// Random numbers come from R's generator, so set.seed() makes it reproducible.
// =====================================================================
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

static inline double mix_cdf(double y, double p, double mu_xi, double s2_xi, double eh) {
  double c0 = R::pnorm(y, 0.0, std::sqrt(eh), 1, 0);
  if (p <= 0.0) return c0;
  double c1 = R::pnorm(y, mu_xi, std::sqrt(s2_xi + eh), 1, 0);
  return p * c1 + (1.0 - p) * c0;
}

// [[Rcpp::export]]
List pf_svjx(NumericVector y, NumericVector p, double m, double mu, double phi, double sigma,
             double mu_xi, double sigma_xi, int n_particles, NumericVector probs,
             int n_quantile_particles = 400) {
  RNGScope scope;
  const int n = y.size(), N = n_particles, K = probs.size();
  if (p.size() != n) stop("length(p) must equal length(y)");
  if (!(std::fabs(phi) < 1.0) || !(sigma > 0.0)) stop("invalid phi or sigma");
  const int M = std::min(n_quantile_particles, N);
  const double s2_xi = sigma_xi * sigma_xi;
  const double LOG2PI = std::log(2.0 * M_PI);

  std::vector<double> h(N), hn(N), lw(N), w(N), cw(N);
  NumericVector logscore(n), pit(n), h_filt(n);
  NumericMatrix q(n, K);

  // initial particles from the stationary distribution of h
  const double sd0 = sigma / std::sqrt(1.0 - phi * phi);
  for (int i = 0; i < N; ++i) h[i] = mu + sd0 * R::norm_rand();

  for (int t = 0; t < n; ++t) {
    // propagate (h currently holds draws from p(h_{t-1} | y_1:t-1), except at t = 0)
    if (t > 0) for (int i = 0; i < N; ++i) h[i] = mu + phi * (h[i] - mu) + sigma * R::norm_rand();
    const double pt = p[t], yt = y[t] - m;          // centred return

    // predictive CDF at y_t, and predictive quantiles from a subset of particles
    double cdf = 0.0;
    for (int i = 0; i < N; ++i) cdf += mix_cdf(yt, pt, mu_xi, s2_xi, std::exp(h[i]));
    pit[t] = cdf / N;

    for (int k = 0; k < K; ++k) {
      double lo = -50.0, hi = 50.0;
      for (int it = 0; it < 60; ++it) {
        double mid = 0.5 * (lo + hi), c = 0.0;
        for (int i = 0; i < M; ++i) c += mix_cdf(mid, pt, mu_xi, s2_xi, std::exp(h[i]));
        if (c / M < probs[k]) lo = mid; else hi = mid;
        if (hi - lo < 1e-6) break;
      }
      q(t, k) = m + 0.5 * (lo + hi);
    }

    // weights = predictive density of y_t for each particle (log scale)
    double mx = -INFINITY;
    for (int i = 0; i < N; ++i) {
      double eh = std::exp(h[i]);
      double l0 = -0.5 * (LOG2PI + h[i] + yt * yt / eh);
      double lwi;
      if (pt > 0.0) {
        double v1 = s2_xi + eh, d1 = yt - mu_xi;
        double l1 = -0.5 * (LOG2PI + std::log(v1) + d1 * d1 / v1);
        double a = std::log(pt) + l1, b = std::log1p(-pt) + l0;
        double m = std::max(a, b);
        lwi = m + std::log(std::exp(a - m) + std::exp(b - m));
      } else {
        lwi = l0;
      }
      lw[i] = lwi;
      if (lwi > mx) mx = lwi;
    }
    double sw = 0.0;
    for (int i = 0; i < N; ++i) { w[i] = std::exp(lw[i] - mx); sw += w[i]; }
    logscore[t] = mx + std::log(sw / N);

    double hm = 0.0;
    for (int i = 0; i < N; ++i) { w[i] /= sw; hm += w[i] * h[i]; }
    h_filt[t] = hm;

    // systematic resampling
    cw[0] = w[0];
    for (int i = 1; i < N; ++i) cw[i] = cw[i - 1] + w[i];
    double u = R::unif_rand() / N;
    int j = 0;
    for (int i = 0; i < N; ++i) {
      double target = u + (double)i / N;
      while (j < N - 1 && cw[j] < target) ++j;
      hn[i] = h[j];
    }
    h.swap(hn);
  }
  return List::create(_["logscore"] = logscore, _["pit"] = pit, _["q"] = q, _["h_filt"] = h_filt);
}
