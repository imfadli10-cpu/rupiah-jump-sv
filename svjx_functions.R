# =====================================================================
# svjx_functions.R
# Simulation and Gibbs sampling for the stochastic volatility model with
# covariate-driven jumps (SVJ-X). SV (no jumps) and SVJ (constant jump
# probability) are special cases.
#
#   y_t   = m + J_t * xi_t + exp(h_t / 2) * eps_t,        eps_t ~ N(0, 1)
#           (m: constant daily drift)
#   h_t   = mu + phi * (h_{t-1} - mu) + sigma * eta_t,    eta_t ~ N(0, 1)
#   h_1   ~ N(mu, sigma^2 / (1 - phi^2))
#   J_t   ~ Bernoulli(p_t),  p_t = Phi(x_t' gamma)        (probit link)
#   xi_t  ~ N(mu_xi, sigma_xi^2)                          (jump size)
#
# Sampler blocks (one Gibbs sweep):
#   1. (J_t, xi_t) | rest, with xi_t integrated out when drawing J_t
#   2. mixture indicators s_t of the 10-component approximation to
#      log chi^2_1 (Omori et al. 2007), then h | rest in one block
#      (tridiagonal precision, C++ routine draw_h_tridiag)
#   3. mu, phi (Metropolis-Hastings, independence proposal), sigma^2;
#      drift m (conjugate normal)
#   4. mu_xi, sigma_xi^2 (conjugate)
#   5. gamma | J by independence Metropolis-Hastings with a multivariate-t
#      proposal at the posterior mode (Laplace approximation of the probit
#      posterior); acceptance is usually above 90%
#
# Returns in % units if y is given in %.
# =====================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

if (!exists("draw_h_tridiag")) {
  stop("Compile the C++ part first: Rcpp::sourceCpp(\"svjx_sampler.cpp\")")
}

# 10-component mixture approximation to log(chi^2_1), Omori et al. (2007)
OMORI <- list(
  p = c(0.00609, 0.04775, 0.13057, 0.20674, 0.22715, 0.18842, 0.12047, 0.05591, 0.01575, 0.00115),
  m = c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000),
  v = c(0.11265, 0.17788, 0.26768, 0.40611, 0.62699, 0.98583, 1.57469, 2.54498, 4.16591, 7.33342))

default_priors <- function(k) {
  list(m_mean = 0, m_sd = 1,                   # drift m ~ N(0, 1) (% per day)
       mu_mean = 0, mu_sd = 10,               # mu ~ N(0, 10^2)
       phi_a = 20, phi_b = 1.5,               # (phi + 1) / 2 ~ Beta(20, 1.5)
       s2_shape = 2.5, s2_rate = 0.025,       # sigma^2 ~ InvGamma(2.5, 0.025)
       muxi_mean = 0, muxi_sd = 5,            # mu_xi ~ N(0, 5^2)
       sxi_shape = 3, sxi_rate = 8,           # sigma_xi^2 ~ InvGamma(3, 8): prior mode 2 (%^2)
       gamma_mean = c(-2.5, rep(0, k - 1)),   # probit intercept: jumps rare a priori
       gamma_sd = c(1, rep(1, k - 1)))
}

# ---- Simulation ---------------------------------------------------------
sim_svjx <- function(n, mu, phi, sigma, gamma, mu_xi, sigma_xi, X = NULL, jumps = TRUE, m = 0) {
  if (is.null(X)) X <- matrix(1, n, 1)
  stopifnot(nrow(X) == n, ncol(X) == length(gamma))
  h <- numeric(n)
  h[1] <- rnorm(1, mu, sigma / sqrt(1 - phi^2))
  for (t in 2:n) h[t] <- mu + phi * (h[t - 1] - mu) + sigma * rnorm(1)
  p <- pnorm(drop(X %*% gamma))
  J <- if (jumps) rbinom(n, 1, p) else rep(0L, n)
  xi <- rnorm(n, mu_xi, sigma_xi)
  y <- m + J * xi + exp(h / 2) * rnorm(n)
  list(y = y, h = h, J = J, xi = xi, p = p, X = X)
}

# ---- Helpers --------------------------------------------------------------
rinvgamma <- function(shape, rate) 1 / rgamma(1, shape = shape, rate = rate)

# Truncated standard normal draws for probit augmentation, stable in the tails.
# For J = 1: w > -eta ; for J = 0: w <= -eta ; returns z = eta + w.
draw_z_probit <- function(eta, J) {
  n <- length(eta)
  z <- numeric(n)
  i1 <- J == 1
  if (any(i1)) {
    pu <- pnorm(eta[i1])                       # P(w > -eta)
    u <- pmax(runif(sum(i1)) * pu, 1e-300)
    z[i1] <- eta[i1] + qnorm(u, lower.tail = FALSE)
  }
  i0 <- !i1
  if (any(i0)) {
    pl <- pnorm(-eta[i0])                      # P(w <= -eta)
    u <- pmax(runif(sum(i0)) * pl, 1e-300)
    z[i0] <- eta[i0] + qnorm(u)
  }
  z
}

# ---- Probit regression helpers (jump probability) -----------------------------
# Log posterior of gamma given jump indicators J, with N(b0, B) prior (B_inv given).
probit_logpost <- function(g, J, X, b0, B_inv) {
  eta <- drop(X %*% g)
  i1 <- J == 1
  sum(pnorm(eta[i1], log.p = TRUE)) + sum(pnorm(eta[!i1], lower.tail = FALSE, log.p = TRUE)) -
    0.5 * drop(crossprod(g - b0, B_inv %*% (g - b0)))
}
# Posterior mode and Cholesky factor of the inverse negative Hessian (Newton steps).
probit_laplace <- function(J, X, b0, B_inv, start, max_it = 25, tol = 1e-8) {
  g <- start
  for (i in seq_len(max_it)) {
    eta <- drop(X %*% g)
    l1 <- exp(dnorm(eta, log = TRUE) - pnorm(eta, log.p = TRUE))                     # phi / Phi(eta)
    l0 <- exp(dnorm(eta, log = TRUE) - pnorm(eta, lower.tail = FALSE, log.p = TRUE)) # phi / Phi(-eta)
    score <- ifelse(J == 1, l1, -l0)
    w <- ifelse(J == 1, l1 * (l1 + eta), l0 * (l0 - eta))                            # > 0
    grad <- crossprod(X, score) - B_inv %*% (g - b0)
    H <- crossprod(X * w, X) + B_inv                                                 # negative Hessian
    step <- solve(H, grad)
    g_new <- drop(g + step)
    if (probit_logpost(g_new, J, X, b0, B_inv) < probit_logpost(g, J, X, b0, B_inv)) {
      g_new <- drop(g + 0.5 * step)                                                  # damped step
    }
    done <- max(abs(g_new - g)) < tol
    g <- g_new
    if (done) break
  }
  eta <- drop(X %*% g)
  l1 <- exp(dnorm(eta, log = TRUE) - pnorm(eta, log.p = TRUE))
  l0 <- exp(dnorm(eta, log = TRUE) - pnorm(eta, lower.tail = FALSE, log.p = TRUE))
  w <- ifelse(J == 1, l1 * (l1 + eta), l0 * (l0 - eta))
  H <- crossprod(X * w, X) + B_inv
  list(mode = g, cov_chol = t(chol(chol2inv(chol(H)))))   # lower-triangular L, cov = L L'
}
# Multivariate t: one draw, and log density up to a constant that cancels in MH.
rmvt1 <- function(m, L, df) drop(m + L %*% rnorm(length(m)) / sqrt(rchisq(1, df) / df))
dmvt_log <- function(x, m, L, df) {
  u <- forwardsolve(L, x - m)
  -0.5 * (df + length(m)) * log1p(sum(u^2) / df)
}

# ---- Gibbs sampler ----------------------------------------------------------
# model: "SV", "SVJ" (X = intercept only) or "SVJX" (X supplied)
fit_svjx <- function(y, X = NULL, model = c("SVJX", "SVJ", "SV"),
                     n_iter = 20000, burn = 5000, thin = 5,
                     priors = NULL, offset = 1e-4, seed = 1, h_draws_keep = 500,
                     verbose = TRUE) {
  model <- match.arg(model)
  set.seed(seed)
  n <- length(y)
  if (anyNA(y)) stop("y contains missing values")
  if (model == "SVJX") {
    if (is.null(X)) stop("model SVJX needs a covariate matrix X (first column = 1)")
  } else {
    X <- matrix(1, n, 1, dimnames = list(NULL, "(Intercept)"))
  }
  X <- as.matrix(X)
  if (nrow(X) != n) stop("nrow(X) must equal length(y)")
  if (anyNA(X)) stop("X contains missing values")
  if (!all(X[, 1] == 1)) stop("first column of X must be the intercept (all ones)")
  k <- ncol(X)
  pr <- default_priors(k)
  if (!is.null(priors)) pr[names(priors)] <- priors
  jumps_on <- model != "SV"

  # Initial values
  h <- rep(log(var(y)), n)
  mu <- log(var(y)); phi <- 0.95; s2 <- 0.05; m <- mean(y)
  J <- integer(n); xi <- numeric(n)
  mu_xi <- 0; s2_xi <- 4
  gamma <- pr$gamma_mean
  B_inv <- diag(1 / pr$gamma_sd^2, k)

  if (burn >= n_iter) stop("burn must be smaller than n_iter")
  keep <- seq(burn + thin, n_iter, by = thin)
  nk <- length(keep)
  is_keep <- logical(n_iter); is_keep[keep] <- TRUE
  par_names <- c("m", "mu", "phi", "sigma", if (jumps_on) c("mu_xi", "sigma_xi", paste0("gamma_", colnames(X) %||% seq_len(k))))
  draws <- matrix(NA_real_, nk, length(par_names), dimnames = list(NULL, par_names))
  h_keep_idx <- unique(round(seq(1, nk, length.out = min(h_draws_keep, nk))))
  h_draws <- matrix(NA_real_, length(h_keep_idx), n)
  pj_sum <- numeric(n); h_sum <- numeric(n); nj_draws <- integer(nk)
  acc_phi <- 0; acc_gamma <- 0
  g_mode <- pr$gamma_mean
  j <- 0; jh <- 0
  t0 <- proc.time()[3]

  for (it in seq_len(n_iter)) {
    # 1. Jumps (J, xi) given h, jump parameters, gamma
    if (jumps_on) {
      vh <- exp(h)
      eta <- drop(X %*% gamma)
      yc <- y - m
      lp1 <- pnorm(eta, log.p = TRUE) + dnorm(yc, mu_xi, sqrt(s2_xi + vh), log = TRUE)
      lp0 <- pnorm(eta, lower.tail = FALSE, log.p = TRUE) + dnorm(yc, 0, sqrt(vh), log = TRUE)
      pj <- 1 / (1 + exp(lp0 - lp1))
      J <- as.integer(runif(n) < pj)
      prec <- 1 / s2_xi + 1 / vh
      mxi <- (mu_xi / s2_xi + yc / vh) / prec
      xi <- ifelse(J == 1, rnorm(n, mxi, sqrt(1 / prec)), 0)
    } else {
      pj <- numeric(n)
    }

    # 2. Volatility path h given jump-adjusted returns
    ys <- y - m - J * xi
    yt <- log(ys^2 + offset)
    # mixture indicators
    lw <- sapply(1:10, function(c) log(OMORI$p[c]) + dnorm(yt, h + OMORI$m[c], sqrt(OMORI$v[c]), log = TRUE))
    lw <- lw - apply(lw, 1, max)
    w <- exp(lw); w <- w / rowSums(w)
    cw <- w %*% upper.tri(diag(10), diag = TRUE)          # row-wise cumulative sums
    s <- 1L + rowSums(cw < runif(n))
    s[s > 10L] <- 10L
    h <- draw_h_tridiag(yt, OMORI$m[s], OMORI$v[s], mu, phi, s2, rnorm(n))

    # 3. mu, phi, sigma^2 given h
    hc <- h - mu
    x0 <- hc[-n]; x1 <- hc[-1]
    sxx <- sum(x0^2)
    phi_hat <- sum(x0 * x1) / sxx
    phi_sd <- sqrt(s2 / sxx)
    repeat { phi_prop <- rnorm(1, phi_hat, phi_sd); if (abs(phi_prop) < 1) break }
    log_target <- function(ph) {
      dbeta((ph + 1) / 2, pr$phi_a, pr$phi_b, log = TRUE) +
        dnorm(hc[1], 0, sqrt(s2 / (1 - ph^2)), log = TRUE)
    }
    if (log(runif(1)) < log_target(phi_prop) - log_target(phi)) { phi <- phi_prop; acc_phi <- acc_phi + 1 }
    res <- c(sqrt(1 - phi^2) * hc[1], x1 - phi * x0)
    s2 <- rinvgamma(pr$s2_shape + n / 2, pr$s2_rate + sum(res^2) / 2)
    prec_mu <- 1 / pr$mu_sd^2 + (1 - phi^2) / s2 + (n - 1) * (1 - phi)^2 / s2
    num_mu <- pr$mu_mean / pr$mu_sd^2 + (1 - phi^2) * h[1] / s2 +
      (1 - phi) * sum(h[-1] - phi * h[-n]) / s2
    mu <- rnorm(1, num_mu / prec_mu, sqrt(1 / prec_mu))

    # 3b. Drift m given h and jumps (conjugate normal)
    iv <- exp(-h)
    prec_m <- 1 / pr$m_sd^2 + sum(iv)
    m <- rnorm(1, (pr$m_mean / pr$m_sd^2 + sum((y - J * xi) * iv)) / prec_m, sqrt(1 / prec_m))

    # 4-5. Jump size distribution and jump probability coefficients
    if (jumps_on) {
      nj <- sum(J)
      xj <- xi[J == 1]
      prec_m <- 1 / pr$muxi_sd^2 + nj / s2_xi
      mu_xi <- rnorm(1, (pr$muxi_mean / pr$muxi_sd^2 + sum(xj) / s2_xi) / prec_m, sqrt(1 / prec_m))
      s2_xi <- rinvgamma(pr$sxi_shape + nj / 2, pr$sxi_rate + sum((xj - mu_xi)^2) / 2)
      # gamma | J: independence Metropolis-Hastings with a multivariate-t
      # proposal centred at the posterior mode (Laplace approximation).
      lap <- probit_laplace(J, X, pr$gamma_mean, B_inv, start = g_mode)
      g_mode <- lap$mode                         # warm start for the next sweep
      g_prop <- rmvt1(lap$mode, lap$cov_chol, df = 6)
      log_r <- probit_logpost(g_prop, J, X, pr$gamma_mean, B_inv) -
        probit_logpost(gamma, J, X, pr$gamma_mean, B_inv) +
        dmvt_log(gamma, lap$mode, lap$cov_chol, df = 6) -
        dmvt_log(g_prop, lap$mode, lap$cov_chol, df = 6)
      if (log(runif(1)) < log_r) { gamma <- g_prop; acc_gamma <- acc_gamma + 1 }
    }

    # Store
    if (is_keep[it]) {
      j <- j + 1
      draws[j, ] <- c(m, mu, phi, sqrt(s2), if (jumps_on) c(mu_xi, sqrt(s2_xi), gamma))
      pj_sum <- pj_sum + pj
      h_sum <- h_sum + h
      nj_draws[j] <- sum(J)
      if (j %in% h_keep_idx) { jh <- jh + 1; h_draws[jh, ] <- h }
    }
    if (verbose && it %% 2000 == 0) {
      message(sprintf("  %s: iteration %d / %d (%.0f s)", model, it, n_iter, proc.time()[3] - t0))
    }
  }

  list(model = model, draws = draws, h_mean = h_sum / nk, h_draws = h_draws,
       jump_prob = pj_sum / nk, n_jumps = nj_draws, acc_phi = acc_phi / n_iter,
       acc_gamma = if (jumps_on) acc_gamma / n_iter else NA,
       X = X, y = y, priors = pr, offset = offset,
       settings = list(n_iter = n_iter, burn = burn, thin = thin, seed = seed),
       seconds = unname(proc.time()[3] - t0))
}
