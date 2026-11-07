# =====================================================================
# eval_functions.R
# Helpers for model evaluation:
#   - GARCH(1,1) with standardised Student-t errors (benchmark), by ML
#   - Value-at-Risk backtests (Kupiec; Christoffersen)
#   - Diebold-Mariano-type test for differences in log scores (HAC)
# =====================================================================

# ---- GARCH(1,1)-t -------------------------------------------------------------
# y_t = m + s_t e_t,  e_t ~ standardised t(nu),  s_t^2 = w + a (y_{t-1} - m)^2 + b s_{t-1}^2
garch_t_filter <- function(par, y) {
  m <- par[1]; w <- par[2]; a <- par[3]; b <- par[4]
  n <- length(y); s2 <- numeric(n)
  s2[1] <- var(y)
  for (t in 2:n) s2[t] <- w + a * (y[t - 1] - m)^2 + b * s2[t - 1]
  s2
}
dstd_t_log <- function(e, nu) {          # log density of a standardised t (variance 1)
  sc <- sqrt((nu - 2) / nu)
  dt(e / sc, df = nu, log = TRUE) - log(sc)
}
fit_garch_t <- function(y) {
  # unconstrained parameters: m, log w, logit-type for (a, b) with a + b < 1, log(nu - 2)
  unpack <- function(th) {
    ab <- exp(th[3:4]) / (1 + sum(exp(th[3:4])))      # a, b > 0, a + b < 1
    c(m = th[1], w = exp(th[2]), a = ab[1], b = ab[2], nu = 2 + exp(th[5]))
  }
  nll <- function(th) {
    p <- unpack(th)
    s2 <- garch_t_filter(p[1:4], y)
    -sum(dstd_t_log((y - p["m"]) / sqrt(s2), p["nu"]) - 0.5 * log(s2))
  }
  starts <- list(c(0, log(0.01 * var(y)), log(0.08 / 0.02), log(0.90 / 0.02), log(4)),
                 c(0, log(0.05 * var(y)), log(0.15 / 0.05), log(0.80 / 0.05), log(2)))
  best <- NULL
  for (s in starts) {
    o <- tryCatch(optim(s, nll, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12)),
                  error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) stop("GARCH-t estimation failed")
  p <- unpack(best$par)
  list(par = p, loglik = -best$value, convergence = best$convergence)
}
# One-step-ahead predictive quantities (log score, PIT, quantiles) for GARCH-t
predict_garch_t <- function(fit, y, probs) {
  p <- fit$par
  s <- sqrt(garch_t_filter(p[1:4], y))
  sc <- sqrt((p["nu"] - 2) / p["nu"])
  e <- (y - p["m"]) / s
  list(logscore = dstd_t_log(e, p["nu"]) - log(s),
       pit = pt(e / sc, df = p["nu"]),
       q = sapply(probs, function(a) p["m"] + s * sc * qt(a, df = p["nu"])))
}

# ---- VaR backtests --------------------------------------------------------------
# hits: 0/1 exceedance indicators; alpha: nominal exceedance probability
kupiec_test <- function(hits, alpha) {
  n <- length(hits); x <- sum(hits); pi_hat <- x / n
  ll0 <- (n - x) * log(1 - alpha) + x * log(alpha)
  ll1 <- if (x == 0 || x == n) 0 else (n - x) * log(1 - pi_hat) + x * log(pi_hat)
  lr <- -2 * (ll0 - ll1)
  c(rate = pi_hat, LR_uc = lr, p_uc = pchisq(lr, 1, lower.tail = FALSE))
}
christoffersen_test <- function(hits, alpha) {
  h0 <- hits[-length(hits)]; h1 <- hits[-1]
  n00 <- sum(h0 == 0 & h1 == 0); n01 <- sum(h0 == 0 & h1 == 1)
  n10 <- sum(h0 == 1 & h1 == 0); n11 <- sum(h0 == 1 & h1 == 1)
  p01 <- n01 / max(n00 + n01, 1); p11 <- n11 / max(n10 + n11, 1)
  p <- (n01 + n11) / max(n00 + n01 + n10 + n11, 1)
  xlx <- function(k, q) if (k == 0) 0 else k * log(q)
  ll_ind <- xlx(n00, 1 - p01) + xlx(n01, p01) + xlx(n10, 1 - p11) + xlx(n11, p11)
  ll_0   <- xlx(n00 + n10, 1 - p) + xlx(n01 + n11, p)
  lr_ind <- -2 * (ll_0 - ll_ind)
  lr_uc <- kupiec_test(hits, alpha)["LR_uc"]
  lr_cc <- unname(lr_uc + lr_ind)
  c(LR_ind = lr_ind, p_ind = pchisq(lr_ind, 1, lower.tail = FALSE),
    LR_cc = lr_cc, p_cc = pchisq(lr_cc, 2, lower.tail = FALSE))
}

# ---- Diebold-Mariano-type test on log-score differences ---------------------------
# d_t = logscore_A - logscore_B ; HAC (Newey-West, Bartlett) standard error
dm_test <- function(d, lag = floor(4 * (length(d) / 100)^(2 / 9))) {
  n <- length(d); dc <- d - mean(d)
  g0 <- sum(dc^2) / n
  lrv <- g0
  for (k in seq_len(lag)) lrv <- lrv + 2 * (1 - k / (lag + 1)) * sum(dc[-(1:k)] * dc[1:(n - k)]) / n
  stat <- mean(d) / sqrt(lrv / n)
  c(mean_diff = mean(d), total_diff = sum(d), z = stat, p = 2 * pnorm(-abs(stat)))
}
