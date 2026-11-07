# =====================================================================
# 09_simulation_study.R
# Simulation study (ADEMP; Morris, White and Crowther 2019).
#
# Aims
#   1. bias and interval coverage of the SVJ-X estimates
#   2. false detection of jump triggers when jumps are not covariate-driven
#   3. spurious jumps when returns are fat-tailed but have no jumps
#   4. whether the log predictive score picks the data-generating model
# Data-generating mechanisms (real covariate matrix, n = 5086 days;
# parameters calibrated to the posterior medians of the full-sample fit)
#   A  SVJ-X, calibrated (dollar effect 0.666)
#   B  SVJ-X, weaker dollar effect (0.30)
#   C  SVJ: constant jump probability equal to the average of A, no triggers
#   D  SV with Student-t(5) errors, no jumps
# Methods: SVJ-X and SVJ fitted to every data set (12,000 iterations,
#          2,000 burn-in, thinning 4), log scores from the particle filter
# Performance: bias, RMSE, 95% coverage, rejection rate of "gamma = 0",
#          number of detected jumps, share of data sets where SVJ-X wins
#
# Every replicate is saved in sim/rep_<scenario>_<r>.rds, so the script can
# be stopped and started again; finished replicates are not recomputed.
# Output: output/sim_replicates.csv, output/table7_simulation.csv,
#         figures/fig9_sim_estimates, figures/fig10_sim_coverage
# Run time: several hours (about 9 CPU-minutes per replicate).
# =====================================================================

N_REP   <- 50            # replicates per scenario
N_ITER  <- 12000; BURN <- 2000; THIN <- 4
N_PART  <- 5000

for (p in c("Rcpp", "ggplot2", "patchwork", "scales")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(scales) })
if (basename(getwd()) != "rupiah-jump-sv" && dir.exists("~/rupiah-jump-sv")) setwd("~/rupiah-jump-sv")
need <- c("svjx_sampler.cpp", "svjx_functions.R", "svjx_filter.cpp", file.path("output", "model_data.rds"))
for (f in need) if (!file.exists(f)) stop("File ", f, " not found. Run 06_fit_models.R first.")
dir.create("sim", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

md <- readRDS(file.path("output", "model_data.rds"))
X <- md$X
n <- nrow(X)

# ---- Data-generating mechanisms ---------------------------------------------------
base <- list(m = 0.014, mu = -2.405, phi = 0.952, sigma = 0.383, mu_xi = 0.091, sigma_xi = 0.817,
             gamma = c(-2.554, 0.134, 0.250, 0.666, 0.017))
names(base$gamma) <- colnames(X)
gB <- base$gamma; gB["abs_d_usd"] <- 0.30
p_avg <- mean(pnorm(drop(X %*% base$gamma)))
gC <- c(qnorm(p_avg), rep(0, ncol(X) - 1)); names(gC) <- colnames(X)
scen <- list(A = list(type = "svjx", gamma = base$gamma),
             B = list(type = "svjx", gamma = gB),
             C = list(type = "svjx", gamma = gC),
             D = list(type = "svt", nu = 5))

simulate_one <- function(sc, base, X, seed) {
  set.seed(seed)
  n <- nrow(X)
  if (sc$type == "svjx") {
    s <- sim_svjx(n, base$mu, base$phi, base$sigma, sc$gamma, base$mu_xi, base$sigma_xi, X = X, m = base$m)
    return(list(y = s$y, n_jumps = sum(s$J)))
  }
  h <- numeric(n)
  h[1] <- rnorm(1, base$mu, base$sigma / sqrt(1 - base$phi^2))
  for (t in 2:n) h[t] <- base$mu + base$phi * (h[t - 1] - base$mu) + base$sigma * rnorm(1)
  e <- rt(n, sc$nu) * sqrt((sc$nu - 2) / sc$nu)            # unit-variance t errors
  list(y = base$m + exp(h / 2) * e, n_jumps = 0L)
}

run_rep <- function(task, scen, base, Xmat, settings, wd, simulate_one) {
  setwd(wd)
  out_file <- file.path("sim", sprintf("rep_%s_%03d.rds", task$s, task$r))
  if (file.exists(out_file)) return(out_file)
  if (!exists("draw_h_tridiag")) Rcpp::sourceCpp("svjx_sampler.cpp")
  if (!exists("pf_svjx")) Rcpp::sourceCpp("svjx_filter.cpp")
  source("svjx_functions.R")
  seed <- 1000 * match(task$s, names(scen)) + task$r
  d <- simulate_one(scen[[task$s]], base, Xmat, seed)
  res <- list(scenario = task$s, rep = task$r, seed = seed, true_jumps = d$n_jumps)
  pm <- function(f, nm) if (nm %in% colnames(f$draws)) mean(f$draws[, nm]) else 0
  for (mod in c("SVJX", "SVJ")) {
    f <- fit_svjx(d$y, if (mod == "SVJX") Xmat else NULL, model = mod, n_iter = settings$n_iter,
                  burn = settings$burn, thin = settings$thin, seed = seed + 7, verbose = FALSE)
    q <- apply(f$draws, 2, quantile, c(0.025, 0.5, 0.975))
    g <- f$draws[, grep("^gamma_", colnames(f$draws)), drop = FALSE]
    pt <- colMeans(pnorm(g %*% t(if (mod == "SVJX") Xmat else matrix(1, nrow(Xmat), 1))))
    set.seed(seed + 11)
    pf <- pf_svjx(d$y, pt, pm(f, "m"), pm(f, "mu"), pm(f, "phi"), pm(f, "sigma"), pm(f, "mu_xi"),
                  max(pm(f, "sigma_xi"), 1e-6), settings$n_part, 0.5, 50)
    res[[mod]] <- list(q = q, n_jumps = mean(f$n_jumps), logscore = sum(pf$logscore), seconds = f$seconds)
  }
  saveRDS(res, out_file)
  out_file
}

# ---- Run (resumable, parallel) -------------------------------------------------------
tasks <- unlist(lapply(names(scen), function(s) lapply(seq_len(N_REP), function(r) list(s = s, r = r))),
                recursive = FALSE)
todo <- Filter(function(tk) !file.exists(file.path("sim", sprintf("rep_%s_%03d.rds", tk$s, tk$r))), tasks)
message(length(tasks) - length(todo), " of ", length(tasks), " replicates already done; ",
        length(todo), " to run.")
settings <- list(n_iter = N_ITER, burn = BURN, thin = THIN, n_part = N_PART)
if (length(todo)) {
  n_workers <- max(1, min(length(todo), parallel::detectCores(logical = FALSE), na.rm = TRUE))
  message("Running on ", n_workers, " core(s). Progress is saved after every replicate.")
  t0 <- proc.time()[3]
  if (n_workers > 1) {
    cl <- parallel::makeCluster(n_workers)
    invisible(tryCatch(parallel::parLapplyLB(cl, todo, run_rep, scen = scen, base = base, Xmat = X,
                                             simulate_one = simulate_one,
                                             settings = settings, wd = getwd()),
                       finally = parallel::stopCluster(cl)))
  } else {
    invisible(lapply(todo, run_rep, scen = scen, base = base, Xmat = X, simulate_one = simulate_one, settings = settings, wd = getwd()))
  }
  message("Finished in ", round((proc.time()[3] - t0) / 3600, 2), " hours.")
}

# ---- Collect --------------------------------------------------------------------------
files <- list.files("sim", pattern = "^rep_[A-D]_[0-9]{3}\\.rds$", full.names = TRUE)
reps <- lapply(files, readRDS)
reps <- Filter(function(z) z$scenario %in% names(scen) && z$rep <= N_REP, reps)
truth_of <- function(s) {
  sc <- scen[[s]]
  if (sc$type != "svjx") return(NULL)
  tv <- c(m = base$m, mu = base$mu, phi = base$phi, sigma = base$sigma, mu_xi = base$mu_xi,
          sigma_xi = base$sigma_xi)
  c(tv, setNames(sc$gamma, paste0("gamma_", names(sc$gamma))))
}
long <- do.call(rbind, lapply(reps, function(z) {
  q <- z$SVJX$q; tv <- truth_of(z$scenario)
  data.frame(scenario = z$scenario, rep = z$rep, parameter = colnames(q),
             lower = q[1, ], median = q[2, ], upper = q[3, ],
             true = if (is.null(tv)) NA_real_ else unname(tv[colnames(q)]), row.names = NULL)
}))
per_rep <- do.call(rbind, lapply(reps, function(z) data.frame(
  scenario = z$scenario, rep = z$rep, true_jumps = z$true_jumps,
  jumps_svjx = z$SVJX$n_jumps, jumps_svj = z$SVJ$n_jumps,
  logscore_diff = z$SVJX$logscore - z$SVJ$logscore,
  usd_excludes_0 = z$SVJX$q[1, "gamma_abs_d_usd"] > 0 | z$SVJX$q[3, "gamma_abs_d_usd"] < 0,
  any_slope_excludes_0 = any(apply(z$SVJX$q[, grep("^gamma_(?!\\(Intercept\\))", colnames(z$SVJX$q), perl = TRUE),
                                            drop = FALSE], 2, function(v) v[1] > 0 | v[3] < 0)))))
write.csv(merge(long, per_rep[, c("scenario", "rep", "logscore_diff")]),
          file.path("output", "sim_replicates.csv"), row.names = FALSE)

# Performance measures with Monte Carlo standard errors
perf <- do.call(rbind, lapply(split(long, list(long$scenario, long$parameter), drop = TRUE), function(z) {
  k <- nrow(z)
  cov <- if (all(is.na(z$true))) NA else mean(z$lower <= z$true & z$true <= z$upper)
  rej <- mean(z$lower > 0 | z$upper < 0)
  data.frame(scenario = z$scenario[1], parameter = z$parameter[1], n_rep = k, true = z$true[1],
             mean_est = mean(z$median), bias = mean(z$median - z$true),
             bias_mcse = sd(z$median) / sqrt(k), rmse = sqrt(mean((z$median - z$true)^2)),
             coverage = cov, coverage_mcse = if (is.na(cov)) NA else sqrt(cov * (1 - cov) / k),
             reject_zero = rej, row.names = NULL)
}))
perf <- perf[order(perf$scenario, perf$parameter), ]
num <- sapply(perf, is.numeric); perf[num] <- lapply(perf[num], function(v) signif(v, 4))
write.csv(perf, file.path("output", "table7_simulation.csv"), row.names = FALSE)

agg <- do.call(rbind, lapply(split(per_rep, per_rep$scenario), function(z) data.frame(
  scenario = z$scenario[1], n_rep = nrow(z), true_jumps = mean(z$true_jumps),
  jumps_svjx = mean(z$jumps_svjx), jumps_svj = mean(z$jumps_svj),
  svjx_wins = mean(z$logscore_diff > 0), mean_logscore_diff = mean(z$logscore_diff),
  usd_detected = mean(z$usd_excludes_0), any_trigger_detected = mean(z$any_slope_excludes_0))))
num <- sapply(agg, is.numeric); agg[num] <- lapply(agg[num], function(v) signif(v, 4))

# ---- Figures ---------------------------------------------------------------------------
COL <- c(A = "#1F3A5F", B = "#4F7CAC", C = "#E0A458", D = "#C8102E")
theme_paper <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        plot.tag = element_text(face = "bold"), legend.position = "none", strip.text = element_text(face = "bold"))
slopes <- c("gamma_log_vix", "gamma_abs_d_ust10", "gamma_abs_d_usd", "gamma_long_break")
lab <- c(gamma_log_vix = "VIX", gamma_abs_d_ust10 = "US 10y yield", gamma_abs_d_usd = "US dollar",
         gamma_long_break = "Long break")
est <- long[long$parameter %in% slopes, ]
est$parameter <- factor(lab[est$parameter], levels = lab)
tr <- unique(est[!is.na(est$true), c("scenario", "parameter", "true")])
p9 <- ggplot(est, aes(scenario, median, colour = scenario, fill = scenario)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_violin(alpha = 0.18, colour = NA, scale = "width", width = 0.8) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.5, linewidth = 0.3) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), size = 0.6, alpha = 0.6) +
  geom_point(data = tr, aes(y = true), shape = 23, size = 2.4, fill = "white", colour = "black") +
  facet_wrap(~ parameter, nrow = 1) + scale_colour_manual(values = COL) + scale_fill_manual(values = COL) +
  labs(x = "Scenario", y = "Posterior median of the trigger effect\n(diamond: true value)") + theme_paper
ggsave(file.path("figures", "fig9_sim_estimates.png"), p9, width = 8, height = 3.4, dpi = 300, bg = "white")
ggsave(file.path("figures", "fig9_sim_estimates.pdf"), p9, width = 8, height = 3.4)

cv <- perf[perf$scenario %in% c("A", "B", "C") & !is.na(perf$coverage), ]
cv$lo <- pmax(0, cv$coverage - 1.96 * cv$coverage_mcse); cv$hi <- pmin(1, cv$coverage + 1.96 * cv$coverage_mcse)
cv$parameter <- factor(cv$parameter, levels = rev(unique(cv$parameter[order(cv$parameter)])))
p10 <- ggplot(cv, aes(coverage, parameter, colour = scenario)) +
  geom_vline(xintercept = 0.95, linetype = 2, colour = "grey30", linewidth = 0.3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.25) +
  scale_colour_manual(values = COL, name = "Scenario") +
  scale_x_continuous(labels = percent) + coord_cartesian(xlim = c(min(0.5, min(cv$lo, na.rm = TRUE)), 1)) +
  labs(x = "Coverage of the 95% credible interval", y = NULL) + theme_paper + theme(legend.position = "bottom")
ggsave(file.path("figures", "fig10_sim_coverage.png"), p10, width = 6, height = 5, dpi = 300, bg = "white")
ggsave(file.path("figures", "fig10_sim_coverage.pdf"), p10, width = 6, height = 5)

# ---- Summary ---------------------------------------------------------------------------
cat("\n================ SIMULATION STUDY ================\n")
cat("Replicates found:", paste(names(table(per_rep$scenario)), table(per_rep$scenario), sep = " = ", collapse = ", "), "\n\n")
cat("Per scenario (jumps; share of data sets where SVJ-X has the better log score;\n",
    "share where the dollar effect / any trigger is declared non-zero):\n")
print(agg, row.names = FALSE)
cat("\nTrigger effects in SVJ-X (bias, coverage, rejection of zero):\n")
print(perf[perf$parameter %in% slopes, c("scenario", "parameter", "true", "mean_est", "bias",
                                          "bias_mcse", "rmse", "coverage", "reject_zero")], row.names = FALSE)
cat("\nFull table: output/table7_simulation.csv\n")
