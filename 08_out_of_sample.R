# =====================================================================
# 08_out_of_sample.R
# Genuine out-of-sample evaluation.
#
#   Training period : 2006 - 2018  (all model parameters estimated here;
#                     covariates standardised with training means/sds)
#   Test period     : 2019 - end of data (never used for estimation)
#
# The particle filter runs through the whole sample with the training
# parameters, so each test-day forecast uses only returns up to the
# previous day and parameters from before 2019. Only test days are scored.
#
# Output: output/table5_oos_logscores.csv, output/table6_oos_var.csv,
#         output/oos_fits.rds, figures/fig8_oos_cumulative_logscore
# Run time: roughly 30-60 minutes (three chains in parallel + filters).
# =====================================================================

for (p in c("Rcpp", "ggplot2", "scales")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({ library(ggplot2); library(scales) })
if (basename(getwd()) != "rupiah-jump-sv" && dir.exists("~/rupiah-jump-sv")) setwd("~/rupiah-jump-sv")
need <- c("svjx_sampler.cpp", "svjx_functions.R", "svjx_filter.cpp", "eval_functions.R",
          file.path("output", "model_data.rds"))
for (f in need) if (!file.exists(f)) stop("File ", f, " not found. Run 06_fit_models.R first.")
dir.create("figures", showWarnings = FALSE)

TRAIN_END <- as.Date("2018-12-31")
N_ITER <- 35000; BURN <- 5000; THIN <- 5
PROBS <- c(0.01, 0.05, 0.95, 0.99)
N_PART <- 10000; SEEDS <- c(101, 202, 303)

Rcpp::sourceCpp("svjx_filter.cpp")
source("eval_functions.R")

# ---- Data: re-standardise covariates with training statistics only ------------
md <- readRDS(file.path("output", "model_data.rds"))
y <- md$y; dates <- md$date; n <- length(y)
sc <- md$cov_scale
X_raw <- sweep(sweep(md$X[, -1, drop = FALSE], 2, sc$sd, `*`), 2, sc$mean, `+`)   # undo full-sample scaling
tr <- dates <= TRAIN_END
te <- !tr
if (sum(tr) < 1000 || sum(te) < 250) stop("Training or test period too short.")
mu_tr <- colMeans(X_raw[tr, , drop = FALSE]); sd_tr <- apply(X_raw[tr, , drop = FALSE], 2, sd)
X <- cbind("(Intercept)" = 1, sweep(sweep(X_raw, 2, mu_tr, `-`), 2, sd_tr, `/`))
message("Training: ", format(min(dates[tr])), " to ", format(max(dates[tr])), " (", sum(tr), " days)")
message("Test:     ", format(min(dates[te])), " to ", format(max(dates[te])), " (", sum(te), " days)")

# ---- Estimate on the training period ------------------------------------------
jobs <- list(list(model = "SV", seed = 21, X = NULL),
             list(model = "SVJ", seed = 22, X = NULL),
             list(model = "SVJX", seed = 23, X = X[tr, , drop = FALSE]))
run_job <- function(job, y, n_iter, burn, thin, wd) {
  setwd(wd)
  if (!exists("draw_h_tridiag")) Rcpp::sourceCpp("svjx_sampler.cpp")
  source("svjx_functions.R")
  fit_svjx(y, job$X, model = job$model, n_iter = n_iter, burn = burn, thin = thin,
           seed = job$seed, verbose = FALSE)
}
n_workers <- max(1, min(length(jobs), parallel::detectCores(logical = FALSE) - 1, na.rm = TRUE))
message("Fitting ", length(jobs), " chains on the training period with ", n_workers, " core(s) ...")
t0 <- proc.time()[3]
if (n_workers > 1) {
  cl <- parallel::makeCluster(n_workers)
  fits <- tryCatch(parallel::parLapplyLB(cl, jobs, run_job, y = y[tr], n_iter = N_ITER, burn = BURN,
                                         thin = THIN, wd = getwd()),
                   finally = parallel::stopCluster(cl))
} else {
  fits <- lapply(jobs, run_job, y = y[tr], n_iter = N_ITER, burn = BURN, thin = THIN, wd = getwd())
}
names(fits) <- c("SV", "SVJ", "SVJX")
for (nm in names(fits)) fits[[nm]][c("h_draws")] <- NULL      # keep the saved file small
saveRDS(list(fits = fits, train_end = TRAIN_END, scale = list(mean = mu_tr, sd = sd_tr)),
        file.path("output", "oos_fits.rds"), compress = "xz")
message("Estimation done in ", round((proc.time()[3] - t0) / 60, 1), " minutes.")

# ---- Forecast the test period ---------------------------------------------------
pm <- function(f, nm) if (nm %in% colnames(f$draws)) mean(f$draws[, nm]) else 0
pred <- list()
for (m in names(fits)) {
  f <- fits[[m]]
  g <- f$draws[, grep("^gamma_", colnames(f$draws)), drop = FALSE]
  pt <- if (ncol(g) == 0) rep(0, n) else {
    Xm <- if (m == "SVJX") X else matrix(1, n, 1)
    colMeans(pnorm(g %*% t(Xm)))
  }
  runs <- lapply(SEEDS, function(s) {
    set.seed(s)
    pf_svjx(y, pt, pm(f, "m"), pm(f, "mu"), pm(f, "phi"), pm(f, "sigma"), pm(f, "mu_xi"),
            max(pm(f, "sigma_xi"), 1e-6), N_PART, PROBS)
  })
  pred[[m]] <- list(logscore = rowMeans(sapply(runs, `[[`, "logscore")), q = runs[[1]]$q,
                    pit = runs[[1]]$pit,
                    mc_sd = sd(sapply(runs, function(r) sum(r$logscore[te]))))
  message("  ", m, ": test-period log score ", round(sum(pred[[m]]$logscore[te]), 1))
}
gfit <- fit_garch_t(y[tr])
pg <- predict_garch_t(gfit, y, PROBS)                 # filter the whole sample, training parameters
pred <- c(list(`GARCH-t` = list(logscore = pg$logscore, q = pg$q, pit = pg$pit, mc_sd = 0)), pred)
message("  GARCH-t: test-period log score ", round(sum(pg$logscore[te]), 1))
models <- names(pred)

# ---- Tables -----------------------------------------------------------------------
tab5 <- do.call(rbind, lapply(models, function(m) {
  dm <- if (m == "SVJX") rep(NA, 4) else dm_test(pred$SVJX$logscore[te] - pred[[m]]$logscore[te])
  data.frame(model = m, test_logscore = sum(pred[[m]]$logscore[te]), mc_sd = pred[[m]]$mc_sd,
             svjx_minus_model = dm[2], dm_z = dm[3], dm_p = dm[4],
             pit_ks_p = suppressWarnings(ks.test(pred[[m]]$pit[te], "punif")$p.value), row.names = NULL)
}))
tab5[, -1] <- signif(tab5[, -1], 5)
write.csv(tab5, file.path("output", "table5_oos_logscores.csv"), row.names = FALSE)

yt <- y[te]
tab6 <- do.call(rbind, lapply(models, function(m) {
  q <- pred[[m]]$q[te, , drop = FALSE]
  do.call(rbind, lapply(seq_along(PROBS), function(k) {
    lower <- PROBS[k] < 0.5
    hits <- as.integer(if (lower) yt < q[, k] else yt > q[, k])
    alpha <- if (lower) PROBS[k] else 1 - PROBS[k]
    ku <- kupiec_test(hits, alpha); ch <- christoffersen_test(hits, alpha)
    data.frame(model = m, tail = if (lower) "appreciation" else "depreciation", alpha = alpha,
               hits = sum(hits), rate = ku["rate"], p_uc = ku["p_uc"], p_ind = ch["p_ind"],
               p_cc = ch["p_cc"], row.names = NULL)
  }))
}))
tab6[, 5:8] <- signif(tab6[, 5:8], 3)
write.csv(tab6, file.path("output", "table6_oos_var.csv"), row.names = FALSE)

# ---- Figure 8: cumulative out-of-sample log-score advantage of SVJ-X ----------------
COL <- c("GARCH-t" = "#9AA0A6", SV = "#8FA9C8", SVJ = "#E0A458")
cum <- do.call(rbind, lapply(setdiff(models, "SVJX"), function(m)
  data.frame(date = dates[te], model = paste("SVJ-X vs", m),
             cum = cumsum(pred$SVJX$logscore[te] - pred[[m]]$logscore[te]))))
cum$model <- factor(cum$model, levels = paste("SVJ-X vs", setdiff(models, "SVJX")))
endlab <- cum[cum$date == max(cum$date), ]
# spread end labels so they do not overlap
endlab <- endlab[order(endlab$cum), ]
gap <- 0.06 * diff(range(cum$cum))
endlab$ylab <- endlab$cum
for (i in seq_len(nrow(endlab))[-1]) endlab$ylab[i] <- max(endlab$ylab[i], endlab$ylab[i - 1] + gap)
p8 <- ggplot(cum, aes(date, cum, colour = model)) +
  annotate("rect", xmin = as.Date("2020-02-20"), xmax = as.Date("2020-04-30"), ymin = -Inf, ymax = Inf,
           fill = "#8FA9C8", alpha = 0.18) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_line(linewidth = 0.6) + geom_point(data = endlab, size = 1.8) +
  geom_text(data = endlab, aes(y = ylab, label = sprintf("%+.1f", cum)), hjust = -0.25, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = setNames(unname(COL[setdiff(models, "SVJX")]), levels(cum$model)), name = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = expansion(c(0.01, 0.08))) +
  labs(x = NULL, y = "Cumulative out-of-sample\nlog-score difference") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path("figures", "fig8_oos_cumulative_logscore.png"), p8, width = 7.5, height = 3.6, dpi = 300, bg = "white")
ggsave(file.path("figures", "fig8_oos_cumulative_logscore.pdf"), p8, width = 7.5, height = 3.6)

# ---- Summary -------------------------------------------------------------------------
cat("\n================ OUT-OF-SAMPLE ================\n")
cat("Training:", format(min(dates[tr])), "to", format(max(dates[tr])), "| Test:",
    format(min(dates[te])), "to", format(max(dates[te])), "(", sum(te), "days )\n\n")
cat("Training-period parameters (posterior medians):\n")
for (m in names(fits)) cat(" ", m, ":", paste(colnames(fits[[m]]$draws),
                                               round(apply(fits[[m]]$draws, 2, median), 3), collapse = ", "), "\n")
cat("\nTable 5. Test-period log predictive score (higher is better):\n"); print(tab5, row.names = FALSE)
cat("\nTable 6. Test-period VaR backtests:\n"); print(tab6, row.names = FALSE)
