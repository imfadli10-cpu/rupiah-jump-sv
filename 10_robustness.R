# =====================================================================
# 10_robustness.R
# Robustness of the SVJ-X jump-trigger estimates. Each variant changes one
# thing relative to the main specification and refits the SVJ-X model:
#
#   baseline        main specification (same MCMC length as the variants)
#   drop_20081028   the largest return in the sample (+13.0%) removed
#   post_2009       sample starts 2009-07-01 (global financial crisis excluded)
#   ecb_source      returns from the ECB reference rate instead of the BI mid
#                   rate; ECB levels more than 2% away from BI are treated as
#                   bad ticks and replaced by log-linear interpolation
#   offset_1e-3     offset in log(y^2 + c) raised to 1e-3
#   offset_1e-5     offset lowered to 1e-5
#   wide_gamma      probit prior N(-2, 2^2) for the intercept, N(0, 2^2) slopes
#   vague_jumpsize  mu_xi ~ N(0, 10^2), sigma_xi^2 ~ InvGamma(2, 1)
#   oil_placebo     adds |change in Brent oil price| as a fifth covariate
#
# Covariates are standardised within each variant's own sample.
# Output: output/table8_robustness.csv, output/robustness_draws.rds, robustness/*.rds,
#         figures/fig11_robustness (ridgelines of trigger effects)
# Run time: roughly 1.5-3 hours (nine chains, run in parallel).
# =====================================================================

for (p in c("Rcpp", "ggplot2", "ggridges", "posterior", "scales")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({ library(ggplot2); library(ggridges) })
if (basename(getwd()) != "rupiah-jump-sv" && dir.exists("~/rupiah-jump-sv")) setwd("~/rupiah-jump-sv")
for (f in c("svjx_sampler.cpp", "svjx_functions.R", file.path("data", "panel_bi.rds"))) {
  if (!file.exists(f)) stop("File ", f, " not found in ", getwd(), ".")
}
dir.create("output", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

N_ITER <- 20000; BURN <- 5000; THIN <- 5
START <- as.Date("2006-01-04")                       # first day of the main sample

# ---- Data ---------------------------------------------------------------------
panel <- readRDS(file.path("data", "panel_bi.rds"))
panel <- panel[order(panel$date), ]
panel$abs_d_ust10 <- abs(c(NA, diff(panel$dgs10_lag))) * 100
panel$abs_d_usd   <- abs(c(NA, 100 * diff(log(panel$dtwexbgs_lag))))
panel$log_vix     <- log(panel$vixcls_lag)
panel$long_break  <- as.numeric(panel$gap >= 4)
panel$abs_d_oil   <- abs(c(NA, 100 * diff(log(panel$dcoilbrenteu_lag))))

# ECB-based returns on BI dates (bad ticks and missing days interpolated in logs)
ecb <- panel$ecb_usdidr
bad <- !is.na(panel$level_gap_pct) & abs(panel$level_gap_pct) > 2
ecb[bad] <- NA
ok_e <- which(!is.na(ecb))
lecb <- approx(ok_e, log(ecb[ok_e]), xout = seq_along(ecb), rule = 2)$y
panel$ret_ecb <- c(NA, 100 * diff(lecb))
n_bad <- sum(bad & panel$date >= START)

carry <- function(x) { for (i in seq_along(x)[-1]) if (is.na(x[i])) x[i] <- x[i - 1]; x }
std <- function(x) (x - mean(x)) / sd(x)
make_data <- function(start = START, drop = NULL, ret = "ret",
                      covs = c("log_vix", "abs_d_ust10", "abs_d_usd", "long_break")) {
  d <- panel[panel$date >= start & !(panel$date %in% drop), ]
  for (v in covs) d[[v]] <- carry(d[[v]])
  d <- d[complete.cases(d[, c(ret, covs)]), ]
  list(y = d[[ret]], X = cbind("(Intercept)" = 1, sapply(d[, covs], std)), date = d$date)
}

variants <- list(
  baseline       = list(data = make_data()),
  drop_20081028  = list(data = make_data(drop = as.Date("2008-10-28"))),
  post_2009      = list(data = make_data(start = as.Date("2009-07-01"))),
  ecb_source     = list(data = make_data(ret = "ret_ecb")),
  `offset_1e-3`  = list(data = make_data(), offset = 1e-3),
  `offset_1e-5`  = list(data = make_data(), offset = 1e-5),
  wide_gamma     = list(data = make_data(), priors = "wide_gamma"),
  vague_jumpsize = list(data = make_data(),
                        priors = list(muxi_sd = 10, sxi_shape = 2, sxi_rate = 1)),
  oil_placebo    = list(data = make_data(covs = c("log_vix", "abs_d_ust10", "abs_d_usd",
                                                  "long_break", "abs_d_oil"))))
for (v in names(variants)) {
  dv <- variants[[v]]$data
  if (length(dv$y) < 2000 || anyNA(dv$y) || anyNA(dv$X)) stop("Variant ", v, ": bad data.")
  if (identical(variants[[v]]$priors, "wide_gamma")) {
    k <- ncol(dv$X)
    variants[[v]]$priors <- list(gamma_mean = c(-2, rep(0, k - 1)), gamma_sd = rep(2, k))
  }
}
message("Variants: ", paste(names(variants), collapse = ", "))
message("ECB bad ticks replaced: ", n_bad, " days")

# ---- Fit (in parallel where possible; each variant is saved, so a rerun resumes) -----
dir.create("robustness", showWarnings = FALSE)
for (v in names(variants)) variants[[v]]$file <- file.path("robustness", paste0(v, ".rds"))
run_variant <- function(v, n_iter, burn, thin, wd) {
  setwd(wd)
  if (file.exists(v$file)) {
    old <- readRDS(v$file)
    if (identical(old$settings, c(n_iter, burn, thin))) return(old)
  }
  if (!exists("draw_h_tridiag")) Rcpp::sourceCpp("svjx_sampler.cpp")
  source("svjx_functions.R")
  f <- fit_svjx(v$data$y, v$data$X, model = "SVJX", n_iter = n_iter, burn = burn, thin = thin,
                priors = v$priors, offset = if (is.null(v$offset)) 1e-4 else v$offset, seed = 31,
                h_draws_keep = 0, verbose = FALSE)
  out <- list(draws = f$draws, n_jumps = f$n_jumps, acc_gamma = f$acc_gamma, n = length(v$data$y),
              start = min(v$data$date), seconds = f$seconds, settings = c(n_iter, burn, thin))
  saveRDS(out, v$file)
  out
}
n_done <- sum(file.exists(sapply(variants, `[[`, "file")))
n_workers <- max(1, min(length(variants), parallel::detectCores(logical = FALSE) - 1, na.rm = TRUE))
message("Fitting ", length(variants), " variants on ", n_workers, " core(s) (", n_done,
        " saved from an earlier run). Please wait ...")
t0 <- proc.time()[3]
if (n_workers > 1) {
  cl <- parallel::makeCluster(n_workers)
  fits <- tryCatch(parallel::parLapplyLB(cl, variants, run_variant, n_iter = N_ITER, burn = BURN,
                                         thin = THIN, wd = getwd()),
                   finally = parallel::stopCluster(cl))
} else {
  fits <- lapply(variants, run_variant, n_iter = N_ITER, burn = BURN, thin = THIN, wd = getwd())
}
names(fits) <- names(variants)
saveRDS(fits, file.path("output", "robustness_draws.rds"), compress = "xz")
message("Done in ", round((proc.time()[3] - t0) / 60, 1), " minutes.")

# ---- Table 8 --------------------------------------------------------------------------
pretty <- c(gamma_log_vix = "VIX (log)", gamma_abs_d_ust10 = "|Change in US 10y yield|",
            gamma_abs_d_usd = "|Change in US dollar index|", gamma_long_break = "Long market break",
            gamma_abs_d_oil = "|Change in Brent oil|")
tab8 <- do.call(rbind, lapply(names(fits), function(v) {
  dr <- fits[[v]]$draws
  g <- grep("^gamma_", colnames(dr), value = TRUE)
  g <- g[g != "gamma_(Intercept)"]
  do.call(rbind, lapply(g, function(p) {
    x <- dr[, p]
    data.frame(variant = v, n = fits[[v]]$n, covariate = pretty[p], median = median(x),
               lower95 = quantile(x, 0.025), upper95 = quantile(x, 0.975),
               prob_positive = mean(x > 0), mean_jump_days = mean(fits[[v]]$n_jumps),
               ess = round(posterior::ess_bulk(x)), row.names = NULL)
  }))
}))
tab8[, c("median", "lower95", "upper95", "prob_positive", "mean_jump_days")] <-
  round(tab8[, c("median", "lower95", "upper95", "prob_positive", "mean_jump_days")], 3)
write.csv(tab8, file.path("output", "table8_robustness.csv"), row.names = FALSE)

# ---- Figure 11: ridgelines of the trigger effects across variants ----------------------
COL <- c(main = "#1F3A5F", accent = "#C8102E", soft = "#8FA9C8")
rd <- do.call(rbind, lapply(names(fits), function(v) {
  dr <- fits[[v]]$draws
  g <- intersect(names(pretty), colnames(dr))
  data.frame(variant = v, covariate = rep(pretty[g], each = nrow(dr)), value = as.vector(dr[, g]))
}))
rd$variant <- factor(rd$variant, levels = rev(names(fits)))
rd$covariate <- factor(rd$covariate, levels = pretty)
rd$is_base <- rd$variant == "baseline"
sm <- tab8; sm$variant <- factor(sm$variant, levels = levels(rd$variant))
sm$covariate <- factor(sm$covariate, levels = pretty)
base_med <- sm[sm$variant == "baseline", c("covariate", "median")]
# trim each panel to the central 99% of all draws so that a few extreme draws of a
# rare binary trigger do not squash the panel; interval bars are clipped to the same range
lim <- tapply(rd$value, rd$covariate, quantile, c(0.005, 0.995))
lo <- sapply(lim, function(z) z[[1]])[as.character(rd$covariate)]; hi <- sapply(lim, function(z) z[[2]])[as.character(rd$covariate)]
rd <- rd[rd$value >= lo & rd$value <= hi, ]
slo <- sapply(lim, function(z) z[[1]])[as.character(sm$covariate)]; shi <- sapply(lim, function(z) z[[2]])[as.character(sm$covariate)]
sm$lower95 <- pmax(sm$lower95, slo); sm$upper95 <- pmin(sm$upper95, shi)

p11 <- ggplot(rd, aes(value, variant)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
  geom_vline(data = base_med, aes(xintercept = median), colour = COL["accent"],
             linetype = 2, linewidth = 0.35) +
  geom_density_ridges(aes(fill = is_base), scale = 0.9, rel_min_height = 0.01,
                      colour = "white", linewidth = 0.25, alpha = 0.85) +
  geom_segment(data = sm, aes(x = lower95, xend = upper95, y = variant, yend = variant),
               colour = "grey15", linewidth = 0.45) +
  geom_point(data = sm, aes(x = median, y = variant), colour = "grey15", size = 1.1) +
  scale_fill_manual(values = c(`TRUE` = COL[["accent"]], `FALSE` = COL[["main"]]), guide = "none") +
  scale_x_continuous(breaks = scales::breaks_extended(4)) +
  facet_wrap(~ covariate, nrow = 1, scales = "free_x", drop = TRUE) +
  labs(x = "Probit coefficient on the jump probability (per standard deviation)", y = NULL,
       caption = "Bars: 95% posterior intervals; dots: medians; dashed red line: baseline median.\nEach panel shows the central 99% of the posterior draws.") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        panel.grid.major.y = element_blank(), panel.spacing.x = unit(1.2, "lines"))
suppressMessages({
  ggsave(file.path("figures", "fig11_robustness.png"), p11, width = 10, height = 4.6, dpi = 300, bg = "white")
  ggsave(file.path("figures", "fig11_robustness.pdf"), p11, width = 10, height = 4.6)
})

# ---- Summary ------------------------------------------------------------------------------
cat("\n================ ROBUSTNESS ================\n")
cat("Iterations:", N_ITER, "| burn-in:", BURN, "| thin:", THIN, "| ECB bad ticks replaced:", n_bad, "\n")
cat("Minutes per variant:", paste(names(fits), round(sapply(fits, `[[`, "seconds") / 60, 1),
                                  collapse = ", "), "\n\n")
print(tab8[, c("variant", "covariate", "median", "lower95", "upper95", "prob_positive", "ess")],
      row.names = FALSE)
cat("\nPosterior mean number of jump days:\n")
print(unique(tab8[, c("variant", "n", "mean_jump_days")]), row.names = FALSE)
