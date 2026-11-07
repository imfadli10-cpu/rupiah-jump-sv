# =====================================================================
# 07_compare_models.R
# One-step-ahead predictive comparison of GARCH-t, SV, SVJ and SVJ-X.
#
#   - log predictive scores from a particle filter (10,000 particles) at
#     posterior-mean parameters (including the daily drift m); for SVJ-X the daily jump probability is
#     the posterior mean of Phi(x_t' gamma)
#   - Diebold-Mariano-type tests of log-score differences against SVJ-X
#   - Value-at-Risk backtests at 1% and 5% in both tails
#     (upper tail = rupiah depreciation)
#   - PIT calibration
#
# Output: output/table3_logscores.csv, output/table4_var_backtests.csv,
#         output/predictive_daily.rds
#         figures/fig5_cumulative_logscore, fig6_pit, fig7_var_barcode
# Requires: output/fit_*.rds and output/model_data.rds from 06_fit_models.R
# Run time: about 5-15 minutes.
# =====================================================================

for (p in c("Rcpp", "ggplot2", "patchwork", "scales")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(scales) })
if (basename(getwd()) != "rupiah-jump-sv" && dir.exists("~/rupiah-jump-sv")) setwd("~/rupiah-jump-sv")
need <- c("svjx_filter.cpp", "eval_functions.R", file.path("output", "model_data.rds"),
          file.path("output", paste0("fit_", c("SV", "SVJ", "SVJX_chain1", "SVJX_chain2"), ".rds")))
for (f in need) if (!file.exists(f)) stop("File ", f, " not found. Run 06_fit_models.R first.")
dir.create("figures", showWarnings = FALSE)

message("Compiling the particle filter ...")
Rcpp::sourceCpp("svjx_filter.cpp")
source("eval_functions.R")

md <- readRDS(file.path("output", "model_data.rds"))
y <- md$y; X <- md$X; dates <- md$date; n <- length(y)
fit <- list(SV = readRDS(file.path("output", "fit_SV.rds")),
            SVJ = readRDS(file.path("output", "fit_SVJ.rds")))
fit$SVJX <- readRDS(file.path("output", "fit_SVJX_chain1.rds"))
fit$SVJX$draws <- rbind(fit$SVJX$draws, readRDS(file.path("output", "fit_SVJX_chain2.rds"))$draws)

PROBS <- c(0.01, 0.05, 0.95, 0.99)
N_PART <- 10000
SEEDS <- c(101, 202, 303)

# ---- Particle-filter evaluation of the Bayesian models -------------------------
pm <- function(f, nm) if (nm %in% colnames(f$draws)) mean(f$draws[, nm]) else 0
jump_prob <- function(f, Xm) {
  g <- f$draws[, grep("^gamma_", colnames(f$draws)), drop = FALSE]
  if (ncol(g) == 0) return(rep(0, n))
  colMeans(pnorm(g %*% t(Xm)))                     # posterior mean of Phi(x_t' gamma)
}
pred <- list(); mc_sd <- c()
for (m in c("SV", "SVJ", "SVJX")) {
  f <- fit[[m]]
  Xm <- if (m == "SVJX") X else matrix(1, n, 1)
  pt <- jump_prob(f, Xm)
  runs <- lapply(SEEDS, function(s) {
    set.seed(s)
    pf_svjx(y, pt, pm(f, "m"), pm(f, "mu"), pm(f, "phi"), pm(f, "sigma"), pm(f, "mu_xi"),
            max(pm(f, "sigma_xi"), 1e-6), N_PART, PROBS)
  })
  tot <- sapply(runs, function(r) sum(r$logscore))
  mc_sd[m] <- sd(tot)
  pred[[m]] <- list(logscore = rowMeans(sapply(runs, `[[`, "logscore")),
                    pit = runs[[1]]$pit, q = runs[[1]]$q, h_filt = runs[[1]]$h_filt, jump_prob = pt)
  message("  ", m, ": total log score ", round(mean(tot), 1), " (MC sd ", round(sd(tot), 2), ")")
}

# ---- GARCH-t benchmark -----------------------------------------------------------
g <- fit_garch_t(y)
pg <- predict_garch_t(g, y, PROBS)
pred <- c(list(`GARCH-t` = list(logscore = pg$logscore, pit = pg$pit, q = pg$q)), pred)
mc_sd["GARCH-t"] <- 0
message("  GARCH-t: total log score ", round(sum(pg$logscore), 1))
models <- names(pred)                                   # GARCH-t, SV, SVJ, SVJX
saveRDS(list(date = dates, y = y, pred = pred, garch = g), file.path("output", "predictive_daily.rds"))

# ---- Table 3: log scores ----------------------------------------------------------
tab3 <- do.call(rbind, lapply(models, function(m) {
  dm <- if (m == "SVJX") c(NA, NA, NA, NA) else dm_test(pred$SVJX$logscore - pred[[m]]$logscore)
  data.frame(model = m, total_logscore = sum(pred[[m]]$logscore), mc_sd = mc_sd[m],
             svjx_minus_model = dm[2], dm_z = dm[3], dm_p = dm[4], row.names = NULL)
}))
tab3[, 2:6] <- signif(tab3[, 2:6], 5)
write.csv(tab3, file.path("output", "table3_logscores.csv"), row.names = FALSE)

# ---- Table 4: VaR backtests ---------------------------------------------------------
tab4 <- do.call(rbind, lapply(models, function(m) {
  q <- pred[[m]]$q
  do.call(rbind, lapply(seq_along(PROBS), function(k) {
    lower <- PROBS[k] < 0.5
    hits <- as.integer(if (lower) y < q[, k] else y > q[, k])
    alpha <- if (lower) PROBS[k] else 1 - PROBS[k]
    ku <- kupiec_test(hits, alpha); ch <- christoffersen_test(hits, alpha)
    data.frame(model = m, tail = if (lower) "appreciation" else "depreciation", alpha = alpha,
               hits = sum(hits), rate = ku["rate"], p_uc = ku["p_uc"], p_ind = ch["p_ind"],
               p_cc = ch["p_cc"], row.names = NULL)
  }))
}))
tab4[, 5:8] <- signif(tab4[, 5:8], 3)
write.csv(tab4, file.path("output", "table4_var_backtests.csv"), row.names = FALSE)

# ---- Figures ----------------------------------------------------------------------
COL <- c("GARCH-t" = "#9AA0A6", SV = "#8FA9C8", SVJ = "#E0A458", SVJX = "#C8102E")
theme_paper <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        plot.tag = element_text(face = "bold"), legend.position = "bottom", strip.text = element_text(face = "bold"))
episodes <- data.frame(start = as.Date(c("2008-09-15", "2013-05-22", "2015-08-01", "2020-02-20")),
                       end = as.Date(c("2009-03-31", "2013-12-31", "2015-10-15", "2020-04-30")))
save_fig <- function(p, name, w, h) {
  ggsave(file.path("figures", paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path("figures", paste0(name, ".pdf")), p, width = w, height = h)
  message("  saved figures/", name)
}

# Fig 5: cumulative log-score advantage of SVJ-X over each competitor
cum <- do.call(rbind, lapply(setdiff(models, "SVJX"), function(m)
  data.frame(date = dates, model = paste("SVJ-X vs", m), cum = cumsum(pred$SVJX$logscore - pred[[m]]$logscore))))
cum$model <- factor(cum$model, levels = paste("SVJ-X vs", setdiff(models, "SVJX")))
endlab <- cum[cum$date == max(cum$date), ]
p5 <- ggplot(cum, aes(date, cum, colour = model)) +
  geom_rect(data = episodes, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
            fill = "#8FA9C8", alpha = 0.18) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  geom_point(data = endlab, size = 1.8) +
  geom_text(data = endlab, aes(label = sprintf("%+.0f", cum)), hjust = -0.25, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = setNames(COL[setdiff(models, "SVJX")], levels(cum$model)), name = NULL) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = expansion(c(0.01, 0.06))) +
  labs(x = NULL, y = "Cumulative log-score difference\n(positive = SVJ-X predicts better)") + theme_paper
save_fig(p5, "fig5_cumulative_logscore", 7.5, 3.8)

# Fig 6: PIT histograms with a 95% band for a uniform histogram
nb <- 20
pit_df <- do.call(rbind, lapply(models, function(m) {
  h <- hist(pred[[m]]$pit, breaks = seq(0, 1, length.out = nb + 1), plot = FALSE)
  data.frame(model = m, mid = h$mids, dens = h$density)
}))
pit_df$model <- factor(pit_df$model, levels = models)
band <- qbinom(c(0.025, 0.975), n, 1 / nb) / (n / nb)
p6 <- ggplot(pit_df, aes(mid, dens, fill = model)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = band[1], ymax = band[2], fill = "grey80", alpha = 0.5) +
  geom_col(width = 1 / nb * 0.92, alpha = 0.9) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey30", linewidth = 0.3) +
  facet_wrap(~ model, nrow = 1) + scale_fill_manual(values = COL, guide = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "PIT value", y = "Density") + theme_paper
save_fig(p6, "fig6_pit", 7.5, 2.6)

# Fig 7: barcode of 1% depreciation-VaR exceedances per model, plus rolling 1-year hit rate
k99 <- which(PROBS == 0.99)
bar <- do.call(rbind, lapply(models, function(m) data.frame(date = dates, model = m, hit = y > pred[[m]]$q[, k99])))
bar$model <- factor(bar$model, levels = rev(models))
p7a <- ggplot(bar[bar$hit, ], aes(date, model, colour = model)) +
  geom_rect(data = episodes, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
            fill = "#8FA9C8", alpha = 0.18) +
  geom_point(shape = "|", size = 4) + scale_colour_manual(values = COL, guide = "none") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(x = NULL, y = NULL, subtitle = "Days on which the depreciation exceeded the 1% VaR") + theme_paper
roll <- do.call(rbind, lapply(models, function(m) {
  hit <- as.numeric(y > pred[[m]]$q[, k99])
  data.frame(date = dates, model = m, rate = 100 * as.numeric(stats::filter(hit, rep(1 / 250, 250), sides = 1)))
}))
roll$model <- factor(roll$model, levels = models)
p7b <- ggplot(roll[!is.na(roll$rate), ], aes(date, rate, colour = model)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey30", linewidth = 0.3) +
  geom_line(linewidth = 0.45) + scale_colour_manual(values = COL, name = NULL) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(x = NULL, y = "Exceedance rate over\nprevious 250 days (%)") + theme_paper
p7 <- p7a / p7b + plot_layout(heights = c(1, 1.2)) + plot_annotation(tag_levels = "a")
save_fig(p7, "fig7_var_barcode", 7.5, 5)

# ---- Summary ---------------------------------------------------------------------
cat("\n================ MODEL COMPARISON ================\n")
cat("Days evaluated:", n, "|", format(min(dates)), "to", format(max(dates)), "\n")
cat("GARCH-t parameters:", paste(names(g$par), round(g$par, 4), collapse = ", "), "\n\n")
cat("Table 3. Total log predictive score (higher is better):\n"); print(tab3, row.names = FALSE)
cat("\nTable 4. VaR backtests (p-values below 0.05 indicate rejection):\n"); print(tab4, row.names = FALSE)
cat("\nPIT uniformity (Kolmogorov-Smirnov p-values):\n")
print(sapply(pred, function(r) signif(suppressWarnings(ks.test(r$pit, "punif")$p.value), 3)))
