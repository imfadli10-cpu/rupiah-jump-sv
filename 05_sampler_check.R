# =====================================================================
# 05_sampler_check.R
# Checks that the SVJ-X Gibbs sampler recovers known parameters from
# simulated data (n = 5000, one covariate). Nothing here uses real data.
#
#   Output: output/sampler_check_recovery.csv
#           figures/figA1_sampler_recovery (posterior vs true value,
#           latent volatility, jump probabilities)
# Run time: about 3-6 minutes (10,000 iterations).
# =====================================================================

for (p in c("Rcpp", "ggplot2", "patchwork", "scales")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(scales) })

if (basename(getwd()) != "rupiah-jump-sv" && dir.exists("~/rupiah-jump-sv")) setwd("~/rupiah-jump-sv")
for (f in c("svjx_sampler.cpp", "svjx_functions.R")) {
  if (!file.exists(f)) stop("File ", f, " not found in ", getwd(), ". Put it in the rupiah-jump-sv folder.")
}
dir.create("output", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

message("Compiling the C++ sampler (first time can take ~30 s) ...")
Rcpp::sourceCpp("svjx_sampler.cpp")
source("svjx_functions.R")

# ---- Simulate ------------------------------------------------------------
set.seed(2026)
n <- 5000
x1 <- as.numeric(arima.sim(list(ar = 0.98), n))       # persistent covariate (like VIX)
x1 <- (x1 - mean(x1)) / sd(x1)
X <- cbind("(Intercept)" = 1, x = x1)
truth <- c(m = 0.02, mu = -1.6, phi = 0.97, sigma = 0.2, mu_xi = 0.5, sigma_xi = 2,
           "gamma_(Intercept)" = -2.4, gamma_x = 0.5)
sim <- sim_svjx(n, truth["mu"], truth["phi"], truth["sigma"],
                truth[c("gamma_(Intercept)", "gamma_x")], truth["mu_xi"], truth["sigma_xi"], X = X,
                m = truth["m"])
message("Simulated ", n, " days with ", sum(sim$J), " jumps.")

# ---- Fit -------------------------------------------------------------------
message("Running the sampler ...")
fit <- fit_svjx(sim$y, X, model = "SVJX", n_iter = 10000, burn = 3000, thin = 5, seed = 1)

q <- t(apply(fit$draws, 2, quantile, c(0.025, 0.5, 0.975)))
rec <- data.frame(parameter = rownames(q), true = unname(truth[rownames(q)]),
                  lower95 = q[, 1], median = q[, 2], upper95 = q[, 3], row.names = NULL)
rec$inside95 <- rec$true >= rec$lower95 & rec$true <= rec$upper95
rec[, 2:5] <- round(rec[, 2:5], 3)
write.csv(rec, file.path("output", "sampler_check_recovery.csv"), row.names = FALSE)

# ---- Figure A1 ---------------------------------------------------------------
COL <- c(main = "#1F3A5F", accent = "#C8102E", soft = "#8FA9C8")
theme_paper <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        plot.tag = element_text(face = "bold"), strip.text = element_text(face = "bold"))

pd <- data.frame(parameter = rep(colnames(fit$draws), each = nrow(fit$draws)),
                 value = as.vector(fit$draws))
pa <- ggplot(pd, aes(value)) +
  geom_density(fill = COL["main"], colour = COL["main"], alpha = 0.3, linewidth = 0.3) +
  geom_vline(data = rec, aes(xintercept = true), colour = COL["accent"], linetype = 2) +
  facet_wrap(~ parameter, scales = "free", nrow = 2) +
  labs(x = "Posterior draws (dashed: true value)", y = NULL) + theme_paper +
  theme(axis.text.y = element_blank())

win <- 1:600
hq <- apply(fit$h_draws[, win], 2, quantile, c(0.05, 0.25, 0.75, 0.95))
hd <- data.frame(t = win, true = sim$h[win], q05 = hq[1, ], q25 = hq[2, ], q75 = hq[3, ], q95 = hq[4, ])
pb <- ggplot(hd, aes(t)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = COL["main"], alpha = 0.2) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = COL["main"], alpha = 0.35) +
  geom_line(aes(y = true), colour = COL["accent"], linewidth = 0.4) +
  labs(x = "Day", y = "Log-variance h") + theme_paper

jd <- data.frame(t = win, prob = fit$jump_prob[win], true_jump = sim$J[win] == 1)
pc <- ggplot(jd, aes(t, prob)) +
  geom_segment(aes(xend = t, y = 0, yend = prob), colour = COL["soft"], linewidth = 0.3) +
  geom_point(data = jd[jd$true_jump, ], colour = COL["accent"], size = 1.6) +
  labs(x = "Day", y = "Posterior jump probability") + theme_paper

figA1 <- pa / (pb | pc) + plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "a",
                  caption = "b: true log-variance (red) with 50% and 90% posterior bands; c: red points are true jump days.")
ggsave(file.path("figures", "figA1_sampler_recovery.png"), figA1, width = 8, height = 6.5, dpi = 300, bg = "white")
ggsave(file.path("figures", "figA1_sampler_recovery.pdf"), figA1, width = 8, height = 6.5)

# ---- Summary -------------------------------------------------------------------
cat("\n================ SAMPLER CHECK ================\n")
cat("Run time:", round(fit$seconds), "seconds | phi acceptance rate:", round(fit$acc_phi, 2), "\n")
cat("True jumps:", sum(sim$J), "| posterior mean number of jumps:", round(mean(fit$n_jumps), 1), "\n\n")
print(rec, row.names = FALSE)
cat("\nParameters inside their 95% interval:", sum(rec$inside95), "of", nrow(rec), "\n")
cat("Correlation of true and estimated log-variance:", round(cor(fit$h_mean, sim$h), 3), "\n")
