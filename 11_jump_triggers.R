# =====================================================================
# 11_jump_triggers.R
# What the SVJ-X jump triggers mean on the probability scale.
#
#   a) posterior jump probability as a function of the size of the
#      previous day's move in the broad US dollar index and in the US
#      10-year yield (other covariates at their sample means, no long
#      break), with 50% and 90% credible bands and the observed values
#      shown as a rug
#   b) daily posterior probability that a jump occurred, 2006-2026
#
# Needs the fits from 06_fit_models.R. Run time: under a minute.
# Output: output/table9_jump_prob_response.csv,
#         output/jump_probability_daily.csv, figures/fig12_jump_triggers
# =====================================================================

for (p in c("ggplot2", "patchwork", "scales")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
if (basename(getwd()) != "rupiah-jump-sv" && dir.exists("~/rupiah-jump-sv")) setwd("~/rupiah-jump-sv")
need <- file.path("output", c("model_data.rds", "fit_SVJX_chain1.rds", "fit_SVJX_chain2.rds"))
for (f in need) if (!file.exists(f)) stop("File ", f, " not found. Run 06_fit_models.R first.")
dir.create("figures", showWarnings = FALSE)

md <- readRDS(need[1])
f1 <- readRDS(need[2]); f2 <- readRDS(need[3])
g <- rbind(f1$draws, f2$draws)
g <- g[, paste0("gamma_", colnames(md$X)), drop = FALSE]
sc <- md$cov_scale
z <- function(v, x) (x - sc$mean[sc$covariate == v]) / sc$sd[sc$covariate == v]
raw <- function(v) md$X[, v] * sc$sd[sc$covariate == v] + sc$mean[sc$covariate == v]

# ---- (a) response curves ------------------------------------------------------
base <- c("(Intercept)" = 1, log_vix = 0, abs_d_ust10 = 0, abs_d_usd = 0,
          long_break = z("long_break", 0))
curve <- function(v, grid, label) {
  out <- lapply(grid, function(x) {
    xv <- base; xv[v] <- z(v, x)
    p <- pnorm(g %*% xv[colnames(md$X)])
    q <- quantile(p, c(0.05, 0.25, 0.5, 0.75, 0.95))
    data.frame(trigger = label, x = x, q05 = q[1], q25 = q[2], q50 = q[3], q75 = q[4], q95 = q[5])
  })
  do.call(rbind, out)
}
lab_usd <- "|Change in US dollar index| (%)"
lab_ust <- "|Change in US 10-year yield| (bp)"
cu <- rbind(curve("abs_d_usd", seq(0, 1.6, by = 0.02), lab_usd),
            curve("abs_d_ust10", seq(0, 25, by = 0.25), lab_ust))
cu$trigger <- factor(cu$trigger, levels = c(lab_usd, lab_ust))
set.seed(1)   # the yield is recorded in whole basis points: jitter the rug slightly
rug <- rbind(data.frame(trigger = lab_usd, x = raw("abs_d_usd")),
             data.frame(trigger = lab_ust, x = raw("abs_d_ust10") + runif(nrow(md$X), -0.45, 0.45)))
rug$trigger <- factor(rug$trigger, levels = levels(cu$trigger))
rug <- rug[(rug$trigger == lab_usd & rug$x <= 1.6) | (rug$trigger == lab_ust & rug$x >= 0 & rug$x <= 25), ]

tab9 <- rbind(curve("abs_d_usd", c(0.1, 0.25, 0.5, 0.75, 1, 1.5), lab_usd),
              curve("abs_d_ust10", c(1, 4, 8, 12, 16, 20), lab_ust))
tab9 <- data.frame(trigger = tab9$trigger, value = tab9$x, median = round(tab9$q50, 4),
                   lower90 = round(tab9$q05, 4), upper90 = round(tab9$q95, 4))
write.csv(tab9, file.path("output", "table9_jump_prob_response.csv"), row.names = FALSE)

COL <- c(main = "#1F3A5F", accent = "#C8102E", soft = "#8FA9C8")
theme_paper <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.tag = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

pa <- ggplot(cu, aes(x)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = COL["main"], alpha = 0.18) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = COL["main"], alpha = 0.35) +
  geom_line(aes(y = q50), colour = COL["main"], linewidth = 0.7) +
  geom_rug(data = rug, aes(x = x), inherit.aes = FALSE, sides = "b", alpha = 0.06,
           colour = COL["accent"], length = unit(0.05, "npc")) +
  facet_wrap(~ trigger, scales = "free_x", strip.position = "bottom") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = NULL, y = "Jump probability") + theme_paper +
  theme(strip.placement = "outside", strip.text = element_text(face = "plain"))

# ---- (b) daily posterior jump probabilities -----------------------------------------
jp <- (f1$jump_prob + f2$jump_prob) / 2
dd <- data.frame(date = md$date, ret = md$y, prob = jp,
                 abs_d_usd = round(raw("abs_d_usd"), 3))
write.csv(dd, file.path("output", "jump_probability_daily.csv"), row.names = FALSE)
dd$dir <- ifelse(dd$ret >= 0, "Rupiah depreciation", "Rupiah appreciation")
hi <- dd[dd$prob > 0.5, ]
episodes <- data.frame(start = as.Date(c("2008-09-15", "2013-05-22", "2020-02-20")),
                       end   = as.Date(c("2009-03-31", "2013-09-30", "2020-04-30")),
                       label = c("Global financial crisis", "Taper tantrum", "COVID-19"))
pb <- ggplot(dd, aes(date, prob)) +
  geom_rect(data = episodes, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = COL["soft"], alpha = 0.25) +
  geom_text(data = episodes, aes(x = start, y = 1.08, label = label), inherit.aes = FALSE,
            hjust = 0, size = 2.7, colour = "grey30") +
  geom_segment(aes(xend = date, y = 0, yend = prob), colour = "grey70", linewidth = 0.15) +
  geom_point(data = hi, aes(colour = dir, size = abs(ret)), alpha = 0.85) +
  scale_colour_manual(values = c("Rupiah depreciation" = COL[["accent"]],
                                 "Rupiah appreciation" = COL[["main"]]), name = NULL) +
  scale_size_continuous(range = c(0.6, 3.2), name = "|Return| (%)", breaks = c(0.5, 1, 2, 3)) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1)) +
  coord_cartesian(ylim = c(0, 1.12), clip = "off") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(x = NULL, y = "Posterior jump probability") + theme_paper +
  theme(legend.position = "bottom")

fig <- pa / pb + plot_layout(heights = c(1, 1.15)) + plot_annotation(tag_levels = "a")
ggsave(file.path("figures", "fig12_jump_triggers.png"), fig, width = 8, height = 7, dpi = 300, bg = "white")
ggsave(file.path("figures", "fig12_jump_triggers.pdf"), fig, width = 8, height = 7)

# ---- Summary ------------------------------------------------------------------------
cat("\n================ JUMP TRIGGERS ================\n")
cat("Jump probability at given trigger values (other covariates at their means):\n")
print(tab9, row.names = FALSE)
cat("\nDays with posterior jump probability above 0.5:", nrow(hi),
    "| of which zero returns:", sum(hi$ret == 0), "\n")
print(table(format(hi$date, "%Y")))
