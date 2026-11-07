# =====================================================================
# 04_stylized_facts.R
# Descriptive statistics, stylized facts and a non-parametric jump test
# for daily USD/IDR returns (Bank Indonesia mid rate).
#
#   Table 1  : descriptive statistics and tests  -> output/table1_descriptives.csv
#   Jumps    : Lee-Mykland (2008) test on daily data -> output/lm_jumps.csv
#   Fig 1    : level, returns with detected jumps, rolling volatility
#   Fig 2    : tails - log-density vs Gaussian, and QQ plot
#   Fig 3    : calendar heatmap of absolute returns (year x month)
#   Fig 4    : autocorrelation of r, |r| and r^2 (lollipop)
#   Figures are saved in figures/ as PNG (300 dpi) and PDF.
#
# Requires: data/panel_bi.rds from 03_build_panel.R
# =====================================================================

for (p in c("ggplot2", "patchwork", "scales")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(scales)
})

if (basename(getwd()) != "rupiah-jump-sv" && dir.exists("~/rupiah-jump-sv")) {
  setwd("~/rupiah-jump-sv")
}
if (!file.exists(file.path("data", "panel_bi.rds"))) {
  stop("data/panel_bi.rds not found. Run 03_build_panel.R first.")
}
dir.create("output", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

panel <- readRDS(file.path("data", "panel_bi.rds"))
d <- panel[!is.na(panel$ret), c("date", "usdidr", "ret", "gap")]
r <- d$ret
n <- length(r)

# ---- Common figure style ----------------------------------------------
COL <- c(main = "#1F3A5F", accent = "#C8102E", soft = "#8FA9C8", grey = "#9AA0A6")
theme_paper <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        axis.title = element_text(colour = "grey25"),
        plot.title = element_text(face = "bold", size = 10),
        plot.tag = element_text(face = "bold"),
        legend.position = "bottom")
save_fig <- function(p, name, w, h) {
  ggsave(file.path("figures", paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path("figures", paste0(name, ".pdf")), p, width = w, height = h)
  message("  saved figures/", name)
}

# Approximate stress episodes, used only for shading in Fig 1
episodes <- data.frame(
  start = as.Date(c("2005-08-01", "2008-09-15", "2013-05-22", "2015-08-01", "2020-02-20")),
  end   = as.Date(c("2005-09-30", "2009-03-31", "2013-12-31", "2015-10-15", "2020-04-30")),
  label = c("2005", "GFC", "Taper\ntantrum", "2015", "COVID-19"))

# ---- 1. Descriptive statistics and tests -------------------------------
skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurt <- function(x) mean((x - mean(x))^4) / sd(x)^4
jarque_bera <- function(x) {
  s <- skew(x); k <- kurt(x); m <- length(x)
  stat <- m / 6 * (s^2 + (k - 3)^2 / 4)
  c(stat = stat, p = pchisq(stat, df = 2, lower.tail = FALSE))
}
arch_lm <- function(x, lags = 5) {
  e2 <- (x - mean(x))^2
  X <- embed(e2, lags + 1)
  fit <- lm(X[, 1] ~ X[, -1])
  stat <- nrow(X) * summary(fit)$r.squared
  c(stat = stat, p = pchisq(stat, df = lags, lower.tail = FALSE))
}
lb_r  <- Box.test(r,   lag = 10, type = "Ljung-Box")
lb_r2 <- Box.test(r^2, lag = 10, type = "Ljung-Box")
jb <- jarque_bera(r)
al <- arch_lm(r, 5)

table1 <- data.frame(
  statistic = c("Observations", "Mean (%)", "Std. dev. (%)", "Minimum (%)", "Maximum (%)",
                "Skewness", "Excess kurtosis", "Share of zero returns",
                "Jarque-Bera", "Ljung-Box Q(10), r", "Ljung-Box Q(10), r^2", "ARCH-LM(5)"),
  value = c(n, mean(r), sd(r), min(r), max(r), skew(r), kurt(r) - 3, mean(r == 0),
            jb["stat"], lb_r$statistic, lb_r2$statistic, al["stat"]),
  p_value = c(rep(NA, 8), jb["p"], lb_r$p.value, lb_r2$p.value, al["p"]))
table1$value <- signif(table1$value, 4)
table1$p_value <- signif(table1$p_value, 3)
write.csv(table1, file.path("output", "table1_descriptives.csv"), row.names = FALSE)

# ---- 2. Lee-Mykland jump test on daily returns ---------------------------
# Local volatility from bipower variation over the previous K returns
# (window excludes day t). Critical value from the Gumbel limit.
lee_mykland <- function(x, K = 16, alpha = 0.01) {
  m <- length(x)
  bv <- rep(NA_real_, m)
  ax <- abs(x)
  for (t in (K + 1):m) {
    idx <- (t - K + 1):(t - 1)                      # K-1 returns before t
    bv[t] <- sum(ax[idx] * ax[idx - 1]) / (K - 2)
  }
  L <- x / sqrt(bv)
  cc <- sqrt(2 / pi)
  Cn <- sqrt(2 * log(m)) / cc - (log(pi) + log(log(m))) / (2 * cc * sqrt(2 * log(m)))
  Sn <- 1 / (cc * sqrt(2 * log(m)))
  beta_star <- -log(-log(1 - alpha))
  stat <- (abs(L) - Cn) / Sn
  list(L = L, stat = stat, jump = !is.na(stat) & stat > beta_star, K = K, alpha = alpha)
}
lm_test <- lee_mykland(r, K = 16, alpha = 0.01)
d$jump <- lm_test$jump
jumps <- d[d$jump, c("date", "usdidr", "ret")]
jumps$ret <- round(jumps$ret, 3)
jumps$lm_stat <- round(lm_test$stat[d$jump], 2)
write.csv(jumps, file.path("output", "lm_jumps.csv"), row.names = FALSE)

# Rolling 60-day volatility (annualised, %)
roll_sd <- function(x, w) {
  out <- rep(NA_real_, length(x))
  for (t in w:length(x)) out[t] <- sd(x[(t - w + 1):t])
  out
}
d$vol60 <- roll_sd(r, 60) * sqrt(252)

# ---- 3. Figures -----------------------------------------------------------
message("\nMaking figures ...")

# Fig 1: three aligned panels sharing the time axis
shade <- geom_rect(data = episodes, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                   inherit.aes = FALSE, fill = COL["soft"], alpha = 0.18)
xscale <- scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = expansion(0.01))
lab_y <- max(d$usdidr) * 1.02
shade_log <- geom_rect(data = episodes, aes(xmin = start, xmax = end,
                                           ymin = min(d$usdidr) * 0.97, ymax = lab_y * 1.12),
                       inherit.aes = FALSE, fill = COL["soft"], alpha = 0.18)
p1a <- ggplot(d, aes(date, usdidr)) + shade_log +
  geom_line(colour = COL["main"], linewidth = 0.35) +
  geom_text(data = episodes, aes(x = start + (end - start) / 2, y = lab_y, label = label),
            inherit.aes = FALSE, size = 2.3, colour = "grey35", vjust = 0, lineheight = 0.85) +
  scale_y_log10(labels = label_comma(), expand = expansion(c(0.02, 0.12))) + xscale +
  labs(x = NULL, y = "IDR per USD (log scale)") + theme_paper +
  theme(axis.text.x = element_blank())
p1b <- ggplot(d, aes(date, ret)) + shade +
  geom_segment(aes(xend = date, y = 0, yend = ret), colour = COL["grey"], linewidth = 0.2) +
  geom_segment(data = d[d$jump, ], aes(xend = date, y = 0, yend = ret),
               colour = COL["accent"], linewidth = 0.4) +
  geom_point(data = d[d$jump, ], aes(size = abs(ret)), colour = COL["accent"], alpha = 0.85) +
  scale_size_area(max_size = 3, guide = "none") + xscale +
  labs(x = NULL, y = "Daily return (%)") + theme_paper +
  theme(axis.text.x = element_blank())
p1c <- ggplot(d, aes(date, vol60)) + shade +
  geom_area(fill = COL["main"], alpha = 0.25, na.rm = TRUE) +
  geom_line(colour = COL["main"], linewidth = 0.35, na.rm = TRUE) + xscale +
  labs(x = NULL, y = "60-day volatility\n(annualised, %)") + theme_paper
fig1 <- p1a / p1b / p1c + plot_layout(heights = c(1.1, 1, 0.8)) +
  plot_annotation(tag_levels = "a",
                  caption = paste0("Red: days flagged by the Lee-Mykland test (K = ", lm_test$K,
                                   ", alpha = ", lm_test$alpha, "). Shaded: approximate stress episodes."))
save_fig(fig1, "fig1_series_jumps_volatility", 7.5, 7)

# Fig 2: tails
z <- (r - mean(r)) / sd(r)
# Empirical density from a fine histogram (bins with at least one observation)
brk <- seq(floor(min(z)) - 0.25, ceiling(max(z)) + 0.25, by = 0.25)
h <- hist(z, breaks = brk, plot = FALSE)
emp <- data.frame(z = h$mids, dens = h$density, count = h$counts)
emp <- emp[emp$count > 0, ]
gx <- seq(-5, 5, by = 0.02)
gau <- data.frame(z = gx, dens = dnorm(gx))
p2a <- ggplot() +
  geom_line(data = gau, aes(z, dens, colour = "Gaussian"), linewidth = 0.5, linetype = 2) +
  geom_point(data = emp, aes(z, dens, colour = "Empirical", size = count), alpha = 0.8) +
  scale_size_area(max_size = 2.2, guide = "none") +
  scale_y_log10(labels = label_number(drop0trailing = TRUE)) +
  coord_cartesian(ylim = c(min(emp$dens) * 0.7, 1)) +
  scale_colour_manual(values = c(Empirical = unname(COL["main"]), Gaussian = unname(COL["accent"])),
                      name = NULL) +
  labs(x = "Standardised return", y = "Density (log scale)") + theme_paper
qq <- data.frame(theo = qnorm(ppoints(n)), samp = sort(z))
# 95% pointwise band for order statistics of a Gaussian sample
qq$lo <- qnorm(qbeta(0.025, seq_len(n), n + 1 - seq_len(n)))
qq$hi <- qnorm(qbeta(0.975, seq_len(n), n + 1 - seq_len(n)))
p2b <- ggplot(qq, aes(theo, samp)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = COL["accent"], alpha = 0.15) +
  geom_abline(slope = 1, intercept = 0, colour = COL["accent"], linetype = 2, linewidth = 0.4) +
  geom_point(colour = COL["main"], size = 0.5, alpha = 0.6) +
  labs(x = "Gaussian quantile", y = "Sample quantile") + theme_paper
fig2 <- p2a + p2b + plot_annotation(tag_levels = "a")
save_fig(fig2, "fig2_tails", 7.5, 3.4)

# Fig 3: calendar heatmap of mean absolute return
d$year <- as.integer(format(d$date, "%Y"))
d$month <- factor(month.abb[as.integer(format(d$date, "%m"))], levels = month.abb)
cal <- aggregate(abs(ret) ~ year + month, data = d, FUN = mean)
names(cal)[3] <- "mabs"
p3 <- ggplot(cal, aes(month, factor(year, levels = rev(sort(unique(year)))), fill = mabs)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  scale_fill_gradientn(colours = c("#F4F6F8", "#8FA9C8", "#1F3A5F", "#C8102E"),
                       trans = "sqrt", name = "Mean |return| (%)") +
  labs(x = NULL, y = NULL) + theme_paper +
  theme(panel.grid = element_blank(), legend.key.width = unit(1.4, "cm"))
save_fig(p3, "fig3_calendar_heatmap", 6, 6.5)

# Fig 4: autocorrelations (lollipop)
acf_df <- function(x, name, lag.max = 40) {
  a <- acf(x, lag.max = lag.max, plot = FALSE)$acf[-1]
  data.frame(lag = seq_along(a), acf = a, series = name)
}
labs_ac <- c("Returns, r", "Absolute returns, |r|", "Squared returns, r^2")
ac <- rbind(acf_df(r, labs_ac[1]), acf_df(abs(r), labs_ac[2]), acf_df(r^2, labs_ac[3]))
ac$series <- factor(ac$series, levels = labs_ac)
band <- 1.96 / sqrt(n)
p4 <- ggplot(ac, aes(lag, acf)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -band, ymax = band, fill = COL["soft"], alpha = 0.3) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_segment(aes(xend = lag, y = 0, yend = acf), colour = COL["main"], linewidth = 0.4) +
  geom_point(colour = COL["main"], size = 1) +
  facet_wrap(~ series, nrow = 1) +
  labs(x = "Lag (days)", y = "Autocorrelation") + theme_paper
save_fig(p4, "fig4_autocorrelation", 7.5, 2.8)

# ---- Summary -------------------------------------------------------------
cat("\n================ STYLIZED FACTS ================\n")
print(table1, row.names = FALSE)
cat("\nLee-Mykland test (K =", lm_test$K, ", alpha =", lm_test$alpha, "): ",
    nrow(jumps), "jump days out of", sum(!is.na(lm_test$stat)), "tested\n")
cat("Jump days per year:\n"); print(table(format(jumps$date, "%Y")))
cat("\nLargest flagged jumps:\n")
print(head(jumps[order(-abs(jumps$ret)), ], 10), row.names = FALSE)
cat("\nTables in output/, figures in figures/\n")
