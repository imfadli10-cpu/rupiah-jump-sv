# =====================================================================
# 06_fit_models.R
# Fits the three Bayesian models to daily USD/IDR returns (BI mid rate):
#
#   SV    : stochastic volatility, no jumps
#   SVJ   : SV + jumps with a constant jump probability
#   SVJX  : SV + jumps whose probability depends on lagged covariates
#           (two independent chains, used to check convergence)
#
# Jump-probability covariates (all known before the BI fixing, standardised):
#   log_vix      log of the VIX on the previous US trading day
#   abs_d_ust10  absolute change in the US 10-year yield (basis points)
#   abs_d_usd    absolute % change in the broad US dollar index
#   long_break   1 if the previous fixing was 4 or more calendar days ago
#                (holiday longer than a normal weekend)
#
# Sample: from the first day on which all covariates exist (Jan 2006).
# Output: output/fit_<model>.rds, output/table2_parameters.csv,
#         output/convergence_svjx.csv
# Run time: roughly 15-45 minutes depending on the computer (models run
# in parallel on several CPU cores when available).
# =====================================================================

for (p in c("Rcpp", "posterior")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
if (basename(getwd()) != "rupiah-jump-sv" && dir.exists("~/rupiah-jump-sv")) setwd("~/rupiah-jump-sv")
for (f in c("svjx_sampler.cpp", "svjx_functions.R", file.path("data", "panel_bi.rds"))) {
  if (!file.exists(f)) stop("File ", f, " not found in ", getwd(), ".")
}
dir.create("output", showWarnings = FALSE)

N_ITER <- 35000; BURN <- 5000; THIN <- 5

# ---- Data and covariates ----------------------------------------------------
panel <- readRDS(file.path("data", "panel_bi.rds"))
panel <- panel[order(panel$date), ]
panel$abs_d_ust10 <- abs(c(NA, diff(panel$dgs10_lag))) * 100            # basis points
panel$abs_d_usd   <- abs(c(NA, 100 * diff(log(panel$dtwexbgs_lag))))    # %
panel$log_vix     <- log(panel$vixcls_lag)
panel$long_break  <- as.numeric(panel$gap >= 4)

covs <- c("log_vix", "abs_d_ust10", "abs_d_usd", "long_break")
ok <- complete.cases(panel[, c("ret", covs)])
first_ok <- which(ok)[1]
d <- panel[first_ok:nrow(panel), ]
if (any(!complete.cases(d[, c("ret", covs)]))) {
  # isolated gaps after the start (e.g. a missing US value) are carried forward
  for (v in covs) {
    x <- d[[v]]
    for (i in seq_along(x)[-1]) if (is.na(x[i])) x[i] <- x[i - 1]
    d[[v]] <- x
  }
}
stopifnot(!anyNA(d[, c("ret", covs)]))

std <- function(x) (x - mean(x)) / sd(x)
X <- cbind("(Intercept)" = 1, sapply(d[, covs], std))
y <- d$ret
cov_scale <- data.frame(covariate = covs, mean = sapply(d[, covs], mean), sd = sapply(d[, covs], sd))
saveRDS(list(date = d$date, y = y, X = X, cov_scale = cov_scale), file.path("output", "model_data.rds"))
message("Sample: ", format(min(d$date)), " to ", format(max(d$date)), ", n = ", length(y))

# ---- Fit (in parallel where possible) ----------------------------------------
jobs <- list(list(model = "SV",   seed = 11, X = NULL),
             list(model = "SVJ",  seed = 12, X = NULL),
             list(model = "SVJX", seed = 13, X = X),
             list(model = "SVJX", seed = 14, X = X))

run_job <- function(job, y, n_iter, burn, thin, wd) {
  setwd(wd)
  if (!exists("draw_h_tridiag")) Rcpp::sourceCpp("svjx_sampler.cpp")
  source("svjx_functions.R")
  fit_svjx(y, job$X, model = job$model, n_iter = n_iter, burn = burn, thin = thin,
           seed = job$seed, verbose = FALSE)
}

n_workers <- max(1, min(length(jobs), parallel::detectCores(logical = FALSE) - 1, na.rm = TRUE))
message("Fitting ", length(jobs), " chains on ", n_workers, " CPU core(s). Please wait ...")
t0 <- proc.time()[3]
if (n_workers > 1) {
  cl <- parallel::makeCluster(n_workers)
  fits <- tryCatch(
    parallel::parLapplyLB(cl, jobs, run_job, y = y, n_iter = N_ITER, burn = BURN,
                          thin = THIN, wd = getwd()),
    finally = parallel::stopCluster(cl))
} else {
  fits <- lapply(jobs, run_job, y = y, n_iter = N_ITER, burn = BURN, thin = THIN, wd = getwd())
}
message("Done in ", round((proc.time()[3] - t0) / 60, 1), " minutes.")

names(fits) <- c("SV", "SVJ", "SVJX_chain1", "SVJX_chain2")
for (nm in names(fits)) {
  fits[[nm]]$date <- d$date
  saveRDS(fits[[nm]], file.path("output", paste0("fit_", nm, ".rds")), compress = "xz")
}

# ---- Convergence of SVJX (two chains) -----------------------------------------
suppressPackageStartupMessages(library(posterior))
c1 <- fits$SVJX_chain1$draws; c2 <- fits$SVJX_chain2$draws
arr <- array(NA_real_, dim = c(nrow(c1), 2, ncol(c1)), dimnames = list(NULL, NULL, colnames(c1)))
arr[, 1, ] <- c1                                   # [iteration, chain, parameter]
arr[, 2, ] <- c2
dr <- posterior::as_draws_array(arr)
conv <- posterior::summarise_draws(dr, "mean", "sd", ~quantile(.x, c(0.025, 0.975)),
                                   "rhat", "ess_bulk", "ess_tail")
write.csv(conv, file.path("output", "convergence_svjx.csv"), row.names = FALSE)

# ---- Parameter table ------------------------------------------------------------
summ <- function(dm, model) {
  q <- t(apply(dm, 2, quantile, c(0.5, 0.025, 0.975)))
  data.frame(model = model, parameter = rownames(q), median = q[, 1],
             lower95 = q[, 2], upper95 = q[, 3], row.names = NULL)
}
svjx_all <- rbind(fits$SVJX_chain1$draws, fits$SVJX_chain2$draws)
tab2 <- rbind(summ(fits$SV$draws, "SV"), summ(fits$SVJ$draws, "SVJ"), summ(svjx_all, "SVJX"))
tab2[, 3:5] <- round(tab2[, 3:5], 3)
write.csv(tab2, file.path("output", "table2_parameters.csv"), row.names = FALSE)

# ---- Summary ----------------------------------------------------------------------
cat("\n================ MODEL FITS ================\n")
cat("Sample:", format(min(d$date)), "to", format(max(d$date)), "| n =", length(y), "\n")
cat("Seconds per chain:", paste(names(fits), round(sapply(fits, `[[`, "seconds")), collapse = ", "), "\n\n")
print(tab2, row.names = FALSE)
cat("\nPosterior mean number of jump days: SVJ =", round(mean(fits$SVJ$n_jumps), 1),
    "| SVJX =", round(mean(c(fits$SVJX_chain1$n_jumps, fits$SVJX_chain2$n_jumps)), 1), "\n")
cat("\nConvergence of SVJX (R-hat should be below 1.01; ESS above 400):\n")
print(conv[, c("variable", "rhat", "ess_bulk", "ess_tail")], row.names = FALSE)
