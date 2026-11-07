# =====================================================================
# 02_check_data.R
# Data-quality check for the ECB USD/IDR series.
#
# Every day with a large ECB move (|return| > 2%) is compared with a
# second, independent source (Yahoo Finance, ticker IDR=X). Because the
# two sources are recorded at different times of day, the Yahoo returns
# on days t-1, t and t+1 are all considered.
#
# A day is flagged as SUSPECT when
#   (a) no Yahoo return in the window t-1..t+1 has the same sign and at
#       least half the size of the ECB move, or
#   (b) the ECB move is reversed by at least 70% on the next ECB day
#       (a typical pattern of a wrong price).
#
# Output (folder output/): check_large_moves.csv
# Requires: data/panel_daily.rds from 01_download_data.R
# =====================================================================

suppressPackageStartupMessages(library(quantmod))

if (basename(getwd()) != "rupiah-jump-sv" && dir.exists("~/rupiah-jump-sv")) {
  setwd("~/rupiah-jump-sv")
}
if (!file.exists(file.path("data", "panel_daily.rds"))) {
  stop("data/panel_daily.rds not found. Run 01_download_data.R first.")
}
dir.create("output", showWarnings = FALSE)

THRESH <- 2      # large move, in %
panel <- readRDS(file.path("data", "panel_daily.rds"))

# ---- Second source: Yahoo Finance --------------------------------------
yh <- NULL
for (k in 1:3) {
  yh <- tryCatch(getSymbols("IDR=X", src = "yahoo", from = min(panel$date) - 10,
                            auto.assign = FALSE), error = function(e) NULL)
  if (!is.null(yh)) break
  Sys.sleep(2 * k)
}
if (is.null(yh)) stop("Could not download Yahoo data (IDR=X). Run the script again later.")
yh <- na.omit(Cl(yh))
yh <- yh[as.numeric(yh) > 0]
yret <- data.frame(date = as.Date(index(yh))[-1],
                   ret  = 100 * diff(log(as.numeric(yh))))

# ---- Compare every large ECB move --------------------------------------
n <- nrow(panel)
big <- which(abs(panel$ret) > THRESH)
res <- data.frame(date = panel$date[big], usdidr = round(panel$usdidr[big], 1),
                  ecb_ret = round(panel$ret[big], 2),
                  yahoo_tm1 = NA_real_, yahoo_t = NA_real_, yahoo_tp1 = NA_real_,
                  ecb_next = NA_real_, flag = "")
for (j in seq_along(big)) {
  d <- res$date[j]
  # Yahoo returns on the previous, same and next Yahoo trading day (within 4 days)
  prev <- yret[yret$date <  d & yret$date >= d - 4, ]
  same <- yret[yret$date == d, ]
  nxt  <- yret[yret$date >  d & yret$date <= d + 4, ]
  w <- c(if (nrow(prev)) prev$ret[nrow(prev)] else NA_real_,
         if (nrow(same)) same$ret[1] else NA_real_,
         if (nrow(nxt))  nxt$ret[1]  else NA_real_)
  res[j, c("yahoo_tm1", "yahoo_t", "yahoo_tp1")] <- round(w, 2)
  i <- big[j]
  if (i < n) res$ecb_next[j] <- round(panel$ret[i + 1], 2)

  confirmed <- any(!is.na(w) & sign(w) == sign(res$ecb_ret[j]) &
                     abs(w) >= 0.5 * abs(res$ecb_ret[j]))
  reversal <- !is.na(res$ecb_next[j]) && sign(res$ecb_next[j]) == -sign(res$ecb_ret[j]) &&
    abs(res$ecb_next[j]) >= 0.7 * abs(res$ecb_ret[j])
  res$flag[j] <- paste(c(if (!confirmed) "not confirmed by Yahoo",
                         if (reversal) "reversed next day"), collapse = "; ")
}
write.csv(res, file.path("output", "check_large_moves.csv"), row.names = FALSE)

# ---- Summary ----------------------------------------------------------
cat("\n================ DATA CHECK ================\n")
cat("Large ECB moves (|ret| >", THRESH, "%):", nrow(res), "\n")
cat("Confirmed by Yahoo and not reversed :", sum(res$flag == ""), "\n")
cat("SUSPECT                              :", sum(res$flag != ""), "\n\n")
if (any(res$flag != "")) {
  cat("Suspect days (returns in %):\n")
  print(res[res$flag != "", ], row.names = FALSE)
}
cat("\nFull table saved in:", normalizePath(file.path("output", "check_large_moves.csv")), "\n")
