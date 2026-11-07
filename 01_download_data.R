# =====================================================================
# 01_download_data.R
# Downloads the raw daily data for the rupiah jump-SV paper and builds
# one daily panel aligned on ECB fixing dates.
#
#   USD/IDR  : ECB euro foreign exchange reference rates (IDR / USD cross)
#   VIXCLS   : CBOE VIX, close (FRED)
#   DGS10    : US 10-year Treasury yield, % (FRED)
#   DTWEXBGS : Nominal broad US dollar index (FRED, starts 2006)
#   DCOILBRENTEU : Brent crude oil price, USD/barrel (FRED)
#
# Output (folder data/):
#   raw_ecb_usdidr.csv, raw_fred_<series>.csv, panel_daily.csv, panel_daily.rds
# Run time: < 1 minute.
# =====================================================================

suppressPackageStartupMessages({
  library(quantmod)   # also loads xts and zoo
})

# ---- 0. Working folder ------------------------------------------------
if (basename(getwd()) != "rupiah-jump-sv" && dir.exists("~/rupiah-jump-sv")) {
  setwd("~/rupiah-jump-sv")
}
if (basename(getwd()) != "rupiah-jump-sv") {
  stop("Please open the folder 'rupiah-jump-sv' first (Session > Set Working Directory).")
}
dir.create("data", showWarnings = FALSE)
message("Working folder: ", normalizePath(getwd()))

START <- as.Date("2005-01-01")
ECB_URL <- "https://www.ecb.europa.eu/stats/eurofxref/eurofxref-hist.zip"
FRED_SERIES <- c("VIXCLS", "DGS10", "DTWEXBGS", "DCOILBRENTEU")

# Retry helper: tries an expression up to 3 times before giving up
with_retry <- function(expr_fun, what, tries = 3) {
  for (k in seq_len(tries)) {
    out <- tryCatch(expr_fun(), error = function(e) e)
    if (!inherits(out, "error")) return(out)
    message("  attempt ", k, " for ", what, " failed: ", conditionMessage(out))
    if (k < tries) Sys.sleep(2 * k)
  }
  NULL
}

# ---- 1. USD/IDR from the ECB ------------------------------------------
message("\n[1/3] Downloading ECB reference rates ...")
ecb_raw <- with_retry(function() {
  tf <- tempfile(fileext = ".zip")
  download.file(ECB_URL, tf, mode = "wb", quiet = TRUE)
  read.csv(unz(tf, "eurofxref-hist.csv"), na.strings = "N/A")
}, "ECB file")
if (is.null(ecb_raw)) stop("Could not download the ECB file. Check the internet connection and run again.")
if (!all(c("Date", "USD", "IDR") %in% names(ecb_raw))) {
  stop("The ECB file does not contain the columns Date, USD and IDR. The file format may have changed.")
}

ecb <- data.frame(date = as.Date(ecb_raw$Date), eur_usd = ecb_raw$USD, eur_idr = ecb_raw$IDR)
ecb <- ecb[!is.na(ecb$date) & !is.na(ecb$eur_usd) & !is.na(ecb$eur_idr), ]
ecb <- ecb[order(ecb$date), ]
ecb$usdidr <- ecb$eur_idr / ecb$eur_usd          # rupiah per US dollar
ecb <- ecb[ecb$date >= START, ]
if (anyDuplicated(ecb$date)) stop("Duplicated dates in the ECB file.")
if (any(ecb$usdidr <= 0)) stop("Non-positive USD/IDR values in the ECB file.")
write.csv(ecb, file.path("data", "raw_ecb_usdidr.csv"), row.names = FALSE)
message("  USD/IDR: ", nrow(ecb), " days, ", format(min(ecb$date)), " to ", format(max(ecb$date)))

# ---- 2. Covariates from FRED ------------------------------------------
message("\n[2/3] Downloading FRED series ...")
fred <- list()
for (s in FRED_SERIES) {
  x <- with_retry(function() getSymbols(s, src = "FRED", auto.assign = FALSE), s)
  if (is.null(x)) { message("  ! ", s, " could not be downloaded"); next }
  df <- data.frame(date = as.Date(index(x)), value = as.numeric(x[, 1]))
  df <- df[df$date >= START, ]
  write.csv(df, file.path("data", paste0("raw_fred_", s, ".csv")), row.names = FALSE)
  fred[[s]] <- df
  message("  ", s, ": ", sum(!is.na(df$value)), " values, ",
          format(min(df$date)), " to ", format(max(df$date)))
}
missing_series <- setdiff(FRED_SERIES, names(fred))

# ---- 3. Daily panel on ECB dates --------------------------------------
message("\n[3/3] Building the daily panel ...")
panel <- ecb[, c("date", "usdidr")]
panel$ret <- c(NA, 100 * diff(log(panel$usdidr)))  # daily log return in %
for (s in names(fred)) {
  # Last available value on or before each ECB date.
  # NOTE: US data for day t are recorded after the ECB fixing (14:15 CET),
  # so covariates must be lagged by one day when they enter the model.
  v <- fred[[s]][!is.na(fred[[s]]$value), ]
  v <- v[order(v$date), ]
  k <- findInterval(as.numeric(panel$date), as.numeric(v$date))
  panel[[tolower(s)]] <- ifelse(k > 0, v$value[pmax(k, 1)], NA)
}
saveRDS(panel, file.path("data", "panel_daily.rds"))
write.csv(panel, file.path("data", "panel_daily.csv"), row.names = FALSE)

# ---- Summary ----------------------------------------------------------
cat("\n================ SUMMARY ================\n")
cat("Panel rows:", nrow(panel), "| from", format(min(panel$date)), "to", format(max(panel$date)), "\n")
cat("Missing values per column:\n"); print(colSums(is.na(panel)))
cat("\nDays with |return| > 2%:", sum(abs(panel$ret) > 2, na.rm = TRUE), "\n")
top <- panel[order(-abs(panel$ret)), c("date", "usdidr", "ret")][1:10, ]
top$ret <- round(top$ret, 2); top$usdidr <- round(top$usdidr, 1)
cat("\nTen largest daily moves (ret in %):\n"); print(top, row.names = FALSE)
if (length(missing_series)) {
  cat("\nWARNING: these FRED series failed and are not in the panel:",
      paste(missing_series, collapse = ", "), "\nRun the script again later.\n")
}
cat("\nFiles saved in:", normalizePath("data"), "\n")
