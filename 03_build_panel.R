# =====================================================================
# 03_build_panel.R
# Builds the final daily panel used in the analysis.
#
#   Primary FX series : Bank Indonesia transaction rate for USD
#                       (Kurs Transaksi BI), mid rate = (sell + buy) / 2.
#                       Downloaded by hand from
#                       https://www.bi.go.id/id/statistik/informasi-kurs/transaksi-bi/default.aspx
#                       (Time Series, USD, 01/01/2005 - latest) and saved
#                       unchanged in data/ (file name starts with "Kurs" and
#                       contains "USD", extension .xlsx).
#   Cross-check       : ECB reference rate (data/raw_ecb_usdidr.csv, script 01).
#   Covariates        : FRED series (data/raw_fred_*.csv, script 01).
#
# Timing: the BI rate is fixed in the Jakarta morning, before US markets
# open. Each covariate therefore takes its last value strictly BEFORE the
# BI date (previous US trading day), so there is no look-ahead.
#
# Output (folder data/): panel_bi.rds, panel_bi.csv
#        (folder output/): bi_vs_ecb_large_moves.csv
# =====================================================================

suppressPackageStartupMessages(library(readxl))

if (basename(getwd()) != "rupiah-jump-sv" && dir.exists("~/rupiah-jump-sv")) {
  setwd("~/rupiah-jump-sv")
}
if (basename(getwd()) != "rupiah-jump-sv") {
  stop("Please open the folder 'rupiah-jump-sv' first (Session > Set Working Directory).")
}
dir.create("output", showWarnings = FALSE)

# ---- 1. Read the BI file ------------------------------------------------
bi_file <- list.files("data", pattern = "^Kurs.*USD.*\\.xlsx$", full.names = TRUE)
if (length(bi_file) == 0) stop("No BI file found in data/. Expected a file like 'Kurs Transaksi USD.xlsx'.")
if (length(bi_file) > 1) stop("More than one BI file in data/: ", paste(basename(bi_file), collapse = ", "),
                              ". Keep only one.")
message("BI file: ", basename(bi_file))

raw <- suppressMessages(read_excel(bi_file, col_names = FALSE, col_types = "text"))
hdr <- which(trimws(raw[[1]]) == "NO")[1]
if (is.na(hdr)) stop("Header row ('NO', 'Nilai', 'Kurs Jual', 'Kurs Beli', 'Tanggal') not found.")
names(raw) <- trimws(as.character(unlist(raw[hdr, ])))
need <- c("Nilai", "Kurs Jual", "Kurs Beli", "Tanggal")
if (!all(need %in% names(raw))) stop("Missing columns in the BI file: ",
                                     paste(setdiff(need, names(raw)), collapse = ", "))
raw <- raw[(hdr + 1):nrow(raw), ]
raw <- raw[!is.na(raw$Tanggal) & trimws(raw$Tanggal) != "", ]

# Indonesian number format: "18.230,70" -> 18230.70
id_num <- function(s) as.numeric(gsub(",", ".", gsub(".", "", trimws(s), fixed = TRUE), fixed = TRUE))
# Dates come as "11/5/2026 12:00:00 AM" (month/day/year)
bi_date <- function(s) as.Date(sub(" .*$", "", trimws(s)), format = "%m/%d/%Y")

bi <- data.frame(date = bi_date(raw$Tanggal), unit = id_num(raw$Nilai),
                 sell = id_num(raw$`Kurs Jual`), buy = id_num(raw$`Kurs Beli`))
bad <- is.na(bi$date) | is.na(bi$sell) | is.na(bi$buy) | is.na(bi$unit)
if (any(bad)) stop(sum(bad), " rows of the BI file could not be read. First one: row ",
                   which(bad)[1] + hdr)
if (any(bi$unit != 1)) stop("Unexpected 'Nilai' (unit) values in the BI file.")
bi <- bi[order(bi$date), ]
if (anyDuplicated(bi$date)) stop("Duplicated dates in the BI file.")
if (any(bi$sell < bi$buy)) stop("Sell rate below buy rate in the BI file.")

bi$usdidr <- (bi$sell + bi$buy) / 2                  # mid rate, rupiah per USD
bi$ret    <- c(NA, 100 * diff(log(bi$usdidr)))       # daily log return in %
bi$gap    <- c(NA, as.numeric(diff(bi$date)))        # calendar days since previous fixing
message("BI USD/IDR: ", nrow(bi), " days, ", format(min(bi$date)), " to ", format(max(bi$date)))

# ---- 2. Covariates (lagged: last value strictly before the BI date) ------
panel <- bi[, c("date", "usdidr", "ret", "gap")]
fred_files <- list.files("data", pattern = "^raw_fred_.*\\.csv$", full.names = TRUE)
if (length(fred_files) == 0) stop("No FRED files in data/. Run 01_download_data.R first.")
for (f in fred_files) {
  v <- read.csv(f, stringsAsFactors = FALSE)
  v$date <- as.Date(v$date)
  v <- v[!is.na(v$value), ]
  v <- v[order(v$date), ]
  k <- findInterval(as.numeric(panel$date) - 1, as.numeric(v$date))   # strictly before
  col <- paste0(tolower(sub("^raw_fred_(.*)\\.csv$", "\\1", basename(f))), "_lag")
  panel[[col]] <- ifelse(k > 0, v$value[pmax(k, 1)], NA)
}

# ---- 3. Cross-check with the ECB rate -----------------------------------
ecb_file <- file.path("data", "raw_ecb_usdidr.csv")
if (file.exists(ecb_file)) {
  ecb <- read.csv(ecb_file, stringsAsFactors = FALSE)
  ecb$date <- as.Date(ecb$date)
  ecb <- ecb[order(ecb$date), ]
  ecb$ecb_ret <- c(NA, 100 * diff(log(ecb$usdidr)))
  m <- match(panel$date, ecb$date)
  panel$ecb_usdidr <- ecb$usdidr[m]
  panel$level_gap_pct <- 100 * (panel$ecb_usdidr / panel$usdidr - 1)
  big <- which(abs(panel$ret) > 2)
  chk <- data.frame(date = panel$date[big], bi_ret = round(panel$ret[big], 2),
                    ecb_ret_same_day = round(ecb$ecb_ret[m[big]], 2))
  write.csv(chk, file.path("output", "bi_vs_ecb_large_moves.csv"), row.names = FALSE)
} else {
  message("ECB file not found; cross-check skipped.")
}

saveRDS(panel, file.path("data", "panel_bi.rds"))
write.csv(panel, file.path("data", "panel_bi.csv"), row.names = FALSE)

# ---- Summary ------------------------------------------------------------
cat("\n================ PANEL (BI) ================\n")
cat("Rows:", nrow(panel), "| from", format(min(panel$date)), "to", format(max(panel$date)), "\n")
cat("Missing values per column:\n"); print(colSums(is.na(panel)))
cat("\nDays with |return| > 2%:", sum(abs(panel$ret) > 2, na.rm = TRUE), "\n")
cat("Days with zero return   :", sum(panel$ret == 0, na.rm = TRUE), "\n")
cat("Gaps between fixings (calendar days):\n"); print(table(panel$gap))
top <- panel[order(-abs(panel$ret)), c("date", "usdidr", "ret")][1:10, ]
top$ret <- round(top$ret, 2); top$usdidr <- round(top$usdidr, 1)
cat("\nTen largest daily moves (ret in %):\n"); print(top, row.names = FALSE)
if ("level_gap_pct" %in% names(panel)) {
  cat("\nBI vs ECB level difference (%), all common days:\n")
  print(round(quantile(panel$level_gap_pct, c(0, .01, .5, .99, 1), na.rm = TRUE), 2))
}
cat("\nFiles saved in:", normalizePath("data"), "\n")
