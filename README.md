# Covariate-driven jumps in the Indonesian rupiah

R and C++ code for the paper

> Ardiansyah, F. *Covariate-driven jumps in the Indonesian rupiah: a Bayesian stochastic volatility model of currency crash risk.* Submitted to *Computational Statistics & Data Analysis*.

The paper proposes SVJ-X, a stochastic volatility model with jumps in which the daily jump probability is a probit function of observed covariates. It applies the model to the Bank Indonesia USD/IDR reference rate (2006–2026). The code covers the Gibbs sampler, a particle filter for density forecasts, a GARCH-t benchmark, the out-of-sample evaluation, a simulation study and robustness checks.

## Requirements

- R (tested with R 4.5.1) and a C++ compiler (Rtools on Windows)
- R packages: `Rcpp`, `quantmod`, `readxl`, `ggplot2`, `patchwork`, `scales`, `ggridges`, `posterior`

Install `quantmod` and `readxl` before running scripts `01`–`03`:

```r
install.packages(c("quantmod", "readxl"))
```

Scripts `04`–`11` install any other missing packages automatically.

## Data

All data are free:

| Source | Series | How it is obtained |
|---|---|---|
| Bank Indonesia | Buying and selling rate of the US dollar (*Kurs Transaksi BI*) | Download manually as Excel from <https://www.bi.go.id/id/statistik/informasi-kurs/transaksi-bi/default.aspx> and save it in `data/` with a file name starting with `Kurs` and containing `USD` |
| FRED | `VIXCLS`, `DGS10`, `DTWEXBGS`, `DCOILBRENTEU` | Downloaded by `01_download_data.R` |
| European Central Bank | Euro reference rates (for the USD/IDR cross rate) | Downloaded by `01_download_data.R` |

## Files

| File | Purpose | Run time |
|---|---|---|
| `01_download_data.R` | Downloads the FRED and ECB data | < 1 min |
| `02_check_data.R` | Data-quality check of the ECB series | < 1 min |
| `03_build_panel.R` | Builds the daily panel from the Bank Indonesia file and the covariates | < 1 min |
| `04_stylized_facts.R` | Table 1, Lee–Mykland jump test, Figures 1–2 and appendix figures | < 1 min |
| `05_sampler_check.R` | Checks the sampler on simulated data | 3–6 min |
| `06_fit_models.R` | Fits SV, SVJ and SVJ-X (two chains) to the full sample | 15–70 min |
| `07_compare_models.R` | In-sample log scores, PIT and VaR backtests, GARCH-t benchmark | 5–15 min |
| `08_out_of_sample.R` | Estimation on 2006–2018, forecasts for 2019–2026 | 30–60 min |
| `09_simulation_study.R` | Simulation study, 4 scenarios × 50 replicates (resumable) | several hours |
| `10_robustness.R` | Nine robustness variants (resumable) | 0.5–3 h |
| `11_jump_triggers.R` | Jump probability as a function of the triggers (needs `06`) | < 1 min |
| `svjx_functions.R` | Simulation, priors and the Gibbs sampler | |
| `svjx_sampler.cpp` | Block sampler for the log-volatility path (tridiagonal precision) | |
| `svjx_filter.cpp` | Particle filter with the jumps integrated out | |
| `eval_functions.R` | GARCH-t, Diebold–Mariano, Kupiec and Christoffersen tests | |

## How to reproduce

1. Put all files in one folder (the scripts look for a folder named `rupiah-jump-sv`, or use the current working directory).
2. Download the Bank Indonesia Excel file into `data/` (see above).
3. Run the scripts in order, `01` to `11`. Each script writes its tables to `output/` and its figures (PNG and PDF) to `figures/`.

The scripts set their random seeds, so the results are reproducible on the same machine. Scripts `09` and `10` save each replicate or variant as soon as it is finished; if they are interrupted, run them again and they continue where they stopped.

## License

MIT License (see `LICENSE`).
