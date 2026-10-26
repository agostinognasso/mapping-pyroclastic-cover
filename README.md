# Pyroclastic cover thickness in Campania: spatial cross-validation and sampling-requirement curves

R code, data and results for the article

> Gnasso, A., Matano, F., Ebrahimi, P., Scepi, G. *Spatial cross-validation and
> sampling-requirement curves for mapping pyroclastic cover thickness.*
> Submitted to *Spatial Statistics*.

The analysis maps the thickness of loose pyroclastic deposits on the slopes of
Campania (southern Italy) from 7,293 field measurements and 16 covariates. It
asks under which data conditions such a map can be relied on:

- spatial cross-validation designs (block, buffered block, nearest-neighbour
  distance matching, leave-one-municipality-out, leave-one-sector-out), sized on
  the residual variogram;
- a sampling-requirement curve, which links the accuracy of the map to the
  amount of in-fill field data;
- a comparison of random forests (RF, RFsp) with ordinary and regression
  kriging, the empirical thickness models, XGBoost and a spatial GAM, on
  identical folds;
- 90% prediction intervals from quantile regression forests, assessed with the
  Winkler score and with coverage within classes;
- an explainable ensemble tree (e2tree) that summarises the forest as decision
  rules; it is never used for prediction;
- the sampled zone and the area of applicability of the delivered map.

## Repository layout

| Path | Content |
|---|---|
| `data/dataset.csv` | PT-Cam-2 field measurements with the covariates (semicolon-separated) |
| `data/prediction_grid.csv` | 1 km prediction grid with the same covariates |
| `R/` | numbered scripts, `00_setup.R` to `39_coverage_figure.R` |
| `run_all.R` | runs the whole pipeline in the right order |
| `output/tables/` | CSV tables written by the scripts (`T*.csv`) |
| `output/figs/` | figures (PDF and PNG) |
| `output/geotiff/` | georeferenced rasters of the delivered map (EPSG:32633) |
| `output/tex/` | LaTeX tables generated for the article |
| `output/models/` | fitted models and intermediate objects (`.rds`); not tracked, see below |

## Data

`data/dataset.csv` comes from the PT-Cam-2 database (Matano et al., 2026), CC-BY
4.0, available on Zenodo (DOI 10.5281/zenodo.23425763). PT-Cam-2 is an upgraded
version of the database published in 2023 (DOI 10.5281/zenodo.8399487). The
file here holds the 8,109 records joined to the covariates. Column names follow
the published data, including the misspelt `localitiy`, and the `note` and
`type` columns keep their Italian values (`parziale` for partial-thickness
measurements, `fioretto` for probing tests).

The modelling set has 7,293 records: the 22 records of Forio and Procida
(islands) and the 794 partial-thickness records are excluded
(`R/01_load_clean.R`). Coordinates are UTM zone 33N (EPSG:32633), in metres.

## Requirements

R 4.6.0 was used. Packages and the versions used:

`ranger` 0.18.0, `randomForest` 4.7.1.2, `e2tree` 1.2.0.9000, `rsample` 1.3.2,
`gstat` 2.1.6, `sp` 2.2.1, `xgboost` 3.2.1.1, `mgcv` 1.9.4, `terra` 1.9.25,
`sf` 1.1.1, `CAST` 1.1.2, `classInt` 0.4.11, `rpart` 4.1.27,
`rpart.plot` 3.1.4, `ggplot2` 4.0.3, `dplyr` 1.2.1, `tidyr` 1.3.2,
`future` 1.70.0, `viridisLite` 0.4.3.

`e2tree` is installed from GitHub (<https://github.com/massimoaria/e2tree>); the
development version 1.2.0.9000 was used. The others are on CRAN.

```r
install.packages(c("ranger", "randomForest", "rsample", "gstat", "sp", "xgboost",
                   "mgcv", "terra", "sf", "CAST", "classInt", "rpart",
                   "rpart.plot", "ggplot2", "dplyr", "tidyr", "future",
                   "viridisLite", "remotes"))
remotes::install_github("massimoaria/e2tree")
```

## Running the analysis

Run everything from the repository root:

```bash
Rscript run_all.R
```

This takes roughly 95 minutes on a laptop. Scripts can also be run one at a time
(`Rscript R/11_scenarios.R`), as long as the scripts they read from have run
first; the order and the dependencies are commented in `run_all.R`. The random
seed is fixed in `R/00_setup.R`.

By default, models, tables and GeoTIFFs go to `output/`, figures to
`output/figs/` and LaTeX tables to `output/tex/`. To run without overwriting
existing results, point the pipeline elsewhere:

```bash
PTCAM_OUT=output_new PTCAM_PAPER=output_new/tex PTCAM_FIGS=output_new/figs Rscript run_all.R
```

### Model objects

The fitted objects in `output/models/` are large (about 165 MB) and are not
tracked in Git. Every script that needs them reads them from there, so a fresh
clone needs `Rscript run_all.R` once before running single scripts. The tables,
figures, GeoTIFFs and LaTeX tables in `output/` are tracked.

## Scripts

| Script | Purpose |
|---|---|
| `00_setup.R` | paths, seed, libraries, shared helpers (maps, figure export) |
| `01`-`02` | data cleaning, exploratory analysis |
| `03`, `03b`, `03c` | random forest, explainable forest, accuracy-fidelity frontier |
| `04`, `04b`, `08`, `18`, `28`, `32` | explainable ensemble tree: fidelity, local explanations, rules, CART control |
| `05`, `09`, `10` | maps, quantile forest, operational model |
| `07`, `11`, `14`, `15`, `33` | validation designs: block, buffered, NNDM, leave-one-out; their geometry |
| `12`, `16`, `20`, `37` | comparisons: empirical models, recalibrated empirical models, XGBoost and GAM, regression kriging |
| `17`, `39` | coverage and Winkler score of the prediction intervals; coverage within classes |
| `19`, `21`, `36` | robustness: removing z0, measurement type, error in z0 |
| `22`, `30` | operating regimes and the registry that generates the LaTeX tables |
| `31`, `35` | sampling-requirement curve in absolute units; thresholds in cm and block means |
| `34`, `38` | southern Cilento check; sampled zone and area of applicability |
| `13`, `24`, `25`, `26`, `29` | figures, GeoTIFF export, banding and smoothing checks, supplement tables |

Table and figure files are numbered by the script that writes them (`T21_...`,
`F52_...`); the registry in `R/30_scenario_registry.R` fixes the acronyms of the
validation designs (R5, B5, B15, B15b, B40, NNDM, LOMO, LOSO, OP10, OP25, OP50).

## Citation

Please cite the article above and the data:

> Matano, F., Angelino, A., Ebrahimi, P., Grimaldi, G. M. (2026). *Upgraded
> database of pyroclastic cover deposit thickness measurements (PT-Cam-2) in
> peri-volcanic areas of Campania (Italy)*. Zenodo.
> https://doi.org/10.5281/zenodo.23425763

## License

Code: MIT, see `LICENSE`. Data: CC-BY 4.0, as published on Zenodo.
