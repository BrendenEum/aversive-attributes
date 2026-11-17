# Aversive Attributes — Does attention to something bad attract or repel choice?

Brenden Eum (Rotman, University of Toronto) · Stephanie M. Smith (Booth, University of Chicago)

Attention is known to bias choice toward whatever is currently being looked at, but that regularity was established almost entirely with stimuli whose values are positive. This project asks what attention does when the attended thing is *bad*.

Two mouse-tracking studies have participants choose between two options described by two attributes, where each attribute row is labelled in advance as always appetitive (**P**), always aversive (**N**), or possibly either (**B**). A value is visible only while the pointer is inside its box, so the hover record is the attention measure. Between subjects the value range is manipulated (**Narrow** vs. **Wide**); within subjects the trial type varies (**PN**, **PB**, **NB**).

The headline results:

- The pull of attention increases with the value of the attended attribute *within* every attribute type. Attending to a mildly aversive attribute still moves choice toward the option carrying it; attending to a severely aversive attribute moves choice away.
- Where that reversal occurs is not fixed — it shifts with the range of values the participant was trained on.
- Choice sensitivity to value is itself range-normalized within attribute type: the same objective value carries less weight when drawn from a wider range.
- A Bayesian sequential-sampling model over attributes, with a prior whose mean acts as the reference point, reproduces the range-dependent crossover.

After exclusions the analyzed samples are 981 participants / 28,646 trials in Study 1 (PAS3, 30 trials each) and 525 participants / 30,523 trials in Study 2 (PAS4, 60 trials each).

---

## Status

Manuscript drafted end to end. One `\textcolor{red}{...}` placeholder is left — the data/code availability URL in the Transparency statement. [docs/TODO.md](docs/TODO.md) is the worklist, but it was last revised 2026-08-10 and predates the Bayesian-model work — read it against the current `main.tex` before trusting it.

Naming note: the two studies in the paper are **Study 1 = PAS3** and **Study 2 = PAS4**. PAS1 and PAS2 were earlier pilots; their processed data is still in `data/processed/`, but no current notebook reads it.

Wording note: what the mouse-tracking record measures is now called a **dwell** throughout the notebook and the manuscript, not a fixation, and "last-fixated" is now "last-seen". The outputs were renamed with it — `dwell_properties.{pdf,tex}`, `dwelldur_{type,value,range}.tex`, `lastseen_bias.tex`, `netdwell_{dwell,range}.tex`. The old `fixation_properties.*`, `fixdur_*`, `lastfix_bias`, `netfix_*` files are still sitting in `paper/figures|tables|text_output/` but nothing writes or includes them any more; the chunk names inside `Main_Analysis.Rmd` (`fix-first`, `fix-duration`, …) and the column names in the data (`FixROI`, `FixDur`) also still use the old word.

---

## How the project fits together

```
data/processed/PAS{3,4}{,long}.csv
        │   code/preprocessing/merge_PAS3_PAS4_to_main_data.R
        ▼
data/processed/main_data.csv, main_data_long.csv
        │   code/utilities/build_boxlong.R   (build_pas(), build_boxlong())
        ├──────────────────────────────┐
        ▼                              ▼
notebooks/Main_Analysis.Rmd      code/fitting/, code/parameter_recovery/
  (paper figures/tables/text)      (Bayesian model fits → results/cache/bayesian_model/)
                                           │
                                           ▼
                                 notebooks/bayesian_model.Rmd
```

Everything is R (plus one C++ file via Rcpp). Paths resolve through [config/paths.yaml](config/paths.yaml) via `here::here()`, anchored by the `.here` marker at the repo root, so notebooks knit and scripts run from either the root or a subdirectory. `paths.yaml` hardcodes an absolute `root` (`/Users/brenden/Desktop/aversive-attributes`) — update it if you clone elsewhere.

`code/utilities/getAllUtilities.R` sources the whole `code/utilities/` directory and is the standard first line of every script and exploratory notebook: it supplies `paths`, the plot theme and scales (`MyPlotOptions.R`), the brms settings (`MyBrmsOptions.R`), the data builders (`build_boxlong.R`), and small helpers (`SE.R`, `getStatsFromPosterior.R`, `write.text.R`, `mute_summarize_warning.R`, `plot_parameter_recovery.R`).

---

## Preprocessing

```bash
Rscript code/preprocessing/merge_PAS3_PAS4_to_main_data.R
```

Combines the two studies into one trial-level file (`main_data.csv`) and one fixation-level file (`main_data_long.csv`), each with a `Study` column and a study-unique `Subject` id. Everything downstream reads these two files.

---

## Notebooks

| Notebook | Role |
|----------|------|
| `Main_Analysis.Rmd` | **The paper.** Produces every figure in `paper/figures/`, every table in `paper/tables/`, and the in-text statistics in `paper/text_output/` |
| `bayesian_model.Rmd` | Exploratory: the Bayesian attribute-level model — sanity checks, parameter recovery, variant fits and cross-validation, out-of-sample attention weights. Written so its sections can be cut and pasted into `Main_Analysis.Rmd` if they hold up |
| `attention_identifies_reference.Rmd` | Exploratory: is the reference point identified from choices alone, without the dwell record? The algebra says no without loss aversion; this evaluates the project's own likelihood over a grid of reference points under four versions of the dwell data |
| `PAS3_Analysis.Rmd`, `PAS4_Analysis.Rmd` | The original per-study notebooks, superseded by `Main_Analysis.Rmd`. Kept for reference |

```r
rmarkdown::render("notebooks/Main_Analysis.Rmd")
```

`Main_Analysis.Rmd` has a `RUN_REGRESSIONS` switch. `FALSE` is a fast, figures-only pass that draws the attribute-level panels from the cached `brms` fits in `notebooks/brm_fits/` and the cached Bayesian-model predictions. `TRUE` fits every regression, writes the LaTeX tables and in-text statistics, and samples any `brms` model whose cache is missing — the two attention-weight fits take many hours each the first time, so run them from a terminal instead (below). `bayesian_model.Rmd` has the analogous `RUN_FITS`. `brm()` caches through its `file` argument; **delete the `.rds` to force a refit.**

Sections of `Main_Analysis.Rmd`, in the order of the Results: Basic Psychometrics (including the range-sensitivity simulations), Dwell Properties, Option-Level Attention (net-dwell and last-seen bias), Attribute-Level Attention (design matrices, fits, figures, tables), Range Sensitivity, and an *Extra Analyses* section that is switched off by default (the pre-registered linear-value-channel specification, the by-study and range-pooled attention weights, the brms range-sensitivity fit, participant-level out-of-sample log-likelihoods, and convergence checks).

Knitted `.html` is tracked, so the rendered results are readable from GitHub without running anything.

**Requirements:** R with `tidyverse`, `brms` (cmdstanr backend and a working CmdStan), `Rcpp`, `glmmTMB`, `lme4`, `bayestestR`, `posterior`, `modelsummary`, `patchwork`, `yaml`, `here`. `ggplot2` 4.x with `patchwork` 1.3.0 is assumed in the figure code (no `&` operator on patchworks).

---

## The Bayesian attribute-level model

An extension of Frömer, Callaway, Griffiths & Shenhav (2025, *Open Mind*) to attributes: every box is an item whose value is estimated by sampling while it is hovered, starting from a prior mean that acts as the reference point, with loss aversion applied to the posterior relative to that prior.

| File | Role |
|------|------|
| `code/simulation/simulate_bayesian_attributes.R` | The model and its choice probabilities — the probit where it is exact (λ = 1), otherwise the exact grid calculation |
| `code/simulation/exact_choice_prob.cpp` | That grid calculation (convolution on a utility grid), C++ via Rcpp |
| `code/fitting/fit_bayesian_attributes.R` | Maximum-likelihood fitting, the list of variants, and the trial-level cross-validation. Per-subject noise is a log-normal random effect whose mean and SD are fitted per range condition; the prior position `q` and loss aversion `λ` are pooled across subjects |
| `code/fitting/run_bayesian_fits.R` | Fits every variant and runs the cross-validation from a terminal, caching to `results/cache/bayesian_model/`. Takes variant names and `fit` / `loo` / `folds=a-b`, so the 12 folds can be split over processes |
| `code/parameter_recovery/recover_bayesian_attributes.R` | Recovery for the main variant (`cond_lam`) on 20 simulated datasets built from the real stimuli and dwell times; `datasets=a-b` splits them over processes |
| `code/fitting/run_all_bayesian_fits.sh` | Runs every fit, fold and recovery dataset as a queue of jobs, seven at a time (roughly 100 process-hours) |
| `code/fitting/watch_bayesian_fits.sh` | Live view of that queue — running jobs, folds and datasets done, latest log lines, memory |

```bash
caffeinate -i nohup zsh code/fitting/run_all_bayesian_fits.sh > results/logs/bayesian_driver.log 2>&1 &
zsh code/fitting/watch_bayesian_fits.sh
```

Every job caches its own output, so the queue can be re-run after an interruption and only the missing pieces are computed.

### Range-sensitivity simulations

Four models of value normalization, sharing the setup in `code/simulation/simulate_normalized_choices.R` and feeding the appendix figure: `simulate_range_normalization.R` (Padoa-Schioppa; Soltani et al.; Kőszegi & Szeidl), `simulate_rank_normalization.R` (decision by sampling), `simulate_context_divisive_normalization.R` and `simulate_trial_divisive_normalization.R` (divisive normalization by the condition's context and within the trial). `simulate_ddm.R` is a bare drift-diffusion simulator, groundwork for an aDDM section.

---

## Standalone figure scripts

`code/visualization/` holds figure scripts that stand outside the notebook pipeline. **They do not read the data** — the numbers are pasted in from the paper, so they have to be re-pasted whenever the underlying fit changes:

| Script | Output | Numbers come from |
|--------|--------|-------------------|
| `bayesian-model-plots.R` | `results/figures/bayesian-model-{narrow,wide}-r[-with-predictions].png` — attention weight against attribute value, with HDIs and the model overlay | Tables S6/S7 and Figure 5 of the manuscript, hardcoded (the model curve was recovered from the PDF's vector paths) |
| `task_attribute_range_subfigure.R` | `attribute_range_subfigure.{png,pdf}` — the Narrow/Wide value sets by attribute type | The task design, hardcoded |

```bash
Rscript code/visualization/bayesian-model-plots.R [output_dir]
Rscript code/visualization/task_attribute_range_subfigure.R
```

Neither resolves paths through `config/paths.yaml`: `bayesian-model-plots.R` defaults to an absolute `results/figures` path on Brenden's machine (pass an output directory to override), and `task_attribute_range_subfigure.R` writes into the working directory under a different name than the copy kept at `results/figures/task_attribute_range_subfigure.png`. Base graphics and `ggplot2` only — no `getAllUtilities.R`, no `brms`.

`results/figures/presentation/` holds talk figures (`attention_anchor_schematic`, `attention_theta_vs_bayesian_schematic`, `sampling_noise_distributions`) that have no generating script in the repository.

---

## Long sampling runs and the cluster

The hierarchical `brms` attention-weight fits are the expensive part of the paper. Three scripts run them outside RStudio, each by evaluating `Main_Analysis.Rmd`'s own chunks so the caches they leave in `notebooks/brm_fits/` are exactly what the notebook loads afterwards:

| Script | Fits |
|--------|------|
| `code/fitting/sample_pooled_attention_weights.R` | The two main attention-weight models (Narrow, Wide), pooled over studies. Optional range argument runs one at a time |
| `code/fitting/sample_bystudy_attention_weights.R` | The by-study robustness version (every coefficient crossed with study), an Extra Analyses check |
| `code/fitting/compute_canada_attention_weights.R` + `.sh` | The 100 fits of the attention-weight regression to the *model-simulated* choices — 5 variants × 2 range conditions × 10 replicates — as a SLURM array on the Alliance clusters (Nibi / Fir / Rorqual), 8 cores per task. These give the model side of the attention-weight overlays and the OOS-MAE column of `tables/bayesian_model_fits.tex` on the same estimator as the data side |

```bash
caffeinate -i nohup Rscript code/fitting/sample_pooled_attention_weights.R > results/logs/sampling.log 2>&1 &
```

The cluster workflow (draw the replicates locally → `rsync` into `~/scratch/attention_fits` → submit the array → `rsync` the draws back to `results/compute_canada/outputs/`) is documented step by step in the header of `compute_canada_attention_weights.sh`.

`code/power_analysis/power_analysis_for_attn_weights.R` is the original simulation-based power analysis for the attention-weight design (`lme4` + `foreach`).

---

## Data

`data/processed/` holds, per study, a **trial-level** table and a `*long*` **fixation-level** table (one row per hover event within a trial), plus the merged pair that everything current reads:

| Files | Read by |
|-------|---------|
| `main_data.csv`, `main_data_long.csv` | `Main_Analysis.Rmd`, every script under `code/` (via `build_boxlong.R`) |
| `PAS3.csv`, `PAS3long.csv` (+ `.RData`, `PAS3_train.csv`) | the merge script; `PAS3_Analysis.Rmd` |
| `PAS4.csv`, `PAS4long.csv` (+ `.RData`) | the merge script; `PAS4_Analysis.Rmd` |
| `PAS1{,long}.RData`, `PAS2{,long}.RData` | — (earlier pilots) |

The `.RData` files load objects named `padf` (trial level) and `padfLong` (fixation level). Key columns — `WN` (Narrow/Wide), `BPN` (trial type), `AttT`/`AttB` (attribute type per row), `VLT`/`VLB`/`VRT`/`VRB` (the four box values), `Dwell*` (per-box dwell in seconds), `FixROI`/`FixDur` (fixation events), `ChooseL`, `Correct`, `RT` — are documented in [experiments/experiment_and_analyses_overview.txt](experiments/experiment_and_analyses_overview.txt), which also records the task procedure, the Narrow/Wide value sets, the trial-construction constraints, the payment scheme and the exclusion criteria.

Study 2 additionally drops six participants whose P/N labels were inconsistent with the displayed value signs. `data/` is gitignored, so none of these files are distributed with the repository.

---

## Repository Structure

```
aversive-attributes/
├── README.md
├── LICENSE
├── Makefile                     # Task runner — all targets currently commented out
├── requirements.txt             # Empty (all code is R / C++)
├── .here                        # Project root marker for here::here()
│
├── config/paths.yaml            # Central path configuration (root, data, results, paper)
│
├── notebooks/
│   ├── Main_Analysis.Rmd                 # The paper's analysis (Studies 1 and 2)
│   ├── bayesian_model.Rmd                # Exploratory: Bayesian attribute-level model
│   ├── attention_identifies_reference.Rmd# Exploratory: reference point without attention
│   ├── PAS3_Analysis.Rmd, PAS4_Analysis.Rmd   # Superseded per-study notebooks
│   ├── *.html                            # Knitted output (tracked)
│   └── brm_fits/                         # Cached brms fits (.rds, gitignored)
│
├── code/
│   ├── preprocessing/merge_PAS3_PAS4_to_main_data.R
│   ├── utilities/               # Sourced helpers; getAllUtilities.R sources the directory
│   ├── simulation/              # Bayesian model + exact_choice_prob.cpp, normalization sims, DDM
│   ├── fitting/                 # ML fits, terminal samplers, Compute Canada array job
│   ├── parameter_recovery/recover_bayesian_attributes.R
│   ├── power_analysis/power_analysis_for_attn_weights.R
│   ├── visualization/           # Standalone figure scripts (hardcoded numbers, no data)
│   ├── behavior/  empirics/     # (empty — scaffold)
│
├── experiments/                 # Task materials, per-study details, trial generation
│   ├── generate_Trials.R        # Builds TrialsNBPN.csv / TrialsPBPN.csv (seed 1337)
│   ├── experiment_and_analyses_overview.txt
│   └── PAS2/  PAS3/  PAS4/      # Experiment details + Study 2 pre-registration
│
├── data/
│   ├── raw/                     # (empty — raw collection lives with the task platform)
│   └── processed/               # main_data{,_long}.csv and the per-study PAS1–PAS4 files
│
├── paper/                       # Manuscript — separate repo, not tracked here (see below)
├── docs/TODO.md                 # Worklist (stale: last revised 2026-08-10)
├── admin/                       # Conference submissions, collaborators, deadlines, funding
├── presentations/
├── results/
│   ├── cache/bayesian_model/    # Model fits, cross-validation folds, recovery datasets
│   ├── compute_canada/          # inputs/, outputs/, smoke/ for the cluster array
│   ├── figures/                 # paper/, exploratory/, presentation/, plus standalone script output
│   └── tables/  models/  logs/  temp/
├── cluster/
└── archive/                     # Deprecated code, old models, previous versions
```

> **Not tracked in this repo.** `paper/` is a **separate git repository** ([aversive-attributes-paper](https://github.com/BrendenEum/aversive-attributes-paper)) nested inside this one and listed in `.gitignore`, along with `archive/`, `notebooks/archive/`, `notebooks/notebook_figures/`, `results/cache|logs|temp/`, and the five template folders (`admin/`, `data/`, `docs/`, `experiments/`, `presentations/`). Model fits (`*.rds`, `notebooks/brm_fits/`), the knitr caches, and large data formats (`*.csv`, `*.zip`, `*.parquet`, …) are additionally excluded by pattern.

---

## Compiling the Paper

`paper/` is its own repository; clone it into place if it is not already present. `Main_Analysis.Rmd` and `bayesian_model.Rmd` write its figures, tables and in-text statistics, so render them first if anything upstream changed.

```bash
latexmk -pdf paper/main.tex
```

**Reference-dependence and range sensitivity determine whether attention to an aversive attribute attracts or repels choice.**

| | Section | Contains |
|---|---------|----------|
| | Statement of Relevance | |
| 1 | Introduction | |
| 2 | Methods | Task · Participants · Mouse Tracking · Mixed-Effects and Hierarchical Regressions · Computational Model · Model Fitting · Out-of-Sample Simulations |
| 3 | Results | Basic Psychometrics · Dwell Properties · Option-Level Attention · Attribute-Level Attention · Bayesian Sequential Sampling |
| 4 | Discussion | |
| | Transparency | Open practices (both studies preregistered on AsPredicted) · data and code availability |
| A | Appendices | Computing Choice Probability · Parameter Recovery · Model Fitting Details · Regressions · Range Sensitivity Simulations · Range Sensitivity by Attribute Type · Bayesian Model Variants · Attention Weights by Study · Attention Weights from Pre-Registered Model |

The body follows an explore/confirm split after Eum et al. (2023): Study 1 was used to settle on the analyses, and Study 2 — preregistered, 60 trials in three blocks of 20 rather than 30 interleaved — replays them as confirmation.

---

## Notes

### Template setup

This project was created from a personal project template. The first time you copy the template as a project folder:

1. Go to `~/Desktop/Google Drive/brendeneum/.` and make a subfolder with the name of the project.
2. Inside that, click and drag the `data`, `admin`, `docs`, `experiments`, and `presentations` folders.
3. Make a symlink for each (right click → Services → Make Symbolic Link), and drag those symlinks back to the main project folder.
4. Go to Overleaf, make your new paper. Sync it to GitHub as its own repo, clone that repo into this project folder as `paper`, and copy the contents of the template's `paper` folder into that sub-repo.

> In this project those five folders are currently **real directories, not symlinks** — step 3 has not been applied here. Step 4 has been applied: `paper/` is the cloned Overleaf repo.
