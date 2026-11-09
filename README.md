# Aversive Attributes — Does attention to something bad attract or repel choice?

Brenden Eum (Rotman, University of Toronto) · Stephanie M. Smith (Booth, University of Chicago)

Attention is known to bias choice toward whatever is currently being looked at, but that regularity was established almost entirely with stimuli whose values are positive. This project asks what attention does when the attended thing is *bad*.

Two mouse-tracking experiments have participants choose between two options described by two attributes, where each attribute is labelled in advance as always appetitive (**P**), always aversive (**N**), or possibly either (**B**). The headline results:

- The pull of attention increases with the value of the attended attribute *within* every attribute type. Attending to a mildly aversive attribute still moves choice toward the option carrying it; attending to a severely aversive attribute moves choice away.
- Where that reversal occurs is not fixed — it shifts with the range of values the participant was trained on.
- Choice sensitivity to value is itself range-normalized within attribute type: the same objective value carries less weight when drawn from a wider range.
- Value-related attention is expressed almost entirely through **refixation count**, not through how long each fixation lasts.

Between-subjects the value range is manipulated (**Narrow** vs. **Wide**); within-subjects the trial type varies (**PN**, **PB**, **NB**, refined by the realized sign of the B attribute into PB+/PB−/PBo and NB+/NB−/NBo).

---

## Status

Both experiments are collected and analyzed; the manuscript is drafted end to end and currently carries `\textcolor{red}{...}` placeholders where numbers and sample details are still outstanding. [docs/TODO.md](docs/TODO.md) is the live worklist — it tracks the model-fit problems that must be resolved before the paper is defensible (MAP vs. posterior medians, Exp 1 value-sensitivity convergence, the preregistered range-sensitivity restriction, the refixation-count model, and the trials with no dwell record), the analyses still to be added, and the items owed by collaborators.

`paper/figures/` is empty: the manuscript's seven figures are referenced by `main.tex` but are not yet produced by any tracked script. The notebooks render their plots inline instead. `requirements.txt` is empty (everything here is R, not Python) and every target in the `Makefile` is still commented out.

Naming note: the two experiments in the paper are **Experiment 1 = PAS3** and **Experiment 2 = PAS4**. PAS1 and PAS2 were earlier pilots; their processed data is still in `data/processed/`, but no current notebook reads it.

---

## Repository Structure

```
aversive-attributes/
├── README.md
├── LICENSE
├── Makefile                     # Task runner — all targets currently commented out
├── requirements.txt             # Empty (all code is R)
├── .here                        # Project root marker for here::here()
│
├── config/
│   └── paths.yaml               # Central path configuration (root, data, results, paper)
│
├── notebooks/
│   ├── PAS3_Analysis.Rmd        # Experiment 1 — exploratory
│   ├── PAS4_Analysis.Rmd        # Experiment 2 — preregistered / confirmatory
│   ├── PAS{3,4}_Analysis.html   # Knitted output (tracked)
│   ├── brm_fits/                # Cached brms fits (.rds, gitignored)
│   └── PAS{3,4}_Analysis_cache/ # knitr chunk cache (gitignored)
│
├── code/
│   ├── power_analysis/power_analysis_for_attn_weights.R
│   ├── simulation/simulate_ddm.R
│   ├── utilities/               # Sourced helpers (see below)
│   ├── behavior/  empirics/  fitting/  parameter_recovery/
│   ├── preprocessing/  visualization/       # (all empty — scaffold)
│
├── experiments/                 # Task materials, per-study details, trial generation
│   ├── generate_Trials.R        # Builds TrialsNBPN.csv / TrialsPBPN.csv
│   ├── experiment_and_analyses_overview.txt
│   ├── PAS2/  PAS3/  PAS4/      # Experiment details docs + Study 4 pre-registration
│   └── AA Next Steps.pdf, experiment_details.pptx
│
├── data/
│   ├── raw/                     # (empty — raw collection lives with the task platform)
│   └── processed/               # PAS1–PAS4 trial-level and *long* fixation-level files
│
├── paper/                       # Manuscript — separate repo, not tracked here (see below)
├── docs/TODO.md                 # Live worklist
├── admin/                       # Conference submissions, collaborators, deadlines, funding
├── presentations/
├── results/                     # figures/, tables/, models/, cache/, logs/, temp/
├── cluster/
└── archive/                     # Deprecated code, old models, previous versions
```

> **Not tracked in this repo.** `paper/` is a **separate git repository** ([aversive-attributes-paper](https://github.com/BrendenEum/aversive-attributes-paper)) nested inside this one and listed in `.gitignore`, along with `archive/`, `notebooks/archive/`, `notebooks/notebook_figures/`, `results/cache|logs|temp/`, and the five template folders (`admin/`, `data/`, `docs/`, `experiments/`, `presentations/`). Model fits (`*.rds`, `notebooks/brm_fits/`), the knitr caches, and large data formats (`*.csv`, `*.zip`, `*.parquet`, …) are additionally excluded by pattern.

---

## Data

`data/processed/` holds two files per study — a **trial-level** table and a `*long*` **fixation-level** table (one row per fixation event within a trial):

| Study | Files | Read by |
|-------|-------|---------|
| PAS1, PAS2 | `PAS1.RData`, `PAS1long.RData`, `PAS2.RData`, `PAS2long.RData` | — (earlier pilots) |
| PAS3 (Exp 1) | `PAS3.csv`, `PAS3long.csv` (+ `.RData`, `PAS3_train.csv`) | `notebooks/PAS3_Analysis.Rmd` |
| PAS4 (Exp 2) | `PAS4.RData`, `PAS4long.RData` (+ `.csv`) | `notebooks/PAS4_Analysis.Rmd` |

The `.RData` files load objects named `padf` (trial level) and `padfLong` (fixation level). Key columns — `WN` (Narrow/Wide), `BPN` (trial type), `AttT`/`AttB` (attribute type per row), `VLT`/`VLB`/`VRT`/`VRB` (the four box values), `Dwell*` (per-box dwell in seconds), `FixROI`/`FixDur` (fixation events), `ChooseL`, `Correct`, `RT` — are documented in full in [experiments/experiment_and_analyses_overview.txt](experiments/experiment_and_analyses_overview.txt), which also records the task procedure, the Narrow/Wide value sets, the trial-construction constraints, the payment scheme, and the exclusion criteria.

After exclusions the analyzed samples are 981 participants / 28,646 trials in Experiment 1 and 525 participants / 30,523 trials in Experiment 2. Experiment 2 additionally drops six participants whose P/N labels were inconsistent with the displayed value signs; they are hard-coded in the `exclude-subjects` chunk of `PAS4_Analysis.Rmd`.

`data/` is gitignored, so these files are not distributed with the repository.

---

## Running the Analysis

The two notebooks are the analysis. They are self-contained R Markdown documents and are the only entry points that produce results.

**Requirements:** R with `tidyverse`, `brms` (with the **cmdstanr** backend and a working CmdStan install), `glmmTMB`, `lme4`, `bayestestR`, `patchwork`, `yaml`, and `here`.

```r
rmarkdown::render("notebooks/PAS3_Analysis.Rmd")   # Experiment 1
rmarkdown::render("notebooks/PAS4_Analysis.Rmd")   # Experiment 2
```

Paths are resolved through [config/paths.yaml](config/paths.yaml) via `here::here()`, anchored by the `.here` file at the repository root, so the notebooks knit from either the root or `notebooks/`. Note that `paths.yaml` hardcodes an absolute `root` (`/Users/brenden/Desktop/aversive-attributes`) — update it if you clone the project elsewhere.

Fitted `brms` models are cached as `.rds` under `notebooks/brm_fits/` and reloaded on re-render; **delete those files to force a refit**. knitr chunk caching is on (`cache = TRUE, autodep = TRUE`), with `PAS{3,4}_Analysis_cache/` holding the chunk cache. Both caches are gitignored; the knitted `.html` and its `_files/figure-html/` PNGs are tracked, so the rendered results are readable from GitHub without running anything.

### What the notebooks do

Both follow the same sequence; sections marked *(Exp 2 only)* exist only in `PAS4_Analysis.Rmd`.

| Section | Content |
|---------|---------|
| Preamble | Libraries, the shared `plt` plot-style list (Rotman palette, named scales for type / range / trial type / box), brms sampler settings (4 chains, 2000 iters, seed 1234) |
| Data Cleaning | Factor coding with Narrow as reference, refined trial type from the realized B valence, normalized net fixation `DwellDiff / (DwellL + DwellR)`, a per-box long representation, and fixation-level cleaning with per-box refixation counts and a total-fixations exposure offset |
| Basic Psychometrics | Choice as a function of attribute value and of value difference, and RT against difficulty, split by trial type × range |
| Fixation Properties | Dwell by box and condition, fixation duration against value, and refixation counts by type and value — the split that isolates revisiting rather than dwell length |
| Option-Level Attention | `glmmTMB` logistic fit of `ChooseL ~ Diff_z * DwellDiff_norm + (1 + Diff_z \| Subject)` — the classic net-fixation-bias test |
| Attribute-Level Attention | The flagship analysis. A signed attention weight per (attribute type × value) cell: the effect of one second of dwell on a box on the log-odds of choosing that box's option, holding value constant. Fit in `brms` with partial pooling — every subject gets a weight per cell, drawn uncorrelated (`\|\|`) for identifiability, while value sensitivity (`SV_P`/`SV_N`/`SV_B`) keeps correlated random slopes. Partial pooling matters here: omitting subject slopes attenuates the population weight toward zero |
| Attention Weights by Range Condition *(Exp 2 only)* | The same estimation refit separately on Narrow and Wide trials, so the value→weight mapping can be compared across ranges. Experiment 1 has too few trials per participant to support this |
| Range Sensitivity | A value-representation test rather than an attention test: choice on the per-type signed value channels interacted with range. Range normalization predicts each `SV_type:WNWide` interaction is negative — value drives choice less when the range is wide |

### Supporting scripts

Nothing in `code/` is called by the notebooks; these are standalone.

| Script | Purpose |
|--------|---------|
| `code/power_analysis/power_analysis_for_attn_weights.R` | Parallel simulation-based power analysis for the attention-weight design (`lme4` + `foreach`); reads `df_t.RData` for effect-size priors |
| `code/simulation/simulate_ddm.R` | A bare drift-diffusion simulator (`simulate_ddm()`), groundwork for the aDDM section of the paper |
| `code/utilities/*.R` | Helpers — `SE()`, `getStatsFromPosterior()`, `write.text()`, plot/brms option lists; `getAllUtilities.R` sources the whole directory |
| `experiments/generate_Trials.R` | Generates the trial sets (`TrialsNBPN.csv`, `TrialsPBPN.csv`) under the dominance and tie constraints, seeded at 1337 |

---

## Compiling the Paper

`paper/` is its own repository; clone it into place if it is not already present.

```bash
latexmk -pdf paper/main.tex
```

**Reference dependence and range sensitivity determine whether attention to an aversive attribute attracts or repels choice.**

| | Section | Contains |
|---|---------|----------|
| — | Statement of Relevance | |
| 1 | Introduction | |
| 2 | Methods | Task · Participants · Mouse Tracking · Inference Strategy · Hierarchical Regressions · Computational Model · Model Fitting · Out-of-Sample Simulations |
| 3 | Results | Basic Psychometrics · Fixation Process · Option-Level Attention · Attribute-Level Attention · Range Sensitivity of Attribute-Level Attention · aDDM and Reference Dependence |
| 4 | Discussion | |
| — | Transparency | |
| S | Supplementary Material | S1 Sample and Exclusions · S2 Basic Psychometrics · S3 Fixation Process · S4 Option-Level Attention · S5 Attribute-Level Attention · S6 Attention Weights by Range Condition · S7 Range Sensitivity of Value · S8 Computational Model |

The body follows an explore/confirm split after Eum et al. (2023): Experiment 1 was used to settle on the analyses, and Experiment 2 — preregistered, 60 trials in three blocks of 20 rather than 30 interleaved — replays them as confirmation. `main.tex` expects seven figures in `paper/figures/` (`fig1_task` … `fig7_value_sensitivity`); that directory is currently empty.

---

## Notes

### Template setup

This project was created from a personal project template. The first time you copy the template as a project folder:

1. Go to `~/Desktop/Google Drive/brendeneum/.` and make a subfolder with the name of the project.
2. Inside that, click and drag the `data`, `admin`, `docs`, `experiments`, and `presentations` folders.
3. Make a symlink for each (right click → Services → Make Symbolic Link), and drag those symlinks back to the main project folder.
4. Go to Overleaf, make your new paper. Sync it to GitHub as its own repo, clone that repo into this project folder as `paper`, and copy the contents of the template's `paper` folder into that sub-repo.

> In this project those five folders are currently **real directories, not symlinks** — step 3 has not been applied here. Step 4 has been applied: `paper/` is the cloned Overleaf repo.
