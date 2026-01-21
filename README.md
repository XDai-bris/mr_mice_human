
# mr_mice_human

MR workflow for **Diversity Outbred (DO) mice** using kinship-aware SNP effect estimation (LOCO–GLS LMM), Steiger filtering, and IVW Mendelian randomization, plus a **cross-species comparison** against **human MR** using **OpenGWAS / TwoSampleMR**.

This repo is primarily a *research script collection* (not an R package). The analysis follows the workflow:

1. **Find QTL top-hits** (per trait)  
2. **Union** top-hits across traits to define a shared instrument candidate set  
3. Estimate signed SNP–trait **effect sizes** with **LMM–GLS (LOCO)**  
4. **Steiger filtering (SF)** to retain reliable instruments / directionality  
5. **MR-IVW** to estimate causal effects across trait pairs  
6. Compare MR vs **observational regression** (pairwise + grouped summaries/figures)  
7. Run **human MR** from OpenGWAS (best-matched phenotypes) and generate mouse–human comparison figures

> **Note:** Input data (phenotypes, genotypes, scan outputs, kinship) are project-specific and are **not** shipped in this GitHub repo. The scripts assume you have local files prepared from upstream DO QTL scans (e.g., qtl2 outputs).

---

## Repository contents (main scripts)

- `FUN.R`  
  Core helper functions used across the pipeline (data wrangling, utilities).

- `RFUN_clean.R`  
  Additional/cleaned helper functions used by downstream scripts.

- **Top-hit / instrument selection**
  - `plot_tophits.R`  
    Visualisation and/or inspection of QTL top-hits.

- **SNP effect estimation (mouse)**
  - `lmm_effect.R`  
    Computes SNP–trait effect sizes and SEs with kinship-aware LMM, using a fast LOCO–GLS strategy (and/or comparisons vs OLS / exact LMM).
  - `plot_cmp_effect.R`  
    Plots comparisons of effect estimates across models.

- **MR (mouse)**
  - `mr_SF.R`  
    MR with **Steiger filtering** (directionality check), then **IVW** MR.
  - `mr.R`  
    MR utilities and/or MR runs without SF / baseline MR.
  - `demo_1pairMR.R`  
    Minimal worked example for **one exposure–outcome pair**.

- **MR vs observational regression**
  - `plot_mr_obs_cmp.R`  
    Pairwise MR vs observational comparison plots.
  - `plot_mr_obs_cmp_group.R`  
    Grouped summaries (e.g., by phenotype groups) and cleaner overview plots.

- **Human MR + mouse–human comparison**
  - `humanMR.R`  
    Human two-sample MR using **OpenGWAS** (via `TwoSampleMR` / `ieugwasr`), with harmonisation and (when feasible) SD-scaling of effects.
  - `plot_miceHuman.R`  
    Mouse–human effect comparison plots (e.g., scatter/forest style summaries).

---

## Outputs already included in this repo (examples)

The repo includes example figures/tables produced by the scripts, such as:

- `tophits.png`
- `cmp_lmm*.png`, `cmp_ols_lmm*.png`
- `cochran_Q_pre_vs_post.png`
- `MR_OBS.png`
- `Fig_MR_vs_OBS_concordance*.pdf`
- `human_vs_mouse_betas.pdf`
- `mr_forest_pairs_fig/` (folder of forest plots)
- `res_lmm.rda` (example saved results object from LMM effect estimation)

These are useful as references to confirm that your pipeline re-runs are producing similar outputs.

---

## End-to-end workflow (mouse)

### Step 0 — Prepare inputs (local, not in repo)
You need (at minimum):

- Phenotype matrix/table (rows = mice, columns = traits)
- Covariates (sex, age at sacrifice, generation, etc.)
- Genotype dosage matrix for marker panel (0/1/2 coding)
- Kinship matrices (preferably LOCO kinship per chromosome)
- QTL scan results / LOD profiles per trait (from upstream scan, e.g., qtl2)

Make sure all data align by a common mouse ID (e.g., `Mouse.ID`).

---

### Step 1 — Identify trait-specific QTL top-hits
**Goal:** from each trait’s LOD profile, pick approximately independent peaks above a threshold, then define support intervals.

Typical script entry point:
```r
source("FUN.R")
source("RFUN_clean.R")
source("plot_tophits.R")
```

Output: a per-trait list of top-hit loci + a global view (`tophits.png`).

---

### Step 2 — Union top-hits across traits (shared instrument candidates)
**Goal:** build one shared marker set so every trait has effect estimates at the same candidate instruments.

This is usually performed inside the main analysis scripts (MR and LMM effect estimation), once top-hits are computed.

---

### Step 3 — Estimate signed SNP–trait effects with LOCO–GLS LMM
**Goal:** obtain `beta` and `SE` for each (trait, marker) pair, correcting for kinship and covariates.

Entry point:
```r
source("lmm_effect.R")
```

Output: saved effect-size objects (example: `res_lmm.rda`) + diagnostic plots (`cmp_lmm*.png`, etc.).

---

### Step 4 — Steiger filtering (SF) for directionality & reliable instruments
**Goal:** keep instruments that explain more variance in the proposed exposure than outcome (and optionally pass a Steiger p-value threshold).

Entry point:
```r
source("mr_SF.R")
```

Output: MR results per ordered pair, instrument counts pre/post SF, and heterogeneity summaries (e.g., Cochran’s Q plot).

---

### Step 5 — MR-IVW (pairwise across traits)
**Goal:** run systematic ordered-pair MR using IVW, using the harmonised signed effects from Step 3 and the filtered instruments from Step 4.

Entry point:
```r
source("mr.R")      # baseline MR utilities / MR runs
source("mr_SF.R")   # primary MR with SF
```

Output: MR tables for downstream plotting and summaries.

---

### Step 6 — Compare against observational regression (pairs + grouped summaries)
**Goal:** compare MR-derived causal estimates to naive phenotype–phenotype regressions (possibly standardised).

Entry points:
```r
source("plot_mr_obs_cmp.R")
source("plot_mr_obs_cmp_group.R")
```

Outputs: figures such as `MR_OBS.png` and `Fig_MR_vs_OBS_concordance*.pdf`.

---

## Human MR (OpenGWAS) + cross-species comparison

### Step 7a — Human two-sample MR using OpenGWAS
**Goal:** for best-matched human traits (mapped from mouse phenotypes), run:
- instrument extraction (genome-wide significant SNPs)
- LD clumping
- outcome extraction
- harmonisation
- IVW (plus optional sensitivity methods)
- SD scaling when feasible (continuous traits)

Entry point:
```r
source("humanMR.R")
```

> You may need an OpenGWAS token depending on API requirements and the datasets you access.

---

### Step 7b — Mouse–human comparison plots
Entry point:
```r
source("plot_miceHuman.R")
```

Outputs: e.g. `human_vs_mouse_betas.pdf` and forest/scatter summaries for matched exposure–outcome pairs.

---

## Minimal “one-pair” demo

If you want to sanity-check the pipeline on a single pair first:
```r
source("demo_1pairMR.R")
```

This is the recommended starting point before launching high-throughput all-by-all MR.

---

## Suggested R dependencies

Your environment will likely need (depending on which scripts you run):

- Mixed-model utilities used by your LMM implementation
- `TwoSampleMR`, `ieugwasr` (human MR)
- Any project-specific packages used in upstream DO scans (e.g., `qtl2`)

---

## Reproducibility notes

- Keep raw data out of git.
- Record thresholds (LOD cutoff, clumping params, Steiger alpha, min instrument count).
- Save intermediate objects (`res_lmm.rda`, MR result tables) to avoid recomputing expensive steps.
- Consider adding a `config.R` or `config.yml` locally for paths and parameters (gitignored).

---

## Citation / acknowledgement

If you use this pipeline in a manuscript, please cite:
- `TwoSampleMR` / OpenGWAS for human MR

---

## Contact

Maintainer: **Xiaoyang Dai**  
GitHub: https://github.com/XDai-bris
