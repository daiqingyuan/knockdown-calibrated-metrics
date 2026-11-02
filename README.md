# Weighted calibrated metrics for perturbation-response prediction reward predicting the knockdown

Analysis code and results for the manuscript of the same title (Qingyuan Dai, University of Western Australia). The repository also supports a shorter Matters Arising version of the same analysis.

The manuscript re-scores the two unseen-gene CRISPRi benchmarks of Miller et al. (*Nature Biotechnology*, 2026; CellSimBench) and measures how much of the reward given by their weighted calibrated metrics (WMSE and weighted delta R²) comes from the knocked-down target gene. It compares published and retrained models with a knockdown-aware baseline that copies the training mean and lowers only the target gene.

## Contents

| Folder | What it holds |
|---|---|
| `code/` | Analysis scripts (Python). |
| `colab/` | Notebooks and helper scripts used to retrain PRESAGE, scLambda and fMLP-GenePT on Google Colab with the authors' code. |
| `results/` | Per-perturbation outputs and summaries behind the numbers in the manuscript. |
| `plans/` | English translations of the locked analysis plans, with SHA-256 digests (see `plans/README.md`). |
| `predictions_index.example.json` | Layout of the model-prediction files that the scripts read. |

## Data

The processed datasets are those distributed with CellSimBench (github.com/shiftbioscience/Perturbation-Models-Outperform-Baselines, commit `15361aec8392fb45a75c8c9da88ca5e6510fca97`), downloadable from the authors' public bucket:

```
https://public-access-bucket-i9sdj34j.s3.amazonaws.com/cellsimbench-data/data/<dataset>/<dataset>_processed_complete.h5ad
```

| Dataset | SHA-256 |
|---|---|
| `replogle22k562` | `594ddf179a743549ebc01a40cc4b6d25d8bde9d2472e89d2d613e5261690d80e` |
| `nadig25hepg2` | `4faf489ce12a4b9b90fc7d1415be89e167c058c0812d260b17971caeaec68525` |

Place each file at `$WMSE_DATA_DIR/<dataset>/<dataset>_processed_complete.h5ad`.

Model predictions are not stored in this repository. List the files you have in a copy of `predictions_index.example.json` (paths relative to that file, optionally with a SHA-256 to verify) and point `WMSE_PRED_INDEX` at it.

## Setup

Python 3.10 or later with the packages in `requirements.txt`.

| Variable | Meaning |
|---|---|
| `WMSE_DATA_DIR` | Folder with the processed datasets (see above). |
| `WMSE_WORK_DIR` | Output root. Caches go to `cache/` (Replogle) or `cache_nadig25hepg2/`; results go to `results/`, `results/nadig25hepg2/` and `results/target_gene/`. |
| `WMSE_CACHE_DIR` | Optional: reuse an existing cache folder. |
| `WMSE_PRED_INDEX` | Prediction index JSON. |
| `WMSE_DS` | `replogle22k562` (default) or `nadig25hepg2`. |
| `WMSE_PRED_SET` | `scgpt_gears` (default) or `retrained` (PRESAGE, scLambda, fMLP-GenePT). |

## Running

Run the scripts from `code/`, once per dataset (`WMSE_DS`):

1. **Statistics (cache):** `01_stats.py` then `02_tstats.py`. These compute per-condition statistics and the authors' differential-expression t statistics, and check them against the values stored in the data file. `01` is resumable and takes a few hours on a laptop.
2. **Reproduction and weight–truth coupling:** Replogle `03 → 05 → 06 → 07 → 08`; Nadig `13 → 15 → 16 → 17`. `06` and `16` are resumable.
3. **Target-gene analysis:** `20_target_gene.py → 21_kd_baseline_table.py → 22_weighted_r2.py`.
4. **Retrained models:** assemble the Colab exports with `30_assemble_colab_predictions.py <dir> [<dir> ...]`, then set `WMSE_PRED_SET=retrained` and run `20 → 22 → 31_retrained_report.py`. Outputs carry the suffix `_retrained`.
5. **Exploratory follow-up:** `python code/23_followup_exploratory.py`, run from the repository root once the outputs above exist. It needs only the per-perturbation results, not the data or the predictions, and writes `results/target_gene/followup_exploratory.json`: positions between the negative and positive controls of Miller et al., models given the knockdown, win rates, influence counts, Holm-adjusted P values, and how the target's contribution to the dataset mean's WMSE is concentrated in highly expressed targets.

All bootstrap intervals use 10,000 resamples of perturbations with seed 20261102.

## Where the manuscript's numbers come from

Items are numbered as in the preprint, with the Matters Arising item in brackets. Files are in `results/target_gene/` (one per dataset) unless a path is given.

| Manuscript item | Source |
|---|---|
| Table 1 | `<dataset>_summary.json` (`desc_*`, `knockdown_factor_by_fold`), `<dataset>_per_condition.csv`, and `followup_exploratory.json` (`target_share_of_mean_baseline_WMSE`) |
| Table 2 (Extended Data Table 1) | `<dataset>_summary[_retrained].json` (`G_nat`, `E1_share_target`, `E2_G_X_target_excluded`, `E3_gap_closed_by_bKD`) and `<dataset>_weighted_r2_summary[_retrained].json` (`E2prime_*`) |
| Table 3 (Extended Data Table 2) | `<dataset>_kd_baseline_table.json`, `<dataset>_weighted_r2_summary.json` (`baselines`) and `followup_exploratory.json` (`bKD_by_set`) |
| Table 4 and Fig. 2 (Fig. 2) | `<dataset>_summary[_retrained].json` (`E3_adv_vs_bKD`), `<dataset>_weighted_r2_summary[_retrained].json` (`E3prime_model_minus_KD`) and `followup_exploratory.json` (`win_rate_vs_bKD`; the Holm-adjusted P values are under `holm_fig1`, a key that predates the renumbering) |
| Table 5 (Extended Data Table 3) | `<dataset>_retrained_report.json` (`fidelity_gate`) for PRESAGE, scLambda and fMLP-GenePT; `results/native_wmse_summary.json`, `results/nadig25hepg2/native_wmse_summary.json` and `<dataset>_weighted_r2_summary.json` for scGPT and GEARS |
| Fig. 1 (Fig. 1) | `followup_exploratory.json` (`position_vs_positive_control`, `models.*.position`) |
| Fig. 3 | `pred_ratio_t` and `true_ratio_t` in the per-condition files, and `E2prime_*` |
| Per-fold estimates, influence counts and knockdown-supplied models in the text | `followup_exploratory.json` |

Where a source is given as `[_retrained]`, rows for scGPT, GEARS and the linear baseline come from the files without `_retrained`, and rows for PRESAGE, scLambda and fMLP-GenePT from the `_retrained` files.

The checks of our t statistics against the values stored in the data files (maximum absolute difference 6.1 × 10⁻⁶ on Replogle and 8.8 × 10⁻⁶ on Nadig) are in `results/tstat_validation.json` and `results/nadig25hepg2/tstat_validation.json`, written by `02_tstats.py`.

The weight–truth coupling results (Methods) are in `results/samecross_native_summary.json` (Replogle) and `results/nadig25hepg2/samecross_native_summary.json` (Nadig).

## Colab retraining

The notebooks rebuild each model's environment from its Dockerfile in the authors' repository and train with the authors' wrappers and default hyperparameters, on one A100 GPU. They were run in this order:

1. `01_retrain_replogle.ipynb`: fMLP-GenePT, scLambda and PRESAGE on Replogle fold 0.
2. `02_fix_sclambda.py`, executed in the same Colab session: the live Ensembl BioMart was unavailable, so scLambda's gene descriptions are fetched from the Ensembl release 115 archive with the wrapper's own code, and scLambda is re-run.
3. `03_run_more_folds.py`, executed in the same session: PRESAGE on Replogle folds 1–4.
4. `04_retrain_nadig_presage.ipynb`: PRESAGE on all five Nadig folds.

Notebook 04 clears PRESAGE's cached knowledge tensor (`WeightedDeepset.*.embeddings.pkl`) before each fold; notebook 01 and `03_run_more_folds.py` reuse it across folds. The tensor is built from external knowledge sources over the dataset's gene list and does not depend on a fold's training data, so reusing it passes no test information to training.

scLambda needs an OpenAI API key, entered in the notebook's settings cell. Never commit a key.

## Licence

The code in this repository is released under the MIT licence (`LICENSE`). PRESAGE (Genentech Non-Commercial Software License) and scLambda (GPL-3.0) are downloaded from their own repositories at run time and are not redistributed here.

## Citation

Dai, Q. Weighted calibrated metrics for perturbation-response prediction reward predicting the knockdown. Preprint at Zenodo (2026). https://doi.org/10.5281/zenodo.23928733
