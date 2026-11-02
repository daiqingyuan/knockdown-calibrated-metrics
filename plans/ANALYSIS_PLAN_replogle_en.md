# Replogle22 K562 Same-Sample and Cross-Sample Scoring Analysis Plan

The original plan was finalized on 2026-11-02 (Asia/Shanghai), before any same-sample scores, cross-sample scores, or I were calculated. Its SHA-256 digest was subsequently recorded.

This document is an English translation of the original plan. The original digest applies to the Chinese version.

## Information seen before the plan was finalized

- First-half DEG t statistics: our implementation matches the authors' stored `uns/scores_df_dict_gt` (2044 conditions, maximum absolute difference 6.1e-6).
- Reproduction of the original-method WMSE (native): GEARS 0.0335, scGPT 0.0268, Dataset Mean 0.0337, Control Mean 0.0349, Linear 0.0344, Tech. Dup. 0.0108; the paper's values are, in the same order, 0.033, 0.027, 0.034, 0.035, 0.034 and 0.011.
- Paired differences from Dataset Mean: GEARS 0.00027 ± 0.00017, scGPT 0.00872 ± 0.00084 (condition-level SEM).
- No same-sample or cross-sample scores, scores under cross-half weights, or I had been calculated.

## Data and predictions

- Data: `replogle22k562_processed_complete.h5ad` with the authors' five folds; the test conditions of each fold are the non-control conditions whose role in that fold is test for all cells.
- Models: previously frozen GEARS fold 0–4 and scGPT fold 0–4 predictions (paths in `code/common.py`), with no retraining and no selection.
- Mean baseline b: the authors' `uns/split_fold_k_mean_baseline` (the paper's Dataset Mean).
- Two half-pools: A = `first_half` (the authors' GT), B = `second_half` (the authors' Tech. Dup.). The native budget is all cells in each half (6–76, median 61).

## Weights and scoring

- qA: computed within partition A by the authors' rule (scanpy `t-test_overestim_var`, reference=rest, non-control conditions only), |t| → per-condition min–max → square. qB: the same rule applied within partition B.
- Normalized on each model's own scored gene axis (GEARS 8370 genes, scGPT 8369 genes). Baselines are compared on the same axis.
- YA, YB: log-normalized means of the two halves. δ = f − b.
- G(q, Y) = Σ q̃ [(b − Y)² − (f − Y)²]; a positive value means the model is better than the mean baseline.
- Same-sample S = [G(qA, YA) + G(qB, YB)] / 2; cross C = [G(qA, YB) + G(qB, YA)] / 2; I = S − C = Σ δ (q̃A − q̃B)(YA − YB).
- Native gap G_nat = G(qA, YA) (that is, WMSE_b − WMSE_f under the paper's protocol).

## Primary endpoints (reported for each model separately; five folds pooled; conditions weighted equally)

1. Mean of I with its 95% interval: 10,000 bootstrap resamples of conditions, seed 20261102.
2. Mean of C with its 95% interval (the advantage once the weights are decoupled from the target).
3. Means of S and G_nat, for comparison.

## Prespecified interpretation rules

- R1: if the 95% interval of C is entirely > 0, conclude "the advantage holds after decoupling"; if the interval contains 0 or is entirely < 0, conclude "the advantage is no longer supported after decoupling".
- R2: if the 95% interval of I excludes 0, conclude "the coupling bias is detectable".
- R3: if the mean of I is ≥ 0.0005 (half of 0.001, the unit in which the paper's tables report, and enough to change a WMSE given to three decimal places), conclude "the bias is material at the reported precision"; otherwise conclude "negligible at the reported precision".
- Per-fold means and signs are also reported, as a description of consistency, without a separate decision.

## Secondary analyses

- Cell-number gradient: n = 8, 16, 32, 64 cells per half, restricted to conditions with at least n cells in both halves. For each condition, 3 nested random subsamples are drawn (seed 20261102); the target cells' weights and means are recomputed from the subsampled cells, while the reference remains all cells of the other conditions in the same partition. Report the means and intervals of I and C at each n.
- Sensitivity analyses: the conditions shared by the two models; I of Linear relative to Dataset Mean (descriptive).

## Exclusions and limitations

- Conditions without a prediction from a model are excluded for that model, without imputation. Conditions with an invalid weight panel (max = min) are excluded and counted.
- Intervals use conditions as the only resampling unit and are conditional on the fixed reference and the fixed training results; they do not cover training randomness (one seed per fold) or the dependence created by the shared reference.
- This analysis is the discovery stage, and its results are not confirmatory evidence. Confirmation on Nadig requires a separately locked plan.
- Endpoints, thresholds, the cell-number grid and inclusion rules will not be changed after the results are seen; any deviation will be labelled separately.
