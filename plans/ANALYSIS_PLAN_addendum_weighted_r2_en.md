# Analysis Plan Addendum: Weighted R²Δ

The original addendum was finalized on 2026-11-02, before any weighted R² calculations, and its SHA-256 digest was subsequently recorded.

The main plan is `ANALYSIS_PLAN_target_gene.md`. The cited original SHA-256 digest, 5ab9ccac…, refers to the Chinese version. The definition of the trivial knockdown baseline b_KD is unchanged: the primary form is multiplicative and the secondary form is additive.

This document is an English translation of the original addendum.

## Known before the addendum was finalized

- WMSE results of the target-gene analysis: b_KD is significantly better than scGPT, and its WMSE is lower than that of every model reported in the paper.
- The weighted R² values in the paper's tables had been read from the supplementary tables.
- No weighted R² had been calculated.

## Metric definitions (following the authors' code)

- **WR²(Δ Ctrl):** computed as r2(y − ctrl, ŷ − ctrl), where ctrl is `ctrl_baseline`.
- **WR²(Δ Pert):** computed as r2(y − b, ŷ − b), where b is the fold's `split_fold_k_mean_baseline`.
- **r2 computation:** the sklearn r2_score formula with sample_weight set to the native DEG weights. The calculation is restricted to the model's scored gene axis, and the result is clipped to ≥ −1. y is the authors' GT.

## Acceptance gate

For Dataset Mean, Control Mean, Linear, Tech. Dup. and Inter. Dup., the reproduced means of both metrics each differ from the paper's values by ≤ 0.01. If the gate is not passed, stop and diagnose.

## Endpoints

Conditions weighted equally, five folds pooled, 10,000 bootstrap resamples of conditions, seed 20261102.

- **E3′:** for each of the two metrics, the paired difference between the model and b_KD, R²(f) − R²(b_KD). Models: scGPT, GEARS and Linear on Replogle, and Linear on Nadig.
- **Descriptive:** the means of b_KD (both forms) over all test conditions, set beside the model values reported in the paper (unpaired).
- **E2′ (secondary):** remove the target gene from the scored genes, recompute the weights from the min–max without the target gene, and then compute the paired difference between the model and Dataset Mean.

## Decision rules

- **D2′:** for each metric, if the 95% interval of E3′ contains 0 or is entirely < 0, conclude "on this metric, b_KD is not significantly worse than the model"; if the interval is entirely > 0, conclude "the model is better than b_KD on this metric".
- Definitions, endpoints and thresholds will not be changed after the results are seen; any deviation will be labelled separately.
