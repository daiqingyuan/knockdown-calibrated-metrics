# Nadig25 HepG2 Confirmatory Analysis Plan

The original plan was finalized on 2026-11-02 (Asia/Shanghai), before any Nadig statistics, WMSE, or same-sample and cross-sample scores were calculated. Its SHA-256 digest was subsequently recorded.

This document is an English translation of the original plan. The original digest applies to the Chinese version.

## Basis

Results of the Replogle discovery stage (original plan SHA-256 eb14e953…51ae2):
- At the native budget, the coupling bias I as a share of the same-sample advantage S: scGPT 4.9%, Linear 12.4%.
- The fewer the cells, the higher this share: at 8 cells per half, scGPT 29%, Linear 43%.

Nadig has 6–31 cells per half (median 23), so its share is expected to be higher than the native Replogle values.

## Known before the plan was finalized

Only the metadata had been read:
- 109,183 cells, 8,746 genes, 2,322 non-control conditions.
- 6–31 cells per half, with the two halves differing by at most 1.
- 2,263 conditions have ≥ 8 cells in both halves, and 1,747 have ≥ 16.
- The uns keys are the same as in Replogle.

No expression statistics, WMSE or I had been calculated.

## Acceptance gates (if a gate is not passed, stop and diagnose; the hypotheses are not evaluated)

1. Maximum absolute difference between the first-half DEG t statistics and the authors' `uns/scores_df_dict_gt` < 1e-4.
2. The original-method WMSE of Dataset Mean, Control Mean, Linear, Tech. Dup. and Inter. Dup. each differs from the paper's value by ≤ 0.001. The paper's values are, in the same order, 0.048, 0.051, 0.048, 0.022 and 0.029.

## Primary hypotheses

The comparison is Linear versus Dataset Mean, five folds pooled, native budget, conditions weighted equally. Intervals come from 10,000 bootstrap resamples of conditions, seed 20261102.
- **H1 (primary):** the lower bound of the 95% interval of I > 0.
- **H2 (prespecified secondary):** the lower bound of the bootstrap 95% interval of the share mean(I)/mean(S) > 0.124 (the native Linear share on Replogle).
- **H3 (prespecified secondary):** the mean of I ≥ 0.0005, i.e., material at the paper's reported precision.
- **R1:** the lower bound of the 95% interval of C > 0, i.e., the advantage holds after decoupling.

H2 and H3 are not corrected for multiple comparisons; all results are reported in full.

## Secondary (descriptive)

- **scGPT fold 0:** performed only if the previously trained scGPT prediction file is available and its SHA-256 equals `383e0d1307c888807537d86f88297efea949d6e7f390e007c6723a71140b1603`. To be reported:
  - The original-method WMSE, compared with the paper's 0.040. Note that this is a single fold, and that the scGPT gene annotation was reconstructed and has not been confirmed to match the one used for the paper.
  - G_nat, C, I, I/S and R1–R3.

  Because there is only one fold, the results are descriptive only, with no confirmatory decision. If the file cannot be obtained, this item is marked as not completed, and no other model is substituted.
- **Cell-number gradient:** n = 8, 16 cells per half, restricted to conditions with ≥ n cells in both halves; 3 nested random subsamples per condition, seed 20261102. The model is Linear (scGPT fold 0 is computed as well when available).

## Other rules

Carried over from the Replogle plan:
- The authors' weight rule: `t-test_overestim_var`, reference=rest, non-control conditions only, |t| → min–max → square. Each half is computed within its own partition.
- Normalization on the model's scored gene axis.
- Missing predictions are not imputed; invalid weight panels are excluded and counted.
- Intervals are conditional on the fixed reference and the training results.

Hypotheses, thresholds, the cell-number grid and inclusion rules will not be changed after the results are seen; any deviation will be labelled separately.
