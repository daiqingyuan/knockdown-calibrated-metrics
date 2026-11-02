# Contribution of the Target Gene to the WMSE Advantage

The original plan was finalized on 2026-11-02 (Asia/Shanghai), before any target-gene-related weights, expression values, predictions, or WMSE decompositions were calculated. Its SHA-256 digest was subsequently recorded.

This document is an English translation of the original plan. The original digest applies to the Chinese version.

## Known before the plan was finalized

- **Target-gene coverage:** only gene names were compared. Conditions whose target gene (the gene symbol corresponding to the condition name) is on the scored gene axis: Replogle 1855/2044, Nadig 2100/2322.
- **Existing results:** native WMSE gaps relative to Dataset Mean were scGPT 0.00872, GEARS 0.00027 and Linear 0.00133 on Replogle, and Linear 0.00203 on Nadig.
- **Not yet seen:** the target gene's weights, expression and predicted values, and its contribution to WMSE.

## Definitions

- **Target gene t(c):** the gene symbol of condition c, provided it is on the model's scored gene axis; conditions whose symbol is not on the axis are treated as having no target gene.
- **Native score:** the native weights q are computed on the GT half-pool by the authors' rule and normalized on the model axis to q̃; the target y is the authors' GT.
- **Per-gene advantage:** g_g = (b_g − y_g)² − (f_g − y_g)². G_nat = Σ q̃_g g_g = WMSE_b − WMSE_f.
- **Target-gene contribution share:** share = mean_c(q̃_t g_t) / mean_c(G_nat), i.e., a ratio of means.
- **Target-excluded metric WMSE_X:** the min–max of |t| is recomputed with the target gene removed, then squared, and normalized over the remaining genes only. G_X = WMSE_X(b) − WMSE_X(f).
- **Knockdown-aware trivial baseline b_KD:**
  1. For each fold, use only the non-control conditions c' with role = train in that fold, and only those whose target gene is on the axis.
  2. For each c', compute r_c' = ȳ_c'[t(c')] / b[t(c')] from the mean over all cells; their median is f_k.
  3. For a test condition c, b_KD equals b, except that b_KD[t(c)] is replaced by f_k · b[t(c)].

  The secondary variant uses an additive form: d_k = median(ȳ_c'[t'] − b[t']).

## Endpoints

Reported for each model separately. Models: scGPT, GEARS and Linear on Replogle; Linear on Nadig, and scGPT fold 0 if its file is obtained. Conditions weighted equally, five folds pooled, intervals from 10,000 bootstrap resamples of conditions, seed 20261102.

- **E1:** the target-gene share and its 95% interval.
- **E2:** the mean and interval of G_X, i.e., whether the model still beats the mean baseline once the target gene is excluded.
- **E3:** the mean and interval of the model's advantage over b_KD, WMSE(b_KD) − WMSE(f); and the fraction of the gap closed by b_KD, closed = mean(WMSE_b − WMSE_bKD) / mean(WMSE_b − WMSE_f).
- **Descriptive measures:**
  - the proportion of conditions in which the target gene has the largest native weight;
  - the distribution of q̃_t;
  - the target gene's share of the mean baseline's WMSE;
  - the medians of the model's target-gene prediction f_t/b_t and of the truth y_t/b_t.

## Decision rules

The primary model is scGPT on Replogle; Nadig and the other models are reported with the same method. The native advantage of GEARS is not itself significant, so the decisions for GEARS are descriptive only.

- **D1:** if share ≥ 50%, conclude "the target gene is the main source of the WMSE advantage".
- **D2:**
  - if the 95% interval of E3 contains 0 or is entirely < 0, conclude "the trivial knockdown-aware baseline is not significantly worse than the model, which directly challenges the main conclusion of the original paper";
  - if the interval is entirely > 0, conclude "the model's advantage cannot be explained by a trivial knockdown of the target gene".
- **D3:** if the interval of E2 is entirely > 0, conclude "the advantage holds after the target gene is excluded".

Definitions, endpoints and thresholds will not be changed after the results are seen; any deviation will be labelled separately.
