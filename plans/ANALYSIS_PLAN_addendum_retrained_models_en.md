# Analysis Plan Addendum: Retrained PRESAGE, scLambda, and fMLP-GenePT

The original addendum was finalized on 2026-11-02, before any predictions from these three models were opened. Its SHA-256 digest was subsequently recorded.

All definitions in the main plan `ANALYSIS_PLAN_target_gene.md` and the weighted R² addendum `ANALYSIS_PLAN_addendum_weighted_r2.md` carry over unchanged, including b_KD, the weights, the GT and the metric algorithms. The cited original SHA-256 digests, 5ab9ccac… and 61a1bf7c…, refer to the corresponding Chinese files.

This document is an English translation of the original addendum.

## Rationale

Limitation 5 of the main report: in the paper, PRESAGE, scLambda and fMLP-GenePT have a higher weighted R² than the trivial knockdown baseline b_KD, but without their prediction files only an unpaired comparison was possible. Here they are retrained with the authors' code, and once the predictions are available, paired comparisons and the target-gene decomposition are carried out.

## Known before the addendum was finalized

- Training follows the docker wrappers and `configs/model/*.yaml` hyperparameters of CellSimBench commit 15361aec and runs on a Colab A100, starting with Replogle fold 0.
- fMLP-GenePT has been trained and has finished predicting; part of its prediction file has been downloaded as bytes but has not been decoded or opened. PRESAGE has finished training and is predicting.
- scLambda's first run stopped because the BioMart of the main Ensembl site was unavailable. In the rerun, the gene descriptions are fetched from the Ensembl 115 archive (sep2025.archive.ensembl.org), with everything else unchanged; 1969 of the 2044 perturbed genes have a description, and the embeddings of the other 75 are zero. This is a known deviation.
- The paper's five-fold means (Replogle): WMSE PRESAGE 0.028, scLambda 0.029, fMLP-GenePT 0.029; WR²(ΔCtrl / ΔPert) PRESAGE 0.108 / 0.134, scLambda 0.066 / 0.092, fMLP-GenePT 0.047 / 0.076.
- Gap between a single fold and the five-fold mean, estimated from existing results: for scGPT retrained in the earlier analysis, fold 0 WR² exceeded the five-fold mean by 0.026 / 0.025, while WMSE was 5% lower.

## Retraining fidelity gate (judged separately for each model)

Computed on all test conditions of Replogle fold 0, in the same way as the paper's tables (the model's own scored gene axis, conditions weighted equally). A model counts as "faithfully retrained" only if all three of the following hold:

- |WR²(ΔCtrl) − paper value| ≤ 0.03;
- |WR²(ΔPert) − paper value| ≤ 0.03;
- |WMSE − paper value| / paper value ≤ 10%.

A model that does not pass the gate still undergoes the paired comparisons, but the conclusions then concern only "our retrained version" and do not represent the model reported in the paper.

## Endpoints (reported for each model separately, fold 0 first)

Conditions weighted equally, 10,000 bootstrap resamples, seed 20261102. The code is `20_target_gene.py` and `22_weighted_r2.py`, with only these three models added to the prediction dictionary.

- **E3:** WMSE(b_KD) − WMSE(f).
- **E3′:** WR²(f) − WR²(b_KD), one each for ΔCtrl and ΔPert.
- **E1:** the target-gene share.
- **E2, E2′:** the advantage after excluding the target gene, one for WMSE and one for each WR².
- **Descriptive:** the medians of f_t/b_t and y_t/b_t; the number of test conditions missing from the predictions.

## Decision rules

- D1, D2, D2′ and D3 are the same as in the earlier plans.
- **D4 (new):** if the interval of E2 is entirely > 0 and the intervals of both E2′ are also entirely > 0, conclude "the model has signal beyond the target gene".
- **Interpretation for the paper's main conclusion:**
  - If, for a model, D2′ concludes "the model is better than b_KD" and D4 holds: the paper's weighted R² advantage stands for this model, and our challenge is limited to WMSE and to scGPT-like models.
  - If D2′ concludes "b_KD is not significantly worse than the model": the challenge extends to this model.

## Scope

- If all three models pass the fidelity gate on fold 0, Replogle folds 1–4 and Nadig will be added (as compute allows), with the code and decision rules unchanged, and the five folds will be pooled and reported again.
- Definitions, endpoints and thresholds will not be changed after the results are seen; any deviation will be labelled separately.
