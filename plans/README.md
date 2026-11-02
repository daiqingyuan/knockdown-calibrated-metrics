# Analysis plans

Each analysis plan was written in Chinese and locked before the outcomes it covers were computed; its SHA-256 digest was recorded at the time shown (UTC). This folder holds English translations of these plans. The Chinese originals are not included; their digests are listed below. Each translation has its own digest file (`*_en.sha256`), computed when the translation was made.

| Plan | Covers | Original SHA-256 | Locked (UTC) | English translation |
|---|---|---|---|---|
| `ANALYSIS_PLAN_replogle.md` | Weight–truth coupling, Replogle K562 (discovery) | `eb14e9537bdb582bb7ff561cb54349f66b7ac50d4dce8b68a135aa55b9351ae2` | 2026-11-01T16:55:38Z | `ANALYSIS_PLAN_replogle_en.md` |
| `ANALYSIS_PLAN_nadig.md` | Weight–truth coupling, Nadig HepG2 (confirmation) | `f33e1640465dc985be1a426b08869ea071fe377b89734f4037a7a0c6d0d7771c` | 2026-11-01T17:44:09Z | `ANALYSIS_PLAN_nadig_en.md` |
| `ANALYSIS_PLAN_target_gene.md` | Target-gene decomposition and the knockdown-aware baseline (WMSE) | `5ab9ccacc7739f6e1f069502211c2441a35e7628d22d40f239c19941d5b474ea` | 2026-11-01T22:28:29Z | `ANALYSIS_PLAN_target_gene_en.md` |
| `ANALYSIS_PLAN_addendum_weighted_r2.md` | Weighted R² addendum | `61a1bf7ca469547abd1f1c7c6449576ecb8ec0ed10768389f85ea8f16f197c1e` | 2026-11-01T23:00:26Z | `ANALYSIS_PLAN_addendum_weighted_r2_en.md` |
| `ANALYSIS_PLAN_addendum_retrained_models.md` | Retrained PRESAGE, scLambda and fMLP-GenePT | `decdbd9e81329732d38a40f3b6a3a5244773cce89978f21f00dca962e906c598` | 2026-11-02T00:49:31Z | `ANALYSIS_PLAN_addendum_retrained_models_en.md` |

File paths, script names and report references in the plans refer to the working directory at the time. In this repository the prediction files are listed in a prediction index rather than in `code/common.py` (see the main README).

To check a translation, run `sha256sum -c ANALYSIS_PLAN_replogle_en.sha256` in this folder; the same works for the other four.
