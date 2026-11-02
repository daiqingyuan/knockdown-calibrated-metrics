"""Locked addendum plans/ANALYSIS_PLAN_addendum_retrained_models.md (sha256 decdbd9e...): retraining
fidelity gate, D1-D3 / D2' (from 20 and 22 run with WMSE_PRED_SET=retrained) and the new D4, per retrained model."""
import json, os
import numpy as np, h5py, pandas as pd
from common import H5, DS, BASE, PRED, dec

OUT = os.path.join(BASE, 'results/target_gene')
s20 = json.load(open(os.path.join(OUT, f'{DS}_summary_retrained.json')))
s22 = json.load(open(os.path.join(OUT, f'{DS}_weighted_r2_summary_retrained.json')))
pc = pd.read_csv(os.path.join(OUT, f'{DS}_per_condition_retrained.csv'))
PUB = {'replogle22k562': {'PRESAGE': (0.028, 0.108, 0.134), 'scLambda': (0.029, 0.066, 0.092), 'fMLP-GenePT': (0.029, 0.047, 0.076)},
       'nadig25hepg2': {'PRESAGE': (0.039, 0.087, 0.140), 'scLambda': (0.042, 0.062, 0.102), 'fMLP-GenePT': (0.041, 0.051, 0.096)}}[DS]
with h5py.File(H5, 'r') as f:
    cats = dec(f['obs/condition/categories'][:]); codes = f['obs/condition/codes'][:]
    n_test = {}
    for k in sorted({k for (_, k) in PRED}):
        sp = f[f'obs/split_fold_{k}']; scat = dec(sp['categories'][:]); roles = {}
        for c, s in zip(codes, sp['codes'][:]):
            roles.setdefault(c, set()).add(scat[s])
        n_test[k] = sum(1 for c, r in roles.items() if cats[c] != 'control' and r == {'test'})
rep = {'addendum_plan_sha256': 'decdbd9e81329732d38a40f3b6a3a5244773cce89978f21f00dca962e906c598', 'dataset': DS,
       'folds': sorted({k for (_, k) in PRED}), 'n_test_conditions_by_fold': n_test, 'models': {}}
for m in sorted({m for (m, _) in PRED}):
    e20, e22 = s20[m], s22[m]
    x = pc[pc.model == m]
    wmse = float(x.WMSE_f.mean()); r2c, r2p = e22['ctrl']['model_mean'], e22['pert']['model_mean']
    pw, pc_, pp = PUB[m]
    gate = {'WMSE': {'ours': round(wmse, 5), 'published': pw, 'rel_diff': round((wmse - pw) / pw, 4)},
            'WR2_ctrl': {'ours': r2c, 'published': pc_, 'diff': round(r2c - pc_, 4)},
            'WR2_pert': {'ours': r2p, 'published': pp, 'diff': round(r2p - pp, 4)}}
    gate['pass'] = bool(abs(gate['WMSE']['rel_diff']) <= 0.10 and abs(gate['WR2_ctrl']['diff']) <= 0.03 and abs(gate['WR2_pert']['diff']) <= 0.03)
    d4 = bool(e20['E2_G_X_target_excluded']['ci95'][0] > 0 and e22['ctrl']['E2prime_target_excluded_model_minus_mean']['ci95'][0] > 0
              and e22['pert']['E2prime_target_excluded_model_minus_mean']['ci95'][0] > 0)
    mfolds = sorted({k for (mm, k) in PRED if mm == m})   # a model may have fewer folds than others
    rep['models'][m] = {
        'folds': mfolds, 'n_conditions': e20['n_conditions'],
        'n_missing_vs_test_set': int(sum(n_test[k] for k in mfolds) - e20['n_conditions']),
        'fidelity_gate': gate,
        'E3_WMSE_adv_vs_bKD': e20['E3_adv_vs_bKD'], 'D2': e20['D2'],
        'E3prime_WR2ctrl_model_minus_KD': e22['ctrl']['E3prime_model_minus_KD'], 'D2prime_ctrl': e22['ctrl']['D2prime'],
        'E3prime_WR2pert_model_minus_KD': e22['pert']['E3prime_model_minus_KD'], 'D2prime_pert': e22['pert']['D2prime'],
        'E1_share_target': e20['E1_share_target'], 'D1': e20['D1_target_major_source'],
        'E2_G_X': e20['E2_G_X_target_excluded'], 'D3': e20['D3_advantage_survives_target_exclusion'],
        'E2prime_ctrl': e22['ctrl']['E2prime_target_excluded_model_minus_mean'], 'E2prime_pert': e22['pert']['E2prime_target_excluded_model_minus_mean'],
        'D4_signal_beyond_target': d4,
        'target_ratio_median_pred_vs_true': [e20['desc_pred_ratio_t_median'], e20['desc_true_ratio_t_median']]}
json.dump(rep, open(os.path.join(OUT, f'{DS}_retrained_report.json'), 'w'), indent=1)
print(json.dumps(rep, indent=1))
