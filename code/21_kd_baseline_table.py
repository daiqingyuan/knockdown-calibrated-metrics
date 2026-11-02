"""EXPLORATORY follow-up to the locked target-gene plan: native WMSE of the knockdown-aware trivial baseline
b_KD (as defined in the plan, plus an exploratory 'target set to zero' variant) over ALL test conditions on the
full gene axis, for comparison with the published table; plus per-fold scGPT-vs-b_KD from the per-condition file."""
import os, json
import numpy as np, h5py, pandas as pd
from common import H5, CACHE, DS, BASE, dec, read_df_rows, author_weights_from_scores

OUT = os.path.join(BASE, 'results/target_gene')
ax = json.load(open(os.path.join(CACHE, 'axes.json'))); cats, genes, ctrl = ax['conditions'], ax['genes'], ax['control_index']
ci = {c: i for i, c in enumerate(cats)}; gpos = {g: j for j, g in enumerate(genes)}
Tauth = np.load(os.path.join(CACHE, 'T_author_first_half.npy'))
mean = np.load(os.path.join(CACHE, 'mean.npy'), mmap_mode='r'); N = np.load(os.path.join(CACHE, 'N.npy'))
f = h5py.File(H5, 'r'); codes = f['obs/condition/codes'][:]
idx, cols, M = read_df_rows(f['uns']['technical_duplicate_first_half_baseline']); gt = {k: M[i] for i, k in enumerate(idx)}
summ = json.load(open(os.path.join(OUT, f'{DS}_summary.json')))
rows = []
for k in range(5):
    sp = f[f'obs/split_fold_{k}']; scat = dec(sp['categories'][:]); scode = sp['codes'][:]
    roles = {}
    for c_code, s_code in zip(codes, scode):
        roles.setdefault(c_code, set()).add(scat[s_code])
    test = [cats[c] for c, r in roles.items() if c != ctrl and r == {'test'}]
    b = read_df_rows(f['uns'][f'split_fold_{k}_mean_baseline'])[2][0]
    fk = summ['knockdown_factor_by_fold'][str(k)]['f_k_median_ratio']; dk = summ['knockdown_factor_by_fold'][str(k)]['d_k_median_shift']
    for c in test:
        q = author_weights_from_scores(Tauth[ci[c]]); q = q / q.sum(); y = gt[f'{DS}_{c}']
        r = {'fold': k, 'condition': c, 'target_in_axis': c in gpos, 'Dataset Mean': float(np.sum(q * (b - y) ** 2))}
        for name, val in (('KD multiplicative (primary)', None), ('KD additive', None), ('KD zero (exploratory)', None)):
            bb = b.copy()
            if c in gpos:
                t = gpos[c]
                bb[t] = fk * b[t] if name.startswith('KD multiplicative') else (b[t] + dk if name == 'KD additive' else 0.0)
            r[name] = float(np.sum(q * (bb - y) ** 2))
        rows.append(r)
df = pd.DataFrame(rows); df.to_csv(os.path.join(OUT, f'{DS}_kd_baseline_all_conditions.csv'), index=False)
pub = {'replogle22k562': {'best published model': ('scGPT', 0.027), 'Dataset Mean': 0.034, 'all models': {'CellFlow': 0.032, 'fMLP-ESM2': 0.030, 'fMLP-Geneformer': 0.032, 'fMLP-GenePT': 0.029, 'fMLP-scGPT': 0.032, 'GEARS': 0.033, 'PRESAGE': 0.028, 'scGPT': 0.027, 'scLambda': 0.029}},
       'nadig25hepg2': {'best published model': ('PRESAGE', 0.039), 'Dataset Mean': 0.048, 'all models': {'CellFlow': 0.046, 'fMLP-ESM2': 0.043, 'fMLP-Geneformer': 0.045, 'fMLP-GenePT': 0.041, 'fMLP-scGPT': 0.044, 'GEARS': 0.046, 'PRESAGE': 0.039, 'scGPT': 0.040, 'scLambda': 0.042}}}[DS]
out = {'dataset': DS, 'n_test_conditions': int(len(df)), 'n_target_in_axis': int(df.target_in_axis.sum()), 'published_models_WMSE': pub['all models']}
for col in ('Dataset Mean', 'KD multiplicative (primary)', 'KD additive', 'KD zero (exploratory)'):
    v = df[col]; out[col] = {'WMSE_mean': round(float(v.mean()), 5), 'sem': round(float(v.std(ddof=1) / np.sqrt(len(v))), 5),
                             'per_fold': [round(float(x), 5) for x in df.groupby('fold')[col].mean()]}
d = df['Dataset Mean'] - df['KD multiplicative (primary)']
out['KD_primary_minus_DatasetMean_paired'] = {'mean': round(float(-d.mean()), 5), 'sem': round(float(d.std(ddof=1) / np.sqrt(len(d))), 5)}
pc = pd.read_csv(os.path.join(OUT, f'{DS}_per_condition.csv'))
for m in sorted(pc.model.unique()):
    x = pc[pc.model == m]
    out[f'{m}: per-fold mean (WMSE_bKD - WMSE_model)'] = [round(float(v), 5) for v in (x.WMSE_bKD - x.WMSE_f).groupby(x.fold).mean()]
    out[f'{m}: share of conditions where model beats b_KD'] = round(float(((x.WMSE_bKD - x.WMSE_f) > 0).mean()), 3)
json.dump(out, open(os.path.join(OUT, f'{DS}_kd_baseline_table.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
