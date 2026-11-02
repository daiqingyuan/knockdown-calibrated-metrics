"""Target-gene analysis, weighted R2 addendum (locked: plans/ANALYSIS_PLAN_addendum_weighted_r2.md, sha256 61a1bf7c...):
author-faithful weighted R2 on deltas (vs control mean and vs fold dataset mean) for baselines, models and the
knockdown-aware trivial baseline b_KD; reproduction gate, paired model-vs-b_KD differences, target-excluded check."""
import os, json, time
import numpy as np, h5py, pandas as pd
from common import H5, CACHE, DS, PRED, PRED_SET, BASE, dec, read_df_rows, load_prediction, author_weights_from_scores

t0 = time.time()
SUF = '' if PRED_SET == 'scgpt_gears' else '_retrained'   # retrained models get their own output files
FOLDS = range(5) if PRED_SET == 'scgpt_gears' else sorted({kk for (_, kk) in PRED})
OUT = os.path.join(BASE, 'results/target_gene')
ax = json.load(open(os.path.join(CACHE, 'axes.json'))); cats, genes, ctrl = ax['conditions'], ax['genes'], ax['control_index']
ci = {c: i for i, c in enumerate(cats)}; G = len(genes); gpos = {g: j for j, g in enumerate(genes)}
Tauth = np.load(os.path.join(CACHE, 'T_author_first_half.npy'))
f = h5py.File(H5, 'r'); codes = f['obs/condition/codes'][:]
def rows(key):
    idx, cols, M = read_df_rows(f['uns'][key]); assert cols == genes, key
    return {k: M[i] for i, k in enumerate(idx)}
gt, td, idup = rows('technical_duplicate_first_half_baseline'), rows('technical_duplicate_second_half_baseline'), rows('interpolated_duplicate_baseline')
ctrl_b = rows('ctrl_baseline')[DS]
kdf = json.load(open(os.path.join(OUT, f'{DS}_summary.json')))['knockdown_factor_by_fold']

def r2w(y, yhat, w):
    w = np.asarray(w, float)
    if w.sum() == 0:
        w = np.ones_like(w)
    ybar = np.average(y, weights=w); num = np.sum(w * (y - yhat) ** 2); den = np.sum(w * (y - ybar) ** 2)
    r = (1.0 if num == 0 else 0.0) if den == 0 else 1.0 - num / den
    return max(r, -1.0)

def weights_excl(tvec, ti):
    a = np.abs(np.asarray(tvec, float)); m = np.ones(G, bool); m[ti] = False
    w = np.zeros(G); w[m] = ((a[m] - a[m].min()) / (a[m].max() - a[m].min())) ** 2
    return np.nan_to_num(w, nan=0.0)

base_rows, pair_rows = [], []
for k in FOLDS:
    sp = f[f'obs/split_fold_{k}']; scat = dec(sp['categories'][:]); scode = sp['codes'][:]
    roles = {}
    for c_code, s_code in zip(codes, scode):
        roles.setdefault(c_code, set()).add(scat[s_code])
    test = [cats[c] for c, r in roles.items() if c != ctrl and r == {'test'}]
    b = read_df_rows(f['uns'][f'split_fold_{k}_mean_baseline'])[2][0]
    lin = rows(f'split_fold_{k}_linear_baseline')
    fk, dk = kdf[str(k)]['f_k_median_ratio'], kdf[str(k)]['d_k_median_shift']
    models = {}
    for (m, kk) in PRED:
        if kk == k:
            pc, pg, PX = load_prediction(m, k)
            pc = [c.split('+')[0] for c in pc]; pc = [c[len(DS) + 1:] if c.startswith(DS + '_') else c for c in pc]
            models[m] = ({c: i for i, c in enumerate(pc)}, np.arange(G) if pg == genes else np.array([gpos[x] for x in pg]), PX)
    for c in test:
        i = ci[c]; key = f'{DS}_{c}'; w_all = author_weights_from_scores(Tauth[i]); y_all = gt[key]; ti = gpos.get(c)
        bkd = b.copy(); bka = b.copy()
        if ti is not None:
            bkd[ti] = fk * b[ti]; bka[ti] = b[ti] + dk
        preds_full = {'Dataset Mean': b, 'Control Mean': ctrl_b, 'Tech. Dup.': td[key], 'Inter. Dup.': idup[key], 'KD': bkd, 'KD_add': bka}
        if key in lin:
            preds_full['Linear'] = lin[key]
        row = {'fold': k, 'condition': c, 'target_in_axis': ti is not None}
        for name, p in preds_full.items():
            row[f'R2ctrl|{name}'] = r2w(y_all - ctrl_b, p - ctrl_b, w_all); row[f'R2pert|{name}'] = r2w(y_all - b, p - b, w_all)
        base_rows.append(row)
        wx_all = weights_excl(Tauth[i], ti) if ti is not None else w_all
        todo = [(m, models[m][2][models[m][0][c]], models[m][1]) for m in models if c in models[m][0]]
        if key in lin:
            todo.append(('Linear', lin[key], np.arange(G)))
        for m, fv, sel in todo:
            y, w, wx, cb, bb = y_all[sel], w_all[sel], wx_all[sel], ctrl_b[sel], b[sel]
            r = {'fold': k, 'condition': c, 'model': m}
            for tag, ref in (('ctrl', cb), ('pert', bb)):
                r[f'R2{tag}|f'] = r2w(y - ref, fv - ref, w)
                r[f'R2{tag}|KD'] = r2w(y - ref, bkd[sel] - ref, w)
                r[f'R2{tag}|KD_add'] = r2w(y - ref, bka[sel] - ref, w)
                r[f'R2{tag}X|f'] = r2w(y - ref, fv - ref, wx)
                r[f'R2{tag}X|mean'] = r2w(y - ref, bb - ref, wx)
            pair_rows.append(r)
    print(f'fold {k} done t={time.time() - t0:.0f}s models={list(models)}', flush=True)

bd = pd.DataFrame(base_rows); pdd = pd.DataFrame(pair_rows)
bd.to_csv(os.path.join(OUT, f'{DS}_weighted_r2_baselines_per_condition{SUF}.csv'), index=False)
pdd.to_csv(os.path.join(OUT, f'{DS}_weighted_r2_models_per_condition{SUF}.csv'), index=False)
PUB = {'replogle22k562': {'Dataset Mean': (-0.069, -0.035), 'Control Mean': (-0.037, -0.020), 'Linear': (-0.042, -0.020), 'Tech. Dup.': (0.319, 0.369), 'Inter. Dup.': (0.245, 0.289),
                          'models': {'PRESAGE': (0.108, 0.134), 'scLambda': (0.066, 0.092), 'scGPT': (0.068, 0.090), 'fMLP-GenePT': (0.047, 0.076), 'fMLP-ESM2': (0.016, 0.049), 'CellFlow': (0.009, 0.031), 'fMLP-Geneformer': (-0.029, 0.001), 'fMLP-scGPT': (-0.049, -0.016), 'GEARS': (-0.196, -0.149)}},
       'nadig25hepg2': {'Dataset Mean': (-0.091, -0.026), 'Control Mean': (-0.027, -0.007), 'Linear': (-0.062, -0.012), 'Tech. Dup.': (0.187, 0.267), 'Inter. Dup.': (0.070, 0.147),
                        'models': {'PRESAGE': (0.087, 0.140), 'scLambda': (0.062, 0.102), 'fMLP-GenePT': (0.051, 0.096), 'CellFlow': (-0.004, 0.038), 'fMLP-ESM2': (-0.009, 0.049), 'scGPT': (-0.014, 0.037), 'fMLP-Geneformer': (-0.033, 0.017), 'fMLP-scGPT': (-0.038, 0.019), 'GEARS': (-0.157, -0.081)}}}[DS]
out = {'plan_sha256': '61a1bf7ca469547abd1f1c7c6449576ecb8ec0ed10768389f85ea8f16f197c1e', 'dataset': DS, 'published_models_(WR2ctrl, WR2pert)': PUB['models'], 'baselines': {}}
gate = True
for name in ('Dataset Mean', 'Control Mean', 'Linear', 'Tech. Dup.', 'Inter. Dup.', 'KD', 'KD_add'):
    e = {}
    for j, tag in enumerate(('ctrl', 'pert')):
        v = bd[f'R2{tag}|{name}'].dropna()
        e[f'WR2_{tag}'] = {'mean': round(float(v.mean()), 4), 'sem': round(float(v.std(ddof=1) / np.sqrt(len(v))), 4), 'n': int(len(v))}
        if name in PUB:
            e[f'WR2_{tag}']['published'] = PUB[name][j]; d = abs(float(v.mean()) - PUB[name][j]); e[f'WR2_{tag}']['abs_diff'] = round(d, 4)
            gate = gate and d <= 0.01
    out['baselines'][name] = e
out['gate_baselines_within_0.01'] = gate
rng = np.random.default_rng(20261102)
def boot(v, B=10000):
    v = np.asarray(v); o = np.empty(B)
    for s in range(0, B, 1000):
        idx = rng.integers(0, len(v), size=(min(1000, B - s), len(v))); o[s:s + idx.shape[0]] = v[idx].mean(1)
    return [round(float(np.percentile(o, 2.5)), 4), round(float(np.percentile(o, 97.5)), 4)]
for m in sorted(pdd.model.unique()):
    if PRED_SET == 'retrained':
        rng = np.random.default_rng(20261102)
    x = pdd[pdd.model == m]; e = {'n_conditions': int(len(x))}
    for tag in ('ctrl', 'pert'):
        d = x[f'R2{tag}|f'] - x[f'R2{tag}|KD']; da = x[f'R2{tag}|f'] - x[f'R2{tag}|KD_add']; dx = x[f'R2{tag}X|f'] - x[f'R2{tag}X|mean']
        e[tag] = {'model_mean': round(float(x[f'R2{tag}|f'].mean()), 4), 'KD_mean_same_conditions': round(float(x[f'R2{tag}|KD'].mean()), 4),
                  'E3prime_model_minus_KD': {'mean': round(float(d.mean()), 4), 'ci95': boot(d.values)},
                  'model_minus_KD_additive': {'mean': round(float(da.mean()), 4), 'ci95': boot(da.values)},
                  'E2prime_target_excluded_model_minus_mean': {'mean': round(float(dx.mean()), 4), 'ci95': boot(dx.values)}}
        lo = e[tag]['E3prime_model_minus_KD']['ci95'][0]
        e[tag]['D2prime'] = 'model_better_than_KD' if lo > 0 else 'KD_not_significantly_worse'
    out[m] = e
out['seconds'] = round(time.time() - t0, 1)
json.dump(out, open(os.path.join(OUT, f'{DS}_weighted_r2_summary{SUF}.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
