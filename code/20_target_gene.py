"""Target-gene analysis (locked plan plans/ANALYSIS_PLAN_target_gene.md, sha256 5ab9ccac...):
contribution of the perturbed (target) gene to native WMSE advantages, target-excluded metric,
and a trivial knockdown-aware baseline built from training conditions only. Run per WMSE_DS."""
import os, json, time
import numpy as np, h5py, pandas as pd
from common import H5, CACHE, DS, PRED, PRED_SET, BASE, dec, read_df_rows, load_prediction, author_weights_from_scores

t0 = time.time()
SUF = '' if PRED_SET == 'scgpt_gears' else '_retrained'   # retrained models get their own output files
FOLDS = range(5) if PRED_SET == 'scgpt_gears' else sorted({kk for (_, kk) in PRED})
OUT = os.path.join(BASE, 'results/target_gene'); os.makedirs(OUT, exist_ok=True)
ax = json.load(open(os.path.join(CACHE, 'axes.json')))
cats, genes, ctrl = ax['conditions'], ax['genes'], ax['control_index']
ci = {c: i for i, c in enumerate(cats)}; G = len(genes); gpos = {g: j for j, g in enumerate(genes)}
Tauth = np.load(os.path.join(CACHE, 'T_author_first_half.npy'))
mean = np.load(os.path.join(CACHE, 'mean.npy'), mmap_mode='r'); N = np.load(os.path.join(CACHE, 'N.npy'))
f = h5py.File(H5, 'r'); codes = f['obs/condition/codes'][:]
def rows(key):
    idx, cols, M = read_df_rows(f['uns'][key]); assert cols == genes, key
    return {k: M[i] for i, k in enumerate(idx)}
gt = rows('technical_duplicate_first_half_baseline')

def weights_excl(tvec, ti):
    a = np.abs(np.asarray(tvec, float)); m = np.ones(G, bool); m[ti] = False
    w = np.zeros(G); w[m] = ((a[m] - a[m].min()) / (a[m].max() - a[m].min())) ** 2
    return np.nan_to_num(w, nan=0.0)

recs, kd = [], {}
for k in FOLDS:
    sp = f[f'obs/split_fold_{k}']; scat = dec(sp['categories'][:]); scode = sp['codes'][:]
    roles = {}
    for c_code, s_code in zip(codes, scode):
        roles.setdefault(c_code, set()).add(scat[s_code])
    test = [cats[c] for c, r in roles.items() if c != ctrl and r == {'test'}]
    train = [cats[c] for c, r in roles.items() if c != ctrl and r == {'train'}]
    b = read_df_rows(f['uns'][f'split_fold_{k}_mean_baseline'])[2][0]
    # knockdown factor from TRAINING conditions only (all cells of each training condition)
    ratios, diffs = [], []
    for c in train:
        if c in gpos:
            i, ti = ci[c], gpos[c]
            ybar = (mean[i, 0, ti] * N[i, 0] + mean[i, 1, ti] * N[i, 1]) / (N[i, 0] + N[i, 1])
            diffs.append(ybar - b[ti])
            if b[ti] > 0:
                ratios.append(ybar / b[ti])
    fk, dk = float(np.median(ratios)), float(np.median(diffs))
    kd[k] = {'f_k_median_ratio': fk, 'd_k_median_shift': dk, 'n_training_conditions_used': len(ratios)}
    lidx, _, Ml = read_df_rows(f['uns'][f'split_fold_{k}_linear_baseline']); lmap = {x: i for i, x in enumerate(lidx)}
    models = {}
    for (m, kk) in PRED:
        if kk == k:
            pc, pg, PX = load_prediction(m, k)
            pc = [c.split('+')[0] for c in pc]; pc = [c[len(DS) + 1:] if c.startswith(DS + '_') else c for c in pc]
            models[m] = ({c: i for i, c in enumerate(pc)}, np.arange(G) if pg == genes else np.array([gpos[x] for x in pg]), PX)
    for c in test:
        i = ci[c]; key = f'{DS}_{c}'
        q_nat = author_weights_from_scores(Tauth[i]); y_all = gt[key]
        ti = gpos.get(c)
        q_x = weights_excl(Tauth[i], ti) if ti is not None else q_nat
        todo = [(m, models[m][2][models[m][0][c]], models[m][1]) for m in models if c in models[m][0]]
        if key in lmap:
            todo.append(('Linear', Ml[lmap[key]], np.arange(G)))
        for m, fv, sel in todo:
            q = q_nat[sel] / q_nat[sel].sum(); qx = q_x[sel] / q_x[sel].sum()
            y, bb = y_all[sel], b[sel]
            g = (bb - y) ** 2 - (fv - y) ** 2
            r = {'fold': k, 'condition': c, 'model': m,
                 'WMSE_b': float(np.sum(q * (bb - y) ** 2)), 'WMSE_f': float(np.sum(q * (fv - y) ** 2)), 'G_nat': float(np.sum(q * g)),
                 'G_X': float(np.sum(qx * g))}
            pos = None
            if ti is not None:
                hit = np.flatnonzero(sel == ti); pos = int(hit[0]) if len(hit) else None
            if pos is not None:
                bkd = bb.copy(); bkd[pos] = fk * bb[pos]
                bka = bb.copy(); bka[pos] = bb[pos] + dk
                r.update({'target_in_axis': True, 'q_t': float(q[pos]), 'target_rank1': bool(np.argmax(q_nat) == ti),
                          'contrib_t': float(q[pos] * g[pos]), 'share_of_WMSE_b_from_t': float(q[pos] * (bb[pos] - y[pos]) ** 2 / np.sum(q * (bb - y) ** 2)),
                          'WMSE_bKD': float(np.sum(q * (bkd - y) ** 2)), 'WMSE_bKD_add': float(np.sum(q * (bka - y) ** 2)),
                          'true_ratio_t': float(y[pos] / bb[pos]) if bb[pos] > 0 else np.nan,
                          'pred_ratio_t': float(fv[pos] / bb[pos]) if bb[pos] > 0 else np.nan})
            else:
                r.update({'target_in_axis': False, 'q_t': 0.0, 'target_rank1': False, 'contrib_t': 0.0, 'share_of_WMSE_b_from_t': 0.0,
                          'WMSE_bKD': r['WMSE_b'], 'WMSE_bKD_add': r['WMSE_b'], 'true_ratio_t': np.nan, 'pred_ratio_t': np.nan})
            recs.append(r)
    print(f'fold {k}: f_k={fk:.3f} d_k={dk:.3f} models={list(models) + ["Linear"]} t={time.time() - t0:.0f}s', flush=True)

df = pd.DataFrame(recs)
df.to_csv(os.path.join(OUT, f'{DS}_per_condition{SUF}.csv'), index=False)
rng = np.random.default_rng(20261102)
def boot_stats(x, fns, B=10000):
    n = len(x); out = np.empty((B, len(fns)))
    arrs = {c: x[c].values for c in x.columns if x[c].dtype != object}
    for s in range(0, B, 500):
        idx = rng.integers(0, n, size=(min(500, B - s), n))
        for j, fn in enumerate(fns):
            out[s:s + idx.shape[0], j] = fn({c: v[idx] for c, v in arrs.items()})
    return out
summ = {'plan_sha256': '5ab9ccacc7739f6e1f069502211c2441a35e7628d22d40f239c19941d5b474ea', 'dataset': DS, 'knockdown_factor_by_fold': kd}
for m in sorted(df.model.unique()):
    if PRED_SET == 'retrained':
        rng = np.random.default_rng(20261102)
    x = df[df.model == m].reset_index(drop=True)
    x['adv_vs_bKD'] = x['WMSE_bKD'] - x['WMSE_f']; x['adv_vs_bKD_add'] = x['WMSE_bKD_add'] - x['WMSE_f']; x['b_minus_bKD'] = x['WMSE_b'] - x['WMSE_bKD']
    fns = [lambda d: d['contrib_t'].mean(1) / d['G_nat'].mean(1), lambda d: d['G_nat'].mean(1), lambda d: d['G_X'].mean(1),
           lambda d: d['adv_vs_bKD'].mean(1), lambda d: d['adv_vs_bKD_add'].mean(1), lambda d: d['b_minus_bKD'].mean(1) / d['G_nat'].mean(1)]
    bs = boot_stats(x, fns)
    ci95 = lambda j: [float(np.percentile(bs[:, j], 2.5)), float(np.percentile(bs[:, j], 97.5))]
    tin = x[x.target_in_axis]
    e = {'n_conditions': int(len(x)), 'n_target_in_axis': int(len(tin)),
         'E1_share_target': {'value': float(x.contrib_t.mean() / x.G_nat.mean()), 'ci95': ci95(0)},
         'G_nat': {'mean': float(x.G_nat.mean()), 'ci95': ci95(1)},
         'E2_G_X_target_excluded': {'mean': float(x.G_X.mean()), 'ci95': ci95(2)},
         'E3_adv_vs_bKD': {'mean': float(x.adv_vs_bKD.mean()), 'ci95': ci95(3)},
         'E3_adv_vs_bKD_additive': {'mean': float(x.adv_vs_bKD_add.mean()), 'ci95': ci95(4)},
         'E3_gap_closed_by_bKD': {'value': float(x.b_minus_bKD.mean() / x.G_nat.mean()), 'ci95': ci95(5)},
         'desc_target_rank1_frac_among_in_axis': float(tin.target_rank1.mean()),
         'desc_q_t_median_among_in_axis': float(tin.q_t.median()),
         'desc_share_of_WMSE_b_from_target_mean': float(x.share_of_WMSE_b_from_t.mean()),
         'desc_true_ratio_t_median': float(tin.true_ratio_t.median()), 'desc_pred_ratio_t_median': float(tin.pred_ratio_t.median())}
    e['D1_target_major_source'] = bool(e['E1_share_target']['value'] >= 0.5)
    lo3, hi3 = e['E3_adv_vs_bKD']['ci95']
    e['D2'] = 'headline_challenged' if lo3 <= 0 else 'advantage_not_explained_by_trivial_knockdown'
    e['D3_advantage_survives_target_exclusion'] = bool(e['E2_G_X_target_excluded']['ci95'][0] > 0)
    summ[m] = e
summ['seconds'] = round(time.time() - t0, 1)
json.dump(summ, open(os.path.join(OUT, f'{DS}_summary{SUF}.json'), 'w'), indent=1)
print(json.dumps(summ, indent=1))
