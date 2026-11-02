"""Same-sample vs cross-half scoring at the native budget for any dataset (WMSE_DS).
For Nadig this evaluates the locked hypotheses H1-H3 and R1 (plans/ANALYSIS_PLAN_nadig.md)."""
import os, json, time, hashlib
import numpy as np, h5py, pandas as pd
from common import H5, CACHE, RESULTS, DS, PRED, PRED_SHA256, dec, read_df_rows, load_prediction, author_weights_from_scores

t0 = time.time()
ax = json.load(open(os.path.join(CACHE, 'axes.json')))
cats, genes, ctrl = ax['conditions'], ax['genes'], ax['control_index']
ci = {c: i for i, c in enumerate(cats)}; G = len(genes); gpos = {g: j for j, g in enumerate(genes)}
mean = np.load(os.path.join(CACHE, 'mean.npy'), mmap_mode='r'); T = np.load(os.path.join(CACHE, 'T.npy'), mmap_mode='r'); N = np.load(os.path.join(CACHE, 'N.npy'))
f = h5py.File(H5, 'r'); codes = f['obs/condition/codes'][:]
ok_pred = {}
for key, p in PRED.items():
    if PRED_SHA256.get(key):
        ok_pred[key] = hashlib.sha256(open(p, 'rb').read()).hexdigest() == PRED_SHA256[key]
def score_set(fv, b, qa, qb, ya, yb):
    qa = qa / qa.sum(); qb = qb / qb.sum()
    g = lambda q, y: float(np.sum(q * ((b - y) ** 2 - (fv - y) ** 2)))
    gaa, gbb, gab, gba = g(qa, ya), g(qb, yb), g(qa, yb), g(qb, ya)
    return {'G_nat': gaa, 'S': (gaa + gbb) / 2, 'C': (gab + gba) / 2, 'I': (gaa + gbb - gab - gba) / 2,
            'I_identity_err': (gaa + gbb - gab - gba) / 2 - float(np.sum((fv - b) * (qa - qb) * (ya - yb)))}
recs = []
for k in range(5):
    sp = f[f'obs/split_fold_{k}']; scat = dec(sp['categories'][:]); scode = sp['codes'][:]
    roles = {}
    for c_code, s_code in zip(codes, scode):
        roles.setdefault(c_code, set()).add(scat[s_code])
    test = [cats[c] for c, r in roles.items() if c != ctrl and r == {'test'}]
    b_all = read_df_rows(f['uns'][f'split_fold_{k}_mean_baseline'])[2][0]
    lidx, _, Ml = read_df_rows(f['uns'][f'split_fold_{k}_linear_baseline']); lmap = {x: i for i, x in enumerate(lidx)}
    models = {}
    for (m, kk) in PRED:
        if kk == k and ok_pred.get((m, kk), True):
            pc, pg, PX = load_prediction(m, k)
            pc = [c.split('+')[0] for c in pc]; pc = [c[len(DS) + 1:] if c.startswith(DS + '_') else c for c in pc]
            models[m] = ({c: i for i, c in enumerate(pc)}, np.arange(G) if pg == genes else np.array([gpos[x] for x in pg]), PX)
    for c in test:
        i = ci[c]; qa_all = author_weights_from_scores(T[i, 0]); qb_all = author_weights_from_scores(T[i, 1])
        if not (np.isfinite(qa_all).all() and np.isfinite(qb_all).all() and qa_all.sum() > 0 and qb_all.sum() > 0):
            recs.append({'fold': k, 'condition': c, 'model': 'INVALID_PANEL'}); continue
        ya_all, yb_all = np.asarray(mean[i, 0]), np.asarray(mean[i, 1])
        todo = [(m, models[m][2][models[m][0][c]], models[m][1]) for m in models if c in models[m][0]]
        if f'{DS}_{c}' in lmap:
            todo.append(('Linear', Ml[lmap[f'{DS}_{c}']], np.arange(G)))
        for m, fv, sel in todo:
            r = score_set(fv, b_all[sel], qa_all[sel], qb_all[sel], ya_all[sel], yb_all[sel])
            r.update({'fold': k, 'condition': c, 'model': m, 'nA': int(N[i, 0]), 'nB': int(N[i, 1])}); recs.append(r)
    print(f'fold {k} done t={time.time() - t0:.0f}s models={list(models)}', flush=True)
df = pd.DataFrame(recs); df.to_csv(os.path.join(RESULTS, 'samecross_native_per_condition.csv'), index=False)
rng = np.random.default_rng(20261102)
def boot_means(x, cols, B=10000):
    v = x[cols].values; n = len(v); out = np.empty((B, len(cols)))
    for s in range(0, B, 1000):
        idx = rng.integers(0, n, size=(min(1000, B - s), n)); out[s:s + idx.shape[0]] = v[idx].mean(1)
    return out
summary = {'plan_sha256': 'f33e1640465dc985be1a426b08869ea071fe377b89734f4037a7a0c6d0d7771c' if DS == 'nadig25hepg2' else 'eb14e9537bdb582bb7ff561cb54349f66b7ac50d4dce8b68a135aa55b9351ae2',
           'invalid_panels': int((df['model'] == 'INVALID_PANEL').sum()), 'max_abs_identity_error': float(df['I_identity_err'].abs().max()),
           'prediction_sha_ok': {f'{m}|fold{k}': v for (m, k), v in ok_pred.items()}}
d2 = df[df.model != 'INVALID_PANEL']
for m in sorted(d2.model.unique()):
    x = d2[d2.model == m]; cols = ['G_nat', 'S', 'C', 'I']; bm = boot_means(x, cols)
    e = {'n_conditions': int(len(x)), 'folds': sorted(int(v) for v in x.fold.unique())}
    for j, col in enumerate(cols):
        e[col] = {'mean': float(x[col].mean()), 'ci95': [float(np.percentile(bm[:, j], 2.5)), float(np.percentile(bm[:, j], 97.5))]}
    ratio = bm[:, 3] / bm[:, 1]
    e['I_over_S'] = {'value': float(x['I'].mean() / x['S'].mean()), 'ci95': [float(np.percentile(ratio, 2.5)), float(np.percentile(ratio, 97.5))]}
    e['I_per_fold'] = [float(v) for v in x.groupby('fold')['I'].mean()]
    e['share_conditions_I_positive'] = float((x['I'] > 0).mean())
    e['R1_advantage_survives_decoupling'] = bool(e['C']['ci95'][0] > 0)
    e['H1_or_R2_coupling_detectable'] = bool(e['I']['ci95'][0] > 0)
    e['H2_ratio_above_replogle_native_0.124'] = bool(e['I_over_S']['ci95'][0] > 0.124)
    e['H3_or_R3_material_I_ge_0.0005'] = bool(e['I']['mean'] >= 0.0005)
    summary[m] = e
summary['seconds'] = round(time.time() - t0, 1)
json.dump(summary, open(os.path.join(RESULTS, 'samecross_native_summary.json'), 'w'), indent=1)
print(json.dumps(summary, indent=1))
