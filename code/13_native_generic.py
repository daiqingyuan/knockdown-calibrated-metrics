"""Native (paper-protocol) WMSE for any dataset (WMSE_DS), all author baselines plus whatever frozen
model predictions are listed in common.PRED; compared with the published table. Used for Nadig."""
import os, json, time, hashlib
import numpy as np, h5py, pandas as pd
from common import H5, CACHE, RESULTS, DS, PRED, PRED_SHA256, dec, read_df_rows, load_prediction, author_weights_from_scores

PUB = {'replogle22k562': {'GEARS': 0.033, 'scGPT': 0.027, 'Dataset Mean': 0.034, 'Control Mean': 0.035, 'Linear': 0.034, 'Tech. Dup.': 0.011, 'Inter. Dup.': 0.013},
       'nadig25hepg2': {'GEARS': 0.046, 'scGPT': 0.040, 'Dataset Mean': 0.048, 'Control Mean': 0.051, 'Linear': 0.048, 'Tech. Dup.': 0.022, 'Inter. Dup.': 0.029}}[DS]
t0 = time.time()
ax = json.load(open(os.path.join(CACHE, 'axes.json')))
cats, genes, ctrl = ax['conditions'], ax['genes'], ax['control_index']
ci = {c: i for i, c in enumerate(cats)}; G = len(genes); gpos = {g: j for j, g in enumerate(genes)}
Tauth = np.load(os.path.join(CACHE, 'T_author_first_half.npy'))
f = h5py.File(H5, 'r'); codes = f['obs/condition/codes'][:]
def rows(key):
    idx, cols, M = read_df_rows(f['uns'][key]); assert cols == genes, key
    return {k: M[i] for i, k in enumerate(idx)}
gt, td, idup = rows('technical_duplicate_first_half_baseline'), rows('technical_duplicate_second_half_baseline'), rows('interpolated_duplicate_baseline')
ctrl_b = rows('ctrl_baseline')[DS]
sha_ok = {}
for (m, k), p in PRED.items():
    if PRED_SHA256.get((m, k)):
        h = hashlib.sha256(open(p, 'rb').read()).hexdigest()
        sha_ok[(m, k)] = (h == PRED_SHA256[(m, k)]); print('prediction sha256 ok:', m, k, sha_ok[(m, k)])
recs = []
for k in range(5):
    sp = f[f'obs/split_fold_{k}']; scat = dec(sp['categories'][:]); scode = sp['codes'][:]
    roles = {}
    for c_code, s_code in zip(codes, scode):
        roles.setdefault(c_code, set()).add(scat[s_code])
    test = [cats[c] for c, r in roles.items() if c != ctrl and r == {'test'}]
    mean_b = rows(f'split_fold_{k}_mean_baseline')[DS]; lin = rows(f'split_fold_{k}_linear_baseline')
    models = {}
    for (m, kk) in PRED:
        if kk != k or not sha_ok.get((m, kk), True):
            continue
        pc, pg, PX = load_prediction(m, k)
        pc = [c.split('+')[0] for c in pc]; pc = [c[len(DS) + 1:] if c.startswith(DS + '_') else c for c in pc]
        sel = None if pg == genes else np.array([gpos[x] for x in pg])
        models[m] = ({c: i for i, c in enumerate(pc)}, sel, PX)
    for c in test:
        key = f'{DS}_{c}'; w_all = author_weights_from_scores(Tauth[ci[c]]); y = gt[key]
        def wm(pred, sel=None):
            w = w_all if sel is None else w_all[sel]; t = y if sel is None else y[sel]
            w = w / w.sum(); return float(np.sum(w * (pred - t) ** 2))
        row = {'fold': k, 'condition': c, 'Dataset Mean': wm(mean_b), 'Control Mean': wm(ctrl_b),
               'Linear': wm(lin[key]) if key in lin else np.nan, 'Tech. Dup.': wm(td[key]), 'Inter. Dup.': wm(idup[key])}
        for m, (pmap, sel, PX) in models.items():
            if c in pmap:
                row[m] = wm(PX[pmap[c]], sel); row[m + '|DatasetMean_same_genes'] = wm(mean_b if sel is None else mean_b[sel], sel)
        recs.append(row)
    print(f'fold {k}: {len(test)} test conditions; models {list(models)}', flush=True)
df = pd.DataFrame(recs); df.to_csv(os.path.join(RESULTS, 'native_wmse_per_condition.csv'), index=False)
out = {}
for m in [x for x in ['GEARS', 'scGPT', 'Dataset Mean', 'Control Mean', 'Linear', 'Tech. Dup.', 'Inter. Dup.'] if x in df]:
    v = df[m].dropna(); pf = df.groupby('fold')[m].mean()
    out[m] = {'published': PUB[m], 'reproduced_mean': round(float(v.mean()), 5), 'abs_diff': round(abs(float(v.mean()) - PUB[m]), 5),
              'sem_over_conditions': round(float(v.std(ddof=1) / np.sqrt(len(v))), 5), 'n_conditions': int(len(v)), 'per_fold': [round(float(x), 5) for x in pf.dropna()]}
    if m + '|DatasetMean_same_genes' in df:
        d = (df[m + '|DatasetMean_same_genes'] - df[m]).dropna()
        out[m]['paired_gap_vs_DatasetMean'] = {'mean': round(float(d.mean()), 5), 'sem': round(float(d.std(ddof=1) / np.sqrt(len(d))), 5)}
gate = all(out[m]['abs_diff'] <= 0.001 for m in ['Dataset Mean', 'Control Mean', 'Linear', 'Tech. Dup.', 'Inter. Dup.'])
out['_gate2_baselines_within_0.001'] = gate
out['_seconds'] = round(time.time() - t0, 1)
json.dump(out, open(os.path.join(RESULTS, 'native_wmse_summary.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
