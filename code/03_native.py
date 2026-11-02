"""Native (paper-protocol) WMSE on Replogle22 K562 for the frozen GEARS/scGPT predictions and the
author baselines, per fold, compared with the published Supplementary Table values."""
import os, json, time
import numpy as np, h5py, pandas as pd
from common import H5, CACHE, RESULTS, DS, dec, read_df_rows, load_prediction, author_weights_from_scores

t0 = time.time()
ax = json.load(open(os.path.join(CACHE, 'axes.json')))
cats, genes, ctrl = ax['conditions'], ax['genes'], ax['control_index']
ci = {c: i for i, c in enumerate(cats)}; G = len(genes)
Tauth = np.load(os.path.join(CACHE, 'T_author_first_half.npy'))
mine_mean = np.load(os.path.join(CACHE, 'mean.npy'), mmap_mode='r')
f = h5py.File(H5, 'r')

def rows_by_condition(key):
    idx, cols, M = read_df_rows(f['uns'][key])
    assert cols == genes, key
    return {k: M[i] for i, k in enumerate(idx)}

gt = rows_by_condition('technical_duplicate_first_half_baseline')
td = rows_by_condition('technical_duplicate_second_half_baseline')
idup = rows_by_condition('interpolated_duplicate_baseline')
ctrl_b = rows_by_condition('ctrl_baseline')[DS]
# sanity: author GT equals my first-half means
dmax = max(float(np.abs(gt[f'{DS}_{c}'] - mine_mean[ci[c], 0]).max()) for c in cats if c != 'control')
print('max |author GT - my first-half mean| =', dmax)

codes = f['obs/condition/codes'][:]
recs = []
for k in range(5):
    part = os.path.join(CACHE, f'native_fold{k}.csv')
    if os.path.exists(part):
        recs.extend(pd.read_csv(part).to_dict('records')); print(f'fold {k}: cached'); continue
    if time.time() - t0 > 100:
        print('time budget reached; rerun to continue'); raise SystemExit(0)
    frecs = []
    sp = f[f'obs/split_fold_{k}']; scat = dec(sp['categories'][:]); scode = sp['codes'][:]
    first = {}
    for c_code, s_code in zip(codes, scode):
        first.setdefault(c_code, set()).add(scat[s_code])
    test = [cats[c] for c, roles in first.items() if c != ctrl and roles == {'test'}]
    mean_b = rows_by_condition(f'split_fold_{k}_mean_baseline')[DS]
    lin = rows_by_condition(f'split_fold_{k}_linear_baseline')
    preds = {m: load_prediction(m, k) for m in ('GEARS', 'scGPT')}
    pidx = {}
    for m, (pc, pg, PX) in preds.items():
        pc = [c.split('+')[0] if '+' in c else c for c in pc]
        pc = [c[len(DS) + 1:] if c.startswith(DS + '_') else c for c in pc]
        gpos = {g: j for j, g in enumerate(genes)}
        gsel_m = None if pg == genes else np.array([gpos[g] for g in pg])
        pidx[m] = ({c: i for i, c in enumerate(pc)}, gsel_m, PX)
    for c in test:
        key = f'{DS}_{c}'
        w_all = author_weights_from_scores(Tauth[ci[c]])
        y = gt[key]
        row = {'fold': k, 'condition': c, 'n_first_half': None}
        def wm(pred, gsel=None):
            if gsel is None:
                w = w_all; p = pred; t = y
            else:
                w, p, t = w_all[gsel], pred, y[gsel]
            w = w / w.sum()
            return float(np.sum(w * (p - t) ** 2))
        row['Dataset Mean'] = wm(mean_b); row['Control Mean'] = wm(ctrl_b)
        row['Linear'] = wm(lin[key]) if key in lin else np.nan
        row['Tech. Dup.'] = wm(td[key]); row['Inter. Dup.'] = wm(idup[key])
        for m, (pmap, gsel, PX) in pidx.items():
            if c in pmap:
                if gsel is None:
                    row[m] = wm(PX[pmap[c]])
                    row[m + '|DatasetMean_same_genes'] = row['Dataset Mean']
                else:
                    row[m] = wm(PX[pmap[c]], gsel)
                    row[m + '|DatasetMean_same_genes'] = wm(mean_b[gsel], gsel)
            else:
                row[m] = np.nan; row[m + '|DatasetMean_same_genes'] = np.nan
        frecs.append(row)
    pd.DataFrame(frecs).to_csv(part, index=False); recs.extend(frecs)
    print(f'fold {k}: {len(test)} test conditions; GEARS predicted {sum(c in pidx["GEARS"][0] for c in test)}, scGPT predicted {sum(c in pidx["scGPT"][0] for c in test)}; t={time.time()-t0:.0f}s', flush=True)

df = pd.DataFrame(recs)
os.makedirs(RESULTS, exist_ok=True)
df.to_csv(os.path.join(RESULTS, 'native_wmse_per_condition.csv'), index=False)
pub = {'GEARS': 0.033, 'scGPT': 0.027, 'Dataset Mean': 0.034, 'Control Mean': 0.035, 'Linear': 0.034, 'Tech. Dup.': 0.011, 'Inter. Dup.': 0.013}
out = {}
for m in ['GEARS', 'scGPT', 'Dataset Mean', 'Control Mean', 'Linear', 'Tech. Dup.', 'Inter. Dup.']:
    v = df[m].dropna()
    perfold = df.groupby('fold')[m].mean()
    out[m] = {'published': pub[m], 'reproduced_mean': round(float(v.mean()), 5), 'sem_over_conditions': round(float(v.std(ddof=1) / np.sqrt(len(v))), 5),
              'sem_over_folds': round(float(perfold.std(ddof=1) / np.sqrt(5)), 5), 'n_conditions': int(len(v)),
              'per_fold': [round(float(x), 5) for x in perfold]}
for m in ['GEARS', 'scGPT']:
    d = (df[m + '|DatasetMean_same_genes'] - df[m]).dropna()
    out[m]['paired_gap_vs_DatasetMean'] = {'mean': round(float(d.mean()), 5), 'sem': round(float(d.std(ddof=1) / np.sqrt(len(d))), 5),
                                           'share_conditions_model_better': round(float((d > 0).mean()), 3)}
out['_check_max_abs_author_GT_minus_my_mean'] = dmax
out['_seconds'] = round(time.time() - t0, 1)
json.dump(out, open(os.path.join(RESULTS, 'native_wmse_summary.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
