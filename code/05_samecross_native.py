"""Primary analysis (locked plan plans/ANALYSIS_PLAN_replogle.md, sha256 eb14e953...):
same-sample vs cross-half WMSE advantage at the native cell budget."""
import os, json, time
import numpy as np, h5py, pandas as pd
from common import H5, CACHE, RESULTS, DS, dec, read_df_rows, load_prediction, author_weights_from_scores

t0 = time.time()
ax = json.load(open(os.path.join(CACHE, 'axes.json')))
cats, genes, ctrl = ax['conditions'], ax['genes'], ax['control_index']
ci = {c: i for i, c in enumerate(cats)}; G = len(genes); gpos = {g: j for j, g in enumerate(genes)}
mean = np.load(os.path.join(CACHE, 'mean.npy'), mmap_mode='r'); T = np.load(os.path.join(CACHE, 'T.npy'), mmap_mode='r')
N = np.load(os.path.join(CACHE, 'N.npy'))
f = h5py.File(H5, 'r'); codes = f['obs/condition/codes'][:]

def score_set(f_, b, qa, qb, ya, yb):
    qa = qa / qa.sum(); qb = qb / qb.sum()
    g = lambda q, y: float(np.sum(q * ((b - y) ** 2 - (f_ - y) ** 2)))
    gaa, gbb, gab, gba = g(qa, ya), g(qb, yb), g(qa, yb), g(qb, ya)
    S = (gaa + gbb) / 2; C = (gab + gba) / 2
    I_direct = float(np.sum((f_ - b) * (qa - qb) * (ya - yb)))
    return {'G_nat': gaa, 'S': S, 'C': C, 'I': S - C, 'I_identity_err': (S - C) - I_direct}

recs = []
for k in range(5):
    sp = f[f'obs/split_fold_{k}']; scat = dec(sp['categories'][:]); scode = sp['codes'][:]
    roles = {}
    for c_code, s_code in zip(codes, scode):
        roles.setdefault(c_code, set()).add(scat[s_code])
    test = [cats[c] for c, r in roles.items() if c != ctrl and r == {'test'}]
    _, _, Mb = read_df_rows(f['uns'][f'split_fold_{k}_mean_baseline']); b_all = Mb[0]
    lidx, _, Ml = read_df_rows(f['uns'][f'split_fold_{k}_linear_baseline']); lmap = {x: i for i, x in enumerate(lidx)}
    models = {}
    for m in ('GEARS', 'scGPT'):
        pc, pg, PX = load_prediction(m, k)
        pc = [c.split('+')[0] for c in pc]; pc = [c[len(DS) + 1:] if c.startswith(DS + '_') else c for c in pc]
        sel = np.arange(G) if pg == genes else np.array([gpos[x] for x in pg])
        models[m] = ({c: i for i, c in enumerate(pc)}, sel, PX)
    for c in test:
        i = ci[c]
        qa_all = author_weights_from_scores(T[i, 0]); qb_all = author_weights_from_scores(T[i, 1])
        if not (np.isfinite(qa_all).all() and np.isfinite(qb_all).all() and qa_all.sum() > 0 and qb_all.sum() > 0):
            recs.append({'fold': k, 'condition': c, 'model': 'INVALID_PANEL'}); continue
        ya_all, yb_all = np.asarray(mean[i, 0]), np.asarray(mean[i, 1])
        todo = [(m, models[m][2][models[m][0][c]], models[m][1]) for m in models if c in models[m][0]]
        key = f'{DS}_{c}'
        if key in lmap:
            todo.append(('Linear', Ml[lmap[key]], np.arange(G)))
        for m, fv, sel in todo:
            r = score_set(fv, b_all[sel], qa_all[sel], qb_all[sel], ya_all[sel], yb_all[sel])
            r.update({'fold': k, 'condition': c, 'model': m, 'nA': int(N[i, 0]), 'nB': int(N[i, 1])})
            recs.append(r)
    print(f'fold {k} done t={time.time() - t0:.0f}s', flush=True)

df = pd.DataFrame(recs)
df.to_csv(os.path.join(RESULTS, 'samecross_native_per_condition.csv'), index=False)
rng = np.random.default_rng(20261102)

def boot(v, B=10000):
    v = np.asarray(v); n = len(v); out = np.empty(B)
    for s in range(0, B, 1000):
        idx = rng.integers(0, n, size=(min(1000, B - s), n)); out[s:s + idx.shape[0]] = v[idx].mean(1)
    return [float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))]

summary = {'plan_sha256': 'eb14e9537bdb582bb7ff561cb54349f66b7ac50d4dce8b68a135aa55b9351ae2',
           'invalid_panels': int((df['model'] == 'INVALID_PANEL').sum()),
           'max_abs_identity_error': float(df['I_identity_err'].abs().max())}
d2 = df[df['model'] != 'INVALID_PANEL']
common = set(d2[d2.model == 'GEARS'].condition) & set(d2[d2.model == 'scGPT'].condition)
for label, sub in [('all', d2), ('common_conditions', d2[d2.condition.isin(common)])]:
    for m in ('GEARS', 'scGPT', 'Linear'):
        x = sub[sub.model == m]
        if len(x) == 0: continue
        e = {'n_conditions': int(len(x))}
        for col in ('G_nat', 'S', 'C', 'I'):
            e[col] = {'mean': float(x[col].mean()), 'ci95': boot(x[col].values)}
        e['I_per_fold'] = [float(v) for v in x.groupby('fold')['I'].mean()]
        e['C_per_fold'] = [float(v) for v in x.groupby('fold')['C'].mean()]
        e['share_conditions_I_positive'] = float((x['I'] > 0).mean())
        if label == 'all' and m in ('GEARS', 'scGPT'):
            e['R1_advantage_survives_decoupling'] = bool(e['C']['ci95'][0] > 0)
            e['R2_coupling_detectable'] = bool(e['I']['ci95'][0] > 0 or e['I']['ci95'][1] < 0)
            e['R3_material_at_reporting_precision'] = bool(e['I']['mean'] >= 0.0005)
        summary[f'{label}:{m}'] = e
summary['seconds'] = round(time.time() - t0, 1)
json.dump(summary, open(os.path.join(RESULTS, 'samecross_native_summary.json'), 'w'), indent=1)
print(json.dumps(summary, indent=1))
