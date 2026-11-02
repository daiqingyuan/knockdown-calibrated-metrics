"""Generic version of 06 for any dataset (WMSE_DS); n grid 8/16/32/64 for Replogle, 8/16 for Nadig (locked plans).
Secondary analysis (locked plan): cell-budget curve. For n = 8/16/32/64 cells per half, three
nested random subsamples per condition (seed 20261102); target weights and means recomputed from the
subsampled cells, reference = all other non-control conditions' cells in the same half partition.
Resumable by parts of 128 conditions."""
import os, sys, json, time
import numpy as np, h5py, pandas as pd
from scipy import sparse
from common import H5, CACHE, RESULTS, DS, PRED, dec, read_df_rows, load_prediction, author_weights_from_scores, var_names

budget = float(sys.argv[1]) if len(sys.argv) > 1 else 140.0
t0 = time.time()
NS = (8, 16, 32, 64) if DS == 'replogle22k562' else (8, 16)
R, SEED, PART = 3, 20261102, 128
OUT = os.path.join(CACHE, 'budget_parts'); os.makedirs(OUT, exist_ok=True)
ax = json.load(open(os.path.join(CACHE, 'axes.json')))
cats, genes, ctrl = ax['conditions'], ax['genes'], ax['control_index']
G = len(genes); gpos = {g: j for j, g in enumerate(genes)}
f = h5py.File(H5, 'r')
codes = f['obs/condition/codes'][:]; half = f['obs/tech_dup_split/codes'][:]
indptr = f['X/indptr'][:].astype(np.int64); Xd, Xi = f['X/data'], f['X/indices']
chg = np.flatnonzero(np.diff(codes)) + 1; starts = np.r_[0, chg]; ends = np.r_[chg, len(codes)]; bcodes = codes[starts]
todo_parts = [p0 for p0 in range(0, len(starts), PART) if not os.path.exists(os.path.join(OUT, f'part_{p0:05d}.csv'))]
if not todo_parts:
    print('all parts done'); raise SystemExit(0)
# full-half sufficient statistics -> partition totals over non-control conditions
P = os.path.join(CACHE, 'stats_parts'); S = np.zeros((len(cats), 2, G)); Q = np.zeros_like(S); N = np.zeros((len(cats), 2), np.int64)
for p in sorted(x for x in os.listdir(P) if x.endswith('.npz') and '.tmp' not in x):
    z = np.load(os.path.join(P, p)); S[z['codes']] = z['S']; Q[z['codes']] = z['Q']; N[z['codes']] = z['N']
nc = np.array([i != ctrl for i in range(len(cats))])
St = S[nc].sum(0); Qt = Q[nc].sum(0); Nt = N[nc].sum(0)          # (2, G), (2,)
# fold membership, predictions, baselines (test rows only)
fold_of, preds, base, lin = {}, {}, {}, {}
for k in range(5):
    sp = f[f'obs/split_fold_{k}']; scat = dec(sp['categories'][:]); scode = sp['codes'][:]
    roles = {}
    for c_code, s_code in zip(codes, scode):
        roles.setdefault(c_code, set()).add(scat[s_code])
    for c, r in roles.items():
        if c != ctrl and r == {'test'}:
            fold_of[c] = k
    base[k] = read_df_rows(f['uns'][f'split_fold_{k}_mean_baseline'])[2][0]
    lidx, _, Ml = read_df_rows(f['uns'][f'split_fold_{k}_linear_baseline'])
    for j, key in enumerate(lidx):
        c = key[len(DS) + 1:]
        if c in cats and fold_of.get(cats.index(c)) == k:
            lin[cats.index(c)] = Ml[j].copy()
    del Ml
    for m in sorted({mm for (mm, kk) in PRED if kk == k}):
        pc, pg, PX = load_prediction(m, k)
        pc = [c.split('+')[0] for c in pc]; pc = [c[len(DS) + 1:] if c.startswith(DS + '_') else c for c in pc]
        sel = np.arange(G) if pg == genes else np.array([gpos[x] for x in pg])
        for j, c in enumerate(pc):
            preds[(m, cats.index(c))] = (PX[j], sel)
print(f'setup done t={time.time() - t0:.0f}s', flush=True)

def tstat(xs, c, h):
    n = xs.shape[0]; m = xs.mean(0); v = xs.var(0, ddof=1)
    nr = Nt[h] - N[c, h]; mr = (St[h] - S[c, h]) / nr
    vr = np.maximum(((Qt[h] - Q[c, h]) - nr * mr ** 2) / (nr - 1), 0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        t = (m - mr) / np.sqrt(v / n + vr / n)
    t[np.isnan(t)] = 0.0
    return t, m

for p0 in todo_parts:
    if time.time() - t0 > budget:
        print('time budget reached; rerun to continue'); break
    rows = []
    for b in range(p0, min(p0 + PART, len(starts))):
        c = int(bcodes[b])
        if c == ctrl:
            continue
        a, z = indptr[starts[b]], indptr[ends[b]]
        X = sparse.csr_matrix((Xd[a:z], Xi[a:z], indptr[starts[b]:ends[b] + 1] - a), shape=(ends[b] - starts[b], G)).toarray().astype(np.float64)
        h = half[starts[b]:ends[b]]; XA, XB = X[h == 0], X[h == 1]
        k = fold_of[c]
        models = [(m, *preds[(m, c)]) for m in ('GEARS', 'scGPT') if (m, c) in preds]
        if c in lin:
            models.append(('Linear', lin[c], np.arange(G)))
        for r in range(R):
            rng = np.random.default_rng([SEED, c, r])
            pa, pb = rng.permutation(XA.shape[0]), rng.permutation(XB.shape[0])
            for n in NS:
                if XA.shape[0] < n or XB.shape[0] < n:
                    continue
                ta, ya = tstat(XA[pa[:n]], c, 0); tb, yb = tstat(XB[pb[:n]], c, 1)
                qa_all, qb_all = author_weights_from_scores(ta), author_weights_from_scores(tb)
                if not (np.isfinite(qa_all).all() and np.isfinite(qb_all).all()):
                    rows.append({'condition': cats[c], 'fold': k, 'n': n, 'rep': r, 'model': 'INVALID_PANEL'}); continue
                for m, fv, sel in models:
                    qa = qa_all[sel] / qa_all[sel].sum(); qb = qb_all[sel] / qb_all[sel].sum()
                    bb, A, Bv = base[k][sel], ya[sel], yb[sel]
                    g = lambda q, y: float(np.sum(q * ((bb - y) ** 2 - (fv - y) ** 2)))
                    gaa, gbb, gab, gba = g(qa, A), g(qb, Bv), g(qa, Bv), g(qb, A)
                    rows.append({'condition': cats[c], 'fold': k, 'n': n, 'rep': r, 'model': m, 'G_nat': gaa,
                                 'S': (gaa + gbb) / 2, 'C': (gab + gba) / 2, 'I': (gaa + gbb - gab - gba) / 2})
    fn = os.path.join(OUT, f'part_{p0:05d}.csv'); pd.DataFrame(rows).to_csv(fn + '.tmp', index=False); os.replace(fn + '.tmp', fn)
    print(f'part {p0} written t={time.time() - t0:.0f}s', flush=True)
left = [p0 for p0 in range(0, len(starts), PART) if not os.path.exists(os.path.join(OUT, f'part_{p0:05d}.csv'))]
print('parts remaining:', len(left))
