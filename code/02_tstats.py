"""Assemble pass-1 statistics; compute author-style DEG t statistics for BOTH technical halves
(scanpy rank_genes_groups, method='t-test_overestim_var', reference='rest', within each half,
non-control conditions with >=4 cells only), and validate the first-half result against the
author's stored uns['scores_df_dict_gt']."""
import os, time, json
import numpy as np, h5py
from common import H5, CACHE, DS, dec, var_names

t0 = time.time()
P = os.path.join(CACHE, 'stats_parts')
parts = sorted(x for x in os.listdir(P) if x.startswith('part_') and x.endswith('.npz') and '.tmp' not in x)
f = h5py.File(H5, 'r')
cats = dec(f['obs/condition/categories'][:]); genes = var_names(f)
C, G = len(cats), len(genes)
S = np.zeros((C, 2, G)); Q = np.zeros((C, 2, G)); N = np.zeros((C, 2), np.int64)
for p in parts:
    z = np.load(os.path.join(P, p))
    S[z['codes']] = z['S']; Q[z['codes']] = z['Q']; N[z['codes']] = z['N']
assert N.sum() == f['obs/condition/codes'].shape[0]
ctrl = cats.index('control')
valid = np.array([(i != ctrl) for i in range(C)])          # all non-control conditions
mean = S / np.maximum(N[..., None], 1)
var = (Q - N[..., None] * mean ** 2) / np.maximum(N[..., None] - 1, 1)
var = np.maximum(var, 0.0)

T = np.zeros((C, 2, G))
for h in (0, 1):
    ok = valid & (N[:, h] >= 4)
    St, Qt, Nt = S[ok, h].sum(0), Q[ok, h].sum(0), N[ok, h].sum()
    for c in np.flatnonzero(ok):
        nr = Nt - N[c, h]
        mr = (St - S[c, h]) / nr
        vr = np.maximum(((Qt - Q[c, h]) - nr * mr ** 2) / (nr - 1), 0.0)
        ng = N[c, h]
        se = np.sqrt(var[c, h] / ng + vr / ng)          # overestim_var: rest uses n_group
        with np.errstate(divide='ignore', invalid='ignore'):
            t = (mean[c, h] - mr) / se
        t[np.isnan(t)] = 0.0
        T[c, h] = t
np.save(os.path.join(CACHE, 'mean.npy'), mean); np.save(os.path.join(CACHE, 'var.npy'), var)
np.save(os.path.join(CACHE, 'N.npy'), N); np.save(os.path.join(CACHE, 'T.npy'), T)
with open(os.path.join(CACHE, 'axes.json'), 'w') as fh:
    json.dump({'conditions': cats, 'genes': genes, 'control_index': ctrl}, fh)

# ---- validation against the author's stored first-half scores ----
gi = {g: j for j, g in enumerate(genes)}
sg, ng_ = f['uns/scores_df_dict_gt'], f['uns/names_df_dict_gt']
maxabs, maxrel, wmax = 0.0, 0.0, 0.0
worst = None
Aauth = np.zeros((C, G))
for c in np.flatnonzero(valid):
    key = f'{DS}_{cats[c]}'
    sc = sg[key][:]; nm = dec(ng_[key][:])
    a = np.zeros(G); a[[gi[x] for x in nm]] = sc
    Aauth[c] = a
    d = np.abs(a - T[c, 0])
    if d.max() > maxabs:
        maxabs = float(d.max()); worst = (cats[c], genes[int(d.argmax())], float(a[d.argmax()]), float(T[c, 0][d.argmax()]))
    maxrel = max(maxrel, float((d / np.maximum(np.abs(a), 1.0)).max()))
    wa = np.abs(a); wa = ((wa - wa.min()) / (wa.max() - wa.min())) ** 2
    wm = np.abs(T[c, 0]); wm = ((wm - wm.min()) / (wm.max() - wm.min())) ** 2
    wmax = max(wmax, float(np.abs(wa / wa.sum() - wm / wm.sum()).max()))
np.save(os.path.join(CACHE, 'T_author_first_half.npy'), Aauth)
rep = {'n_conditions_checked': int(valid.sum()), 'max_abs_t_diff': maxabs, 'max_rel_t_diff(|t|>=1 scale)': maxrel,
       'worst_case(condition, gene, author_t, mine_t)': worst, 'max_abs_normalized_weight_diff': wmax,
       'half_cell_counts_noncontrol': {'first_min': int(N[valid, 0].min()), 'first_median': float(np.median(N[valid, 0])),
                                       'first_max': int(N[valid, 0].max()), 'second_min': int(N[valid, 1].min()),
                                       'second_max': int(N[valid, 1].max())},
       'seconds': round(time.time() - t0, 1)}
print(json.dumps(rep, indent=1, ensure_ascii=False))
with open(os.path.join(CACHE, 'tstat_validation.json'), 'w') as fh:
    json.dump(rep, fh, indent=1)
