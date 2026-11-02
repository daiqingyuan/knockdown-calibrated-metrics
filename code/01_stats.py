"""Pass 1: per-condition, per-technical-half sufficient statistics of log-normalized X.
Cells are stored sorted by condition, so each condition is one contiguous row block.
Resumable: writes one part file per 128 conditions; rerun continues where it stopped."""
import os, sys, time
import numpy as np, h5py
from scipy import sparse
from common import H5, CACHE, dec, var_names

OUT = os.path.join(CACHE, 'stats_parts'); os.makedirs(OUT, exist_ok=True)
PART = 128
budget = float(sys.argv[1]) if len(sys.argv) > 1 else 150.0
t0 = time.time()
f = h5py.File(H5, 'r')
codes = f['obs/condition/codes'][:]
half = f['obs/tech_dup_split/codes'][:]
assert dec(f['obs/tech_dup_split/categories'][:]) == ['first_half', 'second_half']
indptr = f['X/indptr'][:].astype(np.int64)
G = len(var_names(f))
chg = np.flatnonzero(np.diff(codes)) + 1
starts = np.r_[0, chg]; ends = np.r_[chg, len(codes)]
bcodes = codes[starts]
assert len(np.unique(bcodes)) == len(bcodes), 'conditions are not contiguous'
Xd, Xi = f['X/data'], f['X/indices']
for p0 in range(0, len(starts), PART):
    fn = os.path.join(OUT, f'part_{p0:05d}.npz')
    if os.path.exists(fn):
        continue
    if time.time() - t0 > budget:
        print('time budget reached; rerun to continue'); break
    p1 = min(p0 + PART, len(starts))
    S = np.zeros((p1 - p0, 2, G)); Q = np.zeros((p1 - p0, 2, G)); N = np.zeros((p1 - p0, 2), np.int64)
    for j, b in enumerate(range(p0, p1)):
        for c0 in range(starts[b], ends[b], 1024):
            c1 = min(c0 + 1024, ends[b])
            a, z = indptr[c0], indptr[c1]
            X = sparse.csr_matrix((Xd[a:z], Xi[a:z], indptr[c0:c1 + 1] - a), shape=(c1 - c0, G)).toarray().astype(np.float64)
            h = half[c0:c1]
            for hh in (0, 1):
                sel = X[h == hh]
                S[j, hh] += sel.sum(0); Q[j, hh] += np.einsum('ij,ij->j', sel, sel); N[j, hh] += sel.shape[0]
    tmp = fn[:-4] + '.tmp.npz'
    np.savez(tmp, codes=bcodes[p0:p1], S=S, Q=Q, N=N)
    os.replace(tmp, fn)
    print(f'part {p0}-{p1} written at {time.time() - t0:.1f}s', flush=True)
print('parts present:', len([x for x in os.listdir(OUT) if x.startswith('part_') and x.endswith('.npz') and '.tmp' not in x]), 'of', -(-len(starts) // PART))
