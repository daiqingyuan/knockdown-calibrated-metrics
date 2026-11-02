"""Assemble Colab retrained-model predictions (float32 parts + pred_meta.json, each part SHA-256 checked) into the
predictions.npz layout load_prediction() reads. Usage: python 30_assemble_colab_predictions.py <dir> [<dir> ...]"""
import hashlib, json, os, sys
import numpy as np

for d in sys.argv[1:]:
    meta = json.load(open(os.path.join(d, 'pred_meta.json')))
    parts = []
    for p in meta['parts']:
        fn = os.path.join(d, p['file'])
        h = hashlib.sha256(open(fn, 'rb').read()).hexdigest()
        assert h == p['sha256'], f'{fn}: sha256 {h} != {p["sha256"]}'
        a = np.load(fn, allow_pickle=False)
        assert a.shape[0] == p['rows'] and a.dtype == np.float32, (fn, a.shape, a.dtype)
        parts.append(a)
    X = np.concatenate(parts)
    assert list(X.shape) == meta['shape'] and len(meta['obs']['condition']) == X.shape[0] and len(meta['var_names']) == X.shape[1]
    assert np.isfinite(X).all(), f'{d}: non-finite predictions'
    np.savez(os.path.join(d, 'predictions.npz'), conditions=np.array(meta['obs']['condition']),
             gene_names=np.array(meta['var_names']), predictions=X)
    print(f'{d}: {X.shape[0]} conditions x {X.shape[1]} genes, {len(parts)} parts verified; '
          f'source predictions.h5ad sha256 {meta["predictions_h5ad_sha256"][:16]}')
