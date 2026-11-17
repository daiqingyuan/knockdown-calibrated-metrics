"""Shared configuration, loaders and metric helpers for all analysis scripts.

Configuration comes from environment variables (see README.md):
  WMSE_DATA_DIR    directory holding <dataset>/<dataset>_processed_complete.h5ad (CellSimBench processed data)
  WMSE_WORK_DIR    output root: caches go to cache/ (Replogle) or cache_<dataset>/, results to results/ or
                   results/<dataset>/, and the target-gene analysis writes to results/target_gene/
  WMSE_CACHE_DIR   optional: use an existing cache directory instead of the one under WMSE_WORK_DIR
  WMSE_PRED_INDEX  JSON index of model prediction files (see predictions_index.example.json)
  WMSE_DS          replogle22k562 (default) or nadig25hepg2
  WMSE_PRED_SET    scgpt_gears (default: scGPT and GEARS) or retrained (PRESAGE, scLambda, fMLP-GenePT)
"""
import hashlib
import json
import os
import numpy as np
import h5py

DATA_DIR = os.environ.get('WMSE_DATA_DIR', 'data')
BASE = os.environ.get('WMSE_WORK_DIR', '.')
DS = os.environ.get('WMSE_DS', 'replogle22k562')
assert DS in ('replogle22k562', 'nadig25hepg2'), DS
CACHE = os.environ.get('WMSE_CACHE_DIR') or os.path.join(BASE, 'cache' if DS == 'replogle22k562' else f'cache_{DS}')
RESULTS = os.path.join(BASE, 'results' if DS == 'replogle22k562' else f'results/{DS}')
H5 = os.path.join(DATA_DIR, f'{DS}/{DS}_processed_complete.h5ad')
os.makedirs(CACHE, exist_ok=True); os.makedirs(RESULTS, exist_ok=True)

# Prediction sets. 'scgpt_gears': scGPT and GEARS (outputs without suffix). 'retrained': PRESAGE, scLambda and
# fMLP-GenePT retrained on Colab with the authors' code (outputs with suffix _retrained).
PRED_SET = os.environ.get('WMSE_PRED_SET', 'scgpt_gears')
SET_MODELS = {'scgpt_gears': ['GEARS', 'scGPT'], 'retrained': ['fMLP-GenePT', 'PRESAGE', 'scLambda']}
assert PRED_SET in SET_MODELS, PRED_SET

# PRED maps (model, fold) to a prediction file (.npz from code/30_assemble_colab_predictions.py, or the .h5ad a
# CellSimBench model writes). Index entries are a path, or {"path": ..., "sha256": ...}; a listed SHA-256 is checked
# when the file is loaded. Relative paths are taken relative to the index file. A listed file that is absent stops
# the run, unless WMSE_ALLOW_MISSING=1, which skips it with a warning.
PRED, PRED_SHA256 = {}, {}
_index_file = os.environ.get('WMSE_PRED_INDEX', 'predictions_index.json')
if os.path.exists(_index_file):
    _root = os.path.dirname(os.path.abspath(_index_file))
    _index = json.load(open(_index_file)).get(DS, {})
    _missing = []
    for _m in SET_MODELS[PRED_SET]:
        for _k, _e in sorted(_index.get(_m, {}).items(), key=lambda kv: int(kv[0])):
            _e = {'path': _e} if isinstance(_e, str) else _e
            _p = _e['path'] if os.path.isabs(_e['path']) else os.path.join(_root, _e['path'])
            if not os.path.exists(_p):
                _missing.append(_p)
                continue
            PRED[(_m, int(_k))] = _p
            if _e.get('sha256'):
                PRED_SHA256[(_m, int(_k))] = _e['sha256']
    if _missing:
        if os.environ.get('WMSE_ALLOW_MISSING') != '1':
            raise FileNotFoundError('prediction files listed in the index are missing (set WMSE_ALLOW_MISSING=1 to '
                                    'skip them): ' + ', '.join(_missing))
        print('warning: skipping missing prediction files:', ', '.join(_missing))


def dec(a):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in a]


def h5_col(g, name):
    """Read an AnnData dataframe column that may be stored plainly or as a categorical group."""
    o = g[name]
    if isinstance(o, h5py.Group):
        cats = dec(o['categories'][:])
        return [cats[i] if i >= 0 else None for i in o['codes'][:]]
    return dec(o[:])


def var_names(f):
    return h5_col(f['var'], f['var'].attrs['_index'])


def read_df_rows(g):
    """Read an AnnData-encoded dataframe stored column-wise (columns = genes)."""
    idx = dec(g[g.attrs['_index']][:])
    cols = [c for c in dec(g.attrs['column-order'])]
    M = np.empty((len(idx), len(cols)), dtype=np.float64)
    for j, c in enumerate(cols):
        M[:, j] = g[c][:]
    return idx, cols, M


def load_prediction(model, fold):
    """Return (condition names without dataset prefix, gene names, matrix float64)."""
    p = PRED[(model, fold)]
    if (model, fold) in PRED_SHA256:
        h = hashlib.sha256()
        with open(p, 'rb') as fh:
            for block in iter(lambda: fh.read(1 << 20), b''):
                h.update(block)
        if h.hexdigest() != PRED_SHA256[(model, fold)]:
            raise ValueError(f'{p}: sha256 {h.hexdigest()} does not match the index')
    if p.endswith('.npz'):
        z = np.load(p, allow_pickle=False)
        return [str(c) for c in z['conditions']], [str(g) for g in z['gene_names']], z['predictions'].astype(np.float64)
    with h5py.File(p, 'r') as f:
        X = f['X'][:].astype(np.float64)
        genes = h5_col(f['var'], f['var'].attrs['_index'])
        oc = f['obs/condition']
        if isinstance(oc, h5py.Group):
            cond = [dec(oc['categories'][:])[i] for i in oc['codes'][:]]
        else:
            cond = dec(oc[:])
    return cond, genes, X


def author_weights_from_scores(scores):
    """Author DataManager._precompute_deg_weights: |t| -> min-max -> nan->0 -> square."""
    a = np.abs(scores)
    w = (a - np.min(a)) / (np.max(a) - np.min(a))
    w = np.nan_to_num(w, nan=0.0)
    return np.square(w)
