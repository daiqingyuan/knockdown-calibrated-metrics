# ===== scLambda fix =====
# The live Ensembl BioMart (www.ensembl.org) redirected to jun2026.archive.ensembl.org, which answered with a
# 'Service unavailable' page, so scLambda received no gene descriptions, its embeddings were all zero and training stopped.
# Fix: use the scLambda wrapper's own GeneDescriptionFetcher (same BioMart query, same NCBI gene_summary, same cleaning)
# with the BioMart server switched to the Ensembl release 115 archive; build the description cache, test the OpenAI key,
# then re-run scLambda only.
# Finished predictions are also re-exported as parts of at most 4 MB (parts_4mb/), and later exports use 4 MB parts.
import os, json, shutil, pathlib

BIOMART_HOST = 'sep2025.archive.ensembl.org'   # Ensembl release 115: newest archive whose BioMart still answers (checked 2026-11-02)
PART_BYTES = 4_000_000

BUILDER = r'''"""Pre-build scLambda's gene-description cache with the wrapper's own code; only the BioMart server differs."""
import hashlib, json, os, pathlib, pickle, sys
import anndata as ad, h5py, numpy as np, pandas as pd
sys.path.insert(0, os.environ.get('SCLAMBDA_DIR', '/content/cellsimbench/docker/sclambda'))
import sclambda_wrapper as W
from pybiomart import Dataset
try:
    from anndata.io import read_elem
except ImportError:
    from anndata.experimental import read_elem

ds, host = sys.argv[1], sys.argv[2]
data_path = os.environ.get('DATA_ROOT', '/content/data') + f'/{ds}/{ds}_processed_complete.h5ad'
root = pathlib.Path(os.environ.get('CACHE_ROOT', '/content/sclambda_cache')) / ds   # = gene_embedding_config.cache_dir
desc_dir = root / 'descriptions'
desc_dir.mkdir(parents=True, exist_ok=True)

# 1. drop the all-zero embedding pickle left by the failed attempt (a real one is kept)
pkl = root / 'gene_embeddings.pkl'
if pkl.exists():
    E = pickle.load(open(pkl, 'rb'))
    if all(not np.any(v) for g, v in E.items() if g != 'ctrl'):
        pkl.unlink()
        print('removed all-zero gene_embeddings.pkl', flush=True)

# 2. the same gene list SCLambdaWrapper.train() builds from the full dataset
with h5py.File(data_path, 'r') as fh:
    cond = pd.Series(read_elem(fh['obs/condition'])).astype(str).unique()
mini = ad.AnnData(X=np.zeros((len(cond), 1), dtype=np.float32),
                  obs=pd.DataFrame({'condition': cond}, index=[str(i) for i in range(len(cond))]))
w = W.SCLambdaWrapper.__new__(W.SCLambdaWrapper)
genes = w._get_unique_genes(w._convert_to_sclambda_format(mini))

# 3. GeneDescriptionFetcher exactly as the wrapper uses it, but pointed at a BioMart that still answers
cache_file = desc_dir / 'bulk_descriptions_cache.json'
if cache_file.exists():
    old = json.load(open(cache_file, encoding='utf-8'))
    if sum(1 for g in genes if old.get(g)) < 0.9 * len(genes):
        cache_file.unlink()
f = W.GeneDescriptionFetcher.__new__(W.GeneDescriptionFetcher)
f.cache_dir, f.gene_summary_path = desc_dir, desc_dir / 'gene_summary.gz'
f.dataset = Dataset(name='hsapiens_gene_ensembl', host='http://' + host, use_cache=False)
f.dataset._host, f.dataset._port = 'https://' + host, 443      # pybiomart can only build http:// URLs
f._ensure_gene_summary_file()
f._description_cache, f._cache_loaded = {}, False
desc = f.get_gene_descriptions_batch(genes)
sub = {g: desc.get(g, '') for g in genes}
summary = {'dataset': ds, 'biomart': host, 'n_genes': len(genes),
           'n_with_description': sum(1 for v in sub.values() if v),
           'cache_entries': len(json.load(open(cache_file, encoding='utf-8'))) if cache_file.exists() else 0,
           'gene_summary_sha256': hashlib.sha256(f.gene_summary_path.read_bytes()).hexdigest(),
           'descriptions_sha256': hashlib.sha256(json.dumps(sub, sort_keys=True).encode()).hexdigest()}
json.dump(summary, open(desc_dir / 'fix_summary.json', 'w'))
print(json.dumps(summary), flush=True)
sys.exit(0 if summary['cache_entries'] and summary['n_with_description'] >= 0.9 * len(genes) else 1)
'''

OPENAI_TEST = r'''"""One embeddings call with the notebook's key; prints only the outcome, never the key."""
import os, sys
import openai
out = sys.argv[1]
try:
    r = openai.OpenAI(api_key=os.environ['OPENAI_API_KEY']).embeddings.create(
        input=['DNA repair protein'], model='text-embedding-3-large', dimensions=3072)
    msg, rc = f'OpenAI embeddings OK (dim {len(r.data[0].embedding)})', 0
except Exception as e:
    msg, rc = (f'OpenAI embeddings FAILED: {type(e).__name__} status={getattr(e, "status_code", None)} '
               f'code={getattr(e, "code", None)}'), 1
open(out, 'w').write(msg)
print(msg)
sys.exit(rc)
'''

REEXPORT = r'''"""Re-export a finished predictions.h5ad as float32 parts of at most max_bytes each (same layout as wmse_runner export)."""
import hashlib, json, pathlib, sys
import anndata as ad, numpy as np
run, out, max_bytes = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), int(sys.argv[3])
A = ad.read_h5ad(run / 'predictions.h5ad')
X = A.X.toarray() if hasattr(A.X, 'toarray') else np.asarray(A.X)
X = X.astype(np.float32)
out.mkdir(parents=True, exist_ok=True)
meta = {'shape': list(X.shape), 'obs_names': [str(x) for x in A.obs_names], 'var_names': [str(x) for x in A.var_names],
        'obs': {c: [str(v) for v in A.obs[c]] for c in A.obs.columns}, 'parts': [],
        'predictions_h5ad_sha256': hashlib.sha256((run / 'predictions.h5ad').read_bytes()).hexdigest()}
rows = max(1, max_bytes // (X.shape[1] * 4))
for i in range(0, X.shape[0], rows):
    fn = out / f'pred_part{i // rows:02d}.npy'
    np.save(fn, X[i:i + rows])
    meta['parts'].append({'file': fn.name, 'rows': int(min(rows, X.shape[0] - i)),
                          'sha256': hashlib.sha256(fn.read_bytes()).hexdigest()})
json.dump(meta, open(out / 'pred_meta.json', 'w'))
print('re-exported', X.shape, 'as', len(meta['parts']), 'parts', flush=True)
'''

FIX = pathlib.Path('/content/sclambda_fix'); FIX.mkdir(exist_ok=True)   # pybiomart puts its request cache in the working dir
pathlib.Path('/content/cellsimbench/docker/sclambda/.pybiomart.sqlite').unlink(missing_ok=True)   # may hold the cached error page
(FIX / 'build_sclambda_descriptions.py').write_text(BUILDER)
(FIX / 'test_openai_embeddings.py').write_text(OPENAI_TEST)
fix_env = dict(os.environ, PYTHONUNBUFFERED='1', PYTHONPATH='')
(FIX / 'reexport_small_parts.py').write_text(REEXPORT)

# smaller transfer parts: patch the runner for later exports, and re-export models that already finished
runner = pathlib.Path('/content/wmse_runner.py'); rs = runner.read_text()
if '9_000_000 //' in rs:
    runner.write_text(rs.replace('9_000_000 //', f'{PART_BYTES} //')); log(f'runner export parts now <= {PART_BYTES} bytes')
for ds in DATASETS:
    for k in FOLDS:
        for mm in ('fmlp_genept', 'presage'):
            tag = f'{ds}/{mm}/fold{k}'
            src, small = pathlib.Path(f'/content/runs/{tag}'), OUT / tag / 'parts_4mb'
            if (OUT / tag / 'pred_meta.json').exists() and (src / 'predictions.h5ad').exists() and not (small / 'pred_meta.json').exists():
                rc, _ = sh([PY_MAIN, str(FIX / 'reexport_small_parts.py'), str(src), str(small), str(PART_BYTES)],
                           env=fix_env, logfile=PLOG, check=False)
                log(f'{tag}: ' + ('re-exported to parts_4mb' if rc == 0 else 're-export FAILED, see ' + PLOG.name))

if not os.environ.get('OPENAI_API_KEY', '').startswith('sk-'):
    raise RuntimeError('OPENAI_API_KEY not set: fill it in the settings cell, then re-run the model-assets cell')

for ds in DATASETS:
    rc, _ = sh([PY_MAIN, str(FIX / 'build_sclambda_descriptions.py'), ds, BIOMART_HOST], cwd=str(FIX), env=fix_env,
               logfile=PLOG, check=False)
    summ = pathlib.Path(f'/content/sclambda_cache/{ds}/descriptions/fix_summary.json')
    log('scLambda fix: ' + (summ.read_text() if summ.exists() else 'no summary written'))
    if rc != 0:
        raise RuntimeError(f'{ds}: gene descriptions were not built; see the end of {PLOG.name}')

rc, _ = sh([PY_MAIN, str(FIX / 'test_openai_embeddings.py'), str(FIX / 'openai_test.txt')], cwd=str(FIX), env=fix_env,
           logfile=PLOG, check=False)
log('scLambda fix: ' + (FIX / 'openai_test.txt').read_text())
if rc != 0:
    raise RuntimeError('OpenAI embeddings call failed (401: invalid key; 403: key without model access; '
                       '429 insufficient_quota: no credit). Fix the key, re-run the model-assets cell, then this one.')

# re-run only scLambda, exactly like the main loop (anything already exported to Drive is skipped)
m = 'sclambda'
for ds in DATASETS:
    for k in FOLDS:
        tag = f'{ds}/{m}/fold{k}'
        run = pathlib.Path(f'/content/runs/{tag}'); run.mkdir(parents=True, exist_ok=True)
        dst = OUT / tag; dst.mkdir(parents=True, exist_ok=True)
        if (dst / 'pred_meta.json').exists():
            log(f'{tag}: already done, skipped'); continue
        if (run / 'run.log').exists():   # keep each failed attempt's log under its own name
            n = len(list(run.glob('run_failed_attempt*.log'))) + 1
            os.replace(run / 'run.log', run / f'run_failed_attempt{n}.log')
            shutil.copy(run / f'run_failed_attempt{n}.log', dst / f'run_failed_attempt{n}.log')
        try:
            for stage in ('train', 'predict', 'export'):
                _, dt = sh([PY[m], '/content/wmse_runner.py', '--model', m, '--dataset', ds, '--fold', str(k), '--stage', stage],
                           cwd=CWD[m], env=fix_env, logfile=run / 'run.log')
                log(f'{tag}: {stage} done in {dt / 60:.1f} min')
            for f in ['predictions.h5ad', 'pred_meta.json', 'conditions.json', 'train_config.json', 'predict_config.json']:
                shutil.copy(run / f, dst / f)
            for f in sorted(run.glob('pred_part*.npy')):
                shutil.copy(f, dst / f.name)
        except Exception as e:
            log(f'{tag}: FAILED -- {e}')
        finally:
            if (run / 'run.log').exists():
                shutil.copy(run / 'run.log', dst / 'run.log')
log('== scLambda fix run finished')
