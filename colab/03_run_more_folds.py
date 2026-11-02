# ===== Further folds =====
# Usage: in a new notebook cell set EXTRA_MODELS and EXTRA_FOLDS, then exec this file. Same steps as the main loop;
# runs already exported to Drive are skipped.
import os, json, shutil, pathlib

assert 'EXTRA_MODELS' in globals() and 'EXTRA_FOLDS' in globals(), 'set EXTRA_MODELS and EXTRA_FOLDS first'
assert all(m in PY for m in EXTRA_MODELS), f'unknown model in {EXTRA_MODELS}'
assert '4000000 //' in pathlib.Path('/content/wmse_runner.py').read_text(), 'run 02_fix_sclambda.py first (it switches exports to 4 MB parts)'
log(f'== extra folds: models={EXTRA_MODELS} folds={EXTRA_FOLDS}')
for ds in DATASETS:
    for k in EXTRA_FOLDS:
        for m in EXTRA_MODELS:
            tag = f'{ds}/{m}/fold{k}'
            run = pathlib.Path(f'/content/runs/{tag}'); run.mkdir(parents=True, exist_ok=True)
            dst = OUT / tag; dst.mkdir(parents=True, exist_ok=True)
            if (dst / 'pred_meta.json').exists():
                log(f'{tag}: already done, skipped'); continue
            env = dict(os.environ, PYTHONUNBUFFERED='1', PYTHONPATH='/presage_src' if m == 'presage' else '')
            try:
                for stage in ('train', 'predict', 'export'):
                    _, dt = sh([PY[m], '/content/wmse_runner.py', '--model', m, '--dataset', ds, '--fold', str(k), '--stage', stage],
                               cwd=CWD[m], env=env, logfile=run / 'run.log')
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
log('== extra folds finished')
