"""EXPLORATORY (not in the locked plan): same conditions at every n — only conditions with >=64 cells in
both halves — so the n effect is not confounded with a change in which perturbations are evaluated."""
import os, json
import numpy as np, pandas as pd
from common import RESULTS
df = pd.read_csv(os.path.join(RESULTS, 'budget_curve_per_condition_rep.csv'))
nat = pd.read_csv(os.path.join(RESULTS, 'samecross_native_per_condition.csv'))
cond = df.groupby(['model', 'n', 'condition'])[['S', 'C', 'I']].mean().reset_index()
rng = np.random.default_rng(20261102)
out = []
for m in ('scGPT', 'GEARS', 'Linear'):
    keep = set(cond[(cond.model == m) & (cond.n == 64)].condition)
    for n in (8, 16, 32, 64, 'native'):
        x = nat[(nat.model == m) & nat.condition.isin(keep)] if n == 'native' else cond[(cond.model == m) & (cond.n == n) & cond.condition.isin(keep)]
        bs = [x['I'].values[rng.integers(0, len(x), len(x))].mean() for _ in range(4000)]
        out.append({'model': m, 'n': n, 'k': len(x), 'S': round(x.S.mean(), 5), 'C': round(x.C.mean(), 5), 'I': round(x.I.mean(), 5),
                    'I_ci': [round(float(np.percentile(bs, 2.5)), 5), round(float(np.percentile(bs, 97.5)), 5)],
                    'I/S': round(x.I.mean() / x.S.mean(), 2) if x.S.mean() > 0 else None})
for o in out: print(o)
json.dump(out, open(os.path.join(RESULTS, 'exploratory_common_set_curve.json'), 'w'), indent=1)
