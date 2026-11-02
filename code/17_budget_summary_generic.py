"""Generic version of 07 (any WMSE_DS). Summarize the cell-budget curve and draw the figure (I and C by cells per half, native budget at right)."""
import os, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import CACHE, RESULTS, DS

P = os.path.join(CACHE, 'budget_parts')
df = pd.concat([pd.read_csv(os.path.join(P, x)) for x in sorted(os.listdir(P)) if x.endswith('.csv') and os.path.getsize(os.path.join(P, x)) > 1], ignore_index=True)  # empty part = no eligible model rows
inv = int((df.model == 'INVALID_PANEL').sum()); df = df[df.model != 'INVALID_PANEL']
df.to_csv(os.path.join(RESULTS, 'budget_curve_per_condition_rep.csv'), index=False)
cond = df.groupby(['model', 'n', 'condition'])[['G_nat', 'S', 'C', 'I']].mean().reset_index()
nat = pd.read_csv(os.path.join(RESULTS, 'samecross_native_per_condition.csv'))
nat = nat[nat.model.isin(['GEARS', 'scGPT', 'Linear'])]
rng = np.random.default_rng(20261102)
def boot(v, B=10000):
    v = np.asarray(v); out = np.empty(B)
    for s in range(0, B, 1000):
        idx = rng.integers(0, len(v), size=(min(1000, B - s), len(v))); out[s:s + idx.shape[0]] = v[idx].mean(1)
    return [float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))]
rows = []
NS = (8, 16, 32, 64) if DS == 'replogle22k562' else (8, 16)
for m in [x for x in ('GEARS', 'scGPT', 'Linear') if x in set(df.model)]:
    for n in (*NS, 'native'):
        x = nat[nat.model == m] if n == 'native' else cond[(cond.model == m) & (cond.n == n)]
        e = {'model': m, 'n_per_half': n, 'n_conditions': int(len(x))}
        for col in ('G_nat', 'S', 'C', 'I'):
            e[col] = float(x[col].mean()); e[col + '_ci'] = boot(x[col].values)
        e['I_over_S'] = e['I'] / e['S'] if e['S'] > 0 else None
        rows.append(e)
summ = pd.DataFrame(rows)
summ.to_json(os.path.join(RESULTS, 'budget_curve_summary.json'), orient='records', indent=1)
pd.set_option('display.width', 200)
show = summ[['model', 'n_per_half', 'n_conditions', 'G_nat', 'S', 'C', 'I', 'I_over_S']].copy()
for c in ('G_nat', 'S', 'C', 'I'):
    show[c] = show[c].map(lambda v: f'{v:.5f}')
show['I_ci'] = summ['I_ci'].map(lambda v: f'[{v[0]:.5f}, {v[1]:.5f}]')
show['C_ci'] = summ['C_ci'].map(lambda v: f'[{v[0]:.5f}, {v[1]:.5f}]')
show['I_over_S'] = summ['I_over_S'].map(lambda v: '' if v is None else f'{v:.2f}')
print(show.to_string(index=False)); print('invalid panels:', inv)

ms = [x for x in ('scGPT', 'GEARS', 'Linear') if x in set(summ.model)]
fig, axes = plt.subplots(1, len(ms), figsize=(4.4 * len(ms), 4), sharex=True, squeeze=False); axes = axes[0]
xs = list(range(len(NS))) + [len(NS) + 0.5]
plt.rcParams['font.family'] = ['DejaVu Sans']
for a, m in zip(axes, ms):
    s = summ[summ.model == m]
    for col, colr, name in (('S', '#9a9a9a', 'same-sample S'), ('C', '#1f6fb4', 'cross-half C'), ('I', '#d1495b', 'coupling I = S - C')):
        y = s[col].values; lo = [v[0] for v in s[col + '_ci']]; hi = [v[1] for v in s[col + '_ci']]
        a.plot(xs, y, marker='o', color=colr, label=name, lw=1.8)
        a.fill_between(xs, lo, hi, color=colr, alpha=0.15, lw=0)
    a.axhline(0, color='black', lw=0.8)
    a.set_title(f'{m} vs Dataset Mean ({DS})', fontsize=10)
    a.set_xticks(xs); a.set_xticklabels([str(n) for n in NS] + ['native\n(all, median ' + ('61' if DS == 'replogle22k562' else '23') + ')'])
    a.axvline(len(NS) - 0.25, color='#cccccc', lw=0.8, ls='--')
    a.set_xlabel('cells per technical half (all eligible conditions at each n)'); a.grid(alpha=0.25)
axes[0].set_ylabel('WMSE advantage over mean baseline')
axes[0].legend(fontsize=8, frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(RESULTS, 'budget_curve.png'), dpi=160)
print('figure written')
