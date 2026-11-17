"""Exploratory follow-up analyses (not pre-specified in the locked plans).

Computed only from the per-perturbation outputs in results/target_gene/ and the
native WMSE tables, so no predictions or raw data are needed. All intervals are
paired percentile bootstraps over test perturbations (10,000 resamples, seed
20261102, a fresh generator per comparison); for ratios the numerator and
denominator are resampled together.

1. Position of b_KD between the negative control (dataset mean) and the positive
   control of ref. 4 (interpolated duplicate): mean gain of b_KD over the dataset
   mean divided by the mean gain of the positive control, all test perturbations.
2. b_KD scored on each model's own perturbation set (matched leaderboard).
3. Per model: paired differences from b_KD with bootstrap p values, Holm-adjusted
   over the 27 paired comparisons (Fig. 2 of both manuscripts; the key 'holm_fig1'
   predates the renumbering);
   win rates; per-fold estimates; the number of
   perturbations whose removal reverses b_KD's lead (influence).
4. Knockdown-supplied models: the model's prediction at the targeted gene is
   replaced by b_KD's value. WMSE_aug = WMSE_f + contrib_t - (WMSE_b - WMSE_bKD).
   For one perturbation, weighted R^2 = 1 - k * WMSE with k = sum(w) / SST, the
   same for every predictor on the same gene axis (the errors y - yhat are the
   same for both deltas); k is taken from an unfloored predictor of that
   perturbation and the result is floored at -1.
5. Target-excluded advantages over the dataset mean for every model, the target
   share E1, and E3 'gap closed' (pre-specified in the target-gene plan). The
   pre-specified values reported in the manuscripts come from scripts 20 and 22;
   the intervals here use their own generators and can differ in the last digit.
6. Positions of each model, of the model given the knockdown, and of b_KD on the
   model's own perturbations (Fig. 1 of both manuscripts).
7. Concentration of the target's contribution to the dataset mean's WMSE in
   highly expressed targets: the target's training-mean expression b_t is
   recovered per perturbation as sqrt(contribution / q_t) / |1 - y_t / b_t|, and
   the share of the summed contribution coming from the 10% of perturbations
   with the highest b_t is reported, with the Spearman correlation between b_t
   and the contribution.

Usage: python code/23_followup_exploratory.py  (from the repository root)
"""
import json
import pathlib

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / 'results'
TG = RES / 'target_gene'
B, SEED = 10_000, 20261102
NATIVE = {'replogle22k562': RES / 'native_wmse_per_condition.csv',
          'nadig25hepg2': RES / 'nadig25hepg2' / 'native_wmse_per_condition.csv'}


def boot_means(cols, ratio=False):
    """Bootstrap the mean (or the ratio of means of cols[0] / cols[1]) over rows."""
    rng = np.random.default_rng(SEED)
    cols = [np.asarray(c, float) for c in cols]
    n = len(cols[0])
    out = np.empty(B)
    for s in range(0, B, 500):
        idx = rng.integers(0, n, (min(500, B - s), n))
        if ratio:
            out[s:s + len(idx)] = cols[0][idx].mean(1) / cols[1][idx].mean(1)
        else:
            out[s:s + len(idx)] = cols[0][idx].mean(1)
    return out


def summarise(x, ratio_den=None):
    x = np.asarray(x, float)
    if ratio_den is None:
        est = x.mean()
        bt = boot_means([x])
    else:
        d = np.asarray(ratio_den, float)
        est = x.mean() / d.mean()
        bt = boot_means([x, d], ratio=True)
    lo, hi = np.percentile(bt, [2.5, 97.5])
    p = min(1.0, 2 * min((bt <= 0).mean(), (bt >= 0).mean()))
    return {'est': float(est), 'ci95': [float(lo), float(hi)], 'p_boot': max(p, 1 / B)}


def holm(pvals):
    order = np.argsort(pvals)
    m = len(pvals)
    adj = np.empty(m)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, (m - rank) * pvals[i])
        adj[i] = min(1.0, run)
    return adj


def load(ds):
    kd = pd.read_csv(TG / f'{ds}_kd_baseline_all_conditions.csv')
    kd['target_in_axis'] = kd['target_in_axis'].astype(str) == 'True'
    nat = pd.read_csv(NATIVE[ds])[['fold', 'condition', 'Inter. Dup.', 'Tech. Dup.']]
    wb = pd.read_csv(TG / f'{ds}_weighted_r2_baselines_per_condition.csv')
    pc = pd.concat([pd.read_csv(TG / f'{ds}_per_condition.csv'),
                    pd.read_csv(TG / f'{ds}_per_condition_retrained.csv')]).drop_duplicates(['fold', 'condition', 'model'])
    wm = pd.concat([pd.read_csv(TG / f'{ds}_weighted_r2_models_per_condition.csv'),
                    pd.read_csv(TG / f'{ds}_weighted_r2_models_per_condition_retrained.csv')]).drop_duplicates(['fold', 'condition', 'model'])
    return kd, nat, wb, pc.merge(wm, on=['fold', 'condition', 'model'])


def r2_from_wmse(d, metric, wmse_new):
    """Exact per-perturbation weighted R^2 for a new prediction, via k = (1 - R2) / WMSE."""
    k = np.full(len(d), np.nan)
    for r2c, wc in [(f'{metric}|f', 'WMSE_f'), (f'{metric}|KD', 'WMSE_bKD'), (f'{metric}|KD_add', 'WMSE_bKD_add')]:
        ok = np.isnan(k) & (d[r2c].values > -1 + 1e-9) & (d[wc].values > 0)
        k[ok] = (1 - d[r2c].values[ok]) / d[wc].values[ok]
    return np.maximum(-1.0, 1 - k * wmse_new), int(np.isnan(k).sum())


def main():
    out, fig1 = {}, []
    for ds in ['replogle22k562', 'nadig25hepg2']:
        kd, nat, wb, mm = load(ds)
        o = {}
        m = kd.merge(nat, on=['fold', 'condition']).merge(wb, on=['fold', 'condition', 'target_in_axis'])
        dm, bkd = m['Dataset Mean'], m['KD multiplicative (primary)']
        o['position_vs_positive_control'] = {
            'WMSE': summarise(dm - bkd, ratio_den=dm - m['Inter. Dup.']),
            'WR2_ctrl': summarise(m['R2ctrl|KD'] - m['R2ctrl|Dataset Mean'], ratio_den=m['R2ctrl|Inter. Dup.'] - m['R2ctrl|Dataset Mean']),
            'WR2_pert': summarise(m['R2pert|KD'] - m['R2pert|Dataset Mean'], ratio_den=m['R2pert|Inter. Dup.'] - m['R2pert|Dataset Mean']),
        }
        sc = m[m.target_in_axis]
        o['bKD_by_set'] = {
            'all': {'n': len(m), 'WMSE': float(bkd.mean()), 'WR2_ctrl': float(m['R2ctrl|KD'].mean()), 'WR2_pert': float(m['R2pert|KD'].mean())},
            'target_scored': {'n': len(sc), 'WMSE': float(sc['KD multiplicative (primary)'].mean()),
                              'WR2_ctrl': float(sc['R2ctrl|KD'].mean()), 'WR2_pert': float(sc['R2pert|KD'].mean())},
            'zero_variant_all': float(m['KD zero (exploratory)'].mean()),
            'zero_variant_target_scored': float(sc['KD zero (exploratory)'].mean()),
            'positive_control_interdup_all': {'WMSE': float(m['Inter. Dup.'].mean()), 'WR2_ctrl': float(m['R2ctrl|Inter. Dup.'].mean()),
                                              'WR2_pert': float(m['R2pert|Inter. Dup.'].mean())},
        }
        lin = mm[mm.model == 'Linear']   # rows = perturbations with a scored target
        sh = lin['share_of_WMSE_b_from_t']
        o['target_share_of_mean_baseline_WMSE'] = {'mean_of_ratios': float(sh.mean()),
                                                   'ratio_of_means': float((sh * lin.WMSE_b).mean() / lin.WMSE_b.mean()),
                                                   'median': float(sh.median())}
        contrib_b = sh * lin.WMSE_b                           # q_t (b_t - y_t)^2 of the dataset mean
        ok = (lin.q_t > 0) & ((1 - lin.true_ratio_t).abs() > 1e-9)
        b_t = np.sqrt(contrib_b[ok] / lin.q_t[ok]) / (1 - lin.true_ratio_t[ok]).abs()
        top = np.argsort(-b_t.values)[:len(b_t) // 10]
        o['target_share_of_mean_baseline_WMSE']['by_expression'] = {
            'n': int(ok.sum()),
            'top_decile_share_of_contribution': float(contrib_b[ok].values[top].sum() / contrib_b[ok].sum()),
            'spearman_bt_contribution': float(pd.Series(b_t.values).corr(pd.Series(contrib_b[ok].values), method='spearman'))}
        models = {}
        for mod, d in mm.groupby('model'):
            r = {'n': len(d)}
            diff = d.WMSE_bKD - d.WMSE_f                          # > 0 favours the model
            r['WMSE_model'] = float(d.WMSE_f.mean()); r['WMSE_bKD_same_set'] = float(d.WMSE_bKD.mean())
            r['vs_bKD'] = {'WMSE': summarise(diff),
                           'WR2_ctrl': summarise(d['R2ctrl|f'] - d['R2ctrl|KD']),
                           'WR2_pert': summarise(d['R2pert|f'] - d['R2pert|KD'])}
            r['win_rate_vs_bKD'] = float((d.WMSE_f < d.WMSE_bKD).mean())
            r['median_WMSE_diff_vs_bKD'] = float(diff.median())
            r['per_fold_WMSE_diff_vs_bKD'] = [round(float(v), 5) for v in d.assign(x=diff).groupby('fold').x.mean()]
            if diff.mean() < 0:                                   # how many perturbations carry b_KD's lead
                s = np.sort(diff.values)
                r['n_removed_to_reverse_bKD_lead'] = int(next(k for k in range(1, len(s)) if s[k:].mean() >= 0))
                r['most_bKD_favourable'] = d.assign(x=diff).sort_values('x').condition.head(10).tolist()
            aug = d.WMSE_f + d.contrib_t - (d.WMSE_b - d.WMSE_bKD)
            r2c, nanc = r2_from_wmse(d, 'R2ctrl', aug); r2p, nanp = r2_from_wmse(d, 'R2pert', aug)
            okc, okp = ~np.isnan(r2c), ~np.isnan(r2p)
            r['knockdown_supplied'] = {
                'WMSE': float(aug.mean()), 'WR2_ctrl': float(np.nanmean(r2c)), 'WR2_pert': float(np.nanmean(r2p)),
                'n_missing_k': [nanc, nanp],
                'vs_bKD': {'WMSE': summarise(d.WMSE_bKD - aug),
                           'WR2_ctrl': summarise((r2c - d['R2ctrl|KD'].values)[okc]),
                           'WR2_pert': summarise((r2p - d['R2pert|KD'].values)[okp])}}
            r['target_excluded_vs_mean'] = {'WMSE': summarise(d.G_X),
                                            'WR2_ctrl': summarise(d['R2ctrlX|f'] - d['R2ctrlX|mean']),
                                            'WR2_pert': summarise(d['R2pertX|f'] - d['R2pertX|mean']),
                                            'win_rate_WMSE': float((d.G_X > 0).mean())}
            r['WMSE_advantage_vs_mean'] = summarise(d.G_nat)
            r['E1_target_share'] = summarise(d.contrib_t, ratio_den=d.G_nat)
            r['E3_gap_closed'] = summarise(d.WMSE_b - d.WMSE_bKD, ratio_den=d.WMSE_b - d.WMSE_f)
            r['median_pred_ratio_t'] = float(d.pred_ratio_t.median())
            # position on the calibrated scale, on the model's own perturbations (Fig. 1)
            x = d.merge(m[['fold', 'condition', 'Dataset Mean', 'Inter. Dup.', 'R2ctrl|Dataset Mean', 'R2ctrl|Inter. Dup.',
                           'R2pert|Dataset Mean', 'R2pert|Inter. Dup.']], on=['fold', 'condition'])
            def pos(wmse, rc, rp):
                out = {'WMSE': float((x.WMSE_b - wmse).mean()   # gain on the model's own gene axis
                                     / (x['Dataset Mean'] - x['Inter. Dup.']).mean())}
                for met, col, r in [('WR2_ctrl', 'R2ctrl', rc), ('WR2_pert', 'R2pert', rp)]:
                    r = np.asarray(r, float)
                    ok = ~np.isnan(r)   # a perturbation without an unfloored score is left out of every term
                    dm, dup = x.loc[ok, f'{col}|Dataset Mean'], x.loc[ok, f'{col}|Inter. Dup.']
                    out[met] = float((r[ok].mean() - dm.mean()) / (dup - dm).mean())
                return out
            aug_x = x.WMSE_f + x.contrib_t - (x.WMSE_b - x.WMSE_bKD)
            rc_aug, _ = r2_from_wmse(x, 'R2ctrl', aug_x); rp_aug, _ = r2_from_wmse(x, 'R2pert', aug_x)
            r['position'] = {'model': pos(x.WMSE_f, x['R2ctrl|f'], x['R2pert|f']),
                             'knockdown_supplied': pos(aug_x, rc_aug, rp_aug),
                             'bKD_same_set': pos(x.WMSE_bKD, x['R2ctrl|KD'], x['R2pert|KD'])}
            models[mod] = r
            for met in ['WMSE', 'WR2_ctrl', 'WR2_pert']:
                fig1.append((ds, mod, met, r['vs_bKD'][met]['p_boot']))
        o['models'] = models
        out[ds] = o
    adj = holm(np.array([p for *_, p in fig1]))
    out['holm_fig1'] = [{'dataset': ds, 'model': mod, 'metric': met, 'p_boot': p, 'p_holm': float(a)}
                        for (ds, mod, met, p), a in zip(fig1, adj)]
    dst = TG / 'followup_exploratory.json'
    dst.write_text(json.dumps(out, indent=1))
    print('wrote', dst)


if __name__ == '__main__':
    main()
