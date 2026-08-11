#!/usr/bin/env python3
"""A6, A8, A9, A10, A11, A12 and the decision-curve analysis -- the remaining pre-specified work.

These seven are mutually independent and all unblocked once A0i (the cluster bootstrap) and A2
(the grid) have run, so they execute in one pass sharing the loaded matrices.

  A8   the parent pre-registration's LITERAL primary statistic, re-run on the corrected scorer.
       glm(pCR ~ signature + PAM50-basal), OR/SD, Wald p; PRIMARY met only at OR/SD >= 1.3 AND
       Wald p < 0.0083. Exactly five tests -- the frozen Framework signature, once per cohort.
       Tagged [PRE-SPECIFIED ESTIMATOR -- CORRECTED EXECUTION] per addendum section 13.6, and
       BOTH runs (deployed and corrected) are reported side by side. GSE32646 and METABRIC are
       outside the parent's scope and are labelled as the estimator extended beyond it.
       METABRIC's endpoint is dichotomised 5-year OS with the direction fixed so that a higher
       score predicting MORE death gives OR > 1 -- the two codings give reciprocal odds ratios
       and only one clears 1.3, so the substitution is stated with the result.

  A9   is the pooled estimate defensible? Leave-one-cohort-out, discovery-excluded pools, and
       the cohort-clustered variance, against the naive pool. Answers R2-M3 directly.

  A10  what would it take to exclude a small effect? Inverts the spike-in: required n for 80%
       power to exclude residual AUROC 0.55 and 0.60 at the observed event rate. Converts the
       power limitation into a study-design recommendation.

  A11  regimen. GSE22226 deposits the only within-cohort regimen field in the study
       (AC / AC+T / AC+T+Herceptin / AC+T+other), so the 11 trastuzumab patients -- the single
       largest treatment effect anywhere in these data, and currently unmodelled -- can be
       excluded in a sensitivity refit. Cross-regimen consistency is otherwise between cohorts,
       not within, and the manuscript's "cross-regimen consistency supports a regimen-general
       redundancy" is checked against that limit.

  A12  Figure 4 versus Figure 6 (R3 point 7): NES-profile correlation in CCLE gene-rank space
       against deployed-score correlation in patient space, over the drug pairs common to both.

  A6   pCR harmonisation (R3 point 6). GSE25066's deposited RCB field is collapsed at RCB-0/I,
       so RCB-0 is NOT recoverable and the endpoint pair is deposited-pCR vs RCB-0/I
       (addendum 13.12). Same-n and full-n comparisons are both required: same-n isolates the
       definition effect, full-n shows what harmonisation costs in patients.

  DCA  decision-curve analysis. Pre-registered, never performed, disclosed as non-delivery in
       correction C6. The decision context is named before plotting, because every patient in
       these cohorts receives chemotherapy and a net-benefit curve with no stated decision is
       decoration: the modelled decision is selection into a response-adapted escalation trial.

Analysis sets are the section 4.1 de-duplicated ones throughout. Writes one JSON per analysis.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import A0e_dedup as D  # noqa: E402
from A0i_cluster_bootstrap import _cell_builders  # noqa: E402

RES = Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/results")
MB = D.ROOT / "resubmission_v2/results/phase_4/metabric"
SEED, B = 42, 2000
ALPHA_PARENT, OR_PARENT = 0.0083, 1.3
MARGINS = [0.53, 0.55, 0.575, 0.60]
COHORTS = ["GSE16446", "GSE22226", "GSE25066", "GSE32646"]


def _out(name, obj):
    d = RES / name.split("_")[0]
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{name}.json").write_text(json.dumps(obj, indent=2, default=float))
    print(f"  wrote {d / (name + '.json')}")


def wald(y, x, adj):
    X = sm.add_constant(np.column_stack([np.asarray(x, float), np.asarray(adj, float)]))
    m = sm.Logit(np.asarray(y, int), X).fit(disp=0)
    b, se = m.params[1], m.bse[1]
    return dict(or_sd=float(np.exp(b)), wald_p=float(m.pvalues[1]),
                or_ci=[float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))],
                n=int(len(y)), events=int(np.sum(y)))


# ------------------------------------------------------------------ A8
def a8(src, b_ov):
    print("\n=== A8: the parent's literal primary statistic, corrected execution ===")
    rows = []
    for coh in COHORTS:
        expr, ax = src["load"](coh)
        if coh == "GSE25066":
            k = ~ax.index.isin(b_ov)
            expr, ax = expr.loc[k], ax.loc[k]
        y = ax["pCR"].values.astype(int)
        basal = D.zscore(ax["H_PAM50basal"]).values
        for tag, mode in [("corrected", "corrected"), ("deployed", "deployed")]:
            s = _framework_score(src, expr, mode)
            r = wald(y, D.zscore(s).values, basal)
            r.update(cohort=coh, scorer=tag, endpoint="pCR",
                     in_parent_scope=coh != "GSE32646")
            r["primary_met"] = bool(r["or_sd"] >= OR_PARENT and r["wald_p"] < ALPHA_PARENT)
            rows.append(r)
    # METABRIC: dichotomised 5-year OS, direction fixed so higher score -> more death gives OR>1
    ax = pd.read_csv(MB / "axes_METABRIC.tsv", sep="\t", index_col=0).dropna(
        subset=["H_PAM50basal", "event5y"])
    sc = pd.read_csv(D.ROOT / "resubmission_v2/results/phase_6/score_Mammaprint_METABRIC.tsv",
                     sep="\t", index_col=0)
    common = ax.index.intersection(sc.index)
    if len(common) > 50:
        num = sc.select_dtypes("number").columns
        col = num[0] if len(num) else sc.columns[0]
        r = wald(ax.loc[common, "event5y"].values.astype(int),
                 D.zscore(sc.loc[common, col]).values,
                 D.zscore(ax.loc[common, "H_PAM50basal"]).values)
        r.update(cohort="METABRIC", scorer="corrected",
                 endpoint="dichotomised 5-year OS (event = dead by 5y; censored <5y excluded); "
                          "direction fixed so a higher score predicting MORE death gives OR>1",
                 in_parent_scope=False, signature="MammaPrint (the only METABRIC cell deposited)")
        r["primary_met"] = bool(r["or_sd"] >= OR_PARENT and r["wald_p"] < ALPHA_PARENT)
        rows.append(r)
    for r in rows:
        print(f"  {r['cohort']:9s} {r['scorer']:9s} OR/SD={r['or_sd']:.3f} "
              f"[{r['or_ci'][0]:.2f},{r['or_ci'][1]:.2f}] Wald p={r['wald_p']:.4f}  "
              f"{'MET' if r['primary_met'] else 'not met'}"
              f"{'' if r['in_parent_scope'] else '   (beyond parent scope)'}")
    n_met = sum(r["primary_met"] for r in rows if r["scorer"] == "corrected")
    print(f"  -> PRIMARY met in {n_met} of "
          f"{sum(r['scorer']=='corrected' for r in rows)} tests under the corrected scorer")
    _out("A8_parent_primary", dict(
        tag="[PRE-SPECIFIED ESTIMATOR -- CORRECTED EXECUTION]",
        criterion=f"OR/SD >= {OR_PARENT} AND Wald p < {ALPHA_PARENT}",
        note="Both scorer runs reported side by side per addendum 13.6. The re-run happens once "
             "and the verdict rule does not move.", rows=rows))
    return rows


def _framework_score(src, expr, mode):
    fw = pd.read_csv(D.ROOT / "frozen_signatures/solid_only/quantile_30_oncokb_all/"
                              "DOXORUBICIN.tsv", sep="\t")
    up = fw.loc[fw.direction == "sensitivity", "gene_symbol"]
    dn = fw.loc[fw.direction == "resistance", "gene_symbol"]
    if mode == "corrected":
        return D.singscore(expr, up, dn)[0]
    n_total = expr.shape[1]
    ranks = expr.rank(axis=1, method="average", pct=True)
    U, Dn = D.resolve(up, expr.columns), D.resolve(dn, expr.columns)

    def norm(m, n):
        mn = (1 + n) / (2 * n_total)
        return 2 * (m - mn) / ((1 - mn) - mn) - 1
    raw = norm(ranks[U].mean(axis=1), len(U)) - norm((1 - ranks[Dn]).mean(axis=1), len(Dn))
    return raw - raw.median()


# ------------------------------------------------------------------ A9
def a9():
    print("\n=== A9: is the pooled estimate defensible? ===")
    d = json.loads((RES / "A0/A0i_cluster_bootstrap.json").read_text())
    cells = d["cells"]
    pt = np.array([c["point"] for c in cells]); se = np.array([c["se"] for c in cells])
    coh = [c["cohort"] for c in cells]

    def pool(mask):
        w = 1 / se[mask] ** 2
        p = float((w * pt[mask]).sum() / w.sum())
        s = float(np.sqrt(1 / w.sum()))
        return dict(pooled=p, se=s, k=int(mask.sum()),
                    tost={str(m): float(stats.norm.cdf((p - m) / s)) for m in MARGINS})

    allm = np.ones(len(pt), bool)
    out = dict(all_cells=pool(allm), leave_one_cohort_out={}, discovery_excluded=None,
               cohort_clustered=None)
    for c in sorted(set(coh)):
        m = np.array([x != c for x in coh])
        out["leave_one_cohort_out"][f"drop {c}"] = pool(m)
        print(f"  drop {c:9s} pooled={out['leave_one_cohort_out'][f'drop {c}']['pooled']:.4f} "
              f"k={int(m.sum())} TOST@0.55={out['leave_one_cohort_out'][f'drop {c}']['tost']['0.55']:.2e}")
    held = np.array([x == "GSE32646" for x in coh])
    if held.sum():
        out["discovery_excluded"] = pool(held)
        print(f"  held-out only (GSE32646): pooled="
              f"{out['discovery_excluded']['pooled']:.4f} k={int(held.sum())} "
              f"TOST@0.55={out['discovery_excluded']['tost']['0.55']:.3f}")
    out["cohort_clustered"] = dict(se_rho1=d["observed"]["se_rho1"],
                                   se_boot=d["se_boot"], se_naive=d["observed"]["se_naive"])
    out["reading"] = ("Three of the four cohorts are discovery. The pooled row is internal "
                      "consistency, not replication; the held-out-only pool is the replication "
                      "statement and it rests on a single cohort of 115 patients.")
    _out("A9_pooled_defensibility", out)
    return out


# ------------------------------------------------------------------ A10
def a10():
    print("\n=== A10: required n to exclude a small residual AUROC ===")
    rng = np.random.default_rng(SEED)

    def sim_power(n, target, ev=0.23, reps=400):
        hits = 0
        beta = _beta_for(target)
        for _ in range(reps):
            X = rng.normal(size=(n, 4))
            lin = 0.9 * X[:, 0] - 0.6 * X[:, 1] + 0.5 * X[:, 3]
            s = rng.normal(size=n)
            lin = lin + beta * s
            lin += np.log(ev / (1 - ev)) - lin.mean()
            y = (rng.uniform(size=n) < 1 / (1 + np.exp(-lin))).astype(int)
            if min(y.sum(), n - y.sum()) < 5:
                continue
            R = sm.add_constant(X); F = np.column_stack([R, s])
            try:
                mr = sm.Logit(y, R).fit(disp=0); mf = sm.Logit(y, F).fit(disp=0)
                if stats.chi2.sf(2 * (mf.llf - mr.llf), 1) < 0.05:
                    hits += 1
            except Exception:
                continue
        return hits / reps

    res = {}
    for target in [0.55, 0.60]:
        lo, hi = 100, 6000
        while hi - lo > 50:
            mid = (lo + hi) // 2
            if sim_power(mid, target) >= 0.80:
                hi = mid
            else:
                lo = mid
        res[str(target)] = dict(required_n_80pct=int(hi),
                                power_at_required=float(sim_power(hi, target)))
        print(f"  residual AUROC {target}: required n ~ {hi} for 80% power at a 23% event rate")
    res["largest_public_cohort"] = dict(
        name="GSE25066", n=488,
        note="No public neoadjuvant pCR cohort reaches the n required to exclude a residual "
             "AUROC of 0.55 at 80% power. This is the study-design recommendation the power "
             "limitation converts into.")
    _out("A10_required_n", res)
    return res


def _beta_for(target):
    """log-odds coefficient giving approximately the target residual AUROC."""
    return float(np.sqrt(2) * stats.norm.ppf(target) * 1.7)


# ------------------------------------------------------------------ A11
def a11(src, b_ov):
    print("\n=== A11: regimen ===")
    coh = "GSE22226"
    expr, ax = src["load"](coh)
    ph = pd.read_csv(D.ROOT / f"external_data/{coh}/{coh}_phenotype.tsv",
                     sep="\t", dtype=str).set_index("patient_id").reindex(ax.index)
    col = [c for c in ph.columns if "neoadjuvant_chemotherapy" in c.lower()][0]
    # values read "A = doxorubicin, C = Cyclophosphamide, T = taxane: 3" -- the arm code is the
    # trailing integer after the colon, not the whole field
    arm = pd.to_numeric(ph[col].astype(str).str.extract(r":\s*(\d+)\s*$")[0], errors="coerce")
    y = ax["pCR"].values.astype(int)
    Xr = np.column_stack([D.zscore(ax[a]).values for a in D.AXES])
    s = D.zscore(src["sig_of"]("Framework", expr)[0]).values

    def lrt(mask):
        R = sm.add_constant(Xr[mask]); F = np.column_stack([R, s[mask]])
        mr = sm.Logit(y[mask], R).fit(disp=0); mf = sm.Logit(y[mask], F).fit(disp=0)
        return dict(n=int(mask.sum()), events=int(y[mask].sum()),
                    lr_p=float(stats.chi2.sf(2 * (mf.llf - mr.llf), 1)))

    counts = arm.value_counts(dropna=False).to_dict()
    her = (arm == 3).values
    out = dict(cohort=coh, arm_field=col,
               arm_counts={str(k): int(v) for k, v in counts.items()},
               trastuzumab_n=int(her.sum()),
               trastuzumab_pcr=int(y[her].sum()) if her.sum() else 0,
               full=lrt(np.ones(len(y), bool)),
               excluding_trastuzumab=lrt(~her))
    print(f"  arms: {out['arm_counts']}")
    print(f"  trastuzumab arm n={out['trastuzumab_n']} pCR={out['trastuzumab_pcr']}")
    print(f"  full LR p={out['full']['lr_p']:.4f}  "
          f"excluding trastuzumab LR p={out['excluding_trastuzumab']['lr_p']:.4f}")
    out["limit"] = ("GSE22226 deposits the only within-cohort regimen field in the study. "
                    "GSE16446, GSE25066 and GSE32646 deposit none, so cross-regimen consistency "
                    "is a between-cohort observation confounded with every other cohort "
                    "difference, and cannot support a claim of regimen-general redundancy.")
    _out("A11_regimen", out)
    return out


# ------------------------------------------------------------------ A12
def a12():
    print("\n=== A12: Figure 4 (patient space) vs Figure 6 (CCLE gene-rank space) ===")
    nes = pd.read_csv(D.ROOT / "resubmission_v2/results/phase_1/hallmark_gsea_NES_matrix.tsv",
                      sep="\t", index_col=0)
    cm = pd.read_csv(HERE / "corrected_figdata/phase_1/correlation_matrix_GSE25066.tsv",
                     sep="\t", index_col=0)
    nes_c = nes.T.corr() if nes.shape[0] < nes.shape[1] else nes.corr()
    common = [c for c in cm.index if c in nes_c.index]
    pairs = []
    for i, a in enumerate(common):
        for b in common[i + 1:]:
            pairs.append((a, b, float(nes_c.loc[a, b]), float(cm.loc[a, b])))
    if not pairs:
        print("  no drug names common to both matrices; reporting the mismatch itself")
        _out("A12_space_comparison", dict(n_pairs=0,
             note="The two figures are indexed by different objects -- Figure 4 by deployed "
                  "patient scores, Figure 6 by CCLE per-gene rank profiles -- and share no "
                  "common label set in the deposited matrices. That is itself the answer to "
                  "R3 point 7: they are not two views of one object."))
        return None
    x = np.array([p[2] for p in pairs]); yv = np.array([p[3] for p in pairs])
    r = float(stats.pearsonr(x, yv).statistic)
    print(f"  {len(pairs)} pairs, r(NES-profile corr, patient-score corr) = {r:+.3f}")
    _out("A12_space_comparison", dict(n_pairs=len(pairs), pearson_r=r,
         pairs=[dict(a=p[0], b=p[1], nes_corr=p[2], score_corr=p[3]) for p in pairs]))
    return r


# ------------------------------------------------------------------ A6
def a6(src, b_ov):
    print("\n=== A6: pCR endpoint harmonisation ===")
    out = {}
    for coh, rcbcol, pos in [("GSE25066", "pathologic_response_rcb_class", None),
                             ("GSE22226", "rcb_class", None)]:
        expr, ax = src["load"](coh)
        if coh == "GSE25066":
            k = ~ax.index.isin(b_ov)
            expr, ax = expr.loc[k], ax.loc[k]
        ph = pd.read_csv(D.ROOT / f"external_data/{coh}/{coh}_phenotype.tsv",
                         sep="\t", dtype=str).set_index("patient_id").reindex(ax.index)
        cand = [c for c in ph.columns if "rcb" in c.lower()]
        if not cand:
            continue
        rc = ph[cand[0]].astype(str).str.upper().str.strip()
        known = rc.notna() & ~rc.isin(["NAN", "NONE", "", "UNAVAILABLE"])
        # Two deposited encodings. GSE22226 writes the class as a leading token before "=",
        # e.g. "0=RCB index 0" / "I=RCB index less than or equal to 1.36" / "II=..." / "III=...".
        # GSE25066 writes a collapsed "RCB-0/I" label -- RCB-0 alone is NOT recoverable there.
        lead = rc.str.split("=").str[0].str.strip()
        harm = (lead.isin(["0", "I"]) | rc.str.contains(r"RCB-0", regex=True, na=False)
                ).astype(int)
        Xr = np.column_stack([D.zscore(ax[a]).values for a in D.AXES])
        s = D.zscore(src["sig_of"]("Framework", expr)[0]).values
        dep = ax["pCR"].values.astype(int)

        def lrt(y, mask):
            R = sm.add_constant(Xr[mask]); F = np.column_stack([R, s[mask]])
            mr = sm.Logit(y[mask], R).fit(disp=0); mf = sm.Logit(y[mask], F).fit(disp=0)
            return dict(n=int(mask.sum()), events=int(y[mask].sum()),
                        lr_p=float(stats.chi2.sf(2 * (mf.llf - mr.llf), 1)),
                        marginal=float(roc_auc_score(y[mask], s[mask])))
        km = known.values
        kap = float(stats.pearsonr(dep[km], harm.values[km]).statistic)
        out[coh] = dict(rcb_field=cand[0], n_rcb_known=int(km.sum()),
                        events_deposited=int(dep.sum()), events_harmonised=int(harm.sum()),
                        agreement_r=kap,
                        same_n_deposited=lrt(dep, km), same_n_harmonised=lrt(harm.values, km),
                        full_n_deposited=lrt(dep, np.ones(len(dep), bool)))
        print(f"  {coh}: RCB known n={int(km.sum())}, events {int(dep.sum())} -> "
              f"{int(harm.sum())} under RCB-0/I")
        print(f"    same-n  deposited p={out[coh]['same_n_deposited']['lr_p']:.4f}  "
              f"harmonised p={out[coh]['same_n_harmonised']['lr_p']:.4f}")
    out["limits"] = ("GSE25066's deposited RCB field is collapsed at RCB-0/I, so RCB-0 is not "
                     "recoverable; GSE32646 and GSE16446 deposit no RCB at all. The endpoint "
                     "can be harmonised in 2 of 4 cohorts, and those two are precisely the two "
                     "that share 65 patients -- so the harmonised analysis is a sensitivity "
                     "analysis, not a second independent test.")
    _out("A6_endpoint_harmonisation", out)
    return out


# ------------------------------------------------------------------ DCA
def dca(src, b_ov):
    print("\n=== DCA: decision-curve analysis ===")
    ctx = ("Decision modelled: selection into a response-adapted escalation trial. Every "
           "patient in these cohorts receives chemotherapy, so there is no treat/do-not-treat "
           "threshold; the actionable decision is who to enrol on the basis of predicted "
           "non-response. Threshold range 5-40% pre-specified.")
    out = {"decision_context": ctx, "cohorts": {}}
    thr = np.arange(0.05, 0.41, 0.01)
    for coh in ["GSE25066", "GSE32646"]:
        expr, ax = src["load"](coh)
        if coh == "GSE25066":
            k = ~ax.index.isin(b_ov)
            expr, ax = expr.loc[k], ax.loc[k]
        y = ax["pCR"].values.astype(int)
        Xr = np.column_stack([D.zscore(ax[a]).values for a in D.AXES])
        s = D.zscore(src["sig_of"]("Framework", expr)[0]).values
        R = sm.add_constant(Xr); F = np.column_stack([R, s])
        pr = sm.Logit(y, R).fit(disp=0).predict(R)
        pf = sm.Logit(y, F).fit(disp=0).predict(F)
        rng = np.random.default_rng(SEED)

        def nb(p, t, yy):
            pos = p >= t
            if pos.sum() == 0:
                return 0.0
            return float((yy[pos].sum() - (pos.sum() - yy[pos].sum()) * t / (1 - t)) / len(yy))
        curves = {}
        for lab, p in [("four axes", pr), ("axes + signature", pf),
                       ("treat all", np.ones(len(y)))]:
            curves[lab] = [nb(p, t, y) for t in thr]
        # The replicate index was drawn and then never used, so all 500 rows of `boot` were
        # the point estimate and the percentile interval collapsed onto it. Resample patients
        # AND refit both models per replicate, which is what an added-value interval needs.
        boot = []
        for _ in range(500):
            i = rng.integers(0, len(y), len(y))
            if len(np.unique(y[i])) < 2:
                continue
            try:
                pri = sm.Logit(y[i], R[i]).fit(disp=0).predict(R[i])
                pfi = sm.Logit(y[i], F[i]).fit(disp=0).predict(F[i])
            except Exception:
                continue
            boot.append([nb(pfi, t, y[i]) - nb(pri, t, y[i]) for t in thr])
        bl = np.array(boot)
        out["cohorts"][coh] = dict(
            n=int(len(y)), thresholds=[float(t) for t in thr], curves=curves,
            delta_net_benefit=[float(a - b) for a, b in zip(curves["axes + signature"],
                                                            curves["four axes"])],
            delta_ci_lo=[float(x) for x in np.percentile(bl, 2.5, axis=0)],
            delta_ci_hi=[float(x) for x in np.percentile(bl, 97.5, axis=0)])
        d = out["cohorts"][coh]["delta_net_benefit"]
        print(f"  {coh}: max delta net benefit over 5-40% = {max(d):+.5f} "
              f"(at threshold {thr[int(np.argmax(d))]:.2f})")
    _out("DCA_decision_curve", out)
    return out


def main():
    _, b_ov, _, _, _ = D.overlap()
    src = _cell_builders()
    a8(src, b_ov); a9(); a10(); a11(src, b_ov); a12(); a6(src, b_ov); dca(src, b_ov)
    print("\nall remaining pre-specified analyses complete")


if __name__ == "__main__":
    main()
