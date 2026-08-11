#!/usr/bin/env python3
"""
phase_9_revalidation — test the four external-validation adjustments proposed by
the 4-expert panel, on EXISTING data only. Writes nothing into locked dirs.

(1) corrected TOST decision rule (gate equivalence on UPPER CI, not lower)
(2) nested-model likelihood-ratio test as the primary estimator (Expert B)
(3) fixed-effect inverse-variance meta-analysis pooling + forest (Expert A)
(4) orthogonal spike-in POWER CURVE (Expert D) -- replaces the broken control
"""
import json, glob, os, warnings
import numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm
warnings.filterwarnings("ignore")
np.random.seed(42)

REPO = "/home/holiday01/2026_ISMB_code/revise_bioadv"
R = REPO + "/resubmission_v2/results"
OUT = R + "/phase_9_revalidation"
os.makedirs(OUT, exist_ok=True)
AXES = ["P", "H", "B", "H_PAM50basal"]
rep = {}

# ---------------- helpers ----------------
def hm_se(auc, npos, nneg):
    a = min(max(auc, 1e-6), 1 - 1e-6)
    Q1 = a / (2 - a); Q2 = 2 * a * a / (1 + a)
    v = (a * (1 - a) + (npos - 1) * (Q1 - a * a) + (nneg - 1) * (Q2 - a * a)) / (npos * nneg)
    return np.sqrt(max(v, 0))

def auroc(y, s):
    return roc_auc_score(y, s)

def load_cohort(cohort):
    """return df with AXES + sig + y  (pCR for BC cohorts, event5y for METABRIC)."""
    if cohort in ("GSE16446", "GSE25066", "GSE22226"):
        dox = pd.read_csv(f"{REPO}/external_data/{cohort}/score_DOXORUBICIN.tsv", sep="\t")
        ax = pd.read_csv(f"{R}/phase_0/axes_{cohort}.tsv", sep="\t")
        df = dox.merge(ax, left_on="patient_id", right_on="sample_id", how="inner")
        df = df.rename(columns={"signature_score": "sig", "pcr": "y"})
    elif cohort == "GSE32646":
        sc = pd.read_csv(f"{R}/phase_2/score_DOXORUBICIN_GSE32646.tsv", sep="\t")
        ax = pd.read_csv(f"{R}/phase_2/axes_GSE32646.tsv", sep="\t")
        df = sc.merge(ax, left_on="patient_id", right_on="sample_id", how="inner")
        df = df.rename(columns={"sig_dox": "sig", "pcr": "y"})
    elif cohort == "METABRIC":
        sc = pd.read_csv(f"{R}/phase_4/metabric/score_DOXORUBICIN_METABRIC.tsv", sep="\t")
        ax = pd.read_csv(f"{R}/phase_4/metabric/axes_METABRIC.tsv", sep="\t").drop(columns=["event5y"], errors="ignore")
        df = sc.merge(ax, on="sample_id", how="inner")
        df = df.rename(columns={"sig_dox": "sig", "event5y": "y"})
    keep = AXES + ["sig", "y"]
    df = df[keep].apply(pd.to_numeric, errors="coerce").dropna()
    df = df[df["y"].isin([0, 1])]
    return df.reset_index(drop=True)

def zscore(a):
    a = np.asarray(a, float); return (a - a.mean()) / (a.std(ddof=0) + 1e-12)

# ================= PART 2: NESTED-MODEL LRT (Expert B) =================
print("=" * 78); print("PART 2  NESTED-MODEL LIKELIHOOD-RATIO TEST  (primary estimator)"); print("=" * 78)
print(f"{'cohort':10} {'endpt':6} {'n':>4} {'ev':>4} | {'LR_p':>7} {'OR/SD':>7} {'OR_CI95':>16} {'dMcF':>7} {'dAIC':>6} {'dBIC':>6} | {'AUC_4ax':>7} {'AUC_+sig':>8} {'resid':>6}")
nested = []
for c in ["GSE16446", "GSE25066", "GSE22226", "GSE32646", "METABRIC"]:
    df = load_cohort(c)
    y = df["y"].values.astype(int)
    Xax = np.column_stack([zscore(df[a]) for a in AXES])
    sig = zscore(df["sig"])
    Xr = sm.add_constant(Xax)
    Xf = sm.add_constant(np.column_stack([Xax, sig]))
    mr = sm.Logit(y, Xr).fit(disp=0, maxiter=200)
    mf = sm.Logit(y, Xf).fit(disp=0, maxiter=200)
    ll0 = sm.Logit(y, np.ones((len(y), 1))).fit(disp=0).llf
    LR = 2 * (mf.llf - mr.llf); p = stats.chi2.sf(LR, 1)
    beta = mf.params[-1]; se = mf.bse[-1]
    OR = np.exp(beta); OR_lo = np.exp(beta - 1.96 * se); OR_hi = np.exp(beta + 1.96 * se)
    dMcF = (1 - mf.llf / ll0) - (1 - mr.llf / ll0)
    dAIC = mf.aic - mr.aic; dBIC = mf.bic - mr.bic
    aucr = auroc(y, mr.predict(Xr)); aucf = auroc(y, mf.predict(Xf))
    # residual-AUROC cross-check (the OLD estimator)
    res = sig - sm.OLS(sig, Xr).fit().predict(Xr)
    rauc = auroc(y, res)
    endpt = "OS5y" if c == "METABRIC" else "pCR"
    print(f"{c:10} {endpt:6} {len(y):>4} {int(y.sum()):>4} | {p:7.3f} {OR:7.2f} [{OR_lo:5.2f},{OR_hi:5.2f}] {dMcF:7.4f} {dAIC:+6.1f} {dBIC:+6.1f} | {aucr:7.3f} {aucf:8.3f} {rauc:6.3f}")
    nested.append(dict(cohort=c, endpoint=endpt, n=int(len(y)), events=int(y.sum()),
                       LR_p=float(p), OR_per_SD=float(OR), OR_ci=[float(OR_lo), float(OR_hi)],
                       beta=float(beta), beta_se=float(se), dMcFadden=float(dMcF),
                       dAIC=float(dAIC), dBIC=float(dBIC),
                       auc_4axis=float(aucr), auc_plus_sig=float(aucf), residual_auroc=float(rauc)))
rep["nested_lrt"] = nested
# pooled added-effect (standardized beta), fixed-effect, BC-pCR cohorts only
bc = [r for r in nested if r["endpoint"] == "pCR"]
w = np.array([1 / r["beta_se"] ** 2 for r in bc]); b = np.array([r["beta"] for r in bc])
bpool = (w * b).sum() / w.sum(); bse = np.sqrt(1 / w.sum())
print(f"\nPooled added log-OR/SD (4 pCR cohorts, FE): beta={bpool:+.3f}  SE={bse:.3f}  "
      f"OR/SD={np.exp(bpool):.2f} [{np.exp(bpool-1.96*bse):.2f},{np.exp(bpool+1.96*bse):.2f}]  "
      f"p={2*stats.norm.sf(abs(bpool/bse)):.3f}")
rep["nested_pooled_beta"] = dict(beta=float(bpool), se=float(bse), OR=float(np.exp(bpool)),
                                 p=float(2 * stats.norm.sf(abs(bpool / bse))))

# ================= PART 3 + 1: META-POOL residual AUROC + corrected TOST rule =================
print("\n" + "=" * 78); print("PART 3+1  RESIDUAL-AUROC META-ANALYSIS + CORRECTED TOST RULE"); print("=" * 78)
cells = []
def add_cell(name, cohort, point, lo, hi, marg, n):
    cells.append(dict(name=name, cohort=cohort, point=point, lo=lo, hi=hi, marginal=marg, n=n))

# framework held-outs
j = json.load(open(f"{R}/phase_2/gse32646_residual_test.json"))
pr = j["primary_residual_auroc"]
add_cell("Framework(Dox)", "GSE32646", pr["point"], pr["ci_low"], pr["ci_high"],
         j["raw_sig_dox_auroc"]["point"], j["n_with_pcr"])
mt = pd.read_csv(f"{R}/phase_4/metabric/metabric_residual_summary.tsv", sep="\t")
row = mt.iloc[0]
add_cell("Mammaprint*", "METABRIC", float(row.get("residual_auroc_point", row.get("auroc_point", np.nan))),
         float(row.get("residual_ci_low", row.get("ci_low_95", np.nan))),
         float(row.get("residual_ci_high", row.get("ci_high_95", np.nan))),
         float(row.get("marginal_auroc", row.get("raw_auroc", np.nan))), int(row.get("n", row.get("n_used", 0))))
# 9 published signatures
for f in sorted(glob.glob(f"{R}/phase_6/*_residual.json") + glob.glob(f"{R}/phase_7/*_residual.json") + glob.glob(f"{R}/phase_8/*_residual.json")):
    d = json.load(open(f)); ra = d["residual_auroc"]; ma = d.get("marginal_auroc", {})
    nm = os.path.basename(f).replace("_residual.json", "")
    add_cell(nm, "?", ra["point"], ra["ci_low"], ra["ci_high"],
             ma.get("point", np.nan) if isinstance(ma, dict) else np.nan, d.get("n_with_pcr", 0))

cdf = pd.DataFrame(cells).dropna(subset=["point", "lo", "hi"])
cdf["se"] = (cdf["hi"] - cdf["lo"]) / (2 * 1.96)

def pool(sub):
    w = 1 / sub["se"] ** 2
    m = (w * sub["point"]).sum() / w.sum(); se = np.sqrt(1 / w.sum())
    Q = (w * (sub["point"] - m) ** 2).sum(); dfree = len(sub) - 1
    I2 = max(0, (Q - dfree) / Q) * 100 if Q > 0 else 0
    return m, se, m - 1.96 * se, m + 1.96 * se, Q, dfree, I2

def tost(m, se, margin):  # H1: true < margin
    z = (m - margin) / se; return z, stats.norm.cdf(z)

print(f"{'cell':22} {'cohort':9} {'n':>4} {'marg':>6} {'resid':>6} {'95%CI':>16}  OLD->NEW verdict")
for _, r in cdf.iterrows():
    old = "EQUIV" if (r.point <= 0.55 and r.lo <= 0.55) else ("RESID" if (r.point > 0.60 and r.lo > 0.55) else "GRAY")
    if r.marginal < 0.60 or np.isnan(r.marginal):
        new = "VACUOUS(marg<.60)"
    elif r.hi < 0.60:
        new = "EQUIVALENT"
    elif r.lo > 0.55:
        new = "RESIDUAL!"
    else:
        new = "INCONCLUSIVE"
    print(f"{r['name']:22} {r['cohort']:9} {int(r.n):>4} {r.marginal:6.3f} {r.point:6.3f} [{r.lo:5.3f},{r.hi:5.3f}]  {old:5} -> {new}")

print("\n--- POOLED (fixed-effect inverse-variance) ---")
for label, sub in [("ALL cells", cdf),
                   ("held-out framework only (GSE32646+METABRIC)", cdf[cdf.cohort.isin(["GSE32646", "METABRIC"])]),
                   ("exclude METABRIC (near-null marginal)", cdf[cdf.cohort != "METABRIC"]),
                   ("n>=285 subset", cdf[cdf.n >= 285])]:
    if len(sub) == 0:
        continue
    m, se, lo, hi, Q, dfree, I2 = pool(sub)
    _, p55 = tost(m, se, 0.55); _, p60 = tost(m, se, 0.60)
    print(f"{label:46} k={len(sub):2}  pooled={m:.3f} [{lo:.3f},{hi:.3f}]  "
          f"Q={Q:4.1f}/df{dfree} I2={I2:4.1f}%  TOST<.55 p={p55:.1e}  TOST<.60 p={p60:.1e}")
    if label == "ALL cells":
        rep["meta_all"] = dict(k=int(len(sub)), pooled=float(m), ci=[float(lo), float(hi)],
                               Q=float(Q), I2=float(I2), tost_055_p=float(p55), tost_060_p=float(p60))
cdf.to_csv(f"{OUT}/meta_cells.tsv", sep="\t", index=False)

# ================= PART 4: ORTHOGONAL SPIKE-IN POWER CURVE (Expert D) =================
print("\n" + "=" * 78); print("PART 4  ORTHOGONAL SPIKE-IN POWER CURVE (corrected positive control)"); print("=" * 78)
print("plant a known-strength signal ORTHOGONAL to the 4 axes; detection = corrected rule fires (ci_low>0.55)")
grid = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75]
NSIM = 1000
power = {}
print(f"{'cohort':10} {'n':>4} " + " ".join(f"A={g:.2f}" for g in grid) + "   MDE@80%")
for c in ["GSE16446", "GSE22226", "GSE32646", "GSE25066"]:
    df = load_cohort(c)
    y = df["y"].values.astype(int); npos = int(y.sum()); nneg = len(y) - npos
    X = sm.add_constant(np.column_stack([zscore(df[a]) for a in AXES]))
    H = X @ np.linalg.pinv(X)  # hat matrix; residual maker = I-H
    det = []
    for a in grid:
        mu = np.sqrt(2) * stats.norm.ppf(a)  # normal-shift giving AUROC=a
        hits = 0
        for _ in range(NSIM):
            z = np.random.randn(len(y)) + mu * y           # planted signal correlated w/ y
            r = z - H @ z                                  # project OUT the 4 axes -> pure residual
            au = auroc(y, r)
            ci_lo = au - 1.96 * hm_se(au, npos, nneg)
            if ci_lo > 0.55:
                hits += 1
        det.append(hits / NSIM)
    # MDE@80%: interpolate smallest planted AUROC with power>=0.8
    mde = next((g for g, d in zip(grid, det) if d >= 0.8), None)
    if mde is None:
        mde_s = ">0.75"
    else:
        mde_s = f"{mde:.2f}"
    power[c] = dict(n=int(len(y)), grid=grid, detection=det, mde80=mde_s)
    print(f"{c:10} {len(y):>4} " + " ".join(f"{d:5.2f}" for d in det) + f"   {mde_s}")
rep["power_curve"] = power

# ---- figures ----
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    # forest
    fig, ax = plt.subplots(figsize=(7, 5))
    sub = cdf.sort_values("point").reset_index(drop=True)
    ax.errorbar(sub.point, range(len(sub)),
                xerr=[sub.point - sub.lo, sub.hi - sub.point], fmt="o", color="#333", capsize=2)
    m, se, lo, hi, *_ = pool(cdf)
    ax.axvspan(lo, hi, color="#4C72B0", alpha=.25); ax.axvline(m, color="#4C72B0", lw=2)
    ax.axvline(0.5, color="grey", ls=":"); ax.axvline(0.55, color="orange", ls="--"); ax.axvline(0.60, color="red", ls="--")
    ax.set_yticks(range(len(sub))); ax.set_yticklabels([f"{r['name']}/{r.cohort}" for _, r in sub.iterrows()], fontsize=7)
    ax.set_xlabel("residual AUROC"); ax.set_title(f"Pooled residual AUROC = {m:.3f} [{lo:.3f},{hi:.3f}]  (I2={pool(cdf)[6]:.0f}%)")
    fig.tight_layout(); fig.savefig(f"{OUT}/forest_residual_auroc.png", dpi=150); plt.close(fig)
    # power curve
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for c, d in power.items():
        ax.plot(grid, d["detection"], "o-", label=f"{c} (n={d['n']})")
    ax.axhline(0.8, color="grey", ls="--"); ax.axvline(0.55, color="orange", ls=":")
    ax.set_xlabel("planted residual AUROC (orthogonal to axes)"); ax.set_ylabel("detection probability")
    ax.set_title("Spike-in power curve (corrected positive control)"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{OUT}/power_curve.png", dpi=150); plt.close(fig)
    print("\nfigures: forest_residual_auroc.png , power_curve.png")
except Exception as e:
    print("fig skipped:", e)

json.dump(rep, open(f"{OUT}/revalidation_summary.json", "w"), indent=1)
print(f"\nwrote {OUT}/revalidation_summary.json")
