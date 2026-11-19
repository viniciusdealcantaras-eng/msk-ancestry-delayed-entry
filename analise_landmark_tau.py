"""Landmark with delayed entry and conditional Kendall's tau on all patients, by ancestry
and by ancestry within cancer type.

1. Landmark at 6 and 12 months: every patient alive at the landmark is kept. Patients
   already sequenced enter at the landmark and those sequenced later enter at their
   sequencing date (delayed entry). The earlier restricted version (only patients
   sequenced by the landmark) is reported next to it for comparison.
2. Conditional Kendall's tau (Martin & Betensky 2005) computed on all patients of each
   ancestry group, overall and within each cancer type, with a standard error from the
   projection of the U-statistic (delta method for the ratio). Replaces the earlier
   test in random subsamples of up to 3,000 patients.
3. Reason for the difference between the 19,424 patients analysed from diagnosis and
   the 19,278 in the Cox models.

Input : out/banco_ancestralidade_MSK.csv (from build_banco.py)
Output: out/resultados_landmark_tau.json (git-ignored, never redistributed)
"""
import json, warnings
from math import erfc, sqrt
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
warnings.filterwarnings("ignore")

D = pd.read_csv("out/banco_ancestralidade_MSK.csv", low_memory=False)
LAB = {1: "nonASJ-EUR", 2: "ASJ-EUR", 3: "ADMIX", 4: "AFR", 5: "EAS", 6: "SAS", 7: "NAM"}
GR = ["ASJ-EUR", "ADMIX", "AFR", "EAS", "SAS", "NAM"]
D["anc"] = D.ancestry_cat.map(LAB)
A = D[D.anc.notna()]
CH = A[(A.in_chord == 1) & (A.chord_dx_valid == 1) & A.dx_to_seq_months.notna()].copy()
CH["ct"] = CH.chord_cancer_type
CH["female"] = (CH.sex == "Female").astype(int)
an = []
for g in GR:
    CH["a_" + g.replace("-", "_")] = (CH.anc == g).astype(int); an.append("a_" + g.replace("-", "_"))
S = CH[CH.os_months_dx.notna() & (CH.os_months_dx > CH.dx_to_seq_months - 1e-9)].copy()
S["T"] = S.os_months_dx.clip(lower=S.dx_to_seq_months + 1e-4)
S["E"] = S.chord_os_status.astype(int)
S["L"] = S.dx_to_seq_months
covs = an + ["age_dx", "female"]
R = {"n_surv": int(len(S)), "deaths_surv": int(S.E.sum()),
     "missing_age_dx": int(S.age_dx.isna().sum()), "deaths_missing_age": int(S[S.age_dx.isna()].E.sum())}

def cox(df, entry):
    cols = ["T", "E"] + covs + ["ct", entry]
    d = df[cols].dropna().copy()
    for v in cols:
        if v != "ct": d[v] = d[v].astype("float64")
    c = CoxPHFitter().fit(d, "T", "E", strata=["ct"], entry_col=entry, fit_options={"step_size": 0.5})
    hr = {v: [round(float(np.exp(c.params_[v])), 2)] + [round(float(x), 2) for x in np.exp(c.confidence_intervals_.loc[v].values)] for v in an}
    return {"n": int(len(d)), "deaths": int(d.E.sum()), "hr": hr}

S["Lent"] = S.L
R["corrected_check"] = cox(S, "Lent")

# ---------- 1. landmark ----------
R["landmark"] = {}
for t0 in [6, 12]:
    a = S[S["T"] > t0].copy(); a["en"] = np.maximum(a.L, t0)
    r = S[(S.L <= t0) & (S["T"] > t0)].copy(); r["en"] = float(t0)
    R["landmark"][t0] = {"delayed_entry": cox(a, "en"), "restricted_original": cox(r, "en")}

# ---------- 2. conditional Kendall's tau with U-statistic SE ----------
def ctau(df):
    L, X, Dv = df.L.values, df["T"].values, df.E.values.astype(bool)
    n = len(L); a = np.zeros(n); b = np.zeros(n)
    for i0 in range(0, n, 400):
        i = np.arange(i0, min(i0 + 400, n))[:, None]; j = np.arange(n)[None, :]
        Li, Lj, Xi, Xj, Di, Dj = L[i], L[j], X[i], X[j], Dv[i], Dv[j]
        c = (np.maximum(Li, Lj) <= np.minimum(Xi, Xj)) & (j != i)
        c &= np.where(Xi < Xj, Di, np.where(Xj < Xi, Dj, Di & Dj))
        s = np.sign((Li - Lj) * (Xi - Xj)) * c
        a[i[:, 0]] = s.sum(1); b[i[:, 0]] = c.sum(1)
    if b.sum() == 0: return None
    tau = a.sum() / b.sum()
    h = (a - tau * b) / (n - 1); ub = b.sum() / (n * (n - 1))
    se = 2 * np.std(h, ddof=1) / sqrt(n) / ub
    return {"n": int(n), "pairs": int(b.sum() // 2), "tau": round(float(tau), 3),
            "ci": [round(float(tau - 1.96 * se), 3), round(float(tau + 1.96 * se), 3)],
            "p": float(erfc(abs(tau / se) / sqrt(2)))}

R["tau_by_ancestry"] = {g: ctau(S[S.anc == g]) for g in ["nonASJ-EUR"] + GR}
R["tau_by_cancer"] = {ct: ctau(S[S.ct == ct]) for ct in sorted(S.ct.unique())}
R["tau_by_ancestry_cancer"] = {f"{g} | {ct}": ctau(x) for (g, ct), x in S.groupby(["anc", "ct"]) if len(x) >= 40}

json.dump(R, open("out/resultados_landmark_tau.json", "w"), indent=1)
print(json.dumps(R, indent=1))
