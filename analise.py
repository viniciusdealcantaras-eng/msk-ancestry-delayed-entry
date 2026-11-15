"""All analyses of the paper (MSK-IMPACT 50K + MSK-CHORD).

Input : out/banco_ancestralidade_MSK.csv (from build_banco.py)
Output: out/resultados.json and auxiliary files for figuras.R and verificar.R
        (all in out/, git-ignored, never redistributed)
"""
import os, json, warnings
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
from lifelines import CoxPHFitter, KaplanMeierFitter
warnings.filterwarnings("ignore")

OUT = "out/"
os.makedirs(OUT, exist_ok=True)
D = pd.read_csv(OUT + "banco_ancestralidade_MSK.csv", low_memory=False)
D = D.replace({"": np.nan})
LAB = {1: "nonASJ-EUR", 2: "ASJ-EUR", 3: "ADMIX", 4: "AFR", 5: "EAS", 6: "SAS", 7: "NAM"}
GR = ["ASJ-EUR", "ADMIX", "AFR", "EAS", "SAS", "NAM"]
D["anc"] = D.ancestry_cat.map(LAB)
R = {}

def dummies(df):
    for g in GR: df["a_" + g.replace("-", "_")] = (df.anc == g).astype(int)
    return [ "a_" + g.replace("-", "_") for g in GR]

def cox(df, dur, ev, covs, strata, entry=None):
    cols = [dur, ev] + covs + [strata] + ([entry] if entry else [])
    d = df[cols].dropna().copy()
    for v in cols:
        if v != strata: d[v] = d[v].astype("float64")
    covs = [v for v in covs if d[v].nunique() > 1]
    d = d[[dur, ev] + covs + [strata] + ([entry] if entry else [])]
    c = CoxPHFitter()
    c.fit(d, duration_col=dur, event_col=ev, strata=[strata], entry_col=entry, fit_options={"step_size":0.5})
    out = {}
    for v in covs:
        lo, hi = c.confidence_intervals_.loc[v]
        out[v] = [float(np.exp(c.params_[v])), float(np.exp(lo)), float(np.exp(hi)), float(c.summary.loc[v, "p"])]
    return out, len(d), int(d[ev].sum()), c

# ================= 1. Descricao (50K) =================
D.loc[D.tmb < 0, "tmb"] = np.nan
A = D[D.anc.notna()].copy()
R["n_50k_total"] = int(len(D)); R["n_50k_anc"] = int(len(A)); R["n_unknown_anc"] = int(D.anc.isna().sum())
desc = {}
for g in ["nonASJ-EUR"] + GR:
    x = A[A.anc == g]
    desc[g] = {"n": len(x), "female": int((x.sex == "Female").sum()), "age_med": float(x.age_dx.median()),
               "age_q1": float(x.age_dx.quantile(.25)), "age_q3": float(x.age_dx.quantile(.75)),
               "metastatic_sample": int((x.sample_type == "Metastasis").sum()),
               "tmb_med": float(x.tmb.median()), "tmb_q1": float(x.tmb.quantile(.25)), "tmb_q3": float(x.tmb.quantile(.75)),
               "msi_high": int((x.msi_high == 1).sum()), "msi_eval": int(x.msi_high.notna().sum()),
               "fga_med": float(x.fga.median()), "fga_q1": float(x.fga.quantile(.25)), "fga_q3": float(x.fga.quantile(.75)),
               "deaths": int(x.os_status.sum()), "in_chord": int(x.in_chord.sum())}
R["desc"] = desc
top = A.cancer_type.value_counts().head(15).index.tolist()
R["top15"] = top
R["comp_by_cancer"] = {ct: A[A.cancer_type == ct].anc.value_counts().to_dict() for ct in top}

# ================= 2. Intervalo diagnostico-sequenciamento (CHORD) =================
CH = A[(A.in_chord == 1) & (A.chord_dx_valid == 1) & A.dx_to_seq_months.notna()].copy()
CH["stage"] = CH.chord_stage_dx.fillna("Unknown")
CH["ct"] = CH.chord_cancer_type
CH["female"] = (CH.sex == "Female").astype(int)
an = dummies(CH)
R["n_chord"] = int(len(CH))
iv = {}
for g in ["nonASJ-EUR"] + GR:
    x = CH[CH.anc == g].dx_to_seq_months
    iv[g] = {"n": int(len(x)), "med": float(x.median()), "q1": float(x.quantile(.25)), "q3": float(x.quantile(.75)),
             "gt12": int((x > 12).sum()), "gt24": int((x > 24).sum())}
R["interval"] = iv
fq = "dx_to_seq_months ~ " + " + ".join(an) + " + C(ct) + C(stage) + age_dx + female"
qr = smf.quantreg(fq, CH.dropna(subset=["age_dx"])).fit(q=.5, max_iter=5000)
R["qreg"] = {v: [float(qr.params[v]), float(qr.conf_int().loc[v, 0]), float(qr.conf_int().loc[v, 1]), float(qr.pvalues[v])] for v in an}
qr0 = smf.quantreg("dx_to_seq_months ~ " + " + ".join(an), CH).fit(q=.5, max_iter=5000)
R["qreg_crude"] = {v: [float(qr0.params[v]), float(qr0.conf_int().loc[v, 0]), float(qr0.conf_int().loc[v, 1])] for v in an}
for thr in [12, 24]:
    CH[f"gt{thr}"] = (CH.dx_to_seq_months > thr).astype(int)
    m = smf.glm(f"gt{thr} ~ " + " + ".join(an) + " + C(ct) + C(stage) + age_dx + female", CH.dropna(subset=["age_dx"]),
                family=sm.families.Poisson()).fit(cov_type="HC1")
    R[f"pr_gt{thr}"] = {v: [float(np.exp(m.params[v])), float(np.exp(m.conf_int().loc[v, 0])), float(np.exp(m.conf_int().loc[v, 1])), float(m.pvalues[v])] for v in an}

# ================= 3. Sobrevida desde o diagnostico: ingenuo x corrigido (CHORD) =================
S = CH[CH.os_months_dx.notna() & (CH.os_months_dx > CH.dx_to_seq_months - 1e-9)].copy()
S["T"] = S.os_months_dx.clip(lower=S.dx_to_seq_months + 1e-4)
S["E"] = S.chord_os_status.astype(int)
S["L"] = S.dx_to_seq_months
R["n_surv"] = int(len(S)); R["deaths_surv"] = int(S.E.sum())
covs = an + ["age_dx", "female"]
cs = {}
for nm, ent in [("naive", None), ("corrected", "L")]:
    hr, n, e, _ = cox(S, "T", "E", covs, "ct", ent)
    cs[nm] = {"hr": hr, "n": n, "deaths": e}
# corrigido + estadio
S["stage_IV"] = (S.stage == "IV").astype(int); S["stage_unk"] = (S.stage == "Unknown").astype(int)
hr, n, e, _ = cox(S, "T", "E", covs + ["stage_IV", "stage_unk"], "ct", "L"); cs["corrected_stage"] = {"hr": hr, "n": n, "deaths": e}
hr, n, e, _ = cox(S, "T", "E", covs + ["stage_IV", "stage_unk"], "ct", None); cs["naive_stage"] = {"hr": hr, "n": n, "deaths": e}
# desde o sequenciamento na mesma populacao (referencia)
S["Tseq"] = (S["T"] - S.L).clip(lower=1e-4)
hr, n, e, _ = cox(S, "Tseq", "E", covs, "ct", None); cs["from_seq_chord"] = {"hr": hr, "n": n, "deaths": e}
R["cox_dx"] = cs
# medianas KM por grupo
km = {}
for g in ["nonASJ-EUR"] + GR:
    x = S[S.anc == g]
    a = KaplanMeierFitter().fit(x["T"], x.E); b = KaplanMeierFitter().fit(x["T"], x.E, entry=x.L)
    km[g] = {"naive_med": float(a.median_survival_time_), "corr_med": float(b.median_survival_time_),
             "naive_5y": float(a.predict(60)), "corr_5y": float(b.predict(60))}
R["km_dx"] = km

# ================= 4. Sobrevida apos o sequenciamento (50K), modelos sequenciais =================
F = A[A.os_months_seq.notna() & (A.os_months_seq > 0) & A.os_status.notna()].copy()
F["female"] = (F.sex == "Female").astype(int)
F["metast"] = (F.sample_type == "Metastasis").astype(int)
F.loc[F.tmb < 0, "tmb"] = np.nan
F["log2tmb"] = np.log2(F.tmb + 1)
F["msih"] = F.msi_high
F["wgd_"] = F.wgd
GEN = ["mut_TP53", "mut_KRAS", "mut_PIK3CA", "del_CDKN2A", "amp_MYC"]
an50 = dummies(F)
F["E"] = F.os_status.astype(int)
need = an50 + ["age_dx", "female", "metast", "log2tmb", "fga", "wgd_"] + GEN + ["os_months_seq", "E", "cancer_type"]
Fc = F.dropna(subset=need).copy()
Fc["msih"] = Fc.msih.fillna(0); Fc["msi_na"] = F.loc[Fc.index, "msi_high"].isna().astype(int)
R["n_50k_cc"] = int(len(Fc)); R["deaths_50k_cc"] = int(Fc.E.sum())
m1 = an50 + ["age_dx", "female"]; m2 = m1 + ["metast"]; m3 = m2 + ["log2tmb", "msih", "msi_na", "fga", "wgd_"] + GEN
seq = {}
for nm, cv in [("m1", m1), ("m2", m2), ("m3", m3)]:
    hr, n, e, c = cox(Fc, "os_months_seq", "E", cv, "cancer_type")
    seq[nm] = {"hr": hr, "n": n, "deaths": e}
    if nm == "m3":
        R["m3_other"] = {v: hr[v] for v in cv if v not in an50}
R["cox_seq"] = seq
# Schoenfeld (modelo 3) - teste global e da ancestralidade
from lifelines.statistics import proportional_hazard_test
_, _, _, c3 = cox(Fc, "os_months_seq", "E", m3, "cancer_type")
pt = proportional_hazard_test(c3, Fc[m3 + ["os_months_seq", "E", "cancer_type"]], time_transform="rank")
R["schoenfeld_anc_p"] = {v: float(pt.summary.loc[v, "p"]) for v in an50}

# ================= 5. Concordancia ancestralidade x raca autodeclarada =================
K = A[A.in_chord == 1].copy()
rmap = {"White": "White", "Black Or African American": "Black", "Asian-Far East/Indian Subcont": "Asian",
        "Native American-Am Ind/Alaska": "NativeAm"}
K["race"] = K.chord_race.map(rmap).fillna("Other/Unknown")
R["crosstab"] = pd.crosstab(K.anc, K.race).to_dict()
amap = {"nonASJ-EUR": "White", "ASJ-EUR": "White", "AFR": "Black", "EAS": "Asian", "SAS": "Asian", "NAM": "NativeAm"}
k2 = K[K.anc.isin(amap) & K.race.isin(["White", "Black", "Asian", "NativeAm"])]
a = k2.anc.map(amap); b = k2.race
cats = ["White", "Black", "Asian", "NativeAm"]
po = (a == b).mean(); pe = sum((a == c).mean() * (b == c).mean() for c in cats)
R["kappa"] = {"n": int(len(k2)), "agreement": float(po), "kappa": float((po - pe) / (1 - pe))}
R["afr_black"] = {"afr_n": int((K.anc == "AFR").sum()), "afr_selfblack": int(((K.anc == "AFR") & (K.race == "Black")).sum()),
                  "black_n": int((K.race == "Black").sum()), "black_afr": int(((K.anc == "AFR") & (K.race == "Black")).sum()),
                  "black_admix": int(((K.anc == "ADMIX") & (K.race == "Black")).sum())}

# ================= 6. Sensibilidade =================
sens = {}
for lim in [120, 240]:
    x = S[S.L <= lim]
    hr, n, e, _ = cox(x, "T", "E", covs, "ct", "L"); sens[f"interval_le_{lim}"] = {"hr": hr, "n": n, "deaths": e}
# amostra primaria / metastatica (50K, modelo 3)
for st in ["Primary", "Metastasis"]:
    x = Fc[F.loc[Fc.index, "sample_type"] == st]
    hr, n, e, _ = cox(x, "os_months_seq", "E", [v for v in m3 if v != "metast"], "cancer_type")
    sens[f"m3_{st}"] = {"hr": hr, "n": n, "deaths": e}
# maior estadio registrado em vez do estadio ao diagnostico
S["stH"] = S.chord_stage_highest.fillna("Unknown")
S["stH_IV"] = (S.stH == "IV").astype(int); S["stH_unk"] = (S.stH == "Unknown").astype(int)
hr, n, e, _ = cox(S, "T", "E", covs + ["stH_IV", "stH_unk"], "ct", "L"); sens["corrected_stage_highest"] = {"hr": hr, "n": n, "deaths": e}
# landmark 6 e 12 meses
for t0 in [6, 12]:
    x = S[(S.L <= t0) & (S["T"] > t0)].copy(); x["t0"] = float(t0)
    hr, n, e, _ = cox(x, "T", "E", covs, "ct", "t0"); sens[f"landmark_{t0}"] = {"hr": hr, "n": n, "deaths": e}
R["sens"] = sens
# tau de Kendall condicional por grupo (Martin-Betensky), subamostra ate 3000
def ctau(L, X, d):
    n = len(L); conc = disc = 0
    for i in range(n):
        Li, Xi, di = L[i], X[i], d[i]
        j = np.arange(i + 1, n)
        Lj, Xj, dj = L[j], X[j], d[j]
        mx = np.maximum(Li, Lj); mn = np.minimum(Xi, Xj)
        comp = (mx <= mn) & (((Xi < Xj) & (di == 1)) | ((Xj < Xi) & (dj == 1)))
        s = np.sign((Li - Lj) * (Xi - Xj))[comp]
        conc += (s > 0).sum(); disc += (s < 0).sum()
    m = conc + disc
    return (conc - disc) / m if m else np.nan, int(m)
tau = {}
for g in ["nonASJ-EUR"] + GR:
    x = S[S.anc == g]
    x = x.sample(n=min(3000, len(x)), random_state=20261114)
    t, m = ctau(x.L.values, x["T"].values, x.E.values)
    tau[g] = {"n": int(len(x)), "tau": float(t), "pairs": m}
R["tau"] = tau

# ================= 7. Por tipo de cancer (AFR e EAS x nonASJ-EUR, corrigido) com FDR =================
from statsmodels.stats.multitest import multipletests
pc = []
for ct in S.ct.value_counts().index:
    x = S[(S.ct == ct) & S.anc.isin(["nonASJ-EUR", "AFR", "EAS"])].copy()
    if (x.anc == "AFR").sum() < 40 and (x.anc == "EAS").sum() < 40: continue
    x["afr"] = (x.anc == "AFR").astype(int); x["eas"] = (x.anc == "EAS").astype(int)
    cv = [v for v in ["afr", "eas"] if x[v].sum() >= 40 and x.loc[x[v] == 1, "E"].sum() >= 10]
    if not cv: continue
    for nm, ent in [("naive", None), ("corr", "L")]:
        d = x[["T", "E", "L", "age_dx", "female"] + cv].dropna()
        c = CoxPHFitter().fit(d if ent else d.drop(columns="L"), "T", "E", entry_col=ent)
        for v in cv:
            lo, hi = c.confidence_intervals_.loc[v]
            pc.append({"ct": ct, "grp": v.upper(), "model": nm, "n_grp": int(x[v].sum()), "hr": float(np.exp(c.params_[v])),
                       "lo": float(np.exp(lo)), "hi": float(np.exp(hi)), "p": float(c.summary.loc[v, "p"])})
pc = pd.DataFrame(pc)
for nm in ["naive", "corr"]:
    i = pc.model == nm
    pc.loc[i, "q"] = multipletests(pc.loc[i, "p"], method="fdr_bh")[1]
R["per_cancer"] = pc.to_dict("records")

# exportar bases auxiliares para figuras e verificacao em R
S[["patient_id", "anc", "ct", "T", "E", "L", "age_dx", "female", "stage"]].to_csv(OUT + "surv_dx.csv", index=False)
CH[["anc", "dx_to_seq_months"]].to_csv(OUT + "interval.csv", index=False)
Fc[["anc", "os_months_seq", "E", "cancer_type"] + m3].to_csv(OUT + "surv_seq.csv", index=False)
json.dump(R, open(OUT + "resultados.json", "w"), indent=1, default=float)
print(json.dumps({k: R[k] for k in ["n_chord", "n_surv", "deaths_surv", "n_50k_cc", "deaths_50k_cc", "kappa"]}, default=float))
