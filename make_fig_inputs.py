"""Build the small aggregated tables read by figuras.R from out/resultados.json.

Input : out/resultados.json (from analise.py)
Output: out/fig_comp.csv, out/fig_qreg.csv, out/fig_forest.csv, out/fig_percancer.csv
These files hold only aggregated estimates and counts, no patient-level data.
"""
import json
import pandas as pd

OUT = "out/"
R = json.load(open(OUT + "resultados.json"))
GR = ["ASJ-EUR", "ADMIX", "AFR", "EAS", "SAS", "NAM"]
ORD = ["nonASJ-EUR"] + GR
key = lambda g: "a_" + g.replace("-", "_")

# Figure 1: ancestry composition of the 15 most frequent cancer types
comp = pd.DataFrame(R["comp_by_cancer"]).T.reindex(columns=ORD).fillna(0).astype(int)
comp.to_csv(OUT + "fig_comp.csv")

# Figure 2B: crude and adjusted differences in the median interval
rows = []
for g in GR:
    c, a = R["qreg_crude"][key(g)], R["qreg"][key(g)]
    rows.append({"grp": g, "crude": c[0], "clo": c[1], "chi": c[2], "adj": a[0], "alo": a[1], "ahi": a[2]})
pd.DataFrame(rows).to_csv(OUT + "fig_qreg.csv", index=False)

# Figure 4: hazard ratios across models
models = [("MSK-CHORD, without delayed entry", R["cox_dx"]["naive"]),
          ("MSK-CHORD, with delayed entry", R["cox_dx"]["corrected"]),
          ("MSK-CHORD, delayed entry + stage", R["cox_dx"]["corrected_stage"]),
          ("MSK-IMPACT 50K, model 1", R["cox_seq"]["m1"]),
          ("MSK-IMPACT 50K, model 2", R["cox_seq"]["m2"]),
          ("MSK-IMPACT 50K, model 3", R["cox_seq"]["m3"])]
rows = []
for g in GR:
    for name, m in models:
        hr, lo, hi, _ = m["hr"][key(g)]
        rows.append({"grp": g, "model": name, "hr": hr, "lo": lo, "hi": hi})
pd.DataFrame(rows).to_csv(OUT + "fig_forest.csv", index=False)

# Figure 5 (per cancer type, without and with delayed entry)
pc = pd.DataFrame(R["per_cancer"])
pc[["ct", "grp", "model", "n_grp", "hr", "lo", "hi", "q"]].to_csv(OUT + "fig_percancer.csv", index=False)
print("figure inputs written to out/")
