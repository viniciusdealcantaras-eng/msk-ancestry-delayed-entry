# Genetic ancestry, time to sequencing and survival in MSK-IMPACT 50K and MSK-CHORD: reproducibility code

Code to reproduce the numbers, tables and figures of the manuscript
*"Genetic ancestry, time to tumour sequencing and survival in clinical sequencing cohorts: an integrated analysis of MSK-IMPACT 50K and MSK-CHORD with delayed-entry correction"*.

The analysis links the two public studies by the shared de-identified patient
identifier, compares the diagnosis-to-sequencing interval across genetic
ancestry groups and estimates survival by ancestry with and without delayed
entry (left truncation at the date of sequencing).

## Data and license

Public cBioPortal studies `msk_impact_50k_2026` and `msk_chord_2024`, both
distributed by Memorial Sloan Kettering Cancer Center under the
[Creative Commons BY-NC-ND 4.0 license](https://creativecommons.org/licenses/by-nc-nd/4.0/).
For this reason **no raw or derived patient-level data are distributed in this
repository**. The integrated dataset is rebuilt locally by `build_banco.py`.
Please cite the original studies when using the data. No attempt was made to
re-identify any individual.

Place the raw files downloaded from cBioPortal in:

- `data/msk_impact_50k_2026/`: `data_clinical_patient.txt`,
  `data_clinical_sample_50k.tsv`, `mutations_driver_genes_50k.tsv`,
  `data_cna.txt`, `structural_variants_50k.tsv`
- `data/msk_chord_2024/`: `data_clinical_patient.txt`,
  `data_clinical_sample.txt`, `data_timeline_diagnosis.txt`

`data/` and `out/` are git-ignored.

## Environment

- Python >= 3.10 with the packages in `requirements.txt`
- R >= 4.0 with the `survival` package (bundled with base R)
- Stata 18 (optional, independent reproduction)

## Pipeline (run in order, from this folder)

1. `python build_banco.py`: builds the integrated patient-level dataset
   `out/banco_ancestralidade_MSK.csv`.
2. `python analise.py`: descriptive analysis, interval to sequencing (quantile
   and robust Poisson regression), survival from diagnosis without and with
   delayed entry, sequential models after sequencing, ancestry versus
   self-reported race, sensitivity analyses, conditional Kendall's tau and
   per-cancer models with FDR correction. Writes `out/resultados.json`.
3. `python make_fig_inputs.py`: aggregated tables for the figures.
4. `Rscript figuras.R`: Figures 1 to 5 (TIFF 600 dpi, PNG and PDF) in `out/figs/`.
5. `Rscript verificar.R`: independent check of the main Cox models and
   Kaplan-Meier medians with the `survival` package.
6. `python analise_landmark_tau.py`: landmark analyses at 6 and 12 months with
   delayed entry (every patient alive at the landmark, entering at the landmark
   or at sequencing, whichever is later) next to the version restricted to
   patients already sequenced, and the conditional Kendall's tau on all patients
   by ancestry, by cancer type and by ancestry within cancer type, with standard
   errors from the U-statistic projection (Supplementary Tables S3, S5 and S6).
   Writes `out/resultados_landmark_tau.json`.
7. `analise_ancestralidade.do`: independent reproduction in Stata, including the
   landmark analyses with delayed entry. It reads the CSV and never saves or
   overwrites any data file.

## Delayed entry in one line

Each patient enters the risk set at sequencing (entry time = diagnosis-to-sequencing
interval, `L`), with time measured from diagnosis (`T`):

- Python: `CoxPHFitter().fit(df, duration_col="T", event_col="E", entry_col="L", strata=["ct"])`
- R: `coxph(Surv(L, T, E) ~ anc + age_dx + female + strata(ct))`
- Stata: `stset T, failure(E == 1) enter(time L)`
