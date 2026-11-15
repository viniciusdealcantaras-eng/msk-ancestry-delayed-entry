* =====================================================================
* Ancestralidade genetica, tempo ate o sequenciamento e sobrevida
* MSK-IMPACT 50K + MSK-CHORD (cBioPortal). Stata 18.
* Nao altera o banco: le o CSV e nao salva nem sobrescreve nenhum arquivo de dados.
* Licenca dos dados CC BY-NC-ND 4.0: o banco nao e redistribuido.
* =====================================================================
version 18
clear all
set more off
* Run from the repository root after build_banco.py
import delimited "out/banco_ancestralidade_MSK.csv", clear varnames(1) case(preserve) bindquote(strict)

label define anc 1 "Europeia nao asquenaze" 2 "Europeia asquenaze" 3 "Miscigenada/outra" ///
    4 "Africana" 5 "Leste-asiatica" 6 "Sul-asiatica" 7 "Nativo-americana"
label values ancestry_cat anc
gen byte female = sex == "Female"
replace tmb = . if tmb < 0

* ---------------------------------------------------------------------
* 1. Descricao (MSK-IMPACT 50K, ancestralidade conhecida)
* ---------------------------------------------------------------------
preserve
keep if !missing(ancestry_cat)
tab ancestry_cat
tab ancestry_cat female, row
tabstat age_dx tmb fga, by(ancestry_cat) stat(p50 p25 p75) format(%9.2f)
gen byte metast = sample_type == "Metastasis"
tab ancestry_cat metast, row
tab ancestry_cat msi_high, row
tab ancestry_cat os_status, row
tab cancer_type ancestry_cat if inlist(cancer_type, "Non-Small Cell Lung Cancer", "Breast Cancer", ///
    "Colorectal Cancer", "Prostate Cancer", "Pancreatic Cancer"), row
restore

* ---------------------------------------------------------------------
* 2. Intervalo diagnostico-sequenciamento (pacientes no MSK-CHORD)
* ---------------------------------------------------------------------
preserve
keep if in_chord == 1 & chord_dx_valid == 1 & !missing(ancestry_cat) & !missing(dx_to_seq_months)
gen stage = chord_stage_dx
replace stage = "Unknown" if stage == ""
encode chord_cancer_type, gen(ct)
encode stage, gen(stg)
tabstat dx_to_seq_months, by(ancestry_cat) stat(n p50 p25 p75) format(%9.2f)
gen byte gt12 = dx_to_seq_months > 12
gen byte gt24 = dx_to_seq_months > 24
tab ancestry_cat gt12, row
tab ancestry_cat gt24, row
* mediana: bruta e ajustada
qreg dx_to_seq_months ib1.ancestry_cat, quantile(.5)
qreg dx_to_seq_months ib1.ancestry_cat i.ct i.stg age_dx female, quantile(.5)
* proporcao sequenciada > 12 e > 24 meses (Poisson com variancia robusta)
poisson gt12 ib1.ancestry_cat i.ct i.stg age_dx female, vce(robust) irr
poisson gt24 ib1.ancestry_cat i.ct i.stg age_dx female, vce(robust) irr

* ---------------------------------------------------------------------
* 3. Sobrevida desde o diagnostico: ingenua x entrada tardia
* ---------------------------------------------------------------------
keep if !missing(os_months_dx)
gen double T = max(os_months_dx, dx_to_seq_months + 0.0001)
gen double L = dx_to_seq_months
gen byte E = chord_os_status
gen byte stage_iv = stage == "IV"

* versao ingenua (tempo desde o diagnostico, sem entrada tardia)
stset T, failure(E == 1)
stcox ib1.ancestry_cat age_dx female, strata(ct)
stcox ib1.ancestry_cat age_dx female stage_iv, strata(ct)
sts list, by(ancestry_cat) at(60)
stci, by(ancestry_cat) median

* versao corrigida (entrada no sequenciamento)
stset T, failure(E == 1) enter(time L)
stcox ib1.ancestry_cat age_dx female, strata(ct)
stcox ib1.ancestry_cat age_dx female stage_iv, strata(ct)
estat phtest, detail
sts list, by(ancestry_cat) at(60)
stci, by(ancestry_cat) median

* sensibilidade: intervalos ate 10 e 20 anos
stcox ib1.ancestry_cat age_dx female if L <= 120, strata(ct)
stcox ib1.ancestry_cat age_dx female if L <= 240, strata(ct)
* sensibilidade: maior estadio registrado
gen byte stageH_iv = chord_stage_highest == "IV"
stcox ib1.ancestry_cat age_dx female stageH_iv, strata(ct)

* landmark a 6 e 12 meses
foreach t0 of numlist 6 12 {
    stset T if L <= `t0' & T > `t0', failure(E == 1) enter(time `t0')
    stcox ib1.ancestry_cat age_dx female, strata(ct)
}

* por tipo de cancer (africana e leste-asiatica x europeia nao asquenaze)
levelsof chord_cancer_type, local(cts)
foreach c of local cts {
    di as text _n "=== `c' ==="
    stset T if chord_cancer_type == "`c'" & inlist(ancestry_cat, 1, 4, 5), failure(E == 1)
    stcox ib1.ancestry_cat age_dx female
    stset T if chord_cancer_type == "`c'" & inlist(ancestry_cat, 1, 4, 5), failure(E == 1) enter(time L)
    stcox ib1.ancestry_cat age_dx female
}
restore

* ---------------------------------------------------------------------
* 4. Sobrevida apos o sequenciamento (MSK-IMPACT 50K), modelos sequenciais
* ---------------------------------------------------------------------
preserve
keep if !missing(ancestry_cat) & os_months_seq > 0 & !missing(os_months_seq, os_status)
gen byte metast = sample_type == "Metastasis"
gen double log2tmb = ln(tmb + 1) / ln(2)
gen byte msi_na = missing(msi_high)
gen byte msih = msi_high == 1
encode cancer_type, gen(ctype)
* mesma amostra (casos completos) nos tres modelos
egen byte nmiss = rowmiss(age_dx log2tmb fga wgd mut_TP53 mut_KRAS mut_PIK3CA del_CDKN2A amp_MYC)
keep if nmiss == 0
stset os_months_seq, failure(os_status == 1)
stcox ib1.ancestry_cat age_dx female, strata(ctype)
stcox ib1.ancestry_cat age_dx female metast, strata(ctype)
stcox ib1.ancestry_cat age_dx female metast log2tmb msih msi_na fga wgd ///
    mut_TP53 mut_KRAS mut_PIK3CA del_CDKN2A amp_MYC, strata(ctype)
estat phtest, detail
* sensibilidade: amostra primaria e metastatica
stcox ib1.ancestry_cat age_dx female log2tmb msih msi_na fga wgd ///
    mut_TP53 mut_KRAS mut_PIK3CA del_CDKN2A amp_MYC if sample_type == "Primary", strata(ctype)
stcox ib1.ancestry_cat age_dx female log2tmb msih msi_na fga wgd ///
    mut_TP53 mut_KRAS mut_PIK3CA del_CDKN2A amp_MYC if sample_type == "Metastasis", strata(ctype)
restore

* ---------------------------------------------------------------------
* 5. Concordancia ancestralidade genetica x raca autodeclarada
* ---------------------------------------------------------------------
preserve
keep if in_chord == 1 & !missing(ancestry_cat)
tab ancestry_cat chord_race
gen anc4 = .
replace anc4 = 1 if inlist(ancestry_cat, 1, 2)
replace anc4 = 2 if ancestry_cat == 4
replace anc4 = 3 if inlist(ancestry_cat, 5, 6)
replace anc4 = 4 if ancestry_cat == 7
gen race4 = .
replace race4 = 1 if chord_race == "White"
replace race4 = 2 if chord_race == "Black Or African American"
replace race4 = 3 if chord_race == "Asian-Far East/Indian Subcont"
replace race4 = 4 if chord_race == "Native American-Am Ind/Alaska"
kap anc4 race4 if !missing(anc4, race4)
restore
