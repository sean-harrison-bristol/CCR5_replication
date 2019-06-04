
*Load in phenotypes
cd ""
use "Phenotype_data_all.dta"

*UK Biobank variables for this replication
keep n_eid ts_53_0_0  n_31_0_0  n_54_0_0  n_21003_0_0  n_22009_0_*  n_40007_*  ts_40000_* 

save "Phenotype_data.dta", replace

*Import rs113010081 variable and combine with phenotypes

import delimited "rs.csv", clear
merge 1:1 id_ieu using "linker.dta", keep(3) nogen
save "rs.dta", replace

*Add in phenotypes
use "Phenotype_data.dta", clear
merge 1:1 id_phe using "rs.dta", keep(3) nogen
	
*Remove recommended drops
foreach var in _recommended _highly_related _relateds _non_white_british {
	merge 1:1 id_ieu using "exclusions`var'.dta"
	keep if _merge == 1
	drop _merge
}

*And withdrawals
merge 1:1 id_phe using "withdrawals.dta"
drop if _merge == 3
drop _merge

*Merge age and dates at death
foreach var in n_40007 ts_40000 {
	forvalues i = 1/2 {
		qui replace `var'_0_0 = `var'_`i' if `var'_`i' != .
		drop `var'_`i'
	}
}

rename n_40007 death_age
rename ts_40000 death_date
rename ts_53 reg_date

gen fu = (death_date-reg_date)/365.25 //follow-up
replace fu = (20500-reg_date)/365.25 if fu == .

*Deaths before registration can be dropped (n=2)
drop if fu < 0

gen dead = 1 if death_date != .
replace dead = 0 if dead == .

rename n_31_0_0 sex
rename n_54_ centre
rename n_21003 age

forvalues i = 1/40 {
	rename n_22009_0_`i' pc`i'
}

*Note: rs113010081 was imputed if missing, so some values are not 0,1,2
*rs113010081_tri = 0, 1, 2
gen rs113010081_tri = 0 if rs113010081 < 0.5
replace rs113010081_tri = 1 if rs113010081 >= 0.5 & rs113010081 < 1.5
replace rs113010081_tri = 2 if rs113010081 >= 1.5

*rs113010081_bi = 0, 1 [homozygous mutant allele]
gen rs113010081_bi = 0 if rs113010081 < 1.5
replace rs113010081_bi = 1 if rs113010081 >= 1.5

*Analysis
stset fu dead
sts graph, by(rs113010081_tri) ylabel(0.95(0.05) 1) ci
stcox rs113010081 age sex pc* i.centre 
stcox rs113010081_bi age sex pc* i.centre 
stcox i.rs113010081_tri age sex pc* i.centre 
