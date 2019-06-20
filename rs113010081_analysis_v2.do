cd ""

*Load in phenotypes

forvalues i = 1/51 {
	clear all
	set maxvar 15000
	
	use "y`i'.dta"
	*Excel code from "Phenotypes.xlsx" goes here
	keep n_eid  n_31_0_0  n_54_0_0  n_21003_0_0  n_22000_*  n_22009_0_*  n_40007_*  ts_53_0_0  ts_40000_* 

	save "z`i'.dta", replace	

}

use "z1.dta", clear
forvalues i = 2/51 {
	append using "z`i'.dta"
}

compress
rename n_eid id_phe

save "Phenotype_data.dta", replace

*Import GRS and combine with phenotypes

import delimited "grs.csv", clear
rename id id_ieu

*Merge to phenotypic IDs
merge 1:1 id_ieu using "linker2.dta", keep(3) nogen
save "snps.dta", replace

*Add in phenotypes
use "Phenotype_data.dta", clear
merge 1:1 id_phe using "snps.dta", keep(3) nogen
	
*Remove recommended drops
foreach var in _recommended _non_white_british {
	merge 1:1 id_ieu using "exclusions`var'.dta"
	keep if _merge == 1
	drop _merge
}

*Mark relateds (not drop)
merge 1:1 id_ieu using "exclusions_relateds.dta"
drop if _merge == 2
gen related = 1 if _merge == 3
replace related = 0 if related != 1
drop _merge

*And withdrawals
merge 1:1 id_phe using "withdrawals.dta"
drop if _merge == 3
drop _merge

*Drop missing IEU IDs
drop if id_ieu == ""

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

save "analysis.dta", replace

use "analysis.dta", clear

gen fu = (death_date-reg_date)/365.25
replace fu = (20500-reg_date)/365.25 if fu == .

*Deaths before registration can be dropped
drop if fu < 0

format *date %td

gen dead = 1 if death_date != .
replace dead = 0 if dead == .

rename n_31_0_0 sex
rename n_54_ centre
rename n_21003 age

forvalues i = 1/40 {
	rename n_22009_0_`i' pc`i'
}

*Flip rs113341849
replace rs113341849 = -rs113341849

foreach var of varlist rs* {
	rename `var' `var'_raw

	gen `var'_tri = 0 if `var'_raw < 0.5
	replace `var'_tri = 1 if `var'_raw >= 0.5 & `var'_raw < 1.5
	replace `var'_tri = 2 if `var'_raw >= 1.5

	gen `var' = 0 if `var'_raw < 1.5
	replace `var' = 1 if `var'_raw >= 1.5
}

save "analysis2.dta", replace

***



*Original analysis, time = study time
use "analysis2.dta", clear
stset fu, failure(dead)

*With relateds
foreach var of varlist rs113010081 rs113341849 {
	sts graph, by(`var') ylabel(0.95(0.05) 1) ci
	graph export "analysis 1 `var' [relateds].tif", as(tif) replace
	stcox `var' age sex pc* i.centre 
	local beta1_`var' = _b[`var']
	local se1_`var' = _se[`var']
	local n1_`var' = e(N)
}

*Without relateds
drop if related == 1
foreach var of varlist rs113010081 rs113341849 {
	sts graph, by(`var') ylabel(0.95(0.05) 1) ci
	graph export "analysis 1 `var' [unrelateds].tif", as(tif) replace
	stcox `var' age sex pc* i.centre 
	local beta2_`var' = _b[`var']
	local se2_`var' = _se[`var']
	local n2_`var' = e(N)
}

***

*Second analysis, time = age of participant
use "analysis2.dta", clear
gen age_final = age+fu
stset age_final, failure(dead) enter(age)

*With relateds
foreach var of varlist rs113010081 rs113341849 {
	sts graph, by(`var') ylabel(0.75(0.05) 1) ci noorigin
	graph export "analysis 2 `var' [relateds].tif", as(tif) replace
	stcox `var' age sex pc* i.centre 
	local beta3_`var' = _b[`var']
	local se3_`var' = _se[`var']
	local n3_`var' = e(N)
}

*Without relateds
drop if related == 1
foreach var of varlist rs113010081 rs113341849 {
	sts graph, by(`var') ylabel(0.75(0.05) 1) ci noorigin
	graph export "analysis 2 `var' [unrelateds].tif", as(tif) replace
	stcox `var' age sex pc* i.centre 
	local beta4_`var' = _b[`var']
	local se4_`var' = _se[`var']
	local n4_`var' = e(N)
}

***

*Third analysis, time = age of participant, remove covariables
use "analysis2.dta", clear
gen age_final = age+fu
stset age_final, failure(dead) enter(age)

*With relateds
foreach var of varlist rs113010081 rs113341849 {
	stcox `var'
	local beta5_`var' = _b[`var']
	local se5_`var' = _se[`var']
	local n5_`var' = e(N)

	forvalues i = 0/2 {
		qui count if `var'_tri == `i'
		local rn`i'_`var' = r(N)
		qui count if `var'_tri == `i' & dead == 1
		local rd`i'_`var' = r(N)
	}
}

*Without relateds
drop if related == 1
foreach var of varlist rs113010081 rs113341849 {
	stcox `var' 
	local beta6_`var' = _b[`var']
	local se6_`var' = _se[`var']
	local n6_`var' = e(N)
	
	forvalues i = 0/2 {
		qui count if `var'_tri == `i'
		local un`i'_`var' = r(N)
		qui count if `var'_tri == `i' & dead == 1
		local ud`i'_`var' = r(N)
	}
}

*Is poor imputation at fault?
*X = non-integer imputation
use "analysis2.dta", clear
gen age_final = age+fu
stset age_final, failure(dead) enter(age)

foreach var of varlist rs113010081 rs113341849 {
	gen `var'_x = `var'
	replace `var'_x  = . if `var'_raw != 0 & `var'_raw!=1 & `var'_raw!=2
	
	gen `var'_x_tri = `var'_tri
	replace `var'_x_tri  = . if `var'_raw != 0 & `var'_raw!=1 & `var'_raw!=2
}

*Without relateds
drop if related == 1
foreach var of varlist rs113010081_x rs113341849_x {
	sts graph, by(`var') ylabel(0.75(0.05) 1) ci noorigin
	graph export "analysis 3 `var' [unrelateds].tif", as(tif) replace
	stcox `var' age sex pc* i.centre 
	local beta_`var' = _b[`var']
	local se_`var' = _se[`var']
	local n_`var' = e(N)
	
	forvalues i = 0/2 {
		qui count if `var'_tri == `i'
		local un`i'_`var' = r(N)
		qui count if `var'_tri == `i' & dead == 1
		local ud`i'_`var' = r(N)
	}
}

*Y = remove poor imputation (>0.05 from an integer)
use "analysis2.dta", clear
gen age_final = age+fu
stset age_final, failure(dead) enter(age)

foreach var of varlist rs113010081 rs113341849 {
	gen `var'_y = `var'
	replace `var'_y  = . if (`var'_raw > 0.05 & `var'_raw < 0.95) | (`var'_raw > 1.05 & `var'_raw < 1.95)
	
	gen `var'_y_tri = `var'_tri
	replace `var'_y_tri  = . if (`var'_raw > 0.05 & `var'_raw < 0.95) | (`var'_raw > 1.05 & `var'_raw < 1.95)
}

*Without relateds
drop if related == 1
foreach var of varlist rs113010081_y rs113341849_y {
	sts graph, by(`var') ylabel(0.75(0.05) 1) ci noorigin
	graph export "analysis 4 `var' [unrelateds].tif", as(tif) replace
	stcox `var' age sex pc* i.centre 
	local beta_`var' = _b[`var']
	local se_`var' = _se[`var']
	local n_`var' = e(N)
	
	forvalues i = 0/2 {
		qui count if `var'_tri == `i'
		local un`i'_`var' = r(N)
		qui count if `var'_tri == `i' & dead == 1
		local ud`i'_`var' = r(N)
	}
}

***

*Results table
clear
set obs 16
gen snp = "rs113010081"
replace snp = "rs113341849" in 7/12
replace snp = "rs113010081_x" in 13
replace snp = "rs113341849_x" in 14
replace snp = "rs113010081_y" in 15
replace snp = "rs113341849_y" in 16
gen analysis = "Time = follow-up" in 1/2
replace analysis = "Time = age" in 3/4
replace analysis = "Time = age, no covariables" in 5/6
replace analysis = "Time = follow-up" in 7/8
replace analysis = "Time = age" in 9/10
replace analysis = "Time = age, no covariables" in 11/12
replace analysis = "Time = age, only integer imputation" in 13/14
replace analysis = "Time = age, only <=0.05 from integer imputation" in 15/16

gen unrelated = "No"
replace unrelated = "Yes" if mod(_n,2) == 0
replace unrelated = "Yes" in 13/16

gen beta = .
gen se = .
gen n = .

local j = 0
foreach snp in rs113010081 rs113341849 {
	forvalues i = 1/6 {
		local j = `j' + 1
		foreach var of varlist beta se n {
			qui replace `var' = ``var'`i'_`snp'' in `j'
		}
	}
}

local j = 12
foreach snp in rs113010081_x rs113341849_x rs113010081_y rs113341849_y {
	local j = `j' + 1
	foreach var of varlist beta se n {
		qui replace `var' = ``var'_`snp'' in `j'
	}
}


gen hr = exp(beta)
gen l95 = exp(beta-1.96*se)
gen u95 = exp(beta+1.96*se)
gen p = 2*normal(-abs(beta/se))

forvalues i = 0/2 {
	qui gen n_`i' = .
	qui gen d_`i' = .
}

foreach snp in rs113010081 rs113341849 {
	forvalues i = 0/2 {
		qui replace n_`i' = `rn`i'_`snp'' if snp == "`snp'"
		qui replace n_`i' = `un`i'_`snp'' if unrelated == "Yes" & snp == "`snp'"
		qui replace d_`i' = `rd`i'_`snp'' if snp == "`snp'"
		qui replace d_`i' = `ud`i'_`snp'' if unrelated == "Yes" & snp == "`snp'"
	}
}

foreach snp in rs113010081_x rs113341849_x rs113010081_y rs113341849_y  {
	forvalues i = 0/2 {
		qui replace n_`i' = `un`i'_`snp'' if snp == "`snp'"
		qui replace d_`i' = `ud`i'_`snp'' if snp == "`snp'"
	}
}

save "results.dta", replace
export delim "results.csv", replace

*****

use "results.dta", clear
