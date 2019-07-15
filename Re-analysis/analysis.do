cd "M:\projects\ieu2\_working\IEU2_P6_005\data\ccr5\Gib"

/*
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

import delimited "snps.raw", clear delim(" ")
rename fid id_ieu
drop iid-phenotype
foreach var of varlist aff* {
	drop `var'
}

foreach var of varlist rs* {
	capture replace `var' = "" if `var' == "NA"
	capture destring `var', replace force
}

*Merge to phenotypic IDs
merge 1:1 id_ieu using "linker2.dta", keep(3) nogen
save "snps.dta", replace
*/

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

gen age_final = age+fu

save "analysis2.dta", replace

***



*Analyses
use "analysis2.dta", clear

gen snp = ""
gen ea = ""
gen eaf = .
forvalues i = 1/6 {
	gen a`i'_beta = .
	gen a`i'_se = .
	gen a`i'_p = .
	gen a`i'_hr = .
	gen a`i'_l95 = .
	gen a`i'_u95 = .
	gen a`i'_n = .
	forvalues j = 0/2 {
		gen a`i'_n`j' = .
		gen a`i'_d`j' = .
	}
}

*Add in binary versions of the SNPs
foreach snp of varlist rs* {
	qui su `snp'
	if r(mean) > 1 {
		qui replace `snp' = 2-`snp'
		local snp_name = "`snp'"
		local snp_name = subinstr("`snp_name'","_a","_T",.)
		local snp_name = subinstr("`snp_name'","_t","_A",.)
		local snp_name = subinstr("`snp_name'","_g","_C",.)
		local snp_name = subinstr("`snp_name'","_c","_G",.)
		local snp_name = lower("`snp_name'")
		rename `snp' `snp_name'
		local snp = "`snp_name'"
	}
	qui recode `snp' (1 = 0) (2=1), gen("`snp'_bi")
}	
	
*Time = study time
stset fu, failure(dead)

local i = 1
foreach var of varlist rs* {
	dis "Analysing SNP `i'/1010, first analysis"
	qui replace snp = "`var'" in `i'
	qui replace ea = substr("`var'",-1,1) in `i'
	qui su `var'
	qui replace eaf = r(mean)/2 in `i'
	
	*Relateds included
	qui stcox `var' age sex pc* i.centre 
	qui replace a1_beta = _b[`var'] in `i'
	qui replace a1_se = _se[`var'] in `i'
	qui replace a1_n = e(N) in `i'
	
	*Relateds excluded
	qui stcox `var' age sex pc* i.centre if related == 0
	qui replace a2_beta = _b[`var'] in `i'
	qui replace a2_se = _se[`var'] in `i'
	qui replace a2_n = e(N) in `i'
	
	local i = `i' + 1
}

*Analyses with time = age of participant
stset age_final, failure(dead) enter(age)

local i = 1
foreach var of varlist rs* {
	dis "Analysing SNP `i'/1010, second/third analysis"
	*Relateds included
	qui stcox `var' age sex pc* i.centre 
	qui replace a3_beta = _b[`var'] in `i'
	qui replace a3_se = _se[`var'] in `i'
	qui replace a3_n = e(N) in `i'
	
	*Relateds excluded
	qui stcox `var' age sex pc* i.centre if related == 0
	qui replace a4_beta = _b[`var'] in `i'
	qui replace a4_se = _se[`var'] in `i'
	qui replace a4_n = e(N) in `i'
	
	*No covariables included in analysis
	*Relateds included
	qui stcox `var'
	qui replace a5_beta = _b[`var'] in `i'
	qui replace a5_se = _se[`var'] in `i'
	qui replace a5_n = e(N) in `i'
	
	*Relateds excluded
	qui stcox `var' if related == 0
	qui replace a6_beta = _b[`var'] in `i'
	qui replace a6_se = _se[`var'] in `i'
	qui replace a6_n = e(N) in `i'
	
	*Balance of alleles overall & in those that died
	forvalues k = 0/2 {
		qui count if `var' == `k'
		local rn`k' = r(N)
		
		qui count if `var' == `k' & dead == 1
		local rd`k' = r(N)
		
		foreach j in 1 3 5 {
			qui replace a`j'_n`k' = `rn`k'' in `i'
			qui replace a`j'_d`k' = `rd`k'' in `i'
		}
	
		qui count if `var' == `k' & related == 0
		local un`k' = r(N)
		qui count if `var' == `k' & dead == 1 & related == 0
		local ud`k' = r(N)
		
		foreach j in 2 4 6 {
			qui replace a`j'_n`k' = `un`k'' in `i'
			qui replace a`j'_d`k' = `ud`k'' in `i'
		}
	}
	
	local i = `i' + 1
}

***

keep snp - a6_d2
keep if snp != ""

forvalues j = 1/6 {
	qui replace a`j'_hr = exp(a`j'_beta)
	qui replace a`j'_l95 = exp(a`j'_beta-1.96*a`j'_se)
	qui replace a`j'_u95 = exp(a`j'_beta+1.96*a`j'_se)
	qui replace a`j'_p = 2*normal(-abs(a`j'_beta/a`j'_se))
}

foreach x in _a _g _c _t {
	qui replace snp = subinstr(snp,"`x'","",.)
}	

gen binary = 0
replace binary = 1 if strpos(snp,"_bi")>0
qui replace snp = subinstr(snp,"_bi","",.)
sort snp

save "results.dta", replace

*****

import delim "ukb_snp_chr3_v2.bim", clear
keep v2 v4
rename v2 snp
rename v4 pos
merge 1:m snp using "results.dta", nogen keep(3)

qui replace ea = "" if ea == "i"
qui gen x = 1 if ea == "a"
qui replace x = 2 if ea == "c"
qui replace x = 3 if ea == "g"
qui replace x = 4 if ea == "t"

bysort snp: egen y = max(x)
qui replace ea = "a" if y == 1
qui replace ea = "c" if y == 2
qui replace ea = "g" if y == 3
qui replace ea = "t" if y == 4

drop x y

order binary, a(snp)
sort snp binary

save "results_combined.dta", replace

***

use "results_combined.dta", clear

forvalues i = 1/6 {
	gen p`i' = -log10(a`i'_p)
}

twoway (scatter p4 pos if binary == 0, msize(tiny)) , ///
xtitle("Chromosome Position (bp)") ytitle("-log10(p)") xline(46414975)
graph export "Manhattan (best analysis).png", as(png) replace

twoway (scatter p4 pos if binary == 1, msize(tiny)) , ///
xtitle("Chromosome Position (bp)") ytitle("-log10(p)") xline(46414975)
graph export "Manhattan (best analysis, binary).png", as(png) replace

twoway (scatter p1 pos if binary == 0, msize(tiny)) (scatter p2 pos, msize(tiny)) (scatter p3 pos, msize(tiny)) /// 
(scatter p4 pos, msize(tiny)) (scatter p5 pos, msize(tiny)) (scatter p6 pos, msize(tiny)), ///
xtitle("Chromosome Position (bp)") ytitle("-log10(p)") xline(46414975) legend(off)
graph export "Manhattan (all analyses).png", as(png) replace

twoway (scatter p1 pos if binary == 1, msize(tiny)) (scatter p2 pos, msize(tiny)) (scatter p3 pos, msize(tiny)) /// 
(scatter p4 pos, msize(tiny)) (scatter p5 pos, msize(tiny)) (scatter p6 pos, msize(tiny)), ///
xtitle("Chromosome Position (bp)") ytitle("-log10(p)") xline(46414975) legend(off)
graph export "Manhattan (all analyses, binary).png", as(png) replace



****
*rs62625034 analysis
use "analysis2.dta", clear

*Seems like this could be a missing data problem
*Missingness is related to genotype batch
*Genotype is related to death

*Drop SNPs with only 1 allele
foreach var of varlist rs* {
	qui count if `var' != .
	local N = r(N)
	local drop = 0
	forvalues i = 0/2 {
		qui count if `var' == `i'
		if r(N) == `N' {
			local drop = 1
		}
	}
	if `drop' == 1 {
		drop `var'
	}
}
	
mi set flong
mi register imputed rs*
mi impute chained (ologit, augment) rs* = sex i.centre age n_22 pc* dead, add(5) rseed(100) dots 

save "Imputation.dta", replace


*Analyses with time = age of participant
mi stset age_final, failure(dead) enter(age)

mi estimate: stcox imp_rs62625034 age sex pc* i.centre if related == 0

sts graph, by(rs62625034_t) ylabel(0.7(0.05) 1) ci noorigin
graph export "rs62625034.png", as(png) replace

sts graph, by(rs62625034_tri) ylabel(0.7(0.05) 1) ci noorigin
graph export "rs62625034_tri.png", as(png) replace

