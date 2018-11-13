* aa 6aug2018
clear
cap cd "C:\Users\Achim.Ahrens\Dropbox\StataLasso\Current"
cap cd "/home/achim/Dropbox/StataLasso/Current/"
adopath + "`c(pwd)'"

global outputpath "/home/achim/Dropbox/StataLasso/Current/Do_files/mc6"
//global outputpath "C:\Users\Achim.Ahrens\Dropbox\statalasso\current\do_files\mc5"

** load all MC files
local files : dir "$outputpath" files "*.dta"

local fix = 1
foreach file in `files' { 
	if (`fix'==1) {
		use "$outputpath/`file'"
	}
	else {
		append using "$outputpath/`file'"
	}
	local fix = `fix' +1
}
*

** collapse
gen iter=1
//replace delta = . if design == 1
//replace s = . if design == 2
collapse (mean) *_* (sum) iter, by(delta theta sigma s design p obs)


foreach var of varlist shat* {
	local newvar = subinstr("`var'","shat_","",.)
	rename `var' `newvar'_shat
}
foreach var of varlist rmse_oos* {
	local newvar = subinstr("`var'","rmse_oos_","",.)
	rename `var' `newvar'_rmseoos
}
foreach var of varlist rmse_* {
	local newvar = subinstr("`var'","rmse_","",.)
	rename `var' `newvar'_rmse
}
foreach var of varlist l1norm* {
	local newvar = subinstr("`var'","l1norm_","",.)
	rename `var' `newvar'_l1norm
}
foreach var of varlist fpos* {
	local newvar = subinstr("`var'","fpos_","",.)
	rename `var' `newvar'_fpos
}
foreach var of varlist fneg* {
	local newvar = subinstr("`var'","fneg_","",.)
	rename `var' `newvar'_fneg
}
*

reshape long aic_ aicc_ cvlasso_ cvelastic_ bic_ ebic_ ///
		oracle_ rlasso_ sqrt_ cvridge_ ///
		aic_ols_ aicc_ols_ cvlasso_ols_ cvelastic_ols_ bic_ols_ ebic_ols_ ///
		rlasso_ols_ sqrt_ols_ cvridge_ols_, ///
		i(delta theta sigma s design p obs iter ) string

rename _j dim

drop cvridge_ols_
drop if dim =="fneg" & design==2 
drop if dim =="fpos" & design==2 


local plotvars cvlasso_ aic_ ebic_ bic_ aicc_  
//local plotvars cvlasso_ cvelastic_ cvridge_
//local plotvars cvlasso_ rlasso_ sqrt_

egen max = rowmax(`plotvars')

foreach var of varlist `plotvars' {

	replace `var'=`var'/max

}

replace dim = "False pos." if dim=="fpos"
replace dim = "False neg." if dim=="fneg"
replace dim = "Sparsity" if dim=="shat"
replace dim = "L1-norm" if dim=="l1norm"
replace dim = "RMSE (in-sample)" if dim=="rmse"
replace dim = "RMSE (out-of-sample)" if dim=="rmseoos"

keep if design==1 & sigma==3 

radar dim `plotvars'  , aspect(1) lw(thick thick thick thick thick)


