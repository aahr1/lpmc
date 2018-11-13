* aa 13sep2018
clear
cap cd "C:\Users\Achim.Ahrens\Dropbox\StataLasso\Lassopack_MC"
adopath + "`c(pwd)'"

global outputpath "`c(pwd)'\out11\"

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
di `fix'-1
*

//keep if abs(theta-.7)<0.01

** collapse
gen iter=1
replace delta = . if design == 1
replace s = . if design == 2
collapse (mean) *_* (sum) iter, by(delta theta sigma s design p obs alternate)


** reshape
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
		swise_ ///
		aic_ols_ aicc_ols_ cvlasso_ols_ cvelastic_ols_ bic_ols_ ebic_ols_ ///
		sxdep_ sxdep_ols_ rxdep_ rxdep_ols_ ///
		rlasso_ols_ sqrt_ols_ cvridge_ols_, ///
		i(delta theta sigma s design p obs iter alter) string

rename _j dim

** there is no cvridge_ols_
drop cvridge_ols_

sort design theta delta s dim sigma 
order design alternate theta delta s dim sigma p obs iter dim /// 
					aic_* bic_* aicc_* ebic_* ///
					rlasso_* sqrt_* cvlasso_* ///
					cvelastic_* ///
					rxdep_* sxdep_* ///
					oracle cvridg
	
** convert to string
foreach var of varlist aic_ bic_ ebic_ aicc_ /// 
						rlasso_ sqrt_ cvlasso_ ///
						cvelastic_ ///
						rxdep_ sxdep_ ///
						{

	cap gen `var'str = ""
	replace `var'str = strofreal(round(`var',0.001))+"\newline (\emph{"+strofreal(round(`var'ols_,0.001))+"})"  ///
				if dim == "l1norm" | dim =="rmse" | dim=="rmseoos"
	replace `var'str = strofreal(round(`var',0.01)) ///
					if dim != "l1norm" & dim !="rmse" & dim!="rmseoos"
}
*

*** omit fpos & fneg for design 2
drop if dim =="fneg" & design==2 
drop if dim =="fpos" & design==2 

** oracle, stepwise, ridge: l1norm, rmse and rmspe
gen oracle_str = strofreal(round(oracle_,0.001)) + "\newline (--)" ///
						if dim =="l1norm"|dim=="rmse"|dim=="rmseoos"
gen swise_str = strofreal(round(swise_,0.001)) + "\newline (--)" ///
						if dim =="l1norm"|dim=="rmse"|dim=="rmseoos"
gen ridge_str = strofreal(round(cvridge_,0.001)) + "\newline (--)" ///
							if dim =="l1norm"|dim=="rmse"|dim=="rmseoos"

** stepwise: sparsity, fpos & negative
replace swise_str = strofreal(round(swise_,0.01)) ///
					if dim != "l1norm" & dim !="rmse" & dim!="rmseoos"

keep design alter theta delta s dim sigma p obs iter *_str

replace ridge_str = "--" if dim != "l1norm" & dim !="rmse" & dim!="rmseoos"
replace oracle_str = "--" if dim != "l1norm" & dim !="rmse" & dim!="rmseoos"
//replace swise_str = "--" if dim != "l1norm" & dim !="rmse" & dim!="rmseoos"

gen dimix = 1 if dim == "shat"
replace dimix = 2 if dim == "fpos"
replace dimix = 3 if dim == "fneg"
replace dimix = 4 if dim == "l1norm"
replace dimix = 5 if dim == "rmse"
replace dimix = 6 if dim == "rmseoos"
replace dim = "Bias" if dim=="l1norm"
replace dim = "RMSE" if dim == "rmse"
replace dim = "RMSPE" if dim == "rmseoos"
replace dim = "$\hat{s}$" if dim == "shat"
replace dim = "False pos." if dim == "fpos"
replace dim = "False neg." if dim == "fneg"

** rename columns
rename aic_str AIC
rename bic_str BIC
rename aicc_str AICC
rename ebic_str EBIC
rename rlasso_str rlasso
rename sqrt_str	sqrt
rename cvlasso_str cvlasso
rename cvelastic_str cvelastic
rename ridge_str cvridge
rename oracle_str Oracle
rename sxdep_str sqrtxdep
rename rxdep_str rlassoxdep
rename swise_str Stepwise

sort design alternate theta dimix sigma

replace dim ="" if sigma!=.5

ds AIC AICC BIC EBIC cvlasso cvridge cvelastic rlasso* sqrt* Stepwise Ora 

texsave dim sigma AIC AICC BIC EBIC cvlasso cvridge cvelastic rlasso* sqrt* Stepwise Ora ///
					using "tex11\MC_results_alternate1_design1.tex" ///
					if design == 1 & alter==1, ///
					replace hlines(5(5)25) size(tiny) nofix
texify "tex11\MC_results_alternate1_design1.tex"
				
texsave dim sigma AIC AICC BIC EBIC cvlasso cvridge cvelastic rlasso* sqrt* Stepwise   ///
					using "tex11\MC_results_alternate1_design2.tex" ///
						if design == 2 & alter==1, ///
						replace hlines(5(5)15) size(tiny) nofix
texify "tex11\MC_results_alternate1_design2.tex"
						
texsave dim sigma AIC AICC BIC EBIC cvlasso cvridge cvelastic rlasso* sqrt* Stepwise Ora ///
					using "tex11\MC_results_alternate-1_design1.tex" ///
					if design == 1 & alter==-1 & 0.01>abs(theta-.9), ///
					replace hlines(5(5)25) size(tiny) nofix
texify "tex11\MC_results_alternate-1_design1.tex"
				
texsave dim sigma AIC AICC BIC EBIC cvlasso cvridge cvelastic rlasso* sqrt* Stepwise   ///
					using "tex11\MC_results_alternate-1_design2.tex" ///
						if design == 2 & alter==-1 & 0.01>abs(theta-.9), ///
						replace hlines(5(5)15) size(tiny) nofix
texify "tex11\MC_results_alternate-1_design2.tex"
					
texsave dim sigma AIC AICC BIC EBIC cvlasso cvridge cvelastic rlasso* sqrt* Stepwise Ora ///
					using "tex11\MC_results_alternate-1_design1_smalltheta.tex" ///
					if design == 1 & alter==-1 & 0.01>abs(theta-.3), ///
					replace hlines(5(5)25) size(tiny) nofix
texify "tex11\MC_results_alternate-1_design1_smalltheta.tex"
				
texsave dim sigma AIC AICC BIC EBIC cvlasso cvridge cvelastic rlasso* sqrt* Stepwise   ///
					using "tex11\MC_results_alternate-1_design2_smalltheta.tex" ///
						if design == 2 & alter==-1 & 0.01>abs(theta-.3), ///
						replace hlines(5(5)15) size(tiny) nofix
texify "tex11\MC_results_alternate-1_design2_smalltheta.tex"

texsave dim sigma AIC AICC BIC EBIC cvlasso cvridge cvelastic rlasso* sqrt* Stepwise Ora ///
					using "tex11\MC_results_alternate-1_design1_zerotheta.tex" ///
					if design == 1 & alter==-1 & 0.01>abs(theta-0), ///
					replace hlines(5(5)25) size(tiny) nofix
texify "tex11\MC_results_alternate-1_design1_zerotheta.tex"
				
texsave dim sigma AIC AICC BIC EBIC cvlasso cvridge cvelastic rlasso* sqrt* Stepwise   ///
					using "tex11\MC_results_alternate-1_design2_zerotheta.tex" ///
						if design == 2 & alter==-1 & 0.01>abs(theta-0), ///
						replace hlines(5(5)15) size(tiny) nofix
texify "tex11\MC_results_alternate-1_design2_zerotheta.tex"
						
