* aa 
* date of last run started: 13sep2018
clear
cap cd "C:\Users\Achim.Ahrens\Dropbox\statalasso\Lassopack_MC"
cap cd "C:\Users\Achim.Ahrens\Lassopack_MC"

cap log close
log using "`c(pwd)'\out08\mc_log", text replace

cd

which lasso2
which cvlasso
which rlasso
which lassoutils

set seed 666 // the number of the beast

* program to calculate L1 error
cap program drop comparevec
program define comparevec , rclass
	syntax anything [,  NOCONStant]
	
	local truematname	: word 1 of `0'
	local estmatname 	: word 2 of `0'

	tempname truemat estmat 
	mat `truemat' 	= `truematname'
	mat `estmat' 	= `estmatname'

	// compare dim
	local ptrue = colsof(`truemat')
	local pest  = colsof(`estmat')
	assert `ptrue'==`pest'
	
	// don't include constant in L1-norm 
	// constant is always in last column
	if ("`noconstant'"=="") {
		local pcheck = `ptrue' -1
	}
	*
	
	local l1norm = 0
	forvalues i = 1(1)`pcheck' {
		local l1norm = `l1norm' + abs(el(`truemat',1,`i')-el(`estmat',1,`i'))
	}
	return scalar l1norm = `l1norm' // l1 norm
	
end

* program to calculate rmse
cap program drop getrmse
program define getrmse , rclass
	syntax varlist(min=2 max=2) , ///
			esample(varlist min=1 max=1) // estimation smpl
	
	local xb		: word 1 of `varlist'
	local xbhat 	: word 2 of `varlist'
	
	// calculate in-sample rmse
	tempvar pe 
	gen double `pe' = (`xb'-`xbhat')^2 if `esample'==1
	sum `pe', meanonly
	local rmse = sqrt(r(mean))
	
	// calculate out-of-sample rmse
	tempvar pe_oos
	gen double `pe_oos' = (`xb'-`xbhat')^2 if `esample'==0
	sum `pe_oos', meanonly
	local rmse_oos = sqrt(r(mean))
	
	return scalar rmse = `rmse' 			 
	return scalar rmse_oos = `rmse_oos'
end


cap program drop mysim
program define mysim, rclass
     syntax [, obs(integer 100) 	/// total number of obs
				tobs(integer 50)	/// size of estimation sample
				p(integer 30) 		/// ex constant
				design(integer 1) 	/// 0 or 1
				s(integer 20) 		/// sparsity (only relevant if design=1)
				delta(real 0.9)		/// determines rate at which beta decreases
									/// (only relevant if design!=1)
				theta(real 0.5) 	/// determines corr of X
				sigma(real 1) ]

	drop _all
	
	/*
	local obs = 500 
	local p = 100
	local s = 10
	local sigma =1
	local design 1
	local delta = 0.5
	local theta = 0.5
	*/
	
	set obs `obs'
	
	** create beta
	tempname BETA
	if (`design'==1) {
		// sparse beta
		mat `BETA' = J(1,`p'+1,0)
		forvalues j=1(1)`s' {
			mat `BETA'[1,`j'] = 1
		}
		*
		mat `BETA'[1,`p'+1] = 1 // constant
	}
	else {
		// dense beta
		// with beta(j)=(delta)^j and delta<1
		mat `BETA' = J(1,`p'+1,0)
		forvalues j=1(1)`p' {
			mat `BETA'[1,`j'] = (`delta')^(`j')
		}
		*
		mat `BETA'[1,`p'+1] = 1 // constant	
	}
	*
	
	// correlation matrix
	// theta controls the degree of correlation in X
	tempname CORR
	tempname cvec Cmat
	mata: `cvec'=J(1,`p',0)
	mata: for (i=1; i<=`p'; i++) `cvec'[1,i]=(`theta')^(i-1)
	mata: `Cmat'=Toeplitz((`cvec')',`cvec')
	mata: st_matrix("`CORR'",`Cmat')	
	*
	
	************* generate data ************************************************
	// gen e
	gen double e = rnormal()*`sigma'
	// create x
	drawnorm x1-x`p' , double corr(`CORR')
	// standardize x!
	//tempvar xsq
	//gen double `xsq'=.
	foreach var of varlist x1-x`p' {
		qui sum `var' //, meanonly 
		replace `var' = (`var'-r(mean))/r(sd)
	}
	// create xb
	gen double f = 1 // 1 = the intercept
	forvalues j=1(1)`p' {
		replace f = f + x`j'*el(`BETA',1,`j')
	}
	// generate y
	gen y = f + e 
	drop e 
	// estimation and validation sample
	tempvar esample vsample
	gen `esample' = 0
	replace `esample' = 1 if _n <= `tobs'
	gen `vsample' = 0 
	replace `vsample' = 1 if _n > `tobs'
	************* generate data end ********************************************
	
	// oracle estimator
	reg y x1-x`s' if `esample'
	tempvar xb_oracle
	predict double `xb_oracle', xb
	mat beta_oracle = e(b)
	// add zeros to create full/large beta vector
	mat zeros = J(1,`p'-`s',0) 
	mat beta_oracle = (beta_oracle[1,1..`s'],zeros,beta_oracle[1,`s'+1])
	local shat_oracle = `s'
	ds x1-x`s'
	local sel_oracle = r(varlist)
	
	// rlasso
	ereturn clear
	timer on 1
	rlasso y x1-x`p' if `esample'
	timer off 1
	local shat_rlasso = e(s)
	local sel_rlasso = e(selected)
	tempvar xb_rlasso xb_rlasso_ols
	predict double `xb_rlasso', xb
	predict double `xb_rlasso_ols', xb ols
	mat beta_rlasso 	= e(betaAll)
	mat beta_rlasso_ols	= e(betaAllOLS)
	
	// sqrt rlasso
	ereturn clear
	timer on 2
	rlasso y x1-x`p' if `esample', sqrt
	timer off 2
	local shat_sqrt = e(s)
	local sel_sqrt = e(selected)
	tempvar xb_sqrt xb_sqrt_ols
	predict double `xb_sqrt', xb
	predict double `xb_sqrt_ols', xb ols
	mat beta_sqrt 		= e(betaAll)
	mat beta_sqrt_ols 	= e(betaAllOLS)
	
	// rlasso with xdep
	ereturn clear
	timer on 3
	rlasso y x1-x`p' if `esample', xdep
	timer off 3
	local shat_rxdep = e(s)
	local sel_rxdep = e(selected)
	tempvar xb_rxdep xb_rxdep_ols
	predict double `xb_rxdep', xb
	predict double `xb_rxdep_ols', xb ols
	mat beta_rxdep 	= e(betaAll)
	mat beta_rxdep_ols	= e(betaAllOLS)
	
	// sqrt rlasso with xdep
	ereturn clear
	timer on 4
	rlasso y x1-x`p' if `esample', sqrt xdep
	timer off 4
	local shat_sxdep = e(s)
	local sel_sxdep = e(selected)
	tempvar xb_sxdep xb_sxdep_ols
	predict double `xb_sxdep', xb
	predict double `xb_sxdep_ols', xb ols
	mat beta_sxdep 		= e(betaAll)
	mat beta_sxdep_ols 	= e(betaAllOLS)

	// cvlasso
	ereturn clear
	timer on 5
	cvlasso y x1-x`p' if `esample', ///
				lopt postest nfolds(5) unitl nopath
	timer off 5
	local shat_cvlasso = e(s)
	local sel_cvlasso = e(selected)
	tempvar xb_cvlasso xb_cvlasso_ols
	predict double `xb_cvlasso', xb
	predict double `xb_cvlasso_ols', xb ols
	mat beta_cvlasso 		= e(betaAll)
	mat beta_cvlasso_ols	= e(betaAllOLS)

	// cv+ridge
	ereturn clear
	timer on 6
	cvlasso y x1-x`p' if `esample', ///
				lopt postest nfolds(5) unitl nopath alpha(0)
	timer off 6
	local shat_cvridge = e(s)
	local sel_cvridge = e(selected)
	tempvar xb_cvridge 
	predict double `xb_cvridge', xb
	mat beta_cvridge	= e(betaAll)
	
	// cv+elastic
	/*
	ereturn clear
	timer on 7
	cvlasso y x1-x`p' if `esample', ///
				lopt postest nfolds(5) unitl nopath alpha(0 0.1 0.5 0.9 1)
	timer off 7
	local shat_cvelastic = e(s)
	local sel_cvelastic = e(selected)
	tempvar xb_cvelastic xb_cvelastic_ols
	predict double `xb_cvelastic', xb
	predict double `xb_cvelastic_ols', xb ols
	mat beta_cvelastic	 	= e(betaAll)
	mat beta_cvelastic_ols 	= e(betaAllOLS)
	*/

	// lasso2 
	tempname lasso2est
	ereturn clear
	timer on 8
	lasso2 y x1-x`p' if `esample', unitl 
	timer off 8
	estimates store `lasso2est'
		// with aic 
		lasso2, lic(aic) postest
		local shat_aic = e(s)
		local sel_aic = e(selected)
		tempvar xb_aic xb_aic_ols
		predict double `xb_aic', xb
		predict double `xb_aic_ols', xb  ols
		mat beta_aic 		= e(betaAll)
		mat beta_aic_ols 	= e(betaAllOLS)
		// with bic
		estimates restore `lasso2est'
		lasso2, lic(bic) postest
		local shat_bic = e(s)
		local sel_bic = e(selected)
		tempvar xb_bic xb_bic_ols
		predict double `xb_bic', xb
		predict double `xb_bic_ols', xb ols
		mat beta_bic 		= e(betaAll)
		mat beta_bic_ols 	= e(betaAllOLS)
		// with ebic
		estimates restore `lasso2est'
		lasso2, lic(ebic) postest
		local shat_ebic = e(s)
		local sel_ebic = e(selected)
		tempvar xb_ebic xb_ebic_ols
		predict double `xb_ebic', xb
		predict double `xb_ebic_ols', xb ols
		mat beta_ebic 		= e(betaAll)
		mat beta_ebic_ols 	= e(betaAllOLS)
		// with aicc
		estimates restore `lasso2est'
		lasso2, lic(aicc) postest
		local shat_aicc = e(s)
		local sel_aicc = e(selected)
		tempvar xb_aicc xb_aicc_ols
		predict double `xb_aicc', xb
		predict double `xb_aicc_ols', xb ols
		mat beta_aicc 		= e(betaAll)
		mat beta_aicc_ols 	= e(betaAllOLS)
	
	ereturn clear // otherwise lasso2/cvlasso/results are returned
	
	** calculate false positive and false negative
	foreach mthd of newlist  ebic aic bic aicc /// 
							cvlasso /// cvelastic ///
							rlasso sqrt ///
							rxdep sxdep ///
							{
		if (`design'==1) {
			local fpos_vars: list sel_`mthd' - sel_oracle				
			local fneg_vars: list sel_oracle - sel_`mthd' 
			local fpos_`mthd'	: word count `fpos_vars'
			local fneg_`mthd'	: word count `fneg_vars'
		}
		else {
			local fpos_`mthd' = -1
			local fneg_`mthd' = -1
		}
	
	}
	*
	*** calculate L1norm and RMSE
	foreach mthd of newlist oracle ///
							rlasso rlasso_ols ///
							rxdep rxdep_ols ///
							sqrt sqrt_ols ///
							sxdep sxdep_ols /// 
							cvlasso cvlasso_ols ///
							cvridge ///
							/// cvelastic cvelastic_ols ///
							aic aic_ols /// 
							bic bic_ols ///
							aicc aicc_ols ///
							ebic ebic_ols {
		
		// compare true beta and beta-hat
		comparevec `BETA' beta_`mthd' 
		local l1norm_`mthd' 	= r(l1norm)
		return scalar l1norm_`mthd'  	=`l1norm_`mthd'' 
		
		// get rmse
		getrmse f `xb_`mthd'', esample(`esample')
		local rmse_`mthd' 		= r(rmse)
		local rmse_oos_`mthd' 	= r(rmse_oos)		
		return scalar rmse_`mthd' 		=`rmse_`mthd'' 
		return scalar rmse_oos_`mthd' 	=`rmse_oos_`mthd'' 
	}
	*
	*** return fpos, fneg and shat
	foreach mthd of newlist  ebic aic bic aicc /// 
							cvlasso /// cvelastic ///
							rlasso sqrt ///
							rxdep sxdep /// 
							{
	
		return scalar fpos_`mthd' = `fpos_`mthd''
		return scalar fneg_`mthd' = `fneg_`mthd''
		return scalar shat_`mthd' = `shat_`mthd''
	}
	*
	// return simulation parameters
	return scalar obs = `obs'
	return scalar p = `p'
	return scalar design = `design'
	return scalar sigma = `sigma'
	return scalar theta = `theta'
	if (`design'==1) {
		return scalar s = `s'
		return scalar delta = -1
	}
	else {
		return scalar s = -1
		return scalar delta = `delta'
	}
	//
	//for debugging
	//return matrix trueb = `BETA'
	//return matrix ols = beta_oracle
	//return matrix rlassools = beta_rlasso_ols
	*
end

***************** end of program definitions ***********************************

timer clear

local mc=1
foreach te in 0.9 {
foreach s in 0.5 1 2 3 5 {
//foreach s in 0.5 {
	
	simulate , reps(100): mysim, obs(400) tobs(200) p(100) design(1) sigma(`s') theta(`te')
	save "`c(pwd)'\out08\mc_`mc'.dta", replace
	local mc = `mc' + 1	
	
	simulate , reps(100): mysim, obs(400) tobs(200) p(100) design(2) sigma(`s') theta(`te')
	save "`c(pwd)'\out08\mc_`mc'.dta", replace
	local mc = `mc' + 1	
	
}
}
*

timer list

log close
