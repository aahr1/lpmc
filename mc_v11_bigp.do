* aa 
* date of last run started: 08oct2018
clear
cap cd "C:\Users\Achim.Ahrens\Dropbox\statalasso\Lassopack_MC"
cap cd "C:\Users\Achim.Ahrens\Lassopack_MC"
cap cd "Z:\home\achim\Dropbox\Statalasso\Lassopack_MC"

local thep = 220
local thedesign = 2
local thealternate = 1

cap log close
log using "`c(pwd)'\out11\mc_log_p`thep'_alter`thealternate'_design`thedesign'_bigp", text replace

cd

which lasso2
which cvlasso
which rlasso
which lassoutils

* p is equal to:
di `thep'
* design:
di `thedesign'
* alternating beta
di `thealternate'

set seed 666 // the number of the beast

* program to calculate L1 error
cap program drop comparevec
program define comparevec , rclass
	syntax ,  TRUEMATname(string) ESTMATname(string) 
	
	tempname truemat estmat 
	mat `truemat' 	= `truematname'
	mat `estmat' 	= `estmatname'

	// compare dim
	local ptrue = colsof(`truemat')
	local pest  = colsof(`estmat')
	assert `ptrue'==`pest'
	
	// don't include constant in L1-norm 
	// constant is always in last column of estimated beta and true beta
	local pcheck = `pest' - 1

	
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

cap program drop makebAll
program define makebAll, rclass

	syntax , best(string) allvars(string)
	
	// number of vars (incl constant)
	local p : word count `allvars'

	tempname b bAll	
	mat `b' = `best'
	
	// selected vars
	local sel : colnames `b'
	
	// create empty full beta vector
	mat `bAll' = J(1,`p',0)
	mat colnames `bAll' = `allvars'
	
	// fill long beta vector
	foreach i of local sel {
		//di "`i'"
		local bix = colnumb(`b',"`i'")
		local bAllix = colnumb(`bAll',"`i'")
		mat `bAll'[1,`bAllix']=`b'[1,`bix']
	}
	//mat list `bAll'
	
	// return
	return matrix bAll = `bAll'
	return local sel `sel'
	//return scalar p = `p'
end 

cap program drop mysim
program define mysim, rclass
     syntax [, obs(integer 100) 	/// total number of obs
				tobs(integer 50)	/// size of estimation sample
				p(integer 30) 		/// ex constant
				design(integer 1) 	/// 1 or 2
				s(integer 20) 		/// sparsity (only relevant if design=1)
				delta(real 0.9)		/// determines rate at which beta decreases
									/// (only relevant if design!=1)
				theta(real 0.9) 	/// determines corr of X
				sigma(real 1)		///
				ALTernate(integer 1) ]

	drop _all
	
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
	// create xb
	gen double f = el(`BETA',1,`p'+1) // = the intercept (in last column)
	forvalues j=1(1)`p' {
		replace f = f + x`j'*el(`BETA',1,`j')*(`alternate')^(`j')
	}
	// generate y
	gen y = f + e 
	drop e 
	drop f
	// estimation and validation sample
	tempvar esample vsample
	gen `esample' = 0
	replace `esample' = 1 if _n <= `tobs'
	gen `vsample' = 0 
	replace `vsample' = 1 if _n > `tobs'
	************* generate data end ********************************************
	
	// oracle estimator
	if (`design'==1) {
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
		local oraclelocal oracle
	}
	else {
		local oraclelocal
	}
	*
	
	// stepwise
	if (`p'<`tobs') { // only if p<n
		 timer on 11
		stepwise, pr(.1): reg y x1-x`p' if `esample'
		 timer off 11
		local shat_swise = e(df_m)
		tempvar xb_swise
		predict double `xb_swise', xb
		ds x1-x`p'
		local allreg = r(varlist) // list of regressors
		mat beta_swise0 = e(b) // "reduced" beta without zeros
		makebAll , best(beta_swise0) allvars("`allreg' _cons") // convert to long beta with zeros
		mat beta_swise = r(bAll)
		local sel_swise = r(sel)
		local swiselocal swise
	}
	else {
		local swiselocal
	}
	*
	
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

	/*// cvlasso
	ereturn clear
	 timer on 5
	cvlasso y x1-x`p' if `esample', ///
				lopt postest nfolds(5) nopath alpha(1)
	 timer off 5
	local shat_cvlasso = e(s)
	local sel_cvlasso = e(selected)
	tempvar xb_cvlasso xb_cvlasso_ols
	predict double `xb_cvlasso', xb
	predict double `xb_cvlasso_ols', xb ols
	mat beta_cvlasso 		= e(betaAll)
	mat beta_cvlasso_ols	= e(betaAllOLS)
	*/
	
	/*// cv+ridge
	ereturn clear
	 timer on 6
	cvlasso y x1-x`p' if `esample', ///
				lopt postest nfolds(5) nopath alpha(0)
	 timer off 6
	local shat_cvridge = e(s)
	local sel_cvridge = e(selected)
	tempvar xb_cvridge 
	predict double `xb_cvridge', xb
	mat beta_cvridge	= e(betaAll)
	local cvridgelocal cvridge
	*/
	
	// cv+elastic
	/*
	ereturn clear
	timer on 7
	cvlasso y x1-x`p' if `esample', ///
				lopt postest nfolds(5) nopath alpha(0 0.1 0.5 0.9 1)
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
	lasso2 y x1-x`p' if `esample', 
	 timer off 8
	estimates store `lasso2est'
	
		// with aic 
		 timer on 9
		lasso2, lic(aic) postest
		 timer off 9
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
	foreach mthd of newlist  `swiselocal' ///
							ebic aic bic aicc /// 
							/// cvlasso /// cvelastic ///
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
			// placeholder
			local fpos_`mthd' = -1
			local fneg_`mthd' = -1
		}
	
	}
	*
	*** calculate L1norm and RMSE
	foreach mthd of newlist `oraclelocal' ///
							`swiselocal' /// 
							rlasso rlasso_ols ///
							rxdep rxdep_ols ///
							sqrt sqrt_ols ///
							sxdep sxdep_ols /// 
							/// cvlasso cvlasso_ols ///
							///`cvridgelocal' ///
							/// cvelastic cvelastic_ols ///
							aic aic_ols /// 
							bic bic_ols ///
							aicc aicc_ols ///
							ebic ebic_ols {
		
		// compare true beta and beta-hat
		di "`mthd'"
		comparevec, truemat(`BETA') estmat(beta_`mthd')
		local l1norm_`mthd' 	= r(l1norm)
		return scalar l1norm_`mthd'  	=`l1norm_`mthd'' 
		
		// get rmse
		getrmse y `xb_`mthd'', esample(`esample')
		local rmse_`mthd' 		= r(rmse)
		local rmse_oos_`mthd' 	= r(rmse_oos)		
		return scalar rmse_`mthd' 		=`rmse_`mthd'' 
		return scalar rmse_oos_`mthd' 	=`rmse_oos_`mthd'' 
	}
	*
	*** return fpos, fneg and shat
	foreach mthd of newlist `swiselocal' ///
							ebic aic bic aicc /// 
							///cvlasso /// cvelastic ///
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
	return scalar alternate = `alternate'
	if (`design'==1) {
		return scalar s = `s'
		return scalar delta = -1
	}
	else {
		return scalar s = -1
		return scalar delta = `delta'
	}
end

***************** end of program definitions ***********************************

timer clear

local mc=1
foreach s in 0.5 1 2 3 5 {
	
	simulate , reps(300): mysim, obs(400) tobs(200) p(`thep') design(`thedesign') sigma(`s') alter(`thealternate')  
	save "`c(pwd)'\out11\mc_p`thep'_alter`thealternate'_design`thedesign'_`mc'_bigp.dta", replace
	local mc = `mc' + 1	
	
}
*

timer list

log close
