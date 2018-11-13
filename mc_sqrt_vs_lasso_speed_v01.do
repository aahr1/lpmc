cap program drop mydata
program define mydata, rclass
     syntax [, obs(integer 100) 	/// total number of obs
				tobs(integer 50)	/// size of estimation sample
				p(integer 30) 		/// ex constant
				design(integer 1) 	/// 1 or 2
				s(integer 20) 		/// sparsity (only relevant if design=1)
				delta(real 0.9)		/// determines rate at which beta decreases
									/// (only relevant if design!=1)
				theta(real 0.9) 	/// determines corr of X
				sigma(real 1) ]

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
end

****************************************************************
******************* TIMINGS ************************************
****************************************************************

// Default theta=0.9. Sqrt lasso 5x slower.
drop _all
timer clear
set seed 123
forvalues i=1/100 {
	qui mydata
	timer on 1
	qui lasso2 y x*, l(66)
	timer off 1
	timer on 2
	qui lasso2 y x*, sqrt l(37)
	timer off 2
}
timer list


// theta=0. Sqrt lasso 50% faster.
drop _all
timer clear
set seed 123
forvalues i=1/100 {
	qui mydata, theta(0)
	timer on 1
	qui lasso2 y x*, l(66)
	timer off 1
	timer on 2
	qui lasso2 y x*, sqrt l(37)
	timer off 2
}
timer list


// theta=0.4. Sqrt lasso and lasso similar speed.
drop _all
timer clear
set seed 123
forvalues i=1/100 {
	qui mydata, theta(0.4)
	timer on 1
	qui lasso2 y x*, l(66)
	timer off 1
	timer on 2
	qui lasso2 y x*, sqrt l(37)
	timer off 2
}
timer list

