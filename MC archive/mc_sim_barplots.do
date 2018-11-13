* aa 6aug2018
clear
cap cd "C:\Users\Achim.Ahrens\Dropbox\StataLasso\Current"
cap cd "/home/achim/Dropbox/StataLasso/Current/"
adopath + "`c(pwd)'"

//global outputpath "/home/achim/Dropbox/StataLasso/Current/Do_files/mc6"
global outputpath "C:\Users\Achim.Ahrens\Dropbox\statalasso\current\do_files\mc7"

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

reshape long shat_ fneg_ fpos_ rmse_ rmse_oos_ l1norm_, ///
		i(delta theta sigma s design p obs iter ) string

rename _j estimator

gen ols = "-"
replace ols = "ols" if regexm(estimator,"ols")
replace estimator = subinstr(estimator,"_ols","",.)

rename rmse_ rmse_is

replace rmse_is=. if estimator=="oracle" & design==2 
replace rmse_oos=. if estimator=="oracle" & design==2 
replace l1norm=. if estimator=="oracle" & design==2 

gen estid = 1 if estimator=="aic"
replace estid = 2 if estimator=="aicc"
replace estid = 3 if estimator=="bic"
replace estid = 4 if estimator=="ebic"
replace estid = 5 if estimator=="cvlasso"
replace estid = 6 if estimator=="cvelastic"
replace estid = 7 if estimator=="cvridge"
replace estid = 8 if estimator=="rlasso"
replace estid = 9 if estimator=="sqrt"
replace estid = 10 if estimator=="oracle"


foreach s of numlist .5 1 2 3 5 {
foreach t of numlist .7 {
foreach d of numlist 1 2 {

if (`d'==1) {				

//local s = 3
//local t = .7
//local d = 1

graph hbar (asis) fpos if design==1 & (abs(sigma-`s')<10e-4) & (abs(theta-`t')<10e-4) ///
					& ols!="ols" , ///
					 over(estimator, label(nolabels) sort(estid)) ///
					nofill ///
					saving(fpos, replace) ///
					blabel(bar, format(%4.0g)) /// 
					ylab(, nolabels) ///
					ytitle("False pos.") nodraw 

graph hbar (asis) fneg if design==1 & ols!="ols" & ///
					& (abs(sigma-`s')<10e-4) ///
					& (abs(theta-`t')<10e-4), ///
					 over(estimator , sort(estid)) ///
					nofill ///
					saving(fneg, replace) ///
					blabel(bar, format(%4.0g)) /// 
					ylab(, nolabels) ///
					ytitle("False neg.") nodraw 
					
graph hbar (asis) shat if design==1 & ols!="ols" & ///
					& (abs(sigma-`s')<10e-4) ///
					& (abs(theta-`t')<10e-4), ///
					over(estimator, label(nolabels) sort(estid)) ///
					nofill  ///
					saving(shat, replace) ///
					blabel(bar, format(%4.0g)) /// 
					ylab(, nolabels) ///
					ytitle("Sparsity") nodraw 
}
*
		
graph hbar (asis) l1norm if design==`d' & (abs(sigma-`s')<10e-4) & (abs(theta-`t')<10e-4) ///
					, ///
					over(ols) over(estimator, sort(estid)) ///
					nofill ///
					saving(l1norm, replace) ///
					blabel(bar, format(%4.0g)) /// 
					ylab(, nolabels)  ///
					ytitle("L1-Norm") nodraw	
		
graph hbar (asis) rmse_is if design==`d' ///
					& (abs(sigma-`s')<10e-4) ///
					& (abs(theta-`t')<10e-4), ///
					over(ols) over(estimator, label(nolabels) sort(estid)) ///
					nofill ///
					saving(rmse_is, replace) ///
					blabel(bar, format(%4.0g)) /// 
					ylab(, nolabels) ///
					ytitle("RMSE (in-sample)") nodraw 
		
graph hbar (asis) rmse_oos if design==`d' ///
					& (abs(sigma-`s')<10e-4) ///
					& (abs(theta-`t')<10e-4), ///
					over(ols) over(estimator, label(nolabels) sort(estid)) ///
					nofill /// 
					saving(rmse_oos, replace) ///
					blabel(bar, format(%4.0g)) /// 
					ylab(, nolabels) ///
					ytitle("RMSE (out-of-sample)") nodraw 
					


	if (`d'==1) {
		graph combine fneg.gph fpos.gph  shat.gph  l1norm.gph rmse_is.gph rmse_oos.gph, rows(1) ysize(2)
		graph export "`c(pwd)'\Do_files\mc7-fig\design1_sigma`s'_theta`t'.pdf", as(pdf) replace 
	}
	else {
		graph combine  l1norm.gph rmse_is.gph rmse_oos.gph, rows(1) 
		graph export "`c(pwd)'\Do_files\mc7-fig\design`d'_sigma`s'_theta`t'.pdf", as(pdf) replace   
	}
}
}
}
