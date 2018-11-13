drop _all

cap cd "/home/achim/Dropbox/StataLasso/Current/"

insheet using https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data, clear tab
global model lpsa lcavol lweight age lbph svi lcp gleason pgg45

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
		di "`i'"
		local bix = colnumb(`b',"`i'")
		local bAllix = colnumb(`bAll',"`i'")
		mat `bAll'[1,`bAllix']=`b'[1,`bix']
	}
	mat list `bAll'
	
	// return
	return matrix bAll = `bAll'
	return local sel `sel'
	return scalar p = `p'
end 

stepwise, pr(.1): reg $model

mat b = e(b)

makebAll , best("b") allvars(lcavol lweight age lbph svi lcp gleason pgg45 _cons)
