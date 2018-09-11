***********************************************************
*** Utility functions for simulation
**
* Developers: Simo Goshev, Zitong Liu
* 
*
*

********************************************************************************
**** Compute correlation matrix of actual data

capture program drop dataCorrMat
program define dataCorrMat

	syntax namelist, saving(string)  // string here is a path
	
	qui {
		noi di in y "Creating correlation matrix... "
		local all_items ""
		foreach name of local namelist {
			unab scale_items: `name'*
			local all_items "`all_items' `scale_items'"
		}
		
		* scaleclean `varlist' // we will include this utility later; 
		
		** Compute Spearman's correlation
		spearman `all_items', stats(rho)
		
		** Save the matrix
		mat dcorr = r(Rho)
			
		*** fudge the correlation matrix (change missing to 0)
		mata: fixed = editmissing(st_matrix("dcorr"), 0)
		mata: _diag(fixed, 1)
		mata: st_matrix("dcorr", fixed)
		
		** Save the correlation matrix
		preserve
		clear
		svmat dcorr
		save `saving', replace
		restore
		
		noi di in y "Done."
	}
end



********************************************************************************
*** Set up a theoretical correlation matrix


*** Within item correlations
capture program drop genwithin
program define genwithin, rclass
	args dim cor
	mat mywithin = J(`dim', `dim', `cor')
	forval i = 1/`dim' {
		forval j = 1 /`dim' {
			if (`i' == `j') {
				mat mywithin[`i', `j'] = 1
			}
		}
	}
	return matrix within = mywithin
end

*** Build the complete correlation matrix
capture program drop gencombined
program define gencombined, rclass
	args dimw corw dimb corb
	mat mybetween = J(`dimb', `dimb', `corb')
	genwithin `dimw' `corw'
	mat within = r(within)
	forval i = 1(`dimw')`dimb' {
		forval j = 1(`dimw')`dimb' {
			if (`i'==`j') {
				mat mybetween[`i', `j'] = within
			}
		}
	}
	return matrix combined = mybetween
end

*** Build a vector of means
capture program drop genmeans
program define genmeans, rclass

	args nmeans

	mat mymeans = J(1, `nmeans', 0)
	return matrix means = mymeans
end


*** TODO: Build distributions of multi-category variables
*** TODO: Build between-scale-item correlations



********************************************************************************
*** Simulation programme



**** Modification by Zitong:
**** It is better to tell Stata the varlist we want to simulate on for the sake of reading dta files in. 
**** I add this as an argument called "storedvars"
**** Use "storedvars" to create "storedscales"(input argument of catsimcore)


capture program drop ifeats
program define ifeats, rclass
	
	syntax namelist, Nobs(numlist) /// 
	       [NWitems(integer 3) NTitems(integer 12) propmiss(real 0) /// 
		    MBlock(integer 0) SIMCorr(integer 0) SIMMarginals(integer 0) ///
			CORRMatrix(string) MARGinals(string)]
	
	*** Total number of items
	local nTotalItems = 0
	foreach scale of local namelist {
		unab `scale'_items: `scale'*
		local nItems: list sizeof `scale'_items
		local nTotalItems = `nTotalItems' + `nItems'
	}
	
	if `simmarginals' == 1 {
	**** Todo. 
	}
	else {
		* capture unab `storedvars':  `storedvars'
		* qui scaleclean `storedvars' 
	
		mat simmat = J(1, `: word count `nobs'', .)
		local rows = 1
		local cols = 1
		foreach lobs of local nobs {
		* di "`nTotalItems'"
			ifeatsCore `lobs' 1 `cols' `rows' `nwitems' 1 `propmiss' `mblock' `simcorr' `simmarginals' `corrmatrix' `marginals' "`namelist'" "`nTotalItems'"
		}
	}
	mat list simmat
	return matrix sims = simmat
	
end



capture program drop ifeatsCore
program define ifeatsCore
 
	args nobs corb cols rows nwitems corw propmiss mblock simcorr simmarginals corrmatrix marginals namelist nTotalItems
	
/*
	local nobs    = `lobs'     // data has 64
	local nwitems = 3         // number of time periods
	local ntitems = 33        // time periods * number of indicators
	local corw    = 0.9       // within item correlation (over time)
	local corb    = `lcorb'   // between item correlation
	local missing = 0.6       // proportion of missing (based on obs)
	local mblock  = 1         // block missing (boolean)
	local simvcov = 1         // simulate vcov or use data vcov (boolean)

*/

	clear	
	qui set obs `nobs'

	
	**** Need change the display command.
	if `simcorr' == 1 {
		noi di _n in y "************************************************************"
		noi di in y    "* Sample size        : " `nobs'
		noi di in y    "* Within correlation : " `corw'
		noi di in y    "* Between correlation: " `corb'
		noi di in y    "************************************************************"
	}
	else {
		noi di _n in y "************************************************************"
		noi di in y    "* Sample size        : " `nobs'
		noi di in y    "* Using empirical correlation matrix "		
		noi di in y	   "* Using empirical marginal distributions"
		noi di in y "************************************************************"
	}
	
	


	********************************************************************************
	*** Define means and correlation
	
	genmeans "`nTotalItems'"
	mat means = r(means)	
	
	noi di "Simulating data..."
	
	if `simcorr' == 1 {
	
		/* TODO
		noi di "(Simulated correlation matrix)"
		gencombined `nwitems' `corw' `ntitems' `corb'
		mat rho = r(combined)
		*/
	}
	else {
		preserve
		use "`corrmatrix'", clear
		mkmat _all, matrix(rho)
		restore
		noi di "(Correlation matrix of the data loaded.)"
	}
	
	
	if `simmarginals' == 1 {
		**** TODO. Create cut points from simulated marginal distributions 
	}
	else { // Use observed marginals 

		*** Simulate random normal data using the (empirical) distributions and corr mat
		drawnorm myvar1-myvar`nTotalItems', means(means) corr(rho) n(`=_N')
	
		*** Read in scale_cdist's as matrices
		local i = 1
		qui foreach scale of local namelist {
			*** Load the marginals, create matrices
			preserve
			
			use "`marginals'/`scale'_cdist", clear
			order _all, alpha
			order cats, first
			
			unab items: _all
			local `scale'_items = subinstr("`items'", "cats", "", .) // stores names
			local nItems: list sizeof `scale'_items 
			* di "`nItems'"					
			
			qui levelsof cats, local(`scale'_cats)     // categories
			local nCats: list sizeof `scale'_cats
			drop cats
			drop if _n == _N
			mata: `scale'_cd = st_data(.,.)
			restore 
	
			*** Generate the items following the marginals
			local matcol = 1
			foreach item of local `scale'_items {
				* di "`item'"
				mata: mat_ph = `scale'_cd[1...,`matcol'..`matcol']
				mata: _transpose(mat_ph)
				mata: st_local("pctiles", invtokens(strofreal(mat_ph,"%12.10g")))
				* di "`pctiles'"
				egen `item' = xtile(myvar`i'), percentiles(`pctiles') // this requires egenmore
				replace `item' = `item' - 1
				drop myvar`i'
				local ++matcol
				local ++i
			}
		}
	}
	* sum
	
	********************************************************************************
	** Introduce missingness
	********************************************************************************

	noi di _n in y "Introducing missingness... "
	
	*** Random missing pattern
	if `mblock' ~= 1 {
		qui foreach scale of local namelist {
			unab scale_items: `scale'*
			local misslist ""
			tempvar missTotal
			foreach item of local scale_items {
				replace `item' = . if runiform() <= `propmiss'
				local misslist "`misslist', `item'"
			}
			gen missTotal_`scale' = 1- missing(`misslist')
			qui sum missTotal_`scale'
			noi di in y "Complete cases in scale `scale': `=`r(mean)' * 100'%"
		}
	}
	
	*** Block missing pattern
	
	else {
		*** sample a nwithin item (i.e. the time period) to determine missing observations
		tempvar obsmiss
		gen `obsmiss' =(runiform() <= `propmiss')
	
		*** sample the time period that would be missing  for each observation 
		*** and create a new variable called nwithin_miss that containes the 
		*** period for which obs is missing

		mata: ssize = st_nobs()
		mata: timePoints = J(1, strtoreal("`nwitems'"), 1/strtoreal("`nwitems'"))
		mata: myperiod = rdiscrete(ssize, 1, timePoints) // sample one period
		/* st_local("missitem", strofreal(myperiod)) */
		mata: mynewvar = st_addvar("int", "missTimePoint")  // add a new var
		mata: st_store((1, rows(myperiod)), mynewvar, myperiod) // fill in the new var
		
		*** change values to missing for observations with obsmiss == 1 for the
		*** selected time period
		qui forval i = 1/`=_N' {
			if (`obsmiss'[`i'] == 1) {
				local toMissing = missTimePoint[`i'] 
				foreach var of varlist _all {
					replace `var' = . if _n == `i' & substr("`var'", `=length("`var'")', `=length("`var'")') == "`toMissing'"
				}
			}
		}
		qui replace `obsmiss' = 1 - `obsmiss' 
		qui sum `obsmiss'
		noi di in y "Complete cases in the dataset: `=`r(mean)' * 100'%"
		drop missTimePoint  // replace with a temp var
	}
	
	
	**************************************
	*** From here on we have tests that we can spin into separate functions
	**************************************
	
	*misstable sum myvar*
	
	/*
	********************************************************************************
	*** Generate patterns var
	********************************************************************************
	
	noi di "Generating response patterns..."
	local mysum ""
	local myvars ""
	local i = 1
	foreach var of varlist myvar* {
		local myvars "`myvars' `var'"
		if "`i'" == "1" {
			local mysum "string(`var')"
		}
		else {
			local mysum "`mysum' + string(`var')"
		}
		local ++i	
	}
	gen myprofile = `mysum'
	noi sum

	********************************************************************************
	*** Compute tests

	preserve

	keep myprofile
	encode myprofile, gen(mytest)
	gen test = 1
	collapse (count)test, by(mytest)

	decode mytest, gen(test1)
	
	*** Diversity of patterns
	
	if `mblock' == 1 {
		gen mymiss = regexm(test1, "[\.]+")
		replace test = 1 if mymiss == 1
		*** non-missing:
		count if mymiss == 0
	}
	else {
		count
	}
	
	noi di _n "****************************************************************************"
	noi di "* Number of diverse patterns (across obs): " `r(N)' " out of " `nobs' " obs"
	noi di "* Percentage diversity of patterns (out of max): " `r(N)'/(2^`=strlen(test1[1])')*100 "%"
	noi di "****************************************************************************"
		
	restore


	*** Measure the entropy of the string (Shannon)

	* gen mymiss = regexm(myprofile, "[\.]+")
	* entropyetc myprofile if mymiss == 0 // this is not a very good test

	egen mymean = rowmean(myvar*)
	local mystrlen = strlen(myprofile)
	gen ones = mymean * `mystrlen'
	gen p1 = ones/`mystrlen'
	gen p0 = 1 - p1
	gen entropy = cond(p0 == 0 | p1 ==0, 0, - (p0 * ln(p0)/ln(2) + p1 * ln(p1)/ln(2)))

	capture drop mymean ones p0 p1

	// within a subject
	noi tabstat entropy, s(mean sd p5 q p95)
	* hist entropy if `obsmiss' ~= 1, kdensity xline(1) xlab(0(0.1)1) title("Density of (Shannon) Entropy (within obs)") name("Entropy")


	*** Distance between strings (Hamming distance)

	noi di _n "Computing distances between patterns across observations..."
	if `mblock' == 1 {
		qui levelsof myprofile if `obsmiss' ~= 1, local(levs) 
	}
	else {
		qui levelsof myprofile, local(levs) 
	}
	local mydim "`: word count `levs''"
	if `mydim' <= 250 {
		mat mynewmat = J(`mydim', `mydim', .)
		nois _dots 0, title(# of patterns processed) reps(`mydim')
		local i = 1
		foreach ilev of local levs {
			local j = 1
			local mystrlen = strlen("`ilev'")
			foreach jlev of local levs {
				if (`i' ~= `j') {
					qui strdist "`ilev'" "`jlev'"
					local myval =`r(d)'/`mystrlen'
					mat mynewmat[`i', `j'] = `myval'
				}
				else {
					mat mynewmat[`i', `j'] = 0 		
				}
				local ++j
			}
			nois _dots `i' 0
			local ++i
		}

		//across subjects
		noi di _n "Plotting distances..."
		* plotmatrix, mat(mynewmat) lower color(red) title("(Hamming) Distances between patterns across obs") name("DistanceMatrix")
	
		preserve

		svmat mynewmat, name(mynewmat)
		stack mynewmat*, into(distances) clear
		noi tabstat distances, s(mean sd p5 q p95)
		
		* hist distances if distances ~=., xline(1) xlab(0(0.1)1) title("(Hamming) Distances between patterns across obs") name("Distance")

		restore
	
	}
	else {
		noi di "Distances not computed."
	}
	*/

	*** Replace the following with pchained!
	
	
	
	/*
	*** Imputation
	mi set flong
	mi register imputed _all

	local mimodels ""
	forval i= 1/`ntitems' {
		local mimodels "`mimodels' (logit, augment) myvar`i'"
	}
	local mimodels "`mimodels' , add(5) burnin(10) chaindots"

	capture mi impute chained `mimodels'
	*/
	
	*** Store results
	mat simmat[`rows', `cols'] = _rc


end






