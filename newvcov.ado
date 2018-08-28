***********************************************************
**** Utility functions for simulation
*
*
*
*
*
*
* Run with
* mysim, simobs(50(50)150) simcorb(0(.1).2) nwitems(3) ntitems(12) corw(0.9) propmiss(0.2) mblock(0) simvcov(1)
* mysim, simobs(50(50)150) simcorb(0(.1).2) nwitems(3) ntitems(12) corw(0.9) propmiss(0.2) mblock(0) simvcov(0) storedvcov("~/Desktop/test.dta")

********************************************************************************
**** Extract vcov from original data

capture program drop myvcov
program define myvcov

	syntax varlist, saving(string)

	scaleclean `varlist'
	
	
	**** Change by Zitong: 
	**** Use Spearman's correlation for ordinal categorical variables' correlation. 
	
	spearman `varlist', stats(rho)
	mat dcorr = r(Rho)
		
	
	*** fudge the correlation matrix a bit (replace missing with 0)
	mata: fixed = editmissing(st_matrix("dcorr"), 0)
	mata: _diag(fixed, 1)
	mata: st_matrix("dcorr", fixed)
	
	preserve

	clear
	svmat dcorr
	save `saving', replace
*	export delimited using "$folder/categorical_corr.csv", replace 
	restore

end
********************************************************************************	
	
	**** This is for future work. Don't waste time reading this part now. 
	**** Value labels list

	/*
	label dir 
	local labname `r(names)'
	foreach lab of local labname {
	* label list `lab'
	* local value`lab' `r(k)'
	qui ds, has(vall `lab')
	local items_`lab' = "`r(varlist)'"
	* di `items'
	label list `lab'
	di `r(k)'
	* lobal cas_`lab' = 
	}
	*/
	

********************************************************************************
*** Define three matrices 
*** - within item correlations
*** - between item correlations
*** - means


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

capture program drop genmeans
program define genmeans, rclass

	args nmeans

	mat mymeans = J(1, `nmeans', 0)
	return matrix means = mymeans
end






**** MOdification by Zitong:
**** It is better to tell Stata the varlist we want to simulate on for the sake of reading dta files in. 
**** I add this as an argument called "storedvars"
**** Use "storedvars" to create "storedscales"(input argument of catsimcore)

capture program drop catsim
program define catsim, rclass
	
		syntax, simobs(numlist) [nwitems(integer 3) ntitems(integer 12) propmiss(real 0) mblock(integer 0) simvcov(integer 1) simmarginal(integer 1) storedvcov(string) storedmarginal(string) storedvars(string)]
	

	if `simmarginal' == 1 {
	**** Todo. 
	}
	else {
		capture unab `storedvars':  `storedvars'
		qui scaleclean `storedvars' 
	
		mat simmat = J(1, `: word count `simobs'', .)
		local rows = 1
		local cols = 1
		foreach lobs of local simobs {
			catesimcore `lobs' 1 `cols' `rows' `nwitems' `ntitems' 1 `propmiss' `mblock' `simvcov' `simmarginal' `storedvcov' `storedmarginal' `r(scales)'
			local ++cols
		}
	}
	mat list simmat
	return matrix sims = simmat
	
end



capture program drop catesimcore
program define catesimcore
 
	args nobs corb cols rows nwitems ntitems corw propmiss mblock simvcov simmarginal storedvcov storedmarginal storedscales
	
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
	set obs `nobs'

	
	**** Need change the display command.
	if `simvcov' == 1 {
		di _n in y "************************************************************"
		di in y    "* Sample size        : " `nobs'
		di in y    "* Within correlation : " `corw'
		di in y    "* Between correlation: " `corb'
		di in y    "************************************************************"
	}
	else {
		di _n in y "************************************************************"
		di in y    "* Sample size        : " `nobs'
		di in y    "* Using empirical vcov "		
		di in y	   "* Using empirical marginal distribution"
		di in y "************************************************************"
	}
	


		********************************************************************************
		*** Define means and vcov

		genmeans `ntitems'
		mat means = r(means)
		
		noi di "Simulating data..."
		
		if `simvcov' == 1 {
			noi di "(Simulated correlation matrix)"
			gencombined `nwitems' `corw' `ntitems' `corb'
			mat rho = r(combined)
		}
		else {
			preserve
			use "`storedvcov'", clear
			mkmat _all, matrix(rho)
			restore
			noi di "(Correlation matrix of the data)"
		}
		/*
		if `simmarginal' == 1 {
		**** todo. In this case, create cut points from simulated marginal distribution. 
		}
		*/
*		else { // In the else case, create cut points from mu matrix, which is imported. 
		
		**** Modification here (08/13) 
		**** Read in empirical marginal distribution from different scales. 				
		
		local scales `storedscales'
		
		**** Change codes to loop for scales. 
		local i = 1	
		foreach sca of local scales {
			preserve
			use "`storedmarginal'/sca`sca'perc.dta", clear
			qui su
			drop if _n == _N
			drop cate_perc
			mkmat _all, matrix(mu_`i')
			local ++i 
			restore 
		}
			
		local --i 
		di "This is my i : `i'"
			
			local k = 1 // This is the column index in the output
			local mcol = 1 // Column index for each matrix
			local mrow = 1 // Row index for each matrix 
			
			forval j = 1/`i' {  // j is the index for matrices
			forval mcol = 1/`ntitems' {
			local cp`k' ""
			forval mrow = 1/`r(N)' {
			local cutp = mu_`j'[`mrow', `mcol']
			local cp`k' = "`cp`k''" + ", `cutp'"	
			}
			local cp`k' = subinstr("`cp`k''", ", ", "", 1)
			local cp_`k' = "`cp`k''"
			local ++k
			}			
* 			local --k
			noi di "(Marginal distribution from the data)"
			local ++j
		}
		
		drawnorm myvar1-myvar`ntitems', means(means) corr(rho) n(`=_N')
		
		forval i = 1/`ntitems'{
		egen cmyvar`i' = xtile(myvar`i'), percentiles(`cp_`i'')
		drop myvar`i' 
		rename cmyvar`i' myvar`i' 
		replace myvar`i' = myvar`i'- 1 // Changed 08/08, we want every catogorical variable to start from 0 instead of 1. 
		}
		sum myvar*
		
				
		********************************************************************************
		** Introduce missingness
		********************************************************************************

		*** Random missing pattern
		tempvar obsmiss
		if `mblock' ~= 1 {
			foreach var of varlist myvar* {
				replace `var' = . if runiform() <= `propmiss'
			}
		}
		*** Block missing pattern!!!
		else {
			** sample a nwithin item (i.e. the time period). this is hard coded! has to be re-written
			
			gen `obsmiss' = 1 if runiform() <= `propmiss'
			

			*** sample the nwithin item (or the time period) that would be missing 
			*** for each observation
			* mata 
				mata: ssize = st_nobs()
				mata: myperiod = rdiscrete(ssize,1,(0.333333, 0.333333, 0.333334))
				/* st_local("missitem", strofreal(myperiod)) */
				mata: mynewvar = st_addvar("int", "nwithin_miss")
				mata: st_store((1,rows(myperiod)),mynewvar, myperiod)
			
			* end
			
			*** change variable values to missing for obsmiss == 1
			forval i = 1/`=_N' {
				if (`obsmiss'[`i'] == 1) {
					local myval = nwithin_miss[`i'] 
					forval j = `myval'(`nwitems')`ntitems' {
						replace myvar`j' = . if _n == `i'
					}
				}
			}
		}	
		
		*misstable sum myvar*
		
		
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
	

	*** Imputation
	mi set flong
	mi register imputed myvar* 

	local mimodels ""
	forval i= 1/`ntitems' {
		local mimodels "`mimodels' (logit, augment) myvar`i'"
	}
	local mimodels "`mimodels' , add(5) burnin(10) chaindots"

	capture mi impute chained `mimodels'
	
	*** Store results
	mat simmat[`rows', `cols'] = _rc


end






