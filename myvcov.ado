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
	
	order `varlist', alpha
	* corr `varlist'
	
	**** Some people say "don't know", it is coded as 99 or .d, I change these answers into missing values to be imputated. 
	foreach var of varlist `varlist'{
	qui replace `var' = . if (`var' == 99 ) | (`var' == .d)
	}
	
************** Don't read this *************************************************	
	**** Drop the two variable labels "99" and ".d".
	**** Rename label values from 0 to max. 
	
	**** TIP: recode here? 
	/*
	foreach var of varlist `varlist' {
	qui sum `var' 
	if !`r(min)' {
	**** Label? 
	replace `var' = `var' - `r(min)'
	}
	}
	*/
************** Don't read this *************************************************	
	
	
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

/* OLDER VERSION OF MYPERCENT

capture program drop mypercent
	program define mypercent
	
	syntax varlist, saving(string)
	
	order `varlist', alpha
	foreach var of varlist `varlist'{
	qui replace `var' = . if (`var' == 99 ) | (`var' == .d)
	}
	
	**** This is calculate the mean of each categories. 
	
	foreach var of varlist `varlist' {
	preserve
	tempfile perc_`var'
	contract `var', cpercent(percentage_`var') nomiss
	drop if _n == _N
	drop _freq
	rename `var' cate_perc
	save "`perc_`var''", replace
	restore
	}
	

	
	local var1: word 1 of `varlist'
	local newlist: list local(varlist) - local(var1)
	
	preserve
	use "`perc_`var1''", clear	
	foreach var of local newlist {
	merge 1:1 cate_perc using "`perc_`var''", nogen
	}
	recode _all (missing=1e-5)
	drop cate_perc
	save `saving', replace
	
	restore
	
* 	export delimited using "$folder/categorical_percent.csv", replace // Originally used for R. 
end
*/

capture program drop mypercent
program define mypercent
	
	syntax varlist, saving(string)

	capture mkdir "`saving'/percents"
	local saving "`saving'/percents"
	
	scaleclean `varlist'
	
	qui{
		local scales `r(scales)'
		foreach sca of local scales {
			local subsca_`sca' `r(subs`sca')'
		} // This is awkward but I don't know elegant way to inherit the r().
	
	**** Run contract for every argument variables. 
	
	local i = 1 // i: Scale index. 
	foreach sca of local scales {
		local j_`i' = 1 // j_k: Variable(subscale) index k of jth scale. 
		foreach var of local subsca_`sca' { 	
			preserve
			tempfile perc_l`i'v`j_`i''
			contract `var', cpercent(percentage_l`i'v`j_`i'') nomiss
			* drop if _n == _N // Remember you changed here!  
			drop _freq
			rename `var' cate_perc
			save "`perc_l`i'v`j_`i'''", replace
			restore
			
			local ++j_`i'
		}
		local --j_`i'
		local ++i
	}
	local --i 
}	
	**** 1st Merge for subscales from the same scale
	**** Results of this step are several dta files for scales with information of subscales. 
	
 	forval ii = 1/`i' {
		qui{
			preserve
			*tempfile sca`ii'perc
			use "`perc_l`ii'v1'", clear
			
			forval jj = 2/`j_`ii'' {
				merge 1:1 cate_perc using "`perc_l`ii'v`jj''"
				drop _merge
			}
			recode _all (missing=1e-5)
			*drop cate_perc
			local datascale : word `ii' of `scales'
		}
		save "`saving'/sca`datascale'perc", replace // I put the scale name in the file name. 
		restore
	}	
	
	
	**** 2st Merge dta files for scales to a whole dta file. 
	**** Don't know how to do this. Leave it here. 
	/*
	preserve
	use "`sca1perc'", clear
	forval ii = 2/`i' {
	merge 1:1 ???
	}
	*/	
* 	export delimited using "$folder/categorical_percent.csv", replace // Originally used for R. 


end





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





********************************************************************************
*** Simulation function

capture program drop mysim
program define mysim, rclass
	syntax, simobs(numlist) [simcorb(numlist) nwitems(integer 3) ntitems(integer 12) corw(real 0) propmiss(real 0) mblock(integer 0) simvcov(integer 1) storedvcov(string)]

	if `simvcov' == 1 {
	
		mat simmat = J(`: word count `simcorb'' , `: word count `simobs'', .)
		
		local cols = 1
		foreach lobs of local simobs {
			local rows = 1
			foreach lcorb of local simcorb {

				mysimcore `lobs' `lcorb' `cols' `rows' `nwitems' `ntitems' `corw' `propmiss' `mblock' `simvcov' `storedvcov'
				local ++rows
			}
			local ++cols
		}
	}
	else {
		mat simmat = J(1, `: word count `simobs'', .)
		local rows = 1
		local cols = 1
		foreach lobs of local simobs {
			mysimcore `lobs' 1 `cols' `rows' `nwitems' `ntitems' 1 `propmiss' `mblock' `simvcov' `storedvcov'
			local ++cols
		}
	}
	mat list simmat
	return matrix sims = simmat
	
end





capture program drop mysimcore
program define mysimcore

	args nobs corb cols rows nwitems ntitems corw propmiss mblock simvcov storedvcov 
	
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
		di in y "************************************************************"
	}
	
	qui {

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
		
		drawnorm myvar1-myvar`ntitems', means(means) corr(rho) n(`=_N')

		
		
		*** Convert normals to BINARY
		foreach var of varlist myvar* {
			replace `var' = cond(`var' <= 0, 0, 1)
		}

		*** Need another one for categorical
		
		**** Zitong's Part:
		**** This part owes its theoritical background to  Kaiser, Träger and Leisch (2011) 
		**** "Generating Correlated Ordinal Random Values". 
		**** In this paper, authors propose two different methods for simulation 
		**** of ordinal multiple variables: the binary conversion method and mean mapping method. 
		**** In our situation we have at most 4 categories and usually more than 10 variables. 
		**** Mean mapping method is more consistent with our situation. 
		**** My codes are based on mean mapping method in KTL(2011). 
		
		
		**** Try Define a program "ordsim" later. 
		
		* Define rho_c
		
		
		
		
		
		
		/*
		foreach var of varlist myvar* {
		
		}
		*/
		
		

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
	* return matrix sims = simmat

end




********************************************************************************
********************************************************************************
*** Simulation categorical vars
********************************************************************************
********************************************************************************



capture program drop catsim
program define catsim, rclass
	
		syntax, simobs(numlist) [nwitems(integer 3) ntitems(integer 12) propmiss(real 0) mblock(integer 0) simvcov(integer 1) simmarginal(integer 1) storedvcov(string) storedmarginal(string)]
		
	if `simmarginal' == 1 {
	**** Todo. 
	}
	else {
		mat simmat = J(1, `: word count `simobs'', .)
		local rows = 1
		local cols = 1
		foreach lobs of local simobs {
			catesimcore `lobs' 1 `cols' `rows' `nwitems' `ntitems' 1 `propmiss' `mblock' `simvcov' `simmarginal' `storedvcov' `storedmarginal'
			local ++cols
		}
	}
	mat list simmat
	return matrix sims = simmat
	
end



capture program drop catesimcore
program define catesimcore, rclass
 
	args nobs corb cols rows nwitems ntitems corw propmiss mblock simvcov simmarginal storedvcov storedmarginal
	
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
	
	qui {

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
		
		if `simmarginal' == 1 {
		**** todo. In this case, create cut points from simulated marginal distribution. 
		}
		else { // In the else case, create cut points from mu matrix, which is imported. 
			preserve
			use "`storedmarginal'", clear
			qui su
			mkmat _all, matrix(mu)
			restore
			
			local mcol = 1
			local mrow = 1
			forval mcol = 1/`ntitems' {
				local cp`mcol' ""
				forval mrow = 1/`r(N)' {
					local cutp = mu[`mrow', `mcol']
					local cp`mcol' = "`cp`mcol''" + ", `cutp'"
				}
				local cp`mcol' = subinstr("`cp`mcol''", ", ", "", 1)
			}	
			noi di "(Marginal distribution from the data)"
		}
				
		drawnorm myvar1-myvar`ntitems', means(means) corr(rho) n(`=_N')
				
		**** By Zitong:  Convert normals to CATEGORICAL
	
		forval i = 1/`ntitems' {
			*** requires egenmore package ***
			noi di "`myvar`i''"
			noi di "`cp`i''"
			
			egen cmyvar`i' = xtile(myvar`i'), percentiles(`cp`i'')
			drop myvar`i'
			rename cmyvar`i' myvar`i'
		}								
						

		
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
	* return matrix sims = simmat


end





