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

syntax, [corw(real 1) corb(real 1) corbs(real 1) ///
			nWithinItems(string) nwitems(integer 3) nTotalItems(integer 12)]

*	args dimw corw dimb corb
	mat fullcorr = J(`nTotalItems', `nTotalItems', `corbs')
	forval k = 1/`: word count `nWithinItems'' { // Loop over scales
	mat mybetween = J(`:word `k' of `nWithinItems'', `:word `k' of `nWithinItems'', `corb')
	genwithin `nwitems' `nwitems'
	mat within = r(within)
	forval i = 1(`nwitems')`: word count `nWithinItems'' {
		forval j = 1(`nwitems')`: word count `nWithinItems'' {
			if (`i'==`j') {
				mat mybetween[`i', `j'] = within
			}
		}
	}
	* By Zitong: Ask Simo. Rewrite the within Scale part properly
	* I don't know how to do it. 
	} // End of loop over scales
	return matrix combined = fullcorr
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


*** Set up a theoretical uniform marginal distribution
*** Similar with catDist

**** Ask this super large question: Should we use temp files? Anyway we need the long variable names list. 

capture program drop genunimarginals
program define genunimarginals, rclass

		syntax namelist, [scales(string)]
	
	noi di "(Generating Uniform Distributions for `namelist')" // Zitong: Start from uniform case

	qui{
	
		local nscales: word count `namelist'

		*** Parse the scales string: scales(kzf=(numlist) hsclg=(numlist))
		*** Retrieve scale name and scale levels from user input
		
		if "`scales'" == ""{
		**** This part does not make much sense because there are 2 choices: use label or use the "empirical" levels from data. 
		**** Label has the problem of missing value label and empirical thing is just "empirical". Difficult to explain. 
		**** Leave to discuss with Simo
		**** Maybe we just change the option of "simmarginal" and "simscales" to one, string option. 
		**** If you want to simmulate, then you need to give me at least how many categories you want to use. 

		
		noi di "No User-defined Scale Patterns. Use Scale Levels from Data. "
		
		foreach scale of local namelist { //iterate over scales
			unab scale_items: `scale'*
			
			*** Retrieve levels from label
			local labname: value label `: word 1 of `scale_items''
			local labs ""
			if "`labname'" ~= "" {
				mata: values = .; text = ""
				mata: st_vlload("`labname'", values, text); _transpose(values)
				mata: st_local("labs", invtokens(strofreal(values)))
			} 
			
				
			*** Retrieve levels from items
			local catsData ""
			foreach item of local scale_items { // iterate over items of a scale
				qui levelsof `item', local(levs) 
				local catsData: list catsData | levs // keep the union of the lists
			} // This is a little bit empirical. I put it here but I didn't use it. 
			
			
			
		} // End of loop over scales
		} // End of no user-option case
		
		if "`scales'" ~= ""{
		**** Better case. The user will give how many categories they want to use. 
		**** And we always follow the user categories. 
		**** It's good to report how the difference here, because it provides the user information about the "empirical" part. 
		
		noi di "User-defied Scale Pattern will be used. "
		parse_syntax "`scales'" // Parsed result can by EACH stubs		
		}
		
		foreach scale of local namelist { //iterate over scales
			unab scale_items: `scale'*
			
			*** Retrieve levels from label
			local labname: value label `: word 1 of `scale_items''
			local labs ""
			if "`labname'" ~= "" {
				mata: values = .; text = ""
				mata: st_vlload("`labname'", values, text); _transpose(values)
				mata: st_local("labs", invtokens(strofreal(values)))
			}
			
				
			*** Retrieve levels from items
			local catsData ""
			foreach item of local scale_items { // iterate over items of a scale
				qui levelsof `item', local(levs) 
				local catsData: list catsData | levs // keep the union of the lists
			}

			noi di _n in y "****************************************" ///
			_n      "** Scale `scale' " ///
			_n      "****************************************"
			
			if "`labs'" ~= "" {
			
				*** Compare LABELs with user input
				local userInput "`s(`scale')'"
				compare_lists "`labs'" "`userInput'"
						
				if "`s(differences)'" == "" {
					noi di in y "LABELS: user input matches label information."
				}
				else {
					noi di in y "LABELS: label values are DIFFERENT from user input. " ///
					_n "    Differences are: `s(differences)'" /// 
					_n "    Using user input"
				}
			}
			
			*** Compare ITEM CATS with user input
			compare_lists "`catsData'" "`userInput'"
			
			if "`s(differences)'" == "" {
				noi di in y "DATA: user input matches categories in items."
			}
			else {
				noi di in y "DATA: item categories are DIFFERENT from user input. " ///
				_n "    Differences are: `s(differences)'" /// 
				_n "    Using user input"
			}
			
			
			* count how many items i have. Create a matrix with this many columns. 
			
			local nCats: list sizeof userInput
			local marg_prob = 1/ `nCats'
			local nScales: list sizeof scale_items
			mata: mat_marg = J(strtoreal(st_local("nCats")), strtoreal(st_local("nScales")),100 * strtoreal(st_local("marg_prob"))- 1e-5)
			mata: mat_cd = mm_colrunsum(mat_marg)
			mata: st_matrix("mat_cd", mat_cd) 
		return matrix `scale'_cd = mat_cd // Return a matrix for every stub in the "namelist"		
	} 
} // End of qui
end









********************************************************************************
*** Simulation programme



**** Modification by Zitong:
**** It is better to tell Stata the varlist we want to simulate on for the sake of reading dta files in. 
**** I add this as an argument called "storedvars"
**** Use "storedvars" to create "storedscales"(input argument of catsimcore)


capture program drop ifeats
program define ifeats, rclass
	
	syntax namelist, Nobs(numlist) /// 
	       [NWitems(integer 3) NTitems(integer 12) propmiss(string) nwavemiss(string) /// 
		    MBlock(namelist) SIMCorr(integer 0) SIMMarginals(integer 0) SIMScales(string) ///
			CORRMatrix(string) MARGinals(string)]
	
	
	* nwavemiss(numlist) --> numlist per scale!
	* MBlock(namelist) --> scales with block missingness
	* propmiss(numlist) --> per scale proportion of missing
	* Zitong: simScales(string) --> the same thing as scales option in catDist
	
	
	*** Total number of items
	local nTotalItems = 0
	local nWithinItems ""
	foreach scale of local namelist {
		unab `scale'_items: `scale'*
		local nItems: list sizeof `scale'_items 
		local nWithinItems `macval(nWithinItems)' `macval(nItems)' // by Zitong: Add a counter for within items
		local nTotalItems = `nTotalItems' + `nItems'
		local itemNames = "`itemNames'" + " ``scale'_items'" 
	}
	

	/* Written by Simo: 
	if (`nwavemiss' >= `nwitems') {
		noi di in r "Number of missing timepoints should be smaller than the total number of timepoints"
		exit 618
	}
	*/ 
	**** Zitong: Let's deal with simulating marginal distributions before clear the dataset. 
	**** Work here: 
	
	**** I need to follow catdist and define the categorical names etc. 
	
	if `simmarginals' == 1 { 	
	*** Create even distribution. 
	
	noi di _n in y "****************************************" ///
	_n "(Use Simulated Marginal Distribution)" ///
	_n "****************************************"
	
	genunimarginals `namelist', scales(`simscales')
	
	foreach scale of local namelist {
		mat `scale'_cd = r(`scale'_cd)
		}
	}

		* capture unab `storedvars':  `storedvars'
		* qui scaleclean `storedvars' 
	
		mat simmat = J(1, `: word count `nobs'', .)
		local rows = 1
		local cols = 1
		foreach lobs of local nobs {
			ifeatsCore,  nobs(`lobs')  corb(1) cols(`cols') rows(`rows') nwitems(`nwitems')  ///
			corw(1) propmiss(`propmiss')  nwavemiss(`nwavemiss') simcorr(`simcorr')  simmarginals(`simmarginals') simscales("`simscales'") ///
			corrmatrix("`corrmatrix'")  marginals("`marginals'") namelist("`namelist'") nTotalItems(`nTotalItems') ///
			mblock(`mblock') nWithinItems(`nWithinItems') corbs(1) itemlist("`itemNames'")

			local ++cols // Needs revision: propmiss, nwavemiss
		}
		
		* Need to revise the simscales part, set some default for it 
		
	
end



capture program drop ifeatsCore
program define ifeatsCore
 
* 	args nobs corb cols rows nwitems corw propmiss mblock simcorr simmarginals corrmatrix marginals namelist nTotalItems nwavemiss

	syntax, [nobs(integer 1)  corb(real 1) cols(integer 1) rows(integer 1) nwitems(integer 3) ///
	corw(real 1) corbs(real 1) propmiss(string) nwavemiss(string) simcorr(integer 0) simmarginals(integer 0) ///
	corrmatrix(string) marginals(string) namelist(string) nTotalItems(integer 0) mblock(string) ///
	nWithinItems(string) simscales(string) itemlist(string)] // Zitong: Tricky! Need revision!!!
/*
	local nobs    = `lobs'     // data has 64
	local nwitems = 3         // number of time periods
	local ntitems = 33        // time periods * number of indicators
	local corw    = 0.9       // within item correlation (over time)
	local corb    = `lcorb'   // between item correlation
	local missing = 0.6       // proportion of missing (based on obs)
	local mblock  = 1         // block missing (boolean)
	local simvcov = 1         // simulate vcov or use data vcov (boolean)
	local propmiss 			 // missing rate for your stubs
	local nwavemiss 		 // How many blocks are missing, if non block missing, put a zero
	local corbs = 			// between scale correlation  (By Zitong: I think this is usually assumed to 0, while I'm not sure)
	local itemlist 			// A very long string passed to ifeatsCore. I'm not sure about the robustness of this arugment. 

*/

	clear	
	qui set obs `nobs'
	
	
	**** Need change the display command.
	if `simcorr' == 1 {
		noi di _n in y "************************************************************"
		noi di in y    "* Sample size        : " `nobs'
		noi di in y    "* Within correlation : " `corw'
		noi di in y    "* Between-Item, Within-Scale correlation: " `corb'
		noi di in y	   "* Between-scale correlation: " `corbs'	
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
		noi di "(Simulated correlation matrix)"
		noi di "my nwithinitems: `nWithinItems'"

		gencombined, corw(`corw') corb(`corb') corbs(`corbs') ///
			nWithinItems("`nWithinItems'") nwitems(`nwitems') nTotalItems(`nTotalItems') 
		mat rho = r(combined)
		
	}
	else {
		preserve
		use "`corrmatrix'", clear
		mkmat _all, matrix(rho)
		restore
		noi di "(Correlation matrix of the data loaded.)"
	}
	
	*** Simulate random normal data using the (empirical) distributions and corr mat
	drawnorm myvar1-myvar`nTotalItems', means(means) corr(rho) n(`=_N') // Zitong: I change the location of this command. 
	

	
	if `simmarginals' == 1 {
		**** by Zitong: set default to even distribution. 

		* Generate the items using r(`scale'_cd).

		local i = 1
		qui foreach scale of local namelist{
		local `scale'_items ""
		
		foreach name of local itemlist {
		if strpos("`name'", "`scale'") {
		local `scale'_items : list  `scale'_items | name
		}
		} // Zitong: Tricky! Need revision!!!

		
		local matcol = 1
			foreach item of local `scale'_items {
				mata: mat_ph = st_matrix("`scale'_cd")[1...,`matcol'..`matcol']
				mata: _transpose(mat_ph)
				mata: st_local("pctiles", invtokens(strofreal(mat_ph,"%12.10g")))
				egen `item' = xtile(myvar`i'), percentiles(`pctiles') // this requires egenmore
				replace `item' = `item' - 1
				drop myvar`i'
				local ++matcol
				local ++i
				
		}
		}
	} // End of if simmarginals == 1
	

	
	else { // Use observed marginals 	
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
		
	*** Parse the syntax of nwavemiss propmiss
	parse_syntax "`nwavemiss'" "_wmiss"  // nwavemiss(sc1=(1 3) sc2=(0 1))
	parse_syntax "`propmiss'" "_pmiss"   // propmiss(sc1=0 sc2=0.2)
	
	qui foreach scale of local namelist {
		unab scale_items: `scale'*           // get scale item list
		local sca_pmiss `s(`scale'_pmiss)'   // rate of missing per item for scale
		local sca_wmiss `s(`scale'_wmiss)'   // waves missing (a list)
		local sca_bmiss : list posof "`scale'" in mblock  // is scale block missing?
		
		if (`sca_bmiss' == 0) { // Random missing pattern
			
			*** if nwavemiss not specified by user --> missingness across all waves
			if ("`nwavemiss'" == "") {
				* local misslist ""
				* tempvar missTotal
				
				foreach item of local scale_items {
					replace `item' = . if runiform() <= `sca_pmiss'
					* local misslist "`misslist', `item'"
				}
				* noi misstable tree `scale_items'
				* gen missTotal_`scale' = 1 - missing(`misslist')
				* noi sum missTotal_`scale'
				* noi di in y "Complete cases in scale `scale': `=round(`r(mean)', .001) * 100'%" // QUESTION: Why do I need this? 
			}
			else {  // wavemiss is specified by user 					
				* local misslist ""
				* tempvar missTotal
				
				*** build change list 
				local varSubset ""
				*foreach wave of local nwavemiss {
				*	local varSubset "`varSubset' tp`wave'
				*}
				foreach item of local scale_items {
					if regexm("`item'", "[0-9]+$") {
						local mywave `=regexs(0)'
						if (`:list posof "`mywave'" in sca_wmiss') {
							noi di "`=regexs(0)'"
							replace `item' = . if runiform() <= `sca_pmiss'
						}
					}
				}
				* noi misstable tree `scale_items'
				*gen missTotal_`scale' = 1 - missing(`misslist')
				*noi sum missTotal_`scale'
				* noi di in y "Complete cases in scale `scale': `=round(`r(mean)', .001) * 100'%" // QUESTION: Why do I need this? 
			}
		}
		else { //Block random missing fpattern
			
			*** TODO
			
			*** sample a nwithin item (i.e. the time period) to determine missing observations
			tempvar obsmiss 
			local obspropmiss = `nwitems' * `sca_pmiss'/`sca_wmiss' // By Zitong: I changed here 
			gen `obsmiss' =(runiform() <= `obspropmiss')
			
			*** sample the time periods that would be missing  for each observation 
			*** and create a set of new variables that contain the periods for which obs is missing

			mata: ssize = st_nobs()
			mata: timePoints = J(1, strtoreal("`nwitems'"), 1/strtoreal("`nwitems'"))
			mata: myperiod = rdiscrete(ssize, strtoreal("`sca_wmiss'"), timePoints) // sample periods. Zitong changed here. 
			/* st_local("missitem", strofreal(myperiod)) */
			forval i = 1/`sca_wmiss' {
				* tempvar missTime`i'Point    // there is a bug if using tempvars
				mata: mynewvar`i' = st_addvar("int", "missTime`i'Point")
				mata: st_store((1, rows(myperiod)), mynewvar`i', myperiod[1..., `i']) // fill in the new var 
		}
		
		*** change values to missing for observations with obsmiss == 1 for the
		*** selected time period
		qui forval i = 1/`=_N' {
			if (`obsmiss'[`i'] == 1) {
				* noi di in g "`i'"
				forval j = 1/`sca_wmiss' { 
					local toMissing = missTime`j'Point[`i'] - 1
					* no di in y "`toMissing'"
					foreach var of varlist _all {
						replace `var' = . if _n == `i' & substr("`var'", `=length("`var'")', `=length("`var'")') == "`toMissing'"
					}
				}
			}
		}		
		qui replace `obsmiss' = 1 - `obsmiss' 
		qui sum `obsmiss'
		noi di in y "Complete cases in the dataset: `=round(`r(mean)', .001) * 100'%"
		drop missTime?Point `obsmiss'		
		} // End of "else: with block missing"		
	} // End of loop through namelist
	

	********************************************************************************
	** Prepare data for imputation
	********************************************************************************
	
	*** Reshape dataset to long and change var names
	local allItemsReshape ""
	foreach var of varlist _all {
		local stub = substr("`var'", 1, `=length("`var'") - 1')
		local allItemsReshape "`allItemsReshape' `stub'"
	}
	local allItemsReshape: list uniq allItemsReshape
		
	qui gen id = _n
	qui reshape long `allItemsReshape', i(id) j(time)
	
	qui foreach var of varlist *_tp {
		ren `var' `=substr("`var'", 1, `=length("`var'") - 3')'
	}
	
	********************************************************************************
	** Impute
	********************************************************************************
	
	*** Impute; run -pchained-
	noi di _n in y "Imputing with pchained..."
	capture pchained `namelist', p(id) t(time) mio("add(1) burnin(10) chaindots ")
	
	if _rc ~= 0 {
		noi di in r "Failed"
	}
	
	*** Store results
	mat simmat[`rows', `cols'] = _rc

end

