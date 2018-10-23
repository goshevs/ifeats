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
		local mynames: colnames dcorr
		
		*** doctor the correlation matrix (change missing to 0)
		mata: fixed = editmissing(st_matrix("dcorr"), 0)
		mata: _diag(fixed, 1)
		mata: st_matrix("dcorr", fixed)
		
		** Save the correlation matrix
		preserve
		clear
		
		svmat dcorr
		
		local counter = 1
		foreach var of varlist _all {
			* noi di "`: word `counter' of `mynames''"
			* noi di "`var'"
			ren `var' `: word `counter' of `mynames''
			local ++counter
		}
		
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
			nWithinItems(string) nwaves(integer 3) nTotalItems(integer 12)]

*	args dimw corw dimb corb
	mat fullcorr = J(`nTotalItems', `nTotalItems', `corbs')
	forval k = 1/`: word count `nWithinItems'' { // Loop over scales
	mat mybetween = J(`:word `k' of `nWithinItems'', `:word `k' of `nWithinItems'', `corb')
	genwithin `nwaves' `nwaves'
	mat within = r(within)
	forval i = 1(`nwaves')`: word count `nWithinItems'' {
		forval j = 1(`nwaves')`: word count `nWithinItems'' {
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


*** Parse the scales string: scales(kzf=(numlist) hsclg=(numlist))
capture program drop _genMarg
program define _genMarg, rclass

	syntax namelist, scales(string) nwaves(string) nItems(string)
	
	qui{
		*** Parse the scales string: scales(kzf=(numlist) hsclg=(numlist))
		*** Retrieve scale name and scale levels from user input
		parse_syntax "`scales'" // Parsed result can by EACH stubs		
		
		*****
		local scalePos = 1
		foreach scale of local namelist { //iterate over scales	
			local userInput "`s(`scale')'"		
			local nCats: list sizeof userInput
			local marg_prob = 1/`nCats'
			local nScales `=`: word `scalePos' of `nItems''*`nwaves''
			mata: mat_marg = J(strtoreal(st_local("nCats")), strtoreal(st_local("nScales")),100 * strtoreal(st_local("marg_prob"))- 1e-5)
			mata: mat_cd = mm_colrunsum(mat_marg)
			mata: st_matrix("mat_cd", mat_cd) 
			return matrix `scale'_cd = mat_cd // Return a matrix for every stub in the "namelist"	
			local ++scalePos
		}
		
	}
end


*** Retrieve number of items and waves
capture program drop _extractInfo
program define _extractInfo, sclass

	args namelist nwaves scale
	
	unab allVars: _all
	
	*** retrieve number of time periods
	if ("`nwaves'" == "-1") {
		local nWaves ""
		foreach var of local allVars {	
			local uTime = regexr("`var'", "^[a-zA-Z0-9]+_tp", "")
			local nWaves "`nWaves' `uTime'"
		}
		local nWaves: list uniq nWaves
		local nwaves "`: list sizeof nWaves'"
	}
	
	*** get the stubs from the marginals matrix	
	local uStubs ""
	foreach var of local allVars {	
		local uItem = regexr("`var'", "[0-9]+_tp[0-9]+$", "")
		local uStubs "`uStubs' `uItem'"
		
	}
	local uStubs: list uniq uStubs
	
	*** check whether stubs in namelist are present in corrmat stubs
	if ("`scale'" == "") {
		local nItems ""
		foreach scale of local namelist {
			if "`: list scale in uStubs'" ~= "" {
				unab `scale'_items: `scale'*
				local nItems "`nItems' `=`: list sizeof `scale'_items'/`nwaves''"
			}
			else {  // assume that all stubs should be present (we do not allow here for simulating some items)
				noi di "Stub not found in correlation matrix"
				exit 1000
			}
		}
	}
	else {
		if "`: list scale in uStubs'" ~= "" {
			unab `scale'_items: `scale'*
			local nItems "`=`: list sizeof `scale'_items'/`nwaves''"
		}
		else {  // assume that all stubs should be present (we do not allow here for simulating some items)
			noi di "Stub not found in correlation matrix"
			exit 1000
		}
	}
	
	sreturn local nItems "`nItems'"
	sreturn local nwaves "`nwaves'"

end

				
*** Sampling periods of missingness	
cap mata : mata drop sample_periods()
mata:
void function sample_periods(string scalar nwaves, string vector wavemiss) {
		real scalar ssize
		real matrix timePoints
		real matrix myperiod
		real matrix mynewvar

		nw = strtoreal(nwaves)
		mwaves = strtoreal(tokens(wavemiss))
		maxmw = rowmax(mwaves)
		
		ssize = st_nobs()
		timePoints = J(1, nw, 1/nw)
		myperiod = rdiscrete(ssize, maxmw, timePoints)

		for (i=1; i<=maxmw; ++i) {
			mynewvar = st_addvar("int", "missTime" + strofreal(i)+ "Point")
			st_store((1, rows(myperiod)), mynewvar, myperiod[1..., i])
		}
		st_local("maxmwave", strofreal(maxmw))
}
end




*** Introduce missingness
capture program drop _introMiss
program define _introMiss

	args scale_items scaleMiss itemMiss i missType
	
	tempvar missCases unifScale
	gen `unifScale' = runiform()
	sort `unifScale'
	gen `missCases' = (_n <= `scaleMiss' * _N)
	foreach item of local scale_items {
		if regexm("`item'", ".+_tp`i'$") {
			if ("`missType'" == "random"){
				tempvar unif
				gen `unif' = runiform() 
				gsort -`missCases' `unif'
				replace `item' = . if _n <= `itemMiss' * _N
			}
			else {
				replace `item' = . if `missCases' == 1
			}
		}
	}
end



********************************************************************************
*** Simulation programme



**** Modification by Zitong:
**** It is better to tell Stata the varlist we want to simulate on for the sake of reading dta files in. 
**** I add this as an argument called "storedvars"
**** Use "storedvars" to create "storedscales"(input argument of catsimcore)


capture program drop ifeats
program define ifeats, rclass
	
	syntax namelist, Nobs(numlist) scales(string) /// 
	       [nwaves(integer -1) nitems(string)  ///
		    propmiss(string) wavemiss(string) /// 
		    CORRMatrix(string) MARGinals(string)]
	
	
	* scales --> scales(sc1=(0(1)4) sc2=(0(1)4)): item levels of scales
	* nwaves(integer) --> number of waves/time periods
	* nitems(string) --> nitems(sc1=12 sc2=36): number of items per scale
	* propmiss(string) --> propmiss(kzf=(0.05 0.2) hsclg=(0.05 0.3)): for random missing (item missigness; scale missingness)
	*                      propmiss(kzf=(0.2) hsclg=(0.3)): for block missing (scale missingness)
	* wavemiss(sring)  --> wavemiss(kzf=(0 1) hsclg=(1 2)): waves missing for every scale
	* corrmatrix(string) --> path and name of file containing correlation matrix for scales in namelist
	* marginals(string) --> path to directory containing files of marginal distributions for scales in namelist
	
	
	if "`marginals'" == "" local simmarg "1"
	if "`corrmatrix'" == "" local simcorr "1"

	**** Check for existing data
	if ("`nitems'" == "") { // if number of items per scale is not specified
		local nameList ""
		foreach var of local namelist {
			local nameList "`nameList' `var'*"
		}
		* noi di "`nameList'"
		capture unab allItems: `nameList'
		if (_rc == 0) {                           // get info from memory
			* noi di in y "Using variables in dataset loaded in memory"
			
			*** retrieve number of time periods if not specified
			if ("`nwaves'" == "-1") {
				local nWaves ""
				foreach item of local allItems {
					local uTime = regexr("`item'", "^[a-zA-Z0-9]+_tp", "")
					local nWaves "`nWaves' `uTime'"
				}
				local nWaves: list uniq nWaves
				local nwaves "`: list sizeof nWaves'"
			}
	
			local nItems ""
			foreach scale of local namelist {
				unab `scale'_items: `scale'*
				local nItems "`nItems' `=`: list sizeof `scale'_items'/`nwaves''"
			}
			
			*** TODO -->
			*** check for matrices, if none, then run catDist and dataCorrMat
			*** <--
	
		}
		else if (_rc ~= 0) {                 // if no such stubs in the dataset
			if ("`corrmatrix'" ~= "") {      // check correlation matrix if given
				* noi di in y "Using the correlation matrix"
				preserve
				use "`corrmatrix'", clear
				_extractInfo "`namelist'" "`nwaves'"
				
				local nItems "`s(nItems)'"
				local nwaves "`s(nwaves)'"
				
				restore
			}
			else if ("`marginals'" ~= "") {   // check marginals if given
				* noi di in y "Using the marginal distribution(s)"
				preserve
				local nItems ""
				foreach scale of local namelist {
					use "`marginals'/`scale'_cdist", clear
					_extractInfo "1" "`nwaves'" "`scale'"
					local nItems "`nItems' `s(nWithinItems)'"
				}
				local nwaves "`s(nwaves)'"
				restore
			}	
			else {
				* noi di in y "No dataset in memory, no correlation matrix or marginals given" 
				if "`nwaves'" == "-1" | "`nitems'" == "" {
					no di in r "Arguments nwaves and nitems are required in this case"
					exit 1000
				}
			}
		}
		* noi di "`nItems'"
		* noi di "`nwaves'"
	}
	else {
		parse_syntax "`nitems'"
		local nItems ""
		foreach scale of local namelist {
			local nItems "`nItems' `s(`scale')'"
		}
	}
	* noi di "`nItems'"
	* noi di "`nwaves'"
	
	if ("`simmarg'" ~= "") { 	
		*** Create even distribution. 
	
		noi di _n in y "****************************************" ///
		_n "Simulating Marginal Distribution(s)" ///
		_n "****************************************"
	
		noi di "(Generating Uniform Distributions for `namelist')"
		_genMarg `namelist', scales(`scales') nwaves(`nwaves') nItems(`nItems')
	
		*** loading the matrices -- may not need to this here
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
		ifeatsCore,  namelist("`namelist'") nobs(`lobs')  corb(1) cols(`cols') rows(`rows') nwaves(`nwaves') nItems("`nItems'") ///
		corw(1)  corbs(1) propmiss(`propmiss') wavemiss(`wavemiss') simcorr(`simcorr')  simmarg(`simmarg')  ///
		corrmatrix("`corrmatrix'")  marginals("`marginals'") scales("`scales'") itemlist("`itemNames'")

		local ++cols // Needs revision: propmiss, wavemiss
	}
		
	* Need to revise the simscales part, set some default for it 
end



capture program drop ifeatsCore
program define ifeatsCore
 
* 	args nobs corb cols rows nwaves corw propmiss mblock simcorr simmarginals corrmatrix marginals namelist nTotalItems wavemiss

	syntax, [namelist(string) nobs(integer 1) corb(real 1) cols(integer 1) rows(integer 1) nwaves(integer -1) nItems(string) ///
			 corw(real 1) corbs(real 1) propmiss(string) wavemiss(string) simcorr(string) simmarg(string) ///
			 corrmatrix(string) marginals(string) scales(string) itemlist(string)]

	
	/*
	These need to be updated and streamlined
	
	local nobs    = `lobs'     // data has 64
	local nwaves = 3         // number of time periods
	local ntitems = 33        // time periods * number of indicators
	local corw    = 0.9       // within item correlation (over time)
	local corb    = `lcorb'   // between item correlation
	local missing = 0.6       // proportion of missing (based on obs)
	local mblock  = 1         // block missing (boolean)
	local simvcov = 1         // simulate vcov or use data vcov (boolean)
	local propmiss 			 // missing rate for your stubs
	local wavemiss 		 // How many blocks are missing, if non block missing, put a zero
	local corbs = 			// between scale correlation  (By Zitong: I think this is usually assumed to 0, while I'm not sure)
	local itemlist 			// A very long string passed to ifeatsCore. I'm not sure about the robustness of this arugment. 

*/

	clear	
	qui set obs `nobs'
	
	
	**** Display basics
	if ("`simcorr'" == "1") {
		noi di _n in y "************************************************************"
		noi di in y    "* Sample size        : " `nobs'
		noi di in y    "* Within-item correlation (over waves): " `corw'
		noi di in y    "* Between-item, Within-scale correlation: " `corb'
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
	*** Simulation
	********************************************************************************
	
	*** Correlation matrix
	
	*  HAVE TO FIX THIS --->
	if ("`simcorr'" ~= "") {	  // TODO
		noi di "Simulating correlation matrix"
		* noi di "my nwithinitems: `nWithinItems'"

		gencombined, corw(`corw') corb(`corb') corbs(`corbs') ///
			nWithinItems("`nWithinItems'") nwaves(`nwaves') nTotalItems(`nTotalItems') 
		mat rho = r(combined)
		
	} 
	* <---- 
	else { // load matrix of empirical correlations
		preserve
		use "`corrmatrix'", clear
		mkmat _all, matrix(rho)
		restore
		noi di "  Empirical correlation matrix loaded"
	}
	
	*** Means
	local nTotalItems `=`=(`=subinstr("`nItems'", " ", "+",.)')' * `nwaves''
	
	genmeans "`nTotalItems'"
	mat means = r(means)	
	
	*** Simulate random normal data using the (empirical) distributions and corr mat
	drawnorm myvar1-myvar`nTotalItems', means(means) corr(rho) n(`=_N')
		
	if ("`simmarg'" ~= "") {  // simulated marginals
		noi di "Simulating marginal distributions"
		
		*** Generate the item names
		local scalePos = 1	
		local j = 1
		qui foreach scale of local namelist { //iterate over scales	
			local `scale'_items ""
			local items: word `scalePos' of `nItems'
			forval i = 1/`items' {
				forval t = 0/`=`nwaves' - 1' {		
					if (`i' < 10) {
						local sItem "`scale'0`i'_tp`t'"
					}
					else {
						local sItem "`scale'`i'_tp`t'"
					}
					local `scale'_items "``scale'_items' `sItem'" 
				}
			}
			*** Apply the distribution to the latent variable to create cat vars
			local matcol = 1
			foreach item of local `scale'_items {
				mata: mat_ph = st_matrix("`scale'_cd")[1...,`matcol'..`matcol']
				mata: _transpose(mat_ph)
				mata: st_local("pctiles", invtokens(strofreal(mat_ph,"%12.10g")))
				egen `item' = xtile(myvar`j'), percentiles(`pctiles') // this requires egenmore
				replace `item' = `item' - 1
				drop myvar`j'
				local ++matcol
				local ++j
			}
			local ++scalePos
		}
	} // end of if
	
	else { // Load and use  observed marginals	
		*** Read in scale_cdist's as matrices
		noi di "  Empirical marginal distributions loaded..." _c
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
		noi di " and applied to simulated data"
	}
	* sum	
	********************************************************************************
	** Introduce missingness
	********************************************************************************

	noi di _n in y "Introducing missingness... "
	
	*** Identify the type of missingness requested by the user
	
	*** Block missing for all scales
	if regexm("`propmiss'", "^[ ]*0\.[0-9]+[ ]*$") {  // feed only 1 value to propmiss == block missing
	
		noi di "   Missingness in dataset is of block random pattern"
		
		if regexm("`wavemiss'", "^minmax\([0-9]+ [0-9]+\)") {  // min and max number of missing waves
			if regexm("`wavemiss'","[0-9]+ [0-9]+") {
				local wavemiss "`=regexs(0)'"
				capture numlist "`wavemiss'"
				if _rc {
					noi di in r "Error: wavemiss not specified correctly"
					exit 1000
				}
				local wavemiss "`r(numlist)'"
				local wavemiss: list uniq wavemiss	
			}
		}
		
		tempvar missCases unif
		gen `unif' = runiform()
		sort `unif'
		gen `missCases' = (_n <= `propmiss' * _N)
	
		mata: sample_periods("`nwaves'", "`wavemiss'")
		
		*** change values to missing for observations with obsmiss == 1 for the
		*** selected time period
		qui forval i = 1/`=_N' {
			if (`missCases'[`i'] == 1) {
				* noi di in g "`i'"
				forval j = 1/`maxmwave' { 
					local toMissing = missTime`j'Point[`i'] - 1
					* no di in y "`toMissing'"
					foreach var of varlist _all {
						replace `var' = . if _n == `i' & substr("`var'", `=length("`var'")', `=length("`var'")') == "`toMissing'"
					}
				}
			}
		}
		drop  missTime?Point  // does not work well with tempvars
	}
	else {
		*** Parse the syntax of wavemiss propmiss and store into locals
		parse_syntax "`wavemiss'" "_wmiss"  // wavemiss(sc1=(1 3) sc2=(0 1))
		parse_syntax "`propmiss'" "_pmiss"   // propmiss(sc1=(0.2 0.6) sc2=(0.2 0.5)) (item scale)
		qui foreach scale of local namelist {
			local `scale'_pmiss `s(`scale'_pmiss)'   // rate of missing for an item and entire scale
			local `scale'_wmiss `s(`scale'_wmiss)'   // waves missing (a list)
			*local `scale'_bmiss : list posof "`scale'" in mblock  // is scale block missing?
		}
		*noi di "``scale'_pmiss'"
		*noi di "`s(block)'"

		qui foreach scale of local namelist { // loop over scales

			unab scale_items: `scale'*
			*** Check how many values entered in sca_pmiss; 
			*** assign values to item and scale missingness
			local nvals: list sizeof `scale'_pmiss
			if (`nvals' == 2) { // random missing case
				local itemMiss: word 1 of ``scale'_pmiss'
				* noi di "`itemMiss'"
				local scaleMiss: word 2 of ``scale'_pmiss'
				* noi di "`scaleMiss'"
			}
			else if (`nvals' == 1) {  // block missing per scale 
				local scaleMiss "``scale'_pmiss'"
				local `scale'_bmiss = 1
			}	
			else { // error
				noi di in r "Error: have to specify at least one value per scale in propmiss"
				exit 1000
			}

			if ("``scale'_bmiss'" == "") { // Block missing is 0: random missing pattern
				noi di "   Missingness in scale `scale' is of random pattern"
				*** if wavemiss not specified by user --> missingness across all waves
				if ("``scale'_wmiss'" == "") {
					if ("`itemMiss'" == "") {
						noi di in r "Error: item missigness should be specified"
						exit 1000
					}
					else {
						*** Loop over waves and create missing observations in items
						forval i = 0/`=`nwaves' - 1' {   // loop over waves
							*** this could be a function
							_introMiss "`scale_items'" "`scaleMiss'" "`itemMiss'" "`i'" "random"
						}
					}
				}
				else {  // wavemiss is specified by user 					
					foreach i of local `scale'_wmiss {
						_introMiss "`scale_items'" "`scaleMiss'" "`itemMiss'" "`i'" "random"
					}
				}			
			}  // end of random missing
			
			else { //Block random pattern of missingness per scale
				noi di "   Missingness in scale `scale' is of block random pattern"
				*** if wavemiss not specified by user --> missingness across all waves
				if ("``scale'_wmiss'" == "") {
					forval i = 0/`=`nwaves' - 1' {   // loop over waves
						_introMiss "`scale_items'" "`scaleMiss'" 0 "`i'" "block"
					}
				}
				else {  // time periods of missingness are specified by user in wavemiss
					foreach i of local `scale'_wmiss {
						_introMiss "`scale_items'" "`scaleMiss'" 0 "`i'" "block"
					}
				}
			} // end of else
		} // End of loop through namelist
	} // end of else
	exit
	
	********************************************************************************
	** Prepare data for imputation
	********************************************************************************
	
	*** Reshape dataset to long and change var names
	local allItemsReshape ""
	foreach var of varlist _all {
		if !regexm("`var'", "^__[0-9]+.$") {  //use this to filter out tempvars
			local stub = regexr("`var'", "_tp[0-9]+$", "_tp")
			local allItemsReshape "`allItemsReshape' `stub'"
		}
	}
	local allItemsReshape: list uniq allItemsReshape
	noi di "`allItemsReshape'"

	
	qui gen id = _n
	qui reshape long `allItemsReshape', i(id) j(time)
	
	qui foreach var of varlist *_tp {
		ren `var' `=substr("`var'", 1, `=length("`var'") - 3')'  // hard coded, need to change!!!
	}
		
	********************************************************************************
	** Impute
	********************************************************************************
	
	*** Impute; run -pchained-
	noi di _n in y "Imputing with pchained..."
	capture noisily pchained `namelist', p(id) t(time) mio(add(1) burnin(10) chaindots)
	
	if _rc ~= 0 {
		noi di in r "Failed"
	}
	
	*** Store results
	mat simmat[`rows', `cols'] = _rc

end

