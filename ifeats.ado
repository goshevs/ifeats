***********************************************************
*** Core functions
**
* Developers: Simo Goshev, Zitong Liu
* 
*
*

********************************************************************************
*** Utility functions
********************************************************************************



*** Compute empirical correlation matrix
capture program drop dataCorrMat
program define dataCorrMat

	syntax namelist, saving(string)  // string here is a path
	
	qui {
		noi di in y "Creating correlation matrix... " _c
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
		
		noi di in y "done."
	}
end


*** Create the marginal distributions of items
capture program drop catDist
program define catDist

	syntax namelist, scales(string) [saving(string)] // string here is a (path and)filename

	qui { 
		local nscales: word count `namelist'

		*** Parse the scales string: scales(kzf=(numlist) hsclg=(numlist))
		*** Retrieve scale name and scale levels from user input
		if "`scales'" ~= "" {
			parse_syntax "`scales'"
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
			
			local userInput "`s(`scale')'"
			
			if "`labs'" ~= "" {
				*** Compare LABELs with user input
				compare_lists "`labs'" "`userInput'"
				
				noi di "A label exist for this scale. Inferring levels from label."
				
				if "`s(differences)'" == "" {
					noi di in y "LABELS: user input matches label information."
				}
				else {
					noi di in y "LABELS: label values are DIFFERENT from user input. " ///
					_n "    Differences are: `s(differences)'" /// 
					_n "    Using user input"
				}
			}
			else {
				*** Compare ITEM CATS with user input
				if "`scales'" ~= "" {
					parse_syntax "`scales'"
				}
				compare_lists "`catsData'" "`userInput'"

				noi di "No labels for this scale. Inferring levels from data."
			
				if "`s(differences)'" == "" {
					noi di in y "DATA: user input matches categories in items."
				}
				else {
					noi di in y "DATA: item categories are DIFFERENT from user input. " ///
					_n "    Differences are: `s(differences)'" /// 
					_n "Using user input"
				}
			}
			
			noi di in y _n "Creating marginal distributions of items of scale `scale'"
			*** Creates the file with item cummulative distributions		
			local counter = 1
			foreach item of local scale_items { // iterate over items of a scale
				qui levelsof `item', local(levs)   // get the levels of the item
				compare_lists "`levs'" "`userInput'"   // using user input
				local missLevs "`s(differences)'"   // get the missing levels
				local nMissLevs: word count of `missLevs'
				
				*** Get the distributions and create observations for missing levels
				preserve
				contract `item', percent(pct_`item') nomiss // get distribution
				local nLevs: word count of `levs'
				replace pct_`item' = pct_`item' - (`nMissLevs'*1e-5)/`nLevs'
				foreach lev of local missLevs { // add the missing levels
					set obs `=_N + 1'
					replace `item' = `lev' if _n == _N
					replace pct_`item' = 1e-5 if _n == _N
				}
				sort `item'
				gen double cum_`item' = sum(pct_`item')
				ren (`item' cum_`item') (cats `item')
				drop _freq pct_`item'
				
				*** Merge item files
				if (`counter' == 1) {
					save `scale'_cdist, replace // USE TEMP FILES
				}
				else {
					merge 1:1 cats using `scale'_cdist, nogen
					save `scale'_cdist, replace
				}
				
				local ++counter
				restore
			
			} // end of items of a scale
			
			
			if "`saving'" ~= "" {
				qui copy `scale'_cdist.dta `saving'/, replace
				di in gr "`scale'_cdist saved."
			}
		} // end of scale
	}
end

	
*** Parser of the user input with multiple arguments of the type
*** (sc1=(1(1)4) sc2=(0 2)) and variations
capture program drop parse_syntax
program define parse_syntax, sclass
	args myinput type

	local nlistex "[0-9]+\.?[0-9]*(\(?[0-9]+\.?[0-9]*\)?)?[ ]?([0-9]+\.?[0-9]*)?"
	local strregex "[a-zA-Z0-9_-]+[ ]*=[ ]*\(?`nlistex'\)?"
	
	*** Parses pretty general syntax
	while regexm("`myinput'", "`strregex'") {
		local scale `=regexs(0)'
		local myinput = trim(subinstr("`myinput'", "`scale'", "", .))
		gettoken sname levs: scale, parse("=")
		gettoken left levs: levs, parse("=")
		local levs = trim("`levs'")
		if regexm("`levs'", "`nlistex'") {
			numlist "`=regexs(0)'"
		}
		sreturn local `sname'`type' `r(numlist)'
	}
	*** Parses name of correlation martix in simulation
	if regexm("`myinput'", "corrmat=[a-zA-Z0-9]+") {
		local matname `=regexs(0)'
		gettoken left matname: matname, parse("=")
		gettoken left matname: matname, parse("=")
		sreturn local corrmat_arg `matname'
	}
end


*** Compare elements of lists and print elements that differ
capture program drop compare_lists
program define compare_lists, sclass
	args list1 list2
	
	local isect: list list1 & list2
	local union: list list1 | list2
	local lDiff: list union - isect // LONGER SHOULD BE FIRST!
	* di "`lDiff'"
	sreturn local differences `lDiff'

end

*** Generate even marginal distributions of items
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
			if regexm("`var'", "^[a-zA-Z0-9]+_tp[0-9]+$") {
				local uTime = regexr("`var'", "^[a-zA-Z0-9]+_tp", "")
				local nWaves "`nWaves' `uTime'"
			}
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



*** Mata utility functions

cap mata: mata drop withinItemCorr()
cap mata: mata drop means()
cap mata: mata drop sample_periods()
cap mata: mata drop parse_syntax()
cap mata: mata drop nearSPD()
cap mata: mata drop bwCorrMat()

mata:
/*** Generating within-item correlation matrix ***/
function withinItemCorr(real scalar corw, real scalar dim) {
	real matrix wiCorrMat
	wiCorr = J(dim, dim, corw)
	_diag(wiCorr, J(1, dim, 1))
	return(wiCorr)
}

/*** Generating means ***/
void function means(real scalar dim) {
	real matrix means
	means = J(1, dim, 0)
	st_matrix("means", means)
}
	
/*** Sampling periods of missingness ***/
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

/* Parser of syntax */
void function parse_syntax(string scalar tcorr, scaleCorr, iCorrMat) {
	
	nlistex = "[0-9]+\.?[0-9]*(\(?[0-9]+\.?[0-9]*\)?)?[ ]?([0-9]+\.?[0-9]*)?"
	strregex = "[a-zA-Z0-9_-]+[ ]*=[ ]*\(?" + nlistex + "\)?"

	/* parse scale info */
	result = J(1, 2, .)
	while (regexm(tcorr, strregex)) {
		t = tokeninit(char(9))
		tokenset(t, regexs(0))
		scale = tokengetall(t)
		tcorr = strtrim(subinstr(tcorr, scale, ""))
		t = tokeninit("=")
		tokenset(t, scale)
		scale = tokengetall(t)
		f=regexm(scale[2], nlistex)
		result = result \ strtoreal(tokens(regexs(0)))
	}
	scaleCorr = select(result, rowmissing(result):==0 )
	/* retrieve between scale correlation matrix */
	if (regexm(tcorr, "corrmat=[a-zA-Z0-9]+")) {
		t = tokeninit("=")
		tokenset(t, regexs(0))
		iCorrMat = st_matrix(tokengetall(t)[2])
	}
	else {
		iCorrMat = 0
	}
}

/*** generating nearSPD ***/
// Adapted from https://www.statalist.org/forums/forum/general-stata-discussion/general/1340469-drawnorm?p=1340544#post1340544
//Project correlation matrix to nearest definite matrix
//Adapted from http://blogs.sas.com/content/iml/2012/11/28/computing-the-nearest-correlation-matrix.html
function nearSPD(real matrix c, real scalar tol) {
	real matrix W, L2, L
	maxd = 1
	Xold= Yold = c
	ds = J(rows(c), cols(c),0)
	printf("   Finding nearest semi definite matrix...")
	while (maxd > tol) {
		R= Yold - ds
		eigensystem(R, W, L)
		L2 = L
		L2 = Re(L2)
		for(k=1;k<=cols(L);k++) {
			L2[1,k] = max((L2[1,k], 0))
		}
		X = W * diag( L2 ) * W'
		ds = X - R
		Y = makesymmetric(X)
		for(k=1;k<=cols(c);k++) {
			Y[k,k] =  1
		}
		dx = max(abs(rowsum(X - Xold)))/ max(abs(rowsum(X)))
		dy = max(abs(rowsum(Y - Yold)))/max(abs(rowsum(Y)))
		dxy = max(abs(rowsum(Y - X)))/max(abs(rowsum(Y)))
		maxdxy = dx, dy, dxy
		maxd = max(maxdxy)
		// maxd
		Xold = makesymmetric(X)
		Yold = makesymmetric(Y)
	}
	printf("  done\n")
	return(makesymmetric(Re(Y)))
}	


/*** Combined correlation matrix ***/
void function bwCorrMat(string scalar tcorr, string scalar nItems, string scalar nwaves, real scalar tol) {

	real scalar dim, nCorbs, nWaves, start
	real vector nCorwi, nCorbi, vnItems
	real matrix bCorrMat, corrMat

	/* parse the syntax of tcorr */
	parse_syntax(tcorr, scaleCorr, iCorr)

	nWaves = strtoreal(nwaves)
	vnItems = strtoreal(tokens(nItems))

	dim = rowsum(vnItems) * nWaves

	/* initialize the big correlation matrix */
	corrMat = J(1, dim + 1, .)
	dims =  vnItems * nWaves // vector of total items per scale
	/* fill in the big correlation matrix */
	for(i=1;i<=cols(vnItems);i++) {
		rowBlock = J(dims[i],1, .)
		for(j=1;j<=cols(vnItems);j++) {
			if (i == j) {
				bCorrMat = J(dims[i], dims[i], scaleCorr[i, 2])
				wiCorrMat = withinItemCorr(scaleCorr[i, 1], nWaves)
				for(k=nWaves;k<=dims[i];k=k+nWaves){
					l = k - nWaves + 1
					bCorrMat[(l..k), (l..k)] = wiCorrMat
				}
			}
			else {
				if (iCorr == 0){
					bCorrMat = J(dims[i], dims[j], 0)
				}
				else {
					bCorrMat = J(dims[i], dims[j], iCorr[i,j])
				}
			}
			rowBlock =(rowBlock, bCorrMat)
		}
		corrMat = corrMat \ rowBlock
	}
	corrMat = corrMat[(2..dim +1), (2..dim +1 )]
 
	if (det(corrMat) <0) {
		corrMat = nearSPD(corrMat, tol)
		
		// ensures all elements of corrMatt are within [-1,1]
		corrMat = mm_cond(corrMat :<-1, -1, corrMat)
		corrMat = mm_cond(corrMat :>1, 1, corrMat)
	}
	st_matrix("corrMat", corrMat)
}
end

/*
cap mata: mata drop test()
mata:
void function test(real matrix corrMat, string scalar nItems, string scalar nwaves,
				   real matrix iCorr, real matrix scaleCorr) {

	real scalar dim, nWaves
	real vector vnItems, dims
	real matrix bCorrMat, rowBlock

	nWaves = strtoreal(nwaves)
	vnItems = strtoreal(tokens(nItems))
	
	dim =  vnItems * nWaves
	
	iCorr = diag(I(cols(vnItems)))
	scaleCorr =J(1, cols(vnItems), .)
	
	for(i=1;i<=cols(vnItems);i++) {  // rows
		if (i == 1) rowstart = 1 
		else rowstart = rowsum(dim[1::(i - 1)]) + 1
		
		for(j=1;j<=cols(vnItems);j++) {
			if (j == 1) colstart = 1 
			else colstart = rowsum(dim[1::(j - 1)]) + 1

			bCorrMat = corrMat[rowstart::(rowstart + dim[i] - 1), colstart::(colstart + dim[j] - 1)]
			
			if (i == j) {
				wiCorrMat = bCorrMat[1..nWaves, 1..nWaves]
				wiCorr = mean(subdiagget(wiCorrMat, 1..nWaves, 0))
				bCorr = mean(subdiagget(bCorrMat, nWaves..dim[i],0))
				scaleCorr = scaleCorr \ wiCorr, bCorr
			}
			else {
				iCorr[i,j] = mean(vec(bCorrMat))
			
			}
		}
	}
	scaleCorr[2..(cols(vnItems) + 1),1...]
	iCorr
}
end

mata: ph1 = ph2 = .
mata: test(st_matrix("corrMat"), "12 25", "3", ph1, ph2)
*/


/* Old program
mata:
void function bwCorrMat(string scalar corwi, string scalar corbi, string scalar corbs, ///
						string scalar nItems, string scalar nwaves ) {

	real scalar dim, nCorbs, nWaves, start
	real vector nCorwi, nCorbi, vnItems
	real matrix bCorrMat, corrMat
	
	nCorwi = strtoreal(tokens(corwi))
	nCorbi = strtoreal(tokens(corbi))
	nCorbs = strtoreal(corbs)
	nWaves = strtoreal(nwaves)
	vnItems = strtoreal(tokens(nItems))

	dim = rowsum(vnItems) * nWaves
	corrMat = J(dim, dim, nCorbs)
	
	start = 1
	for(i=1;i<=cols(vnItems);i++) {
		dim = vnItems[i] * nWaves
		bCorrMat = J(dim, dim, nCorbi[i])
		wiCorrMat = withinItemCorr(nCorwi[i], nWaves)
		for(j=nWaves;j<=dim;j=j+nWaves){
			k = j - nWaves + 1
			bCorrMat[(k..j), (k..j)] = wiCorrMat
		}
		corrMat[(start..(start+dim - 1)), (start..(start+dim - 1))] = bCorrMat
		start = start + dim
	}
	st_matrix("corrMat", corrMat)
}

end

*/




********************************************************************************
*** Core simulation functions
********************************************************************************

capture program drop ifeats
program define ifeats, rclass
	
	syntax namelist, Nobs(numlist) scales(string) /// 
	       [nwaves(integer -1) nitems(string) tcorr(string) ///
		    propmiss(string) wavemiss(string) psdtol(real 1e-08) /// 
		    CORRMatrix(string) MARGinals(string)]
	
	
	* scales --> scales(sc1=(0(1)4) sc2=(0(1)4)): item levels of scales
	* nwaves(integer) --> number of waves/time periods
	* nitems(string) --> nitems(sc1=12 sc2=36): number of items per scale
	* tcorr(string) --> tcorr(scl=(0.9 0.2) sc2=(0.8 0.6) corrmat=matName) : theoretical correlations (within-item between-items); corrmat: inter-scale item correlations
	* propmiss(string) --> propmiss(kzf=(0.05 0.2) hsclg=(0.05 0.3)): for random missing (item missigness; scale missingness)
	*                      propmiss(kzf=(0.2) hsclg=(0.3)): for block missing (scale missingness)
	* wavemiss(sring)  --> wavemiss(kzf=(0 1) hsclg=(1 2)): waves missing for every scale
	* corrmatrix(string) --> path and name of file containing correlation matrix for scales in namelist
	* marginals(string) --> path to directory containing files of marginal distributions for scales in namelist
	* psdtol(real)      --> tolerance in positive semi-definite matrix projection

	**** Check for data in memory
	
	*** Build a list of stub*
	local nameList ""
	foreach var of local namelist {
		local nameList "`nameList' `var'*"
	}
	capture unab allItems: `nameList'
	if (_rc == 0) {                           // get info from data
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
		
		*** retrieve number if items
		if ("`nitems'" == "") { // if number of items per scale is not specified
			*** retrieve number of items
			local nItems ""
			foreach scale of local namelist {
				unab `scale'_items: `scale'*
				local nItems "`nItems' `=`: list sizeof `scale'_items'/`nwaves''"
			}
		}
		* noi di "`nameList'"
		else {
			parse_syntax "`nitems'"
			local nItems ""
			foreach scale of local namelist {
				local nItems "`nItems' `s(`scale')'"
			}
		}
		
		**** obtain the missing empirical component
		if ("`marginals'" == "") {
			catDist `namelist', scales(`scales') saving(`c(tmpdir)')
			local marginals `"`c(tmpdir)'"'
		}
		if ("`corrmatrix'" == "") {
			tempfile corrMat
			dataCorrMat `namelist', saving(`corrMat') 
			local corrmatrix `corrMat'
		}
	}
	else {                 // if no such stubs in the dataset
		if "`marginals'" == "" local simmarg "1"
		if "`corrmatrix'" == "" local simcorr "1"
	
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
				local nItems "`nItems' `s(nItems)'"
			}
			local nwaves "`s(nwaves)'"
			restore
			* noi di "`nItems'"
			* noi di "`nwaves'"
		}	
		else {
			if ("`nitems'" ~= "") {
				parse_syntax "`nitems'"
				local nItems ""
				foreach scale of local namelist {
					local nItems "`nItems' `s(`scale')'"
				}
			}
			* noi di in y "No dataset in memory, no correlation matrix or marginals given" 
			else if ("`nwaves'" == "-1" | "`nitems'" == "") {
				no di in r "Error: Arguments nwaves and nitems are required in this case"
				exit 1000
			}
		}
	}
	* noi di "`nItems'"
	* noi di "`nwaves'"
	
	if ("`simmarg'" ~= "") { 	
		*** Create even distribution. 
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
		ifeatsCore,  namelist("`namelist'") nobs(`lobs') cols(`cols') rows(`rows') nwaves(`nwaves') nItems("`nItems'") ///
		tcorr(`tcorr') propmiss(`propmiss') wavemiss(`wavemiss') simcorr(`simcorr')  simmarg(`simmarg')  ///
		corrmatrix("`corrmatrix'") marginals("`marginals'") scales("`scales'") itemlist("`itemNames'") psdtol(`psdtol')

		local ++cols
	}
		
end



capture program drop ifeatsCore
program define ifeatsCore
 
* 	args nobs corb cols rows nwaves corw propmiss mblock simcorr simmarginals corrmatrix marginals namelist nTotalItems wavemiss

	syntax, [namelist(string) nobs(integer 1) corb(real 1) cols(integer 1) rows(integer 1) nwaves(integer -1) nItems(string) ///
			 tcorr(string) propmiss(string) wavemiss(string) simcorr(string) simmarg(string) ///
			 corrmatrix(string) marginals(string) scales(string) itemlist(string) psdtol(real 1)]

	
	clear	
	qui set obs `nobs'
	
	
	**** Display basics
	if ("`simcorr'" == "1") {
		noi di _n in y "************************************************************"
		noi di in y    "* Sample size        : " `nobs'
		noi di in y    "* Using simulated correlation matrix"	
		if ("`simmarg'" ~= "") {
			noi di in y	   "* Using simulated marginal distributions"
		}
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
	
	*** total number of items across all scales
	
	local nTotalItems `=`=(`=subinstr("`nItems'", " ", "+",.)')' * `nwaves''
	
	*** Correlation matrix
	
	if ("`simcorr'" ~= "") { // simulate a correlation matrix	  
		noi di "Simulating correlation matrix..." 
		noi mata: bwCorrMat("`tcorr'", "`nItems'", "`nwaves'", `psdtol')
	}
	else { // load matrix of empirical correlations
		preserve
		use "`corrmatrix'", clear
		mkmat _all, matrix(corrMat)
		restore
		noi di "Empirical correlation matrix loaded"
	}
	
	*** Means
	mata: means(`nTotalItems')

	*** Simulate random normal data using the (empirical) distributions and corr mat
	drawnorm myvar1-myvar`nTotalItems', means(means) corr(corrMat) n(`=_N')
		
	if ("`simmarg'" ~= "") {  // simulated marginals
		noi di "Simulating marginal distributions... " _c
		
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
		noi di " done"
	} // end of if
	
	else { // Load and use observed marginals	
		*** Read in scale_cdist's as matrices
		noi di "Empirical marginal distributions loaded..." _c
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
		drop missTime?Point  // does not work well with tempvars
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
	capture noisily pchained `namelist', i(id) t(time) mio(add(1) burnin(10) chaindots)
	
	if _rc ~= 0 {
		noi di in r "Failed"
	}
	
	*** Store results
	mat simmat[`rows', `cols'] = _rc

end

