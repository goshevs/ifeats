********************************************************************************
*** Utility functions for diagnosing problems with the data
*
*
*
*
*
*
*


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
