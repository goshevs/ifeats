
********************************************************************************
**** Include 3 programs: scaleclean, mypercent and cregrid. 
**** The third used in simulation(?)
**** Author: Zitong Liu
**** Created: Week 2, 2018
**** cregrid is not created now. 

********************************************************************************
**** Scale cleanings: used in simulation and imputation

* First step in simulation and imputation. 


capture program drop scaleclean
program define scaleclean, rclass

	syntax varlist 
	
	**** First, there are some variables that are "complementary information",
	**** like "hsclg20_a: CG indicated suicidal thoughts; Name of mental health services referral."
	**** When we drop these type of variables from our varlist and only work only on numerical variables. 

	local nonum ""
	foreach var of varlist `varlist' {
	local typevar: type `var'
	if !(("`typevar'" == "byte")|("`typevar'" =="float")) {
	local nonum = "`nonum'" +  " `var'"
	}
	}
	if "`nonum'"!= "" {
	noi di("These variables are not numeric scales and will be ignored: `nonum'")
 	local varlist: list varlist - nonum
	}
		
	***** Order my varlist	
	order `varlist', alpha
	
	**** Some people say "don't know", it is coded as 99 or .d, I change these answers into missing values to be imputated. 
	**** Then I remove the "Don't Know" Label.
	foreach var of varlist `varlist'{
	qui replace `var' = . if (`var' == 99 ) | (`var' == .d)
	}
	foreach var of varlist `varlist' {
	local scaletype_`var': value label `var'
	local lab99: label `scaletype_`var'' 99, strict
	if !("`lab99'"=="") {
	label define `scaletype_`var'' 99 ., modify
	}
	local labd: label `scaletype_`var'' .d, strict
	if !("`labd'"=="") {
	label define `scaletype_`var'' .d ., modify
	}
	}

	**** Recode (temporarily) internally all argument variables to have values from zero. 
	**** Those have "1 2 3 4" will be recoded to "0 1 2 3".
	**** Ask Simo: do we want to keep this? 
	
	foreach var of varlist `varlist' {
	qui{
	sum `var' 
	if !`r(min)' {
	**** Make it starts from zeros 
	replace `var' = `var' - `r(min)'
	}
	}
	}
	
	**** Short Summary over the Scale Type. 

	***** Keyword: grid creation. 

	qui{
	label dir 
	local labname `r(names)'
	local hasscale "" // hasscale stores the scales appears in the arugment varlist. 
	
	foreach lab of local labname {
	* label list `lab'
	* local value`lab' `r(k)'
	**** Search and list all argument variables that has each labels. If some label is not used by the input argument. 
	qui ds `varlist', has(vall `lab')
	local subsca_`lab' = "`r(varlist)'"   
	if !("`subsca_`lab''" == "")  local hasscale : list hasscale | lab 
	}
	}
	
	local numscale: list sizeof hasscale
	noi di "Argument subscales from `numscale' scales. " _newline ///
	"They are: `hasscale'. "
	
	foreach lab of local hasscale {
	noi di "Scale `lab' with subscale variables: `subsca_`lab''. "
	}
	
	
	**** Return
	
	return local scales `"`hasscale'"'
	foreach lab of local hasscale {
	return local subs`lab' `"`subsca_`lab''"'
	}
	**** Then we can focus on using label list `hasscale' instead of `labname'.
	
	**** This part is creating categorical "storage" for each type of scales. 
	/*
	foreach var of varlist `varlist' {
	local scaletype_`var': value label `var'
	}
	*/

end



********************************************************************************
**** PROGRAM mypercent help you to get the percentages for each categories 
**** Instructions:
**** 1. scaleclean is called. 
**** 2. Outputs are stored in `saving'/percents folder. 
**** 3. Works for both binary and categorical variables. 
**** 4. Generate CUMULATIVE PERCENTAGES.  When used in simulation, need to drop the last line. 
********************************************************************************



capture program drop mypercent
	program define mypercent
	
	syntax varlist, saving(string)

	capture mkdir "`saving'/empiricals"
	local saving "`saving'/empiricals"
	
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
* 	drop if _n == _N // Remember you changed here!  
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
*	tempfile sca`ii'perc
	use "`perc_l`ii'v1'", clear
	forval jj = 2/`j_`ii'' {
	merge 1:1 cate_perc using "`perc_l`ii'v`jj''"
	drop _merge
	}
	recode _all (missing=1e-5)
* 	drop cate_perc
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
**** PROGRAM cregrid is a functional program for 


capture program drop cregrid



