********************************************************************************
*** Utility function which is used for pre-processing (cleaning) of the scales
*
*
*
*
*
*


********************************************************************************
*** Utility function for scale cleaning 
** - filter scales from non-numeric items
** - deals with missing values
** - recodes levels to start from 0
** - define items that are part of a scale
** - reports filtered items of a scale and associated label
********************************************************************************


capture program drop scaleclean
program define scaleclean, rclass

	syntax varlist 
	
	**** First, there are some variables that are "complementary information",
	**** like "hsclg20_a: CG indicated suicidal thoughts; Name of mental health services referral."
	**** When we drop these type of variables from our varlist and only work only on numerical variables. 
	
	**** Check for non-numeric variables
	local nonum "" // list of non-numberic items
	foreach var of varlist `varlist' {
		local typevar: type `var'     // get type
		*** If type is not one of byte, float or int, add to list 
		if !(("`typevar'" == "byte")|("`typevar'" =="float")|("`typevar'" =="int")) {
			local nonum = "`nonum' `var'"
		}
	}
	if "`nonum'"!= "" {
		noi di("These items are non-numeric and will be dropped: `nonum'")
		local varlist: list varlist - nonum
	}
	
	***** Order my varlist	
	order `varlist', alpha
	
	
	********************************************
	*** Check for missing
	*** SG: This is tricky! We can perhaps look at the label and tease this out?
	*** There are other ways of doing this as well -- let's not worry about this yet
	
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

	************************************************
	*** SG: Probably do not want to touch the original vars! 
	**** Recode (temporarily) internally all argument variables to have values from zero. 
	**** Those have "1 2 3 4" will be recoded to "0 1 2 3".
	**** Ask Simo: do we want to keep this? 
	
	qui foreach var of varlist `varlist' {
		sum `var' 
		if (`r(min)' ~= 0) {
			**** Make it start from zeros 
			replace `var' = `var' - `r(min)'
		}
	}
	
	*** Produce scale mamber items and associated label
	
	**** Short Summary over the Scale Type. 
	***** Keyword: grid creation. 

	label dir 
	local labname `r(names)'
	local hasscale "" // hasscale stores the scales appears in the argument varlist. 
	
	foreach lab of local labname {
		* label list `lab'
		* local value`lab' `r(k)'
		di "`lab'"
		**** Search and list all argument variables that has each labels. If some label is not used by the input argument. 
		ds `varlist', has(vall `lab')
		local subsca_`lab' = "`r(varlist)'"   
		di "`subsca_`lab''"
		
		if !("`subsca_`lab''" == "")  local hasscale : list hasscale | lab 
	}
	

	local numscale: list sizeof hasscale
	noi di "Argument subscales from `numscale' scales. " _newline ///
	"They are: `hasscale'. "
	
	foreach lab of local hasscale {
		noi di "Scale `lab' with subscale variables: `subsca_`lab''. "
	}

	**** Return scale items and associated label
	
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
