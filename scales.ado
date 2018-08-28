
********************************************************************************
**** Include 3 programs: scaleclean, mypercent and cregrid (not yet created!). 
**** The third used in simulation(?)
**** Developer: Simo Goshev, Zitong Liu
**** Created: Week 2, 2018
**** cregrid is not created now. 

********************************************************************************
**** Scale cleanings: used in simulation and imputation

* First step in simulation and imputation. 



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







********************************************************************************
**** PROGRAM mypercent help you to get the percentages for each categories 
**** Instructions:
**** 1. scaleclean is called. 
**** 2. Outputs are stored in `saving'/percents folder. 
**** 3. Works for both binary and categorical variables. 
**** 4. Generate CUMULATIVE PERCENTAGES.  When used in simulation, need to drop the last line. 
********************************************************************************

*** SG: This is a bit wonky

capture program drop mypercent
program define mypercent
	
	syntax varlist, saving(string)

	*** SG: We can send this to a temp file later on
	capture mkdir "`saving'/empiricals"
	local saving "`saving'/empiricals"
	
	*** Store scale items and associated labels in r()
	qui scaleclean `varlist'

	*** Get the labels of the scales from r()
	local scales `r(scales)'
	*** Get the items of scales from r()
	foreach sca of local scales {
		local subsca_`sca' `r(subs`sca')'
	}

	
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
** Re-write of mypercent

*******
** - check levels of each item and take union
** - check label and compare with union, take union
** - store item distributions

********************************************************************************

capture program drop catdist
program define catdist

	syntax namelist [, saving(string)]
	
	local nscales: word count `namelist'
	foreach scale of local namelist { //iterate over scales
		unab scale_items: `scale'*
		* di "`scale_items'"
		
		*** If a label exists:  // FIX THIS!
		*** Retrieve levels from label
		local labname: value label `: word 1 of `scale_items''
		mata: values = .; text = ""
		mata: st_vlload("`labname'", values, text); _transpose(values)
		mata: st_local("labs", invtokens(strofreal(values)))
		
		*** Retrieve levels from items
		local all_levels ""
		local counter = 1
		foreach item of local scale_items { // iterate over items of a scale
			qui levelsof `item', local(levs)   // get the levels of the item
			local missLevs: list labs - levs   // get the missing levels
			preserve
			contract `item', cpercent(cum_`item') nomiss  // get cum distribution
			foreach lev of local missLevs { // add the missing levels
				set obs `=_N + 1'
				replace `item' = `lev' if _n == _N
				replace cum_`item' = 1e-5 if _n == _N
			}
			ren (`item' cum_`item') (cats `item')
			drop _freq
			
			*** Merge item files
			if (`counter' == 1) {
				save `scale'_complete, replace // USE TEMP FILES
			}
			else {
				merge 1:1 cats using `scale'_complete, nogen
				save `scale'_complete, replace
			}
			local ++counter
			restore
		}
	}
end

	
	



********************************************************************************
**** TODO: PROGRAM cregrid is a functional program for 
capture program drop cregrid



