********************************************************************************
*** Create cummulative distributions for each scale item
*
*
*
*
*

********************************************************************************
** Here we:
** - retrieve the user input for each scale
** - check levels of each item and take union
** - check label and compare with user
** - pick user but alert user for discrapancies
** - store item distributions

********************************************************************************


*** Assume a clean dataset, no problems with variables

capture program drop catDist
program define catDist

	syntax namelist, scales(string) [saving(string)]

	qui { 
		local nscales: word count `namelist'

		*** Parse the scales string: scales(kzf=(numlist) hsclg=(numlist))
		*** Retrieve scale name and scale levels from user input
		if "`scales'" ~= "" {
			parse_scales "`scales'"
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

	
*** Parser of the user input for scale values
capture program drop parse_scales
program define parse_scales, sclass
	args myscales
	
	foreach scale of local myscales {
		gettoken sname 0: scale, parse("=")
		gettoken left levs: 0, parse("=")
		local levs = substr("`levs'", 2, `=length("`levs'") - 2')
		numlist "`levs'"
		sreturn local `sname' `r(numlist)'
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
