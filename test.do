********************************************************************************
*** This is the testing script
**
*
*
*
*
*
*


clear all
set more off

set maxvar 32000

*** Point to directories
local data_folder "~/Desktop/imputations/data"
local script_folder "~/Desktop/gitProjects/ifeats"
local output_folder "~/Desktop/imputations/temp-data"


*** Load data and programmes
use "`data_folder'/test-data", clear
do "`script_folder'/scales.ado"
do "`script_folder'/newvcov.ado"


* scaleclean kzf* hsclg*
* mypercent kzf saving("`output_folder'") 
catdist kzf hsclg


* use kzf_complete, clear
* use hsclg_complete, clear


*order _all, alpha
*order levels, first

exit

********************************************************************************
**** Use Extracted information

**** Test for binary variable. 
/*
local binavars "ka_ecdm*"
**** marginal distribution
**** Create a subfolder called percents under the `saving' folder. 
mypercent `binavars', saving("$folder") // Here the saving option standard input should be a folder.  
**** covariance matrix
myvcov `binavars', saving("$folder/empiricalvcov.dta")
mysim, simobs(50(50)100) nwitems(3) ntitems(33) propmiss(0.2) mblock(1) simvcov(0) storedvcov("$folder/empiricalvcov.dta")
*/




*******************************
**** Test for Categorical variable 

/*
local catevars "kzf*"
qui unab inputvars: `catevars'
mypercent `catevars', saving("`output_folder'") 
myvcov `catevars', saving("`output_folder'/empiricals/empiricalvcov.dta") 
catsim, simobs(50(50)100) nwitems(3) ntitems(36) propmiss(0.2) mblock(1) simvcov(0) ///
        simmarginal(0) storedvcov("`output_folder'/empiricals/empiricalvcov.dta") /// 
		storedmarginal("`output_folder'/empiricals/") storedvars("`inputvars'")

*/



/*		
*** Exracting the empirical vcov of ka_ecdm* and saving it to a file
 myvcov kzf*, saving("`output_folder'/categorivcov.dta")
 mypercent kzf*, saving("`output_folder'")  // , saving("`temp_data'/catepercent.dta")

**** By Zitong: Think about categorical variables. 
* mysim, simobs(50(50)100) nwitems(3) ntitems(33) propmiss(0.2) mblock(1) simvcov(0) storedvcov("`temp_data'/categorivcov.dta")
*/






