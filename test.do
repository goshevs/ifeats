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

*** Point to directories depending on user
if (c(username) == "goshev") {
	local data_folder "~/Desktop/imputations/data"
	local script_folder "~/Desktop/gitProjects/ifeats"
	local output_folder "~/Desktop/imputations/temp-data"
}

else {
	local data_folder "/Users/zitongliu/Dropbox/2018/BTProject/ifeats_old"
	local script_folder "/Users/zitongliu/Dropbox/2018/BTProject/ifeats"
	local output_folder "/Users/zitongliu/Dropbox/2018/BTProject/ifeats_old"
	cd "`output_folder'"
}


*** Load data and programmes
use "`data_folder'/test-data", clear
qui do "`script_folder'/scales.ado"
qui do "`script_folder'/core-programs.ado"


*** Load pchained from GitHub
qui do "https://raw.githubusercontent.com/goshevs/pchained-github/master/pchained.ado"


*** Create the cummulative distributions of all categorical items
* catDist kzf hsclg, scales(kzf=(0(1)4) hsclg=(1(1)4)) saving("`output_folder'")

*** Obtain the empirical correlation matrix of the data
* dataCorrMat kzf hsclg, saving("`output_folder'/empirCorrMat.dta") 



/*
***** Options for ifeats
nobs        : sample sizes
simcorb     : between-item correlation
nwitems     : number of within items (essentially the number of time periods)
ntitems     : number of total items (time periods * number of indicators)
corw        : within-item correlation (over time)
propmiss    : proportion of missing (based on obs)
mblock      : block missing (boolean); all items missing in a period
simcorr     : simulate corrMat or use empirical corrMat (boolean)
storedcorr  :  if simvcov(0) then location of stored empirical corrMat 
simmarginal :
storedmarginal
storedvars

*/


ifeats kzf hsclg, nobs(50(50)100) nwitems(3) ntitems(36) propmiss(kzf=0.2 hsclg=0.3)  simcorr(0) ///
       simmarginals(1)  corrmatrix("`output_folder'/empirCorrMat.dta") /// 
	   marginals("`output_folder'")   simscales(kzf=(0(1)4) hsclg=(1(1)4)) 

exit	

nwavemiss(kzf=(0 1) hsclg=(1 2)) mblock(kzf hsclg)

*order _all, alpha
*order levels, first



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






