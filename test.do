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

*** Using observed distributions
ifeats kzf hsclg, nobs(50(50)100) nwitems(3) ntitems(36) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1))  simcorr(0) ///
        simmarginals(0)  corrmatrix("`output_folder'/empirCorrMat.dta") /// 
		marginals("`output_folder'")  mblock(kzf hsclg) nwavemiss(kzf=(0 1) hsclg=(1 2))

*** Simulating
ifeats kzf hsclg, nobs(50(50)100) nwitems(3) ntitems(36) propmiss(kzf=0.2 hsclg=0.3)  simcorr(0) ///
       simmarginals(1)  corrmatrix("`output_folder'/empirCorrMat.dta") /// 
	   marginals("`output_folder'")   simscales(kzf=(0(1)4) hsclg=(1(1)4)) 		
		
		
exit
	
* (error 900: maxvar issue; no room to add more variables)
		



gen test0 = missing(kzf01_tp0, kzf02_tp0, kzf03_tp0, kzf04_tp0, kzf05_tp0, ///
				   kzf06_tp0, kzf08_tp0, kzf10_tp0, kzf11_tp0, kzf12_tp0)
gen test2 = missing(kzf01_tp2, kzf02_tp2, kzf03_tp2, kzf04_tp2, kzf05_tp2, ///
				   kzf06_tp2, kzf08_tp2, kzf10_tp2, kzf11_tp2, kzf12_tp2)
				   
tab test0, missing
tab test2, missing
tab kzf01_tp0, missing
tab kzf01_tp2, missing
drop test*



gen test0 = missing(hsclg01_tp0, hsclg02_tp0, hsclg03_tp0, hsclg04_tp0, ///
					hsclg05_tp0, hsclg06_tp0, hsclg07_tp0, hsclg08_tp0, ///
					hsclg09_tp0, hsclg10_tp0, hsclg11_tp0, hsclg12_tp0, ///
					hsclg13_tp0, hsclg14_tp0, hsclg15_tp0, hsclg16_tp0, ///
					hsclg17_tp0, hsclg18_tp0, hsclg19_tp0, hsclg20_tp0, ///
					hsclg21_tp0, hsclg22_tp0, hsclg23_tp0, hsclg24_tp0, hsclg25_tp0)
gen test1 = missing(hsclg01_tp1, hsclg02_tp1, hsclg03_tp1, hsclg04_tp1, ///
					hsclg05_tp1, hsclg06_tp1, hsclg07_tp1, hsclg08_tp1, ///
					hsclg09_tp1, hsclg10_tp1, hsclg11_tp1, hsclg12_tp1, ///
					hsclg13_tp1, hsclg14_tp1, hsclg15_tp1, hsclg16_tp1, ///
					hsclg17_tp1, hsclg18_tp1, hsclg19_tp1, hsclg20_tp1, ///
					hsclg21_tp1, hsclg22_tp1, hsclg23_tp1, hsclg24_tp1, hsclg25_tp1)
tab test0, missing
tab test1, missing
tab hsclg01_tp0, missing
tab hsclg01_tp2, missing
drop test*


/*
* mblock(kzf hsclg)

*order _all, alpha
*order levels, first

propmiss(kzf=0.2 hsclg=0.3) does not work with only one scale:
propmiss(hsclg=0.3)
*/


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






