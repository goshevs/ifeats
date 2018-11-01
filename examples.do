********************************************************************************
*** Examples
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

snapshot save

*** Load pchained from GitHub
qui do "https://raw.githubusercontent.com/goshevs/pchained-github/master/pchained.ado"


*** Create the cummulative distributions of all categorical items
catDist kzf hsclg, scales(kzf=(0(1)4) hsclg=(1(1)4)) saving("`output_folder'")

*** Obtain the empirical correlation matrix of the data
dataCorrMat kzf hsclg, saving("`output_folder'/empirCorrMat.dta") 


********************************************************************************
*** Using data in memory

**** Complete syntax
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(kzf=(0 1) hsclg=(1 2))
		
snapshot restore 1
**** Randomly picking missing waves (random pattern)
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4))

snapshot restore 1	
**** Block missing pattern per variable (mixed pattern also possible)	
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=0.3 hsclg=0.1) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4))

snapshot restore 1
**** Block missing pattern over all specified scales 
ifeats kzf hsclg, nobs(50(50)100) propmiss(0.3) wavemiss(minmax(1 2)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4))




********************************************************************************
*** Using data in memory and previously stored corr matrix and marginals

snapshot restore 1
**** Complete syntax
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(kzf=(0 1) hsclg=(1 2)) ///
        corrmatrix("`output_folder'/empirCorrMat.dta") marginals("`output_folder'") 

snapshot restore 1
**** Randomly picking missing waves (random pattern)
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) ///
        corrmatrix("`output_folder'/empirCorrMat.dta") marginals("`output_folder'") 

**** Block missing pattern per variable (mixed pattern also possible)	
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=0.3 hsclg=0.1) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) ///
        corrmatrix("`output_folder'/empirCorrMat.dta") marginals("`output_folder'")

snapshot restore 1
**** Block missing pattern over all specified scales
ifeats kzf hsclg, nobs(50(50)100) propmiss(0.3) wavemiss(minmax(1 2)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) ///
        corrmatrix("`output_folder'/empirCorrMat.dta") marginals("`output_folder'")



********************************************************************************
*** Using previously stored corr matrix/marginals, simulating what has not been provided

*=== corr matrix and marginals ===*

**** Complete syntax
clear
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(kzf=(0 1) hsclg=(1 2)) ///
        corrmatrix("`output_folder'/empirCorrMat.dta") marginals("`output_folder'")

**** Randomly picking missing waves (random pattern)
clear
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) ///
        corrmatrix("`output_folder'/empirCorrMat.dta") marginals("`output_folder'")
		
**** Block missing pattern per variable	(mixed pattern also possible)	
clear
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.3) hsclg=(0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) ///
        corrmatrix("`output_folder'/empirCorrMat.dta") marginals("`output_folder'")
		
**** Block missing pattern over all specified scales
clear
ifeats kzf hsclg, nobs(50(50)100) propmiss(0.3) wavemiss(minmax(1 2)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) ///
        corrmatrix("`output_folder'/empirCorrMat.dta") marginals("`output_folder'")
		

*=== corr matrix only ===*

**** Complete syntax
clear
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(kzf=(0 1) hsclg=(1 2))  ///
        corrmatrix("`output_folder'/empirCorrMat.dta")

**** Randomly picking missing waves (random pattern)
clear
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) ///
        corrmatrix("`output_folder'/empirCorrMat.dta")
		

**** Block missing pattern per variable (mixed pattern also possible) 
clear
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.3) hsclg=(0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4))  ///
        corrmatrix("`output_folder'/empirCorrMat.dta")
		
**** Block missing pattern over all specified scales
clear
ifeats kzf hsclg, nobs(50(50)100) propmiss(0.3) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(minmax(1 2))  ///
        corrmatrix("`output_folder'/empirCorrMat.dta")
		


*=== marginals only ===*
		
**** Complete syntax
clear
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(kzf=(0 1) hsclg=(1 2)) ///
        marginals("`output_folder'") tcorr(kzf=(0.8 0.2) hsclg=(0.8 0.2)) 

		
**** Randomly picking missing waves (random pattern)
clear
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) ///
        marginals("`output_folder'") tcorr(kzf=(0.8 0.2) hsclg=(0.8 0.2)) 

		
**** Block missing pattern per variable	(mixed pattern also possible)
clear
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.3) hsclg=(0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) ///
        marginals("`output_folder'") tcorr(kzf=(0.8 0.2) hsclg=(0.8 0.2)) 

		
**** Block missing pattern over all specified scales
clear
ifeats kzf hsclg, nobs(50(50)100) propmiss(0.3) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(minmax(1 2)) ///
        marginals("`output_folder'") tcorr(kzf=(0.8 0.2) hsclg=(0.8 0.2)) 

		
********************************************************************************
*** Full simulation (simulating everything)

**** Complete syntax
mat myCorrMat = (1, 0.5 \ 0.5, 1) // inter-scale correlation matrix (at item level)

clear
ifeats kzf hsclg, nobs(50(50)100) nwaves(3) nitems(kzf=12 hsclg=25) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
      scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(kzf=(0 1) hsclg=(1 2)) ///
	  tcorr(kzf=(0.9 0.2) hsclg=(0.8 0.6) corrmat=myCorrMat) // if corrmat not specified --> inter-scale corr = 0

**** Randomly picking missing waves (random pattern)
clear
ifeats kzf hsclg, nobs(50(50)100) nwaves(3) nitems(kzf=12 hsclg=25) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
      scales(kzf=(0(1)4) hsclg=(1(1)4)) ///
	  tcorr(kzf=(0.9 0.2) hsclg=(0.8 0.6)) 

**** Block missing pattern per variable	(mixed pattern also possible)
clear
ifeats kzf hsclg, nobs(50(50)100) nwaves(3) nitems(kzf=12 hsclg=25) propmiss(kzf=(0.3) hsclg=(0.1)) ///
      scales(kzf=(0(1)4) hsclg=(1(1)4)) ///
	  tcorr(kzf=(0.9 0.2) hsclg=(0.8 0.6) corrmat=myCorrMat) 
	  
**** Block missing pattern over all specified scales
clear
ifeats kzf hsclg, nobs(50(50)100) nwaves(3) nitems(kzf=12 hsclg=25) propmiss(0.2) ///
      scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(minmax(1 2)) ///   
	  tcorr(kzf=(0.9 0.2) hsclg=(0.8 0.6))
	
