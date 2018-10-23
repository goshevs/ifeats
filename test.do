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
*catDist kzf hsclg, scales(kzf=(0(1)4) hsclg=(1(1)4)) saving("`output_folder'")

*** Obtain the empirical correlation matrix of the data
* dataCorrMat kzf hsclg, saving("`output_folder'/empirCorrMat.dta") 


/*
***** Options for ifeats
	* scales --> scales(sc1=(0(1)4) sc2=(0(1)4)): item levels of scales (required)
	* nwaves(integer) --> number of waves/time periods (required for full simulation)
	* nitems(string) --> nitems(sc1=12 sc2=36): number of items per scale (required for full simulation)
	* propmiss(string) --> propmiss(kzf=(0.05 0.2) hsclg=(0.05 0.3)): for random missing (item missigness; scale missingness)
	*                      propmiss(kzf=(0.2) hsclg=(0.3)): for block missing (scale missingness)
						   propmiss(0.2): block missing for entire dataset (required for full simulation)
	* wavemiss(sring)  --> wavemiss(kzf=(0 1) hsclg=(1 2)): waves missing for every scale
						   wavemiss(minmax(1 2)): min and max number of missing waves (required for full simulation)
	* corrmatrix(string) --> path and name of file containing correlation matrix for scales in namelist
	* marginals(string) --> path to directory containing files of marginal distributions for scales in namelist

*/


********************************************************************************
*** Using data in memory

**** Complete syntax
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(kzf=(0 1) hsclg=(1 2))
        
		
**** Randomly picking missing waves (random pattern)
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4))

**** Block missing pattern per variable (mixed pattern also possible)	
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=0.3 hsclg=0.1) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4))

**** Block missing pattern over all specified scales 
ifeats kzf hsclg, nobs(50(50)100) propmiss(0.3) wavemiss(minmax(1 2)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4))




********************************************************************************
*** Using data in memory and previously stored corr matrix and marginals

**** Complete syntax
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(kzf=(0 1) hsclg=(1 2)) ///
        corrmatrix("`output_folder'/empirCorrMat.dta") marginals("`output_folder'") 
		
**** Randomly picking missing waves (random pattern)
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) ///
        corrmatrix("`output_folder'/empirCorrMat.dta") marginals("`output_folder'") 

**** Block missing pattern per variable (mixed pattern also possible)	
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=0.3 hsclg=0.1) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) ///
        corrmatrix("`output_folder'/empirCorrMat.dta") marginals("`output_folder'")

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
		scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(kzf=(0 1) hsclg=(1 2)) ///
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
		

*=== marginals only ===* DOES NOT WORK AT THIS TIME
		
**** Complete syntax
clear
capture ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(kzf=(0 1) hsclg=(1 2)) ///
        marginals("`output_folder'")
		
**** Randomly picking missing waves (random pattern)
**** Block missing pattern per variable	(mixed pattern also possible)
**** Block missing pattern over all specified scales
		


********************************************************************************
*** Full simulation (simulating everything)

**** Complete syntax
clear
ifeats kzf hsclg, nobs(50(50)100) nwaves(3) nitems(kzf=12 hsclg=25) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
      scales(kzf=(0(1)4) hsclg=(0(1)4)) wavemiss(kzf=(0 1) hsclg=(1 2)) //wavemiss could be skipped

**** Block missing pattern per variable	(mixed pattern also possible)
clear
ifeats kzf hsclg, nobs(50(50)100) nwaves(3) nitems(kzf=12 hsclg=25) propmiss(kzf=(0.3) hsclg=(0.1)) ///
      scales(kzf=(0(1)4) hsclg=(0(1)4)) wavemiss(kzf=(0 1) hsclg=(1 2)) //wavemiss could be skipped
	  
**** Block missing pattern over all specified scales
clear
ifeats kzf hsclg, nobs(50(50)100) nwaves(3) nitems(kzf=12 hsclg=25) propmiss(0.2) ///
      scales(kzf=(0(1)4) hsclg=(0(1)4)) wavemiss(minmax(1 2)) 	    
	  
exit
	  
		

		
		
***************** END OF FILE **************************************************
		
		

		
		
		
		
		
		
		
*** Simulating
ifeats kzf hsclg, nobs(50(50)100) nwitems(3) ntitems(36) propmiss(kzf=(0.05 0.2) hsclg=(0.05 0.3))  simcorr(0) ///
       simmarginals(1)  corrmatrix("`output_folder'/empirCorrMat.dta") /// 
	   marginals("`output_folder'") simscales(kzf=(0(1)4) hsclg=(1(1)4))
	   
*/


*** New syntax

* catDist kzf, scales(kzf=(0(1)4)) saving("`output_folder'")
* dataCorrMat kzf hsclg, saving("`output_folder'/empirCorrMatKZF.dta") 

/*
ifeats, nobs(500(500)1000) nwaves(3) propmiss(kzf=(0.05 0.2)) ///
       corrm("`output_folder'/empirCorrMatKZF.dta") /// 
	   marg("`output_folder'") scales(kzf=(0(1)4) hsclg=(0(1)4))

clear
ifeats kzf hsclg, nobs(100(100)200) nwaves(3) propmiss(kzf=(0.05 0.2) hsclg=(0.3)) ///
	   corrm("`output_folder'/empirCorrMat.dta") /// 
	   scales(kzf=(0(1)4) hsclg=(0(1)4))

clear
ifeats kzf hsclg, nobs(500(500)1000) propmiss(kzf=(0.05 0.2)) ///
       marg("`output_folder'/") scales(kzf=(0(1)4) hsclg=(0(1)4))

clear
ifeats kzf hsclg, nobs(500(500)1000) propmiss(kzf=(0.05 0.2)) ///
       scales(kzf=(0(1)4) hsclg=(0(1)4))

*/
  

	   

	   
exit
	   

	   
exit

	
ifeats kzf hsclg, nobs(500(500)1000) propmiss(kzf=(0.05 0.2)) ///
       corrm("`output_folder'/empirCorrMat.dta") /// 
	   marg("`output_folder'/") scales(kzf=(0(1)4) hsclg=(0(1)4))
	   
exit
	   
ifeats kzf hsclg, nobs(500(500)1000) nwaves(3) propmiss(kzf=(0.05 0.2)) ///
       corrm("`output_folder'/empirCorrMatKZF.dta") /// 
	   marg("`output_folder'") scales(kzf=(0(1)4) hsclg=(0(1)4))
	   
		*/   
	   
	   
clear
ifeats kzf hsclg, nobs(500(500)1000) nwaves(3) nitems(kzf=12 hsclg=25) propmiss(kzf=(0.05 0.2)) ///
       scales(kzf=(0(1)4) hsclg=(0(1)4))

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






