****** 
clear all

set more off

set maxvar 32000

global folder "/Users/zitongliu/Dropbox/2018/BTProject/label_manu_2"

use "$folder/cg_data_pre_imputation_07_10_2018.dta", clear

adopath + "$folder"

do "$folder/scales.ado"
do "$folder/newvcov.ado"

* adopath + "C:/Users/goshev/Desktop/gitProjects/plumpton"


snapshot erase _all

order prelim_id2, first
order ratimepoint, after(prelim_id2)
sort prelim_id2 ratimepoint



********************************************************************************
*** Filling in the missing observations

tsset prelim_id2 ratimepoint
tsfill, full

*** For all cov_fixed
*** Assuming all variables have been encoded
foreach var of varlist rct_arm base_cg_age cgrel base_cg_educat family_structure {
	sort prelim_id2 `var' // this could take a long time too; is there a better way?
	bys prelim_id2: replace `var' = `var'[1]
}
sort prelim_id2 ratimepoint
tsset, clear



snapshot save, label(originalData) //snapshot 1



********************************************************************************
** Limit the dataset to variables that we would need
********************************************************************************

********************************************************************************
*** Create a common stub for the scale to facilitate handling (best practice!)

ren (ka_ecdm3 ka_ecdm4 ka_ecdm5 ka_ecdm11 ka_ecdm12)(ka_o_ecdm3 ka_o_ecdm4 ka_o_ecdm5 ka_o_ecdm11 ka_o_ecdm12)
ren (ka_ecdm3_r ka_ecdm4_r ka_ecdm5_r ka_ecdm11_r ka_ecdm12_r) (ka_ecdm3 ka_ecdm4 ka_ecdm5 ka_ecdm11 ka_ecdm12)


********************************************************************************
*** Create a common stub for the scale to facilitate handling (best practice!)

ren (dmn1 dmn2 dmn3 dmn4 dmn5 dmn6 dmn7) (dmn_o_1 dmn_o_2 dmn_o_3 dmn_o_4 dmn_o_5 dmn_o_6 dmn_o_7)
ren (dmn1_share dmn2_share dmn3_share dmn4_share dmn5_share dmn6_share dmn7_share) (dmn1 dmn2 dmn3 dmn4 dmn5 dmn6 dmn7)


********************************************************************************
*** Keep only relevant vars

keep prelim_id2 ratimepoint ka_ecdm* dmn1 dmn2 dmn3 dmn4 dmn5 dmn6 dmn7 ///
     kzf1 kzf2 kzf3 kzf4 kzf5 kzf6 kzf8 kzf10 kzf11 kzf12 kzf14p kzf17p /// 
     hsclg1-hsclg10 hsclg11-hsclg19 hsclg20 hsclg21-hsclg25


 
foreach var of varlist ka_ecdm? dmn? kzf? hsclg? {
	if regexm("`var'", "[0-9]$") {
		local lastdigit = substr("`var'",`=strlen("`var'")', 1)
		ren `var' `=substr("`var'",1, `=strlen("`var'")-1')'0`lastdigit'
	}
}
	 
 
order _all, alpha
order prelim_id2, first
order ratimepoint, after(prelim_id2)
sort prelim_id2 ratimepoint
	 


	 
********************************************************************************
*** Reshape to wide

foreach var of varlist ka_ecdm* dmn* kzf* hsclg* {
	ren `var' `var'_tp
}

reshape wide ka_ecdm* dmn* kzf* hsclg*, i(prelim_id2) j(ratimepoint)




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

local catevars "kzf*"
qui unab inputvars: `catevars'
mypercent `catevars', saving("$folder") 
myvcov `catevars', saving("$folder/empiricals/empiricalvcov.dta") 
catsim, simobs(50(50)100) nwitems(3) ntitems(36) propmiss(0.2) mblock(1) simvcov(0) simmarginal(0) storedvcov("$folder/empiricals/empiricalvcov.dta") storedmarginal("$folder/empiricals") storedvars("`inputvars'")









