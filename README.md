## Imputation Feasibility Testing for Scales

This repo contains programs that help identify sufficient conditions for successful imputation of scales with mi impute chained.

Currently, the repo is in a messy working state. Clean-up is underway.


## Content

The repo containes Stata ado files that define the following programs:

myvcov -- extract the empirial correlation matrix from the data -- works both for binary and categorical vars

Functions that create theoretical correlation matrices and a mean vector -- works both for binary and categorical vars
genwithin
gencombined
genmeans


mysim -- shell function which switches between empirical and simulated correlation matrix
mysimcore -- meat function that does the simulation -- works with binary variables only


catsim -- shell function which switches between empirical and simulated correlation matrix
catesimcore --  meat function that does the simulation -- works with categorical variables (currently fails on line 

	forval i = 1/`ntitems' {
			*** requires egenmore package ***
			egen cmyvar`i' = xtile(myvar`i'), percentiles(`cp`i'')
			drop myvar`i'
			rename cmyvar`i' myvar`i'



scaleclean -- figures out variables types and returns all levels per label
mypercent -- saves the commulative distribution of a categorical variable to a file (to be used for simulation!)

