## Imputation Feasibility Testing for Scales

This repo contains programs that help identify sufficient conditions for successful imputation of scales with mi impute chained.

The repo is still a bit messy but we are working on cleaning it up.


## Contents of this repo

### Working programs

* core-programs.ado      : core simulation programs that are currently under development

The core program here is ifeats. It has the following syntax:

```
syntax namelist, Nobs(numlist) /// 
	       [NWitems(integer 3) NTitems(integer 12) propmiss(string) nwavemiss(string) /// 
		    MBlock(namelist) SIMCorr(integer 0) SIMMarginals(integer 0) ///
			CORRMatrix(string) MARGinals(string)]
```

The syntax needs cleaning up. A desciptions of its arguments follows:

**Required arguments**


| input       | description            |
|-------------|------------------------|
| *namelist*  | the stubs of the scales that will be simulated and imputed|
| *Nobs*      | a number list of sample sizes to be simulated |


<br>

**Options available to the user**


| option         | description            |
|----------------|------------------------|
| *NWitems*      | number of within items; in other words, the number of time points or waves (we need different name for this) |
| *NTitems*      | total number of items to be simulated (this argument needs chaning) |
| *propmiss*     | proportion missing per period; syntax for random missing patter: propmiss(sc1=(0.1 0.3) sc2=(0.2 0.5)); the first number is missingness per item, 
second value is missingness in the scale as a whole; syntax for block missing: propmiss(sc1=(0.3) sc2=(0.5)): the values is missingness in the scale as a whole |
| *nwavemiss*    | use specified time period with missingness per scale; syntax: nwavemiss(sc1=(1 3) sc2=(0 2)) |
| *mblock*       | block missing pattern; syntax: mblock(sc1 sc2) |
| *SIMCorr*      | a binary flag for whether correlation matrix is simulated or observed (needs chaning)|
| *SIMMarginals* | a binary flag for whether the marginal distributions of items are simulated or observed (needs chaning)|
| *CORRMatrix*   | location of empirical correlation matrix |
| *MARGinals*    | location of empirical marginal distributions of items |

<br>

*Examples*

```
ifeats kzf hsclg, nobs(50(50)100) nwitems(3) ntitems(36) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1))  simcorr(0) ///
        simmarginals(0)  corrmatrix("`output_folder'/empirCorrMat.dta") /// 
		marginals("`output_folder'")  mblock(kzf hsclg) nwavemiss(kzf=(0 1) hsclg=(1 2))


ifeats kzf hsclg, nobs(50(50)100) nwitems(3) ntitems(36) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1))  simcorr(0) ///
        simmarginals(0)  corrmatrix("`output_folder'/empirCorrMat.dta") /// 
		marginals("`output_folder'")  mblock(hsclg) nwavemiss(kzf=(0 1) hsclg=(1 2))

```

* scales.ado             : file for creating cummulative distributions of scale items




### Drafts

* diagnostic-tools.ado   : draft of programs to identify problems with data
* scale-preprocessing.ado: a draft scale item cleaning utility

### Older scripts

* myvcov.ado             : *old script* for simulating binary and multi-category items
* binary.ado             : *old script* for simulating binary variables

## Testing script

* test.do                : file for testing with working examples
