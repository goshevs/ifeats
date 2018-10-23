## Imputation Feasibility Testing for Scales

This repo contains programs that help identify sufficient conditions for successful imputation of scales with mi impute chained.

The repo is still a bit messy but we are working on cleaning it up.


## Contents of this repo

### Working programs

****

`core-programs.ado` -- core simulation programs that are currently under development

The core program here is `ifeats`. It has the following syntax:

```
syntax namelist, Nobs(numlist) scales(string) propmiss(string)  /// 
	       [nwaves(integer) nitems(string) wavemiss(string) /// 
		    CORRMatrix(string) MARGinals(string)]
```

`ifeats` takes the following arguments:

**Required**

| argument    | description            |
|-------------|------------------------|
| *namelist*  | the stubs of the scales that will be simulated and imputed|
| *Nobs*      | a number list of sample sizes to be simulated |
| *scales*    | a list of the item levels of every scale in `namelist`: `scales(sc1=(0(1)4) sc2=(0(1)4))` |
| *propmiss*  | proportion missing observations; see below for mode specific syntax |

<br>

**Options and conditionally required arguments:**


| argument       | description            |
|----------------|------------------------|
| *nwaves*       | number of waves/time periods; required if conducting a full simulation |
| *nitems*       | number of items per scale: `nitems(sc1=12 sc2=36)`; required if conducting a full simulation |
| *wavemiss*     | missing waves, conditinally required; see below for mode specific syntax |
| *CORRMatrix*   | location of empirical correlation matrix |
| *MARGinals*    | location of empirical marginal distributions of items |

<br>


`ifeats` can be used in several different modes:

1. With a dataset loaded in memory

`ifeats` retrieves the empirical correlation matrix and marginal 
distibutions of the levels of every item (if not provided by the user), and then conducts the core simulation.

`propmiss` can be specified as:

- `propmiss(kzf=(0.05 0.2) hsclg=(0.05 0.3))` for random missing (item missigness, scale missingness)
- `propmiss(kzf=(0.2) hsclg=(0.3))` for block missing (scale missingness)
- `propmiss(kzf=(0.2) hsclg=(0.05 0.3)`) for mixed missing pattern
- `propmiss(0.2)` for block missingness across all scales in `namelist`

`wavemiss` can be specified as:

- `wavemiss(kzf=(0 1) hsclg=(1 2))` waves missing for every scale
- can be omitted in which case waves are selected at random
- **must be specified** as `wavemiss(minmax(1 2))` if `propmiss` is selected to be block missingness across all scales in `namelist`;
the user must specify in `minmax()` the minimum and maximum number of waves with missing values


2. With an existing correlation matrix and/or marginal distributions of item levels

Irrespectively of whether data are present in memory, `ifeats` will use
either or both stored (if provided) pieces of data, simulate the one that is not provided, 
and then conduct the core simulate.

`propmiss` and `wavemiss` are specified as above


3. Full simulation

`ifeats` will simulate the correlation matrix and marginal distributions of item levels and
use both in the core simulation.

`propmiss` and `wavemiss` are specified as above; `nwaves` and `nitems` are required arguments
<br>


*Types of missingness simulated*:

- Random pattern by item and scale
- Block missing by scale
- Mixed pattern (a combination of the previous two)
- Block missing across scales in `namelist`



*Some examples (more available in test.do)*

```

********************************************************************************
*** Peviously stored corr matrix and marginal distributions

**** Complete syntax
ifeats kzf hsclg, nobs(50(50)100) propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(kzf=(0 1) hsclg=(1 2)) ///
        corrmatrix("`output_folder'/empirCorrMat.dta") marginals("`output_folder'") 


**** Block missing pattern over all specified scales
ifeats kzf hsclg, nobs(50(50)100) propmiss(0.3) wavemiss(minmax(1 2)) ///
		scales(kzf=(0(1)4) hsclg=(1(1)4)) ///
        corrmatrix("`output_folder'/empirCorrMat.dta") marginals("`output_folder'")


********************************************************************************
*** Full simulation (simulating everything)

**** Complete syntax
clear
ifeats kzf hsclg, nobs(50(50)100) nwaves(3) nitems(kzf=12 hsclg=25) ///
	  propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
      scales(kzf=(0(1)4) hsclg=(0(1)4)) wavemiss(kzf=(0 1) hsclg=(1 2))


**** Block missing pattern over all specified scales
clear
ifeats kzf hsclg, nobs(50(50)100) nwaves(3) nitems(kzf=12 hsclg=25) propmiss(0.2) ///
      scales(kzf=(0(1)4) hsclg=(0(1)4)) wavemiss(minmax(1 2)) 	    
```
<br>
<br>

****

`scales.ado` a file for creating cummulative distributions of scale items

The core programs here are `catDist` and `dataCorrMat`. Both programs are utility programs
used to extract the marginal distributions of item levels and correlation matrix of data.
More about them will be posted soon.

****

### Drafts

* diagnostic-tools.ado   : draft of programs to identify problems with data
* scale-preprocessing.ado: a draft scale item cleaning utility

### Older scripts

* myvcov.ado             : *old script* for simulating binary and multi-category items
* binary.ado             : *old script* for simulating binary variables

### Testing script

* test.do                : file for testing with working examples
