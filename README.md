## Imputation Feasibility Testing for Scales

The repo contains programs that help identify sufficient conditions for successful imputation of scales with mi impute chained.


### Introduction

Working with small samples and scales of multiple items that are subject to moderate missingness is challenging. 
In such settings, imputation could help tackle the missingness problem and avoid loosing observations.
However, in ultra-wide datasets with a small number of observations, chained imputation fails routinely.
Our goal is to develop a package that would help researchers identify sufficient conditions for succesful chained imputation.

### Installation

There is no need of installation. Simply include the following line in your do file to load ifeats
and accompanying utility programs:

```
do https://raw.githubusercontent.com/goshevs/ifeats/master/ifeats.ado
```

### Core simulation program

The core program `ifeats` is located in `ifeats.ado`. It has the following syntax:

```
syntax namelist, Nobs(numlist) scales(string) /// 
	       [nwaves(integer) nitems(string) tcorr(string) ///
		    propmiss(string) wavemiss(string) psdtol(real) /// 
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
| *tcorr*        | theoretical correlations, within-item and between-items as well as inter-scale item corration; required if simulating correlation matrix; see below for mode specific syntax|
| *wavemiss*     | missing waves, conditinally required; see below for mode specific syntax |
| *psdtol*       | tolerance in positive semi-definite matrix projection |
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
and then conduct the core simulation.

`propmiss` and `wavemiss` are specified as above

If simulating the correlation matrix, option `tcorr` is **required**. `tcorr` can be specified as:
- `corr(kzf=(0.9 0.2) hsclg=(0.8 0.6) corrmat=matName)`: per scale within item and between item correlations; `corrmat` is a matrix 
of inter-scale item correlations. Item correlations are assumed to be identical across all items of two scale and hence for 2 scales, `corrmat`
would be a 2x2 matrix. If `corrmat` is not specified, inter-scale correlations are assumed to be 0.

3. Full simulation (simulating both marginal distributions and correlation matrix)

`ifeats` will simulate the correlation matrix and marginal distributions of item levels and
use both in the core simulation.

`propmiss` and `wavemiss` are specified as above; `nwaves`, `nitems` and `tcorr` are required arguments
<br>


*Types of missingness simulated*:

- Random pattern by item and scale
- Block missing by scale
- Mixed pattern (a combination of the previous two)
- Block missing across scales in `namelist`



*Some examples (more available in examples.do)*

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
mat myCorrMat = (1, 0.5 \ 0.5, 1)
ifeats kzf hsclg, nobs(50(50)100) nwaves(3) nitems(kzf=12 hsclg=25) /// 
	  propmiss(kzf=(0.1 0.3) hsclg=(0.05 0.1)) ///
      scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(kzf=(0 1) hsclg=(1 2)) ///
	  tcorr(kzf=(0.9 0.2) hsclg=(0.8 0.6) corrmat=myCorrMat)

**** Block missing pattern over all specified scales
clear
ifeats kzf hsclg, nobs(50(50)100) nwaves(3) nitems(kzf=12 hsclg=25) /// 
	  propmiss(0.2) scales(kzf=(0(1)4) hsclg=(1(1)4)) wavemiss(minmax(1 2)) ///   
	  tcorr(kzf=(0.9 0.2) hsclg=(0.8 0.6))	    
```
<br>
<br>



### Utility programs 

There are two important utility programs in `ifeats.ado`: `catDist` and `dataCorrMat`.
Both programs are called by `ifeats` in some modes but can also be used as stand-alone programs.

<br>

*******

`catDist`

<br>

`catDist` creates the marginal commulative distributions of scale items. It creates
separate files for every scale. The syntax is as follows:

```
syntax namelist, saving(string)
```

`catDist` takes the following arguments:

**Required**

| argument    | description            |
|-------------|------------------------|
| *namelist*  | the stubs of the scales that will be simulated and imputed|
| *saving*    | a path pointing to a directory where the marginal distributions should be stored |

<br>

*******

<br>

`dataCorrMat`

<br>

`dataCorrMat` creates the empirical correlation matrix of the data in memory. It has the following
syntax:

```
syntax namelist, scales(string) [saving(string)]
```

`dataCorrMat` takes the following arguments:

**Required**

| argument    | description            |
|-------------|------------------------|
| *namelist*  | the stubs of the scales that will be simulated and imputed |
| *scales*    | a list of the item levels of every scale in `namelist`: `scales(sc1=(0(1)4) sc2=(0(1)4))` |

<br>

**Optional arguments:**

| argument       | description            |
|----------------|------------------------|
| *saving*       | a directory and filename to be used for the file containing the empirical correlation matrix  |

<br>

*******


### TODO

- pull out outcomes of psd projection and display them
- check security (file persistence)
- plug and play module for fiml
- variance assessment for successful imputation



### Other programs/files

#### Drafts

* diagnostic-tools.ado   : drafts of programs to identify problems in data
* scale-preprocessing.ado: a draft scale item cleaning utility

#### Examples/testing

* examples.do                : file for testing with working examples
