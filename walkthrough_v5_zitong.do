
********************************************************************************
*** Using the commands in myvcov.ado

*** Exracting the empirical vcov of ka_ecdm* and saving it to a file
 myvcov kzf*, saving("`temp_data'/categorivcov.dta")
 mypercent kzf*, saving("`temp_data'")  // , saving("`temp_data'/catepercent.dta")

**** By Zitong: Think about categorical variables. 
* mysim, simobs(50(50)100) nwitems(3) ntitems(33) propmiss(0.2) mblock(1) simvcov(0) storedvcov("`temp_data'/categorivcov.dta")

catsim, simobs(50(50)100) nwitems(3) ntitems(36) propmiss(0.2) mblock(1) simvcov(0) simmarginal(0) ///
        storedvcov("`temp_data'/categorivcov.dta") storedmarginal("`temp_data'/percents/scakzperc.dta")

		
* mat list r(sims)


*** Exracting the empirical vcov of kzf* and saving it to a file


/*
*** Simulating with the empirical vcov of ka_ecdm*
mysim, simobs(50(50)100) nwitems(3) ntitems(33) propmiss(0.2) mblock(1) simvcov(0) storedvcov("$folder/empiricalvcov.dta")
*** See the simulation results:
mat list r(sims)

*** Simulating with theoretical vcov for ka_ecdm*
mysim, simobs(50(50)100) simcorb(0(.1).2) nwitems(3) ntitems(33) corw(0.9) propmiss(0.2) mblock(1) simvcov(1)
*** See the simulation results:
mat list r(sims)

/*
***** Options for mysim
simobs    : sample size
simcorb   : between-item correlation
nwitems   : number of within items (essentially the number of time periods)
ntitems   : number of total items (time periods * number of indicators)
corw      : within-item correlation (over time)
propmiss  : proportion of missing (based on obs)
mblock    : block missing (boolean); all items missing in a period
simvcov   : simulate vcov or use empirical vcov (boolean)
storedvcov: if simvcov(0) then where is the empirical vcov stored

*/





