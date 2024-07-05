*** The weaktsiv function carries out weak-iv robust inference for two-sample IV regressions. ***
*** It reports the classic TS2SLS estimator as well as the TSAR, TSK, and TSCLR proposed in ***
*** Choi, et. al. (2018). 

*** Note that the benchmark TSAR, TSK, and TSCLR methods require homoskedasticity as well as equal***
*** moments of exogeneous regressors, just like the classic TS2SLS estimation method. But the ***
*** robust TSAR, TSK, and TSCLR methods allow for heteroskedasticity as well as unequal moments of ***
*** exogeneous regressors. 

*** Our codes for the benchmark tests and confidence regions are modified from the condivreg codes ***
*** by Mikusheva and Poi (2006), who propose fast computation methods for the classic one-sample *** 
*** AR, K, and CLR methods. Benchmark confidence regions are calculated by analytical solutions. ***

*** Our codes for the robust tests and confidence regions are modified from Finlay and Magnusson (2009), ***
*** who propose minimum distance based tests and confidence regions for classic one-sample iv models ***
*** with potential weak first stage and also heteroskedasticity. Robust confidence regions are calculated ***
*** by grid search.



**********************
*** Data Structure ***
**********************

*** y  x  z  ***
*** y1 .  z1 *** (First sample, we observe y and z.)
*** .  x2 z2 *** (Second sample, we observe x and z.)

*************************************
************* weaktsiv **************
************* main program **********
*************************************
program define weaktsiv, eclass

        version 14.2
                 
        /* First get the list of variables.     */
        gettoken lhs 0 : 0, match(paren)
        IsStop `lhs'
        if `s(stop)' { 
                error 198 
        }
        while `s(stop)'==0 {
                if "`paren'"=="(" {
                        gettoken p lhs : lhs, parse(" =")
                        while "`p'"!="=" {
                                local end1  `p'
                                gettoken p lhs : lhs, parse(" =")
                        }
                        tsunab end1 : `end1'
                        tsunab rinst : `lhs'
                }
                else {
                        local exog `exog' `lhs'
                }
                gettoken lhs 0 : 0, match(paren)
                IsStop `lhs'
        }
        local 0 `"`lhs' `0'"'

        tsunab exog : `exog'
        tokenize `exog'
        local ry1 "`1'"
        local 1 " "
        local rexog `*'
        loc ry2 "`end1'"
        syntax [if] [in], [nocons robust       /* 
		        */         level(integer $S_level) test(real 0)        /*
                */         grid(numlist ascending) points(integer 100) /*
				*/         gridmult(real 2) retmat]

		if "`robust'" == "" & "`grid'" != ""{
             di as error "Benchmark confidence sets are calculated analytically, grid option not comptible."
             exit
        }	
		if "`robust'" == "" & "`points'" != "100"{
             di as error "Benchmark confidence sets are calculated analytically. Points option not comptible."
             exit
        }
		
*************************************
*** Do not exclude missing values ***
*** Create samp = 1 and samp = 2  ***
*************************************

	marksample touse
	markout `touse' `rexog' `rinst'
	    
		qui drop if `ry1'==. & `ry2'==.
		qui drop if `ry1'~=. & `ry2'~=.
		
		tempvar samp
		gen `samp' = 1
		qui replace `samp' = 2 if `ry1'==.


/* check for multicollinearity */
                _rmcoll `rexog' `rinst', `cons'
                local newcnt : word count `r(varlist)'
                local oldcnt : word count `rexog' `rinst'
                if `newcnt' != `oldcnt' {
                        noi di as error "Multicollinearity!"
                        exit 498
                }

loc k : word count `rinst'
**********************************
*** Compute the TS2SLS results ***
**********************************
               tempname beta var rmse rss  ivbetahat ivbetahatse fstat2 dfr2 dfm dfr n1 n2 mss  rsq adjrsq chi2 Vy1het Vx2het seb2shet var1het
			   if "`robust'" != "" {
               myts2sls_robust `ry1' `ry2' "`rexog'" "`rinst'" "`cons'"  `touse'  `samp' 		
			   }
			   else {
               myts2sls `ry1' `ry2' "`rexog'" "`rinst'" "`cons'"  `touse'  `samp' 					   
			   }
               mat `beta' = r(beta)
			   				if "`robust'" != "" {			
			   mat `var1het' = r(var1het)
	           if "`cons'" == "" {
                        mat colnames `beta' = `ry2' `rexog' _cons
                        mat rownames `var1het' = `ry2' `rexog' _cons
                        mat colnames `var1het' = `ry2' `rexog' _cons
                                                                
                }
                else {
                        mat colnames `beta' = `ry2' `rexog'
                        mat rownames `var1het' = `ry2' `rexog'
                        mat colnames `var1het' = `ry2' `rexog'
                }		   
			   }
			   else {
			   mat `var' = r(var)
               if "`cons'" == "" {
                        mat colnames `beta' = `ry2' `rexog' _cons
                        mat rownames `var' = `ry2' `rexog' _cons
                        mat colnames `var' = `ry2' `rexog' _cons
                                                                
                }
                else {
                        mat colnames `beta' = `ry2' `rexog'
                        mat rownames `var' = `ry2' `rexog'
                        mat colnames `var' = `ry2' `rexog'
                }
				}
	           sca `ivbetahat' = `beta'[1,1]
			   if "`robust'" != "" {	
	           sca `ivbetahatse' = sqrt(`var1het'[1,1])
			   }
			   else {
	           sca `ivbetahatse' = sqrt(`var'[1,1])
			   }
		 
			   
			   sca `fstat2' = r(fstat2)
               sca `dfr2' = r(dfr2)
			   sca `rsq' = r(rsq)
               sca `adjrsq' = r(adjrsq)
			   sca `dfm' = r(dfm)
               sca `dfr' = r(dfr)
			   sca `n1' = r(n1)
			   sca `n2' = r(n2)
			   sca `rss'= r(rss)
			   sca `mss'= r(mss)
			   sca `rmse' = sqrt(`rss'/`dfr')
                if "`cons'" != "" {
                        sca `chi2' = .
                }
                else {
                        tempname varslopes dif vinv
                        mat `dif' = `beta'[1,1..`dfm']
						
						if "`robust'" != "" {	
						mat `varslopes'=`var1het'[1..`dfm',1..`dfm']
						}
						else {
						mat `varslopes'=`var'[1..`dfm',1..`dfm']						
						}
                        mat `vinv' = inv(`varslopes')
                        mat `chi2' = `dif'*`vinv'*`dif''
                        sca `chi2' = trace(`chi2')
                }
			   tempname bbb 
               sca `bbb'=`test'
			   
/* Post TS2SLS results */
                tempvar touse2
							   				
                gen byte `touse2' = `touse'   
							
				if "`robust'" != "" {
				eret post `beta' `var1het', esample(`touse2') depname("`ry1'")  
				}
				else {
                eret post `beta' `var', esample(`touse2') depname("`ry1'")  	
				}
								
                eret scalar df_r_first = `dfr2'
				eret scalar numinst=`k'
				eret scalar F_first = `fstat2'	
				eret scalar rmse =`rmse'
				eret scalar rss = `rss'
				eret scalar mss = `mss'
				eret scalar r2_a = `adjrsq'
				eret scalar r2 = `rsq'	
                eret scalar df_r = `dfr'
				eret scalar df_m = `dfm'				
				eret scalar chi2 = `chi2' 				
				eret scalar n2 = `n2'
				eret scalar n1 = `n1'
				eret scalar H0_b = `bbb'   
                eret scalar level = `level'	
                eret local exog "`rexog'"
                eret local insts "`rinst'"
                eret local instd "`ry2'"
				
                if "`cons'" == "" {
                        eret local cons "yes"
                }
                else {
                        eret local cons "no"
                }			
			      
/* preserve currently posted results*/
               tempname myest
              _estimates hold `myest'
			  
               tempvar y1 y2
               qui gen double `y1' = `ry1' if `touse' & `samp' == 1
               qui gen double `y2' = `ry2' if `touse' & `samp' == 2
			   
			   tempvar one
			   gen `one' = 1

			   if "`cons'" == "" {
                        loc exog "`rexog' `one'"
                }
                else {
                        loc exog "`rexog'"
                }	
				
							
****************************************************************
*** Start to calculate weak iv robust inference methods. *******
*** First reconstruct iv and make iv orthogonal to exog var. ***
****************************************************************
                tempname ehold
                loc inst0 = ""		
                qui foreach v in `rinst' {
                        tempvar inst00 inst01 inst02
                        reg `v' `exog' if `touse' & `samp' == 1, nocons
                        predict double `inst01' if `touse' & `samp' == 1, residuals
                        reg `v' `exog' if `touse' & `samp' == 2, nocons
                        predict double `inst02' if `touse' & `samp' == 2, residuals
						gen double `inst00' = `inst02'
						replace `inst00' = `inst01' if `touse' & `samp' == 1
                        loc inst0 "`inst0' `inst00'"
                  }	
				  
***************************************************				
*** Compute Omega following Choi et. al. (2018) ***
***************************************************

                tempname mzy1_rss mzy2_rss n1 n2 df1 df2 omega cross1 cross2  /*
				*/ZPZ1 ZPZ2 sigma1 sigma2 sigma sqrtsigma hateta hatpi ARstat /*
				*/ SIGMAU1 SIGMAE2 z1pz1 z2pz2 SIGMA1 SIGMA2 SIGMA sqrtSIGMA
				tempvar mzy1 mzy2 w1hat

                qui reg `y1' `inst0' `exog' if `touse' & `samp' == 1, nocons 
				scalar `n1' = e(N)
				scalar `mzy1_rss' = e(rss)
				scalar `df1' = e(df_r)
                qui predict double `mzy1' if `touse' & `samp' == 1, residuals

                qui reg `y2' `inst0' `exog' if `touse' & `samp' == 2, nocons
				scalar `n2' = e(N)
				scalar `mzy2_rss' = e(rss)
				scalar `df2' = e(df_r)	
                qui predict double `w1hat' if `touse' & `samp' == 1, xb
                qui predict double `mzy2' if `touse' & `samp' == 2, residuals				
                mat `omega' = J(2,2,0)
				mat `omega'[1,1] = `mzy1_rss' / `df1'
				mat `omega'[2,2] = `mzy2_rss' / `df2' * `n1' / `n2'
				
qui if "`robust'" != ""{
*************************************************************				
*** Compute SIGMA (robust) following  Choi et. al. (2018) ***
*************************************************************
				*** hateta hatpi ***
				reg `y1' `inst0' if `touse' & `samp' == 1, nocon
                matrix `hateta' = e(b)
				reg `y2' `inst0' if `touse' & `samp' == 2, nocon
				matrix `hatpi' = e(b)
			
				*** SIGMAU1 SIGMAE2 ***
				local AAAA = ""
				local BBBB = ""
				foreach v in `inst0' {
				   tempvar tempinst1`v' tempinst2`v'
                   gen double `tempinst1`v'' = `v' * `mzy1' if `touse' & `samp' == 1
			       local AAAA "`AAAA' `tempinst1`v''"   
				   gen double `tempinst2`v'' = `v' * `mzy2' if `touse' & `samp' == 2
			       local BBBB "`BBBB' `tempinst2`v''"  
						}
                mat accum `SIGMAU1' = `AAAA' if `touse' & `samp' == 1, noconstant		
                mat accum `SIGMAE2' = `BBBB' if `touse' & `samp' == 2, noconstant
				
				*** SIGMA1 and SIGMA2 ***
				mat accum `z1pz1' = `inst0' if `touse' & `samp' == 1, noconstant
				mat accum `z2pz2' = `inst0' if `touse' & `samp' == 2, noconstant		
				mat `SIGMA1' = invsym(`z1pz1')*`SIGMAU1'*invsym(`z1pz1')*`n1'^2 / `df1'
				mat `SIGMA2' = invsym(`z2pz2')*`SIGMAE2'*invsym(`z2pz2')*`n2'^2 / `df2'

****************************************************************************
*** Robust TSAR, TSK, TSCLR confidence regions using grid search. **********
*** The next block is modified from rivtest.ado by Finlay and Magnusson. ***
****************************************************************************	
                        local alpha=1-`level'/100			
						local gridinput : length local grid
						if `gridinput'==0 {
							* default grid radius is twice that of the confidence interval from the TS2SLS estimation
								local gridradius = abs(`gridmult') * `ivbetahatse' * invnormal(1-`alpha'/2)
							* create grid for confidence sets
								local gridmin = `ivbetahat' - `gridradius'
								local gridmax = `ivbetahat' + `gridradius'
								local gridinterval = .999999999*(`gridmax'-`gridmin')/(`points'-1)
								local grid "`gridmin'(`gridinterval')`gridmax'"
								local gridbegin : di %8.0g `gridmin'
								local gridend : di %8.0g `gridmax'
						}
						numlist "`grid'"
						local gridlist "`r(numlist)'"
						local points : word count `gridlist'
						if `gridinput'>0 {
							local gridbegin : word 1 of `gridlist'
							local gridend : word `points' of `gridlist'
						}
						local grid_description "[`gridbegin',`gridend']"

				* create macros for storing confidence sets
					local testlist "tsar tsk tsclr"
					foreach testname in `testlist' {
						local `testname'_cset ""
						local `testname'_rbegin=0
						local `testname'_rend=0
						local `testname'_rbegin_null=0
						local `testname'_rend_null=0
					}
				local counter = 0
				_dots `counter' 0, title(Estimating confidence sets over grid points)
				foreach gridnull in `gridlist' {
					local ++counter
					_dots `counter' 0
					tempname rk numk tsar_p tsar_chi2 tsar_df tsk_p tsk_chi2 tsk_df tsclr_p tsclr_stat tsclr_df tsar_r tsk_r j_r tsclr_r 
					
					* calculate test stats
						computeivtests_robust, n1(`n1') n2(`n2') hateta(`hateta') hatpi(`hatpi') sigmaf(`SIGMA1') sigmas(`SIGMA2') null(`gridnull')
						scalar `tsar_chi2'=r(tsar_chi2)
						scalar `tsk_chi2'=r(tsk_chi2)
						scalar `tsclr_stat'=r(tsclr_stat)
						scalar `rk'=r(rk)					   						

                        loc numk : word count `inst0'
						
					* calculate test statistics, p-values, and rejection indicators from above matrices
						compute_pvals, null(`gridnull') rk(`rk') numk(`numk') level(`level')  ///
							tsar_p(`tsar_p') tsar_chi2(`tsar_chi2') tsar_df(`tsar_df')  ///
							tsk_p(`tsk_p') tsk_chi2(`tsk_chi2') tsk_df(`tsk_df') ///
							tsclr_p(`tsclr_p') tsclr_stat(`tsclr_stat')  tsclr_df(`tsclr_df') ///
							tsar_r(`tsar_r') tsk_r(`tsk_r') tsclr_r(`tsclr_r') 							
				   
					* write out confidence sets from rejection indicators
						if `tsclr_stat'==.								local tsclr_cset "."
						foreach testname in `testlist' {
							if "``testname'_cset'"!="." { 
								if ``testname'_r'==0 {
									if ``testname'_rbegin'==0 {
										local `testname'_rbegin=`counter'
										local `testname'_rbegin_null=`gridnull'
									}
									local `testname'_rend=`counter'
									local `testname'_rend_null=`gridnull'
								}
								if ``testname'_r'==1 | (``testname'_r'==0 & `counter'==`points') {
									if ``testname'_rbegin'>0 & ``testname'_rend'>0 & ``testname'_rbegin'==``testname'_rend' {
										local rnull : di %8.0g ``testname'_rbegin_null'
										if length("``testname'_cset'")==0	local `testname'_cset "`rnull'"
										else								local `testname'_cset "``testname'_cset' U `rnull'"
										local `testname'_rbegin=0
										local `testname'_rend=0
									}
									else if ``testname'_rbegin'>0 & ``testname'_rend'>0 & ``testname'_rbegin'<``testname'_rend' {
										local rnull1 : di %8.0g ``testname'_rbegin_null'
										local rnull2 : di %8.0g ``testname'_rend_null'
										if length("``testname'_cset'")==0	local `testname'_cset "[`rnull1',`rnull2']"
										else								local `testname'_cset "``testname'_cset' U [`rnull1',`rnull2']"
										local `testname'_rbegin=0
										local `testname'_rend=0
									}
								}
							}
						}
				
}
				
				foreach testname in `testlist' { 
					if length("``testname'_cset'")==0 		local `testname'_cset "empty"
				}
				
***************************************
*** Robust TSAR, TSK, TSCLR tests.  ***
***************************************
	tempname rk numk tsar_p tsar_chi2 tsar_df tsk_p tsk_chi2 tsk_df tsclr_p tsclr_stat tsclr_df tsar_r tsk_r j_r tsclr_r ///

	* calculate test stats
	computeivtests_robust, n1(`n1') n2(`n2') hateta(`hateta') hatpi(`hatpi') sigmaf(`SIGMA1') sigmas(`SIGMA2') null(`bbb')
						
	scalar `tsar_chi2'=r(tsar_chi2)
	scalar `tsk_chi2'=r(tsk_chi2)
	scalar `tsclr_stat'=r(tsclr_stat)
	scalar `rk'=r(rk)	
	scalar points=`points'

    loc numk : word count `inst0'
						
	* calculate test statistics, p-values, and rejection indicators from above matrices
	compute_pvals, null(`bbb') rk(`rk') numk(`numk') level(`level')  ///
	               tsar_p(`tsar_p') tsar_chi2(`tsar_chi2') tsar_df(`tsar_df')  ///
				   tsk_p(`tsk_p') tsk_chi2(`tsk_chi2') tsk_df(`tsk_df') ///
				   tsclr_p(`tsclr_p') tsclr_stat(`tsclr_stat')  tsclr_df(`tsclr_df') ///
				   tsar_r(`tsar_r') tsk_r(`tsk_r') tsclr_r(`tsclr_r')	
				   
	  		   
}

/* Safe to restore my results */

             _estimates unhold `myest'

if "`robust'" != ""{	
		eret scalar points=`points'
		eret local grid `grid_description'		
		eret scalar p_TSAR = `tsar_p'
        eret scalar p_TSK = `tsk_p'
        eret scalar p_TSCLR = `tsclr_p'
		eret local tsk_cset="`tsk_cset'"
		eret local tsar_cset="`tsar_cset'"
		eret local tsclr_cset="`tsclr_cset'"			   
}



/*------Benchmark TSAR, TSK, TSCLR Confidence sets, Analytical Solution -------------*/
/*
Test            Result type             Interval
-----------------------------------------------------------------------
TSCLR           1                       Empty set
                2                       [x1, x2]
                3                       (-infty, +infty)
                4                   (-infty, x1] U [x2, infty)
                
TSAR            1                       Empty set
                2                       [x1, x2]
                3                       (-infty, +infty)
                4                   (-infty, x1] U [x2, infty)
                
TSK             1                      Not used (not possible)
                2                       [x1, x2]                                
                3                       (-infty, +infty)
                4               (-infty, x1] U [x2, infty)
                5               (-infty, x1] U [x2, x3] U [x4, infty)
                6                   [x1, x2] U [x3, x4]

-----------------------------------------------------------------------
*/
    if "`robust'" == ""{            
                tempname cross zpz sqrtzpzi zpy MM ypz sqrtomegai v     /*
                        */      d M N alpha C A D aa x1 x2 g type											

                qui mat accum `cross' = `inst0' `y1' `w1hat'                /*
                        */      if `touse' & `samp' == 1, noconstant
				
                mat `zpz' = `cross'[1..`k', 1..`k']
                mat `zpy' = `cross'[1..`k', (`k'+1)...]
                mat_inv_sqrt `zpz' `sqrtzpzi'
                mat_inv_sqrt `omega' `sqrtomegai'
                mat `ypz'=`zpy''
                mat `MM' = `sqrtomegai'*`ypz'*invsym(`zpz')*`zpy'*`sqrtomegai'
				
                mat symeigen `v' `d' = `MM'
                sca `M' = `d'[1,1]
                sca `N' =`d'[1,2]
                sca `alpha' = 1-`level'/100
				
/* inversion of TSCLR*/
                inversefun `M' `k' `alpha' `C'

                mat `A' =invsym(`omega')*`ypz'*invsym(`zpz')*`zpy'*invsym(`omega')- /*
                        */      `C'*invsym(`omega')
                sca `D' = -det(`A')
                sca `aa' = `A'[1,1]

                if (`aa'<0) {
                        if (`D' <0) {
                                sca `type'=1
                        }
                        else{
                                sca `type'=2
                                sca `x1'= (-`A'[1,2] + sqrt(`D'))/`aa'
                                sca `x2' = (-`A'[1,2] - sqrt(`D'))/`aa'
                                mat `g'=(`x1'\ `x2')
                        }
                }
                else{
                        if (`D'<0) {
                                sca `type'=3
                        }
                        else {
                                sca `type'=4
                                sca `x1'= (-`A'[1,2]-sqrt(`D'))/`aa'
                                sca `x2'= (-`A'[1,2]+sqrt(`D'))/`aa'
                                mat `g'=(`x1' \ `x2')
                        }
                }
                ereturn local TSCLR_type `"`=`type''"'
                if `type' == 2 | `type' == 4 {
                        eret scalar TSCLR_x1 = `x1'
                        eret scalar TSCLR_x2 = `x2'
                }

/* inversion of TSK */
                tempname tskcv q1 q2 A1 A2 D1 D2 y1 y2 y3 y4 type1
                sca `tskcv' = invchi2tail(1, (1-`level'/100))
                if (`k'==1) {
                        sca `q1' = `M'-`tskcv'
                        mat `A1' =invsym(`omega')*`ypz'*invsym(`zpz')*`zpy'*  /*
                                */      invsym(`omega')-`q1'*invsym(`omega')
                        sca `D1' = -4*det(`A1')
                        sca `y1'= (-2*`A1'[1,2]+ sqrt(`D1'))/2/`A1'[1,1]
                        sca `y2'= (-2*`A1'[1,2]-sqrt(`D1'))/2/`A1'[1,1]
                        if (`A1'[1,1]>0) { 
                                if (`D1'>0) {
                                        sca `type1'=4 
                                                /* two infinite intervals*/
                                        eret scalar TSK_x1 = `y2'
                                        eret scalar TSK_x2 = `y1'
                                }
                                else {
                                        sca `type1'=3
                                }
                        }
                        else{
                                if (`D1'>0) {
                                        sca `type1'=2 /*one interval */
                                        eret scalar TSK_x1 = `y1'
                                        eret scalar TSK_x2 = `y2'

                                }
                                else {
                                        sca `type1'=3
                                }
                        }
                }
                else {
                        if ((`M' +`N' - `tskcv')^2-4*`M'*`N'<0) { 
                                sca `type1' = 3
                        }
                        else {
                            sca `q1' = (`M'+ `N' - `tskcv' -     /*
                                */ sqrt((`M'+`N' - `tskcv')^2 -  /*
                                */      4*`M'*`N'))/2
                            sca `q2' = (`M'+`N' - `tskcv' +              /*
                                */ sqrt((`M'+`N'-`tskcv')^2 -    /*
                                */      4*`M'*`N'))/2
                            if ((`q1' < `N') | (`q2' > `M')) {
                                sca `type1' = 3
                            }
                            else {              
                                mat `A1' = invsym(`omega')*`ypz'*invsym(`zpz')*   /*
                                        */ `zpy'*invsym(`omega')-`q1'*invsym(`omega')
                                mat `A2' = invsym(`omega')*`ypz'*invsym(`zpz')*   /*
                                        */ `zpy'*invsym(`omega')-`q2'*invsym(`omega')
                                sca `D1' = -4*det(`A1')
                                sca `D2' = -4*det(`A2')
                                if (`A1'[1,1]>0) { 
                                        if (`A2'[1,1]>0) { 
                                                sca `type1' = 5
                                                sca `y1' = (-2*`A1'[1,2] +  /*
                                                    */ sqrt(`D1'))/2/`A1'[1,1]
                                                sca `y2' = (-2*`A1'[1,2] -  /*
                                                    */ sqrt(`D1'))/2/`A1'[1,1]
                                                sca `y3' = (-2*`A2'[1,2] +  /*
                                                    */ sqrt(`D2'))/2/`A2'[1,1]
                                                sca `y4' = (-2*`A2'[1,2] -  /*
                                                    */ sqrt(`D2'))/2/`A2'[1,1]
                                                eret scalar TSK_x1 = `y4'
                                                eret scalar TSK_x2 = `y2'
                                                eret scalar TSK_x3 = `y1'
                                                eret scalar TSK_x4 = `y3'
                                        }
                                        else {
                                                sca `type1' = 6
                                                sca `y1' = (-2*`A1'[1,2] +  /*
                                                    */ sqrt(`D1'))/2/`A1'[1,1]
                                                sca `y2' = (-2*`A1'[1,2] -  /*
                                                    */ sqrt(`D1'))/2/`A1'[1,1]
                                                sca `y3' = (-2*`A2'[1,2] +  /*
                                                    */ sqrt(`D2'))/2/`A2'[1,1]
                                                sca `y4' = (-2*`A2'[1,2] -  /*
                                                    */ sqrt(`D2'))/2/`A2'[1,1]
                                                eret scalar TSK_x1 = `y3'
                                                eret scalar TSK_x2 = `y4'
                                                eret scalar TSK_x3 = `y2'
                                                eret scalar TSK_x4 = `y1'
                                        }
                                }
                                if (`A1'[1,1]<=0) {
                                        sca `type1' =5
                                        sca `y1' = (-2*`A1'[1,2] +      /*
                                                */  sqrt(`D1'))/2/`A1'[1,1]
                                        sca `y2' = (-2*`A1'[1,2] -      /*
                                                */  sqrt(`D1'))/2/`A1'[1,1]
                                        sca `y3' = (-2*`A2'[1,2] +      /*
                                                */  sqrt(`D2'))/2/`A2'[1,1]
                                        sca `y4' = (-2*`A2'[1,2] -      /*
                                                */  sqrt(`D2'))/2/`A2'[1,1]
                                        eret scalar TSK_x1 = `y1'
                                        eret scalar TSK_x2 = `y3'
                                        eret scalar TSK_x3 = `y4'
                                        eret scalar TSK_x4 = `y2'
                                }
                            }
                        }
                }
                eret local TSK_type `"`=`type1''"'			

/* inversion of TSAR */
                tempname tsarcv  AAA type2 xx1 xx2 DDD aaa
                sca `tsarcv' = invchi2tail(`k', (1-`level'/100))
                mat `AAA' =`ypz'*invsym(`zpz')*`zpy'-`tsarcv'*`omega'
                sca `DDD' = -det(`AAA')
                sca `aaa' = `AAA'[2,2]
                if (`aaa'<0) {
                        if (`DDD' <0) {
                                sca `type2'=3
                        }
                        else{
                                sca `type2'=4
                                sca `xx1'= (`AAA'[1,2] + sqrt(`DDD'))/`aaa'
                                sca `xx2' = (`AAA'[1,2] - sqrt(`DDD'))/`aaa'
                                eret scalar TSAR_x1 = `xx1'
                                eret scalar TSAR_x2 = `xx2'
                         }
                }
                else {
                        if (`DDD'<0) {
                                sca `type2'=1
                        }
                        else {
                                sca `type2'=2
                                sca `xx1'= (`AAA'[1,2]-sqrt(`DDD'))/`aaa'
                                sca `xx2'= (`AAA'[1,2]+sqrt(`DDD'))/`aaa'
                                eret scalar TSAR_x1 = `xx1'
                                eret scalar TSAR_x2 = `xx2'
                        }
                }
                eret local TSAR_type `"`=`type2''"'					   
			   
/*------Benchmark TSAR, TSK, TSCLR Tests-----*/
        tempname  a b aprime bprime     
        
        mat `a' = (`bbb'\1)
        mat `b' = (1\(-1*`bbb'))
        mat `aprime' = `a''             
        mat `bprime' = `b''             

        tempname oia
        mat `oia' = invsym(`omega')*`a'
        tempname bob aoia
        matrix `bob' = `bprime'*`omega'*`b'
        scalar `bob' = trace(`bob')
        matrix `aoia' = `aprime'*invsym(`omega')*`a'
        scalar `aoia' = trace(`aoia')
        
        tempname  sbar tbar
        mat `sbar' = `sqrtzpzi'*`zpy'*`b'/sqrt(`bob')
        mat `tbar' = `sqrtzpzi'*`zpy'*`oia'/sqrt(`aoia')
 
        tempname tsar tsk tsclr wald qt tsclrpv
        calcstat `sbar' `tbar' `tsar' `tsk' `qt' `tsclr' 
        loc kk=`k'-1

        tsclr_pcal `k' `qt' `tsclr' `tsclrpv'     
		tempname tsarpv tskpv
        sca `tsarpv' = 1-chi2(`k', `tsar')  	
        sca `tskpv' = 1-chi2(1, `tsk')			
        if `tsclrpv'<0	scalar `tsclrpv'=0.0000  
 
/* Post results */   
        eret scalar p_TSAR = `tsarpv'
        eret scalar p_TSK = `tskpv'
        eret scalar p_TSCLR = `tsclrpv'
	}	
	
				if "`robust'" != "" {
                        eret local robust "yes"
                }
                else {
                        eret local robust "no"
                }	
	eret local cmd "weaktsiv"
/*------Print results-----*/	
				if "`robust'" != "" {
                        myprint_robust  
                }
                else {
                        myprint  
                }	
          
		
		
if `k'>1 & "`robust'" == ""{
        myprintregions 	
} 
else if `k'==1 & "`robust'" == ""{
        myprintregions_just 
} 
else if `k'>1 & "`robust'" != ""{
        myprintregions_robust "`tsar_cset'" "`tsk_cset'" "`tsclr_cset'" "`tsar_p'" "`tsk_p'" "`tsclr_p'" 
}			
else {
        myprintregions_robust_just "`tsclr_cset'" "`tsclr_p'"
}
if "`robust'" != ""{
        myprintregions_footnote
} 
end					
		
		
/*----------------------------------------------------------------------*/
/* Computes two-sample test statistics.       */
prog def calcstat

        args sbar tbar tsar tsk qt tsclr 
        tempname sbarp tbarp st
        
        mat `sbarp' = `sbar''
        mat `tbarp' = `tbar''   
        mat `tsar' = `sbarp'*`sbar'
        sca `tsar' = trace(`tsar')
        mat `tsk' = (trace(`sbarp'*`tbar')^2) / trace(`tbarp'*`tbar')
        sca `tsk' = trace(`tsk')
        mat `qt' = `tbarp'*`tbar'
        sca `qt' = trace(`qt')           
        mat `st' = `sbarp'*`tbar'
        sca `st' = trace(`st')
        sca `tsclr' = 0.5*(`tsar' - `qt' + sqrt((`tsar' + `qt')^2 -    /*
                        */      4*(`tsar'*`qt' - (`st')^2)))       
end

/* Computes TSCLR p-values.     */
program define tsclr_pcal

	args k qt tsclrstat tsclrpv

	tempname gamma pval  u s2 qs farg1 farg2 farg wt
	
    sca `gamma' = 2*exp(lngamma(`k'/2) - log(sqrt(_pi)) -  lngamma((`k'-1)/2))
				  
				  
	if("`k'" == "1") {
		sca `pval' = 1 - chi2(`k', `tsclrstat')
	}
	else if ("`k'"== "2") {
		local ni 100
		mat `u' = J(`ni'+1,1,0)
		mat `s2' = J(`ni'+1,1,0)
		mat `qs' = J(`ni'+1,1,0)
		mat `wt' = J(1,`ni'+1,2)
		mat `farg1' = J(`ni'+1,1,0)
		mat `qs'[1,1] = (`qt'+`tsclrstat')
		mat `farg1'[1,1] = `gamma'*chi2(`k',`qs'[1,1])
		forv i =1(1)`ni'{
			mat `u'[`i'+1,1] = `i'*_pi/2/`ni'
			mat `s2'[`i'+1,1] = sin(`u'[`i'+1,1])
			mat `qs'[`i'+1,1] = (`qt'+`tsclrstat') / 	/*
				*/ (1+(`qt'/`tsclrstat')*`s2'[`i'+1,1]*`s2'[`i'+1,1])
			mat `farg1'[`i'+1,1] = `gamma'*chi2(`k',`qs'[`i'+1,1])
		}
		mat `wt'[1,1] = 1
		mat `wt'[1,`ni'+1] = 1
		local ni = `ni'/2
		forv i =1(1)`ni'{
			mat `wt'[1,`i'*2] = 4
		}
		local ni = `ni'*2
		mat `wt' = `wt'*_pi/2/3/`ni'
		mat `pval' = `wt'*`farg1'
		sca `pval' = 1-trace(`pval')
	}
	else if ("`k'"== "3") {
		local ni 100
		mat `s2' = J(`ni'+1,1,0)
		mat `qs' = J(`ni'+1,1,0)
		mat `wt' = J(1,`ni'+1,2)
		mat `farg1' = J(`ni'+1,1,0)
		mat `qs'[1,1] = (`qt'+`tsclrstat')
		mat `farg1'[1,1] = `gamma'*chi2(`k',`qs'[1,1])
		forv i =1(1)`ni'{
			mat `s2'[`i'+1,1] = `i'/`ni'
			mat `qs'[`i'+1,1] = (`qt'+`tsclrstat') / 	/*
				*/ (1+(`qt'/`tsclrstat')*`s2'[`i'+1,1]*`s2'[`i'+1,1])
			mat `farg1'[`i'+1,1] = `gamma'*chi2(`k',`qs'[`i'+1,1])
		}
		mat `wt'[1,1] = 1
		mat `wt'[1,`ni'+1] = 1
		local ni = `ni'/2
		forv i =1(1)`ni'{
			mat `wt'[1,`i'*2] = 4
		}
		local ni = `ni'*2
		mat `wt' = `wt'/3/`ni'
		mat `pval' = `wt'*`farg1'
		sca `pval' = 1-trace(`pval')
	}
	else if ("`k'"== "4") {
		local eps .02
		local ni 100
		mat `s2' = J(`ni'+1,1,0)
		mat `qs' = J(`ni'+1,1,0)
		mat `wt' = J(1,`ni'+1,2)
		mat `farg' = J(`ni'+1,1,0)
		mat `farg1' = J(`ni'+1,1,0)
		mat `farg2' = J(`ni'+1,1,1)
		mat `qs'[1,1] = (`qt'+`tsclrstat')
		mat `farg1'[1,1] = `gamma'*chi2(`k',`qs'[1,1])
		mat `farg'[1,1] = `farg1'[1,1]*`farg2'[1,1]
		forv i = 1(1)`ni'{
			mat `s2'[`i'+1,1] = `i'/`ni'*(1-`eps')
			mat `qs'[`i'+1,1] = (`qt'+`tsclrstat') / 	/*
				*/ (1+(`qt'/`tsclrstat')*`s2'[`i'+1,1]*`s2'[`i'+1,1])
			mat `farg1'[`i'+1,1] = `gamma'*chi2(`k',`qs'[`i'+1,1])
			mat `farg2'[`i'+1,1] = sqrt(1-`s2'[`i'+1,1]*`s2'[`i'+1,1])
			mat `farg'[`i'+1,1] = `farg1'[`i'+1,1]*`farg2'[`i'+1,1]
		}
		mat `wt'[1,1] = 1
		mat `wt'[1,`ni'+1] = 1
		local ni = `ni'/2
		forv i = 1(1)`ni'{
			mat `wt'[1,`i'*2] = 4
		}
		local ni = `ni'*2
		mat `wt' = `wt'/3/`ni'*(1-`eps')
		mat `pval' = `wt'*`farg'
		sca `pval' = 1-trace(`pval')
		sca `s2' = 1-`eps'/2
		sca `qs' = (`qt'+`tsclrstat')/(1+(`qt'/`tsclrstat')*`s2'*`s2')
		sca `farg1' = `gamma'*chi2(`k',`qs')
		sca `farg2' = 0.5*(asin(1)-asin(1-`eps'))-(1-`eps') / /*
			*/	2*sqrt(1-(1-`eps')*(1-`eps'))
		sca `pval' = `pval'-`farg1'*`farg2'
	}
	else {
		local ni 100
		mat `s2' = J(`ni'+1,1,0)
		mat `qs' = J(`ni'+1,1,0)
		mat `wt' = J(1,`ni'+1,2)
		mat `farg' = J(`ni'+1,1,0)
		mat `farg1' = J(`ni'+1,1,0)
		mat `farg2' = J(`ni'+1,1,1)
		mat `qs'[1,1] = (`qt'+`tsclrstat')
		mat `farg1'[1,1] = `gamma'*chi2(`k',`qs'[1,1])
		mat `farg'[1,1] = `farg1'[1,1]*`farg2'[1,1]
		forv i =1(1)`ni'{
			mat `s2'[`i'+1,1] = `i'/`ni'
			mat `qs'[`i'+1,1] = (`qt'+`tsclrstat') / 	/*
				*/ (1+(`qt'/`tsclrstat')*`s2'[`i'+1,1]*`s2'[`i'+1,1])
			mat `farg1'[`i'+1,1] = `gamma'*chi2(`k',`qs'[`i'+1,1])
			if "`i'" == "`ni'"{
				mat `farg2'[`i'+1,1] = 0
			}
			else{
				mat `farg2'[`i'+1,1] = /*
				*/   (1-`s2'[`i'+1,1]*`s2'[`i'+1,1])^((`k'-3)/2)
			}
			mat `farg'[`i'+1,1] = `farg1'[`i'+1,1]*`farg2'[`i'+1,1]
		}
		mat `wt'[1,1] = 1
		mat `wt'[1,`ni'+1] = 1
		local ni = `ni'/2
		forv i = 1(1)`ni'{
			mat `wt'[1,`i'*2] = 4
		}
		local ni = `ni'*2
		mat `wt' = `wt'/3/`ni'
		mat `pval' = `wt'*`farg'
		sca `pval' = 1-trace(`pval')
	}
	
	sca `tsclrpv' = `pval'

end 

prog def inversefun

	args M k alpha C

	tempname eps a  b x fa fb tsclrstat fx

	sca `eps' = 0.000001
	sca `a' = `eps' 
	sca `b' = `M' - `eps'
	sca `tsclrstat'= `M' - `a'

	tsclr_pcal `k' `a' `tsclrstat' `fa'
	sca `tsclrstat' = `M' - `b'
	tsclr_pcal `k' `b' `tsclrstat' `fb'

	if(`fa' > `alpha') {
		sca `C' = `a'
	}
	else  if ( `fb' <`alpha') {
		sca `C' = `b'
	}
	else {
		while (`b'-`a'>`eps') {
			sca `x' = (`b'-`a')/2+`a'
			sca `tsclrstat'= `M'-`x'
			tsclr_pcal `k' `x' `tsclrstat' `fx'
			if (`fx' >`alpha') {
			 	  sca `b' = `x'
			}
			else {
				sca `a' = `x'
			}				 
		}
		sca `C' = `x'
	}
	  
end
	
prog def mat_inv_sqrt

        args in out
        tempname v vpri lam srlam

        loc k = rowsof(`in')
        mat symeigen `v' `lam' = `in'
        mat `vpri' = `v''
        mat `srlam' = diag(`lam')
        forv i = 1/`k' {
                mat `srlam'[`i', `i'] = 1/sqrt(`srlam'[`i', `i'])
        }
        mat `out' = `v'*`srlam'*`vpri'

end

program define myts2sls, rclass

        args ry1 ry2 rexog rinst cons touse samp

        quietly {
                tempname n1 n2 v2sigmasq bts2sls bts2sls11 dfr dfr2 dfm rss mss e1sigmasq var_unadj varinosol fstat2 rsq adjrsq ch ke Vx2het Vy1het var1het seb2shet constant ki
                tempvar w1hat residsq1 residsq2 
                
				reg `ry2' `rinst' `rexog' if `touse' & `samp' == 2, `cons'				
				sca `n2' = e(N)
				sca `dfr2' = e(df_r)	
				sca `v2sigmasq' = e(rss) / `dfr2'
				test `rinst'
				sca `fstat2' = r(F)
                predict double `w1hat' if `touse' & `samp' == 1, xb
				
				reg `ry1' `w1hat' `rexog' if `touse' & `samp' == 1, `cons'
				matrix `bts2sls' = e(b)
				sca `bts2sls11' = `bts2sls'[1,1]
				sca `dfr' = e(df_r)
				sca `dfm' = e(df_m)
				sca `n1' = e(N)
				sca `rsq' = e(r2)
                sca `adjrsq' = e(r2_a)
				sca `mss'=e(mss)
				sca `rss'=e(rss)
                matrix `var_unadj' = e(V)
				scalar `e1sigmasq' = `rss' / `dfr'
				matrix `varinosol' =  `var_unadj'*(1+`n1' / `n2'*`bts2sls11'^2*`v2sigmasq'/`e1sigmasq')
				ret sca rss = `rss'
				ret sca mss = `mss'
				ret sca dfr = `dfr'
				ret sca dfm = `dfm'
				ret sca dfr2 = `dfr2'
				ret sca fstat2 = `fstat2' 
				ret sca rsq = `rsq' 
                ret sca adjrsq =`adjrsq'
				ret sca n1 = `n1'
				ret sca n2 = `n2'
				ret mat beta = `bts2sls'
                ret mat var = `varinosol'			
        }

end

program define myts2sls_robust, rclass

        args ry1 ry2 rexog rinst cons touse samp

        quietly {
                tempname n1 n2 v2sigmasq bts2sls bts2sls11 dfr dfr2 dfm rss mss e1sigmasq var_unadj varinosol fstat2 rsq adjrsq ch ke Vx2het Vy1het var1het seb2shet constant ki
                tempvar w1hat residsq1 residsq2 
                
				reg `ry2' `rinst' `rexog' if `touse' & `samp' == 2, `cons'				
				sca `n2' = e(N)
				sca `dfr2' = e(df_r)	
				sca `v2sigmasq' = e(rss) / `dfr2'
				test `rinst'
				sca `fstat2' = r(F)
                predict double `w1hat' if `touse' & `samp' == 1, xb
				
				reg `ry1' `w1hat' `rexog' if `touse' & `samp' == 1, `cons'
				matrix `bts2sls' = e(b)
				sca `bts2sls11' = `bts2sls'[1,1]
				sca `dfr' = e(df_r)
				sca `dfm' = e(df_m)
				sca `n1' = e(N)
				sca `rsq' = e(r2)
                sca `adjrsq' = e(r2_a)
				sca `mss'=e(mss)
				sca `rss'=e(rss)
                matrix `var_unadj' = e(V)
				scalar `e1sigmasq' = `rss' / `dfr'
				matrix `varinosol' =  `var_unadj'*(1+`n1' / `n2'*`bts2sls11'^2*`v2sigmasq'/`e1sigmasq')
				ret sca rss = `rss'
				ret sca mss = `mss'
				ret sca dfr = `dfr'
				ret sca dfm = `dfm'
				ret sca dfr2 = `dfr2'
				ret sca fstat2 = `fstat2' 
				ret sca rsq = `rsq' 
                ret sca adjrsq =`adjrsq'
				ret sca n1 = `n1'
				ret sca n2 = `n2'
				ret mat beta = `bts2sls'
                ret mat var = `varinosol'
				
				****************************************************************
				/* Robust variance estimate of pix */
				****************************************************************
				gen `constant' = 1
				gmm (`ry2' - {xb: `rinst' `rexog' `constant'}) if `touse' & `samp' == 2, ///
				instruments(`rinst' `rexog') winit(unadjusted,independent) ///
				onestep deriv(1/xb = -1)
				mat `Vx2het' = e(V) 
				****************************************************************
				
				****************************************************************
				/* Robust variance estimate of piy */
				****************************************************************
				reg `ry1' `rinst' `rexog' if `touse' & `samp' == 1, `cons' robust
				mat `Vy1het' = e(V)*e(df_r)/_N
				****************************************************************

				****************************************************************
				/* Constructing C hat */
				****************************************************************

				loc ke : word count `rexog'
				loc ki : word count `rinst'
                
	  
				mat `ch' = J(`ke'+2,1,0)
				
                qui foreach v in `rinst' {
                        reg `v' `w1hat' `rexog' if `touse' & `samp' == 1, `cons'
						mat `ch' = `ch',e(b)'
                }					
				matrix `ch' = `ch'[1..`ke'+2,2..`ki'+1]
				mat `ch' = `ch',(J(1,`ke'+1,0)\I(`ke'+1))
				****************************************************************
				
				****************************************************************
				/* Calculating robust standard errors */	
				****************************************************************				
				mat `var1het' = `ch'*`Vy1het'*`ch'' + (`bts2sls11'' * `ch') *`Vx2het' * (`bts2sls11' * `ch'')
		        mat `seb2shet' = vecdiag(cholesky(diag(vecdiag(`var1het'))))'
				ret mat var1het = `var1het'
				ret mat Vy1het = `Vy1het'
				ret mat Vx2het = `Vx2het'
				ret mat seb2shet = `seb2shet'		
				****************************************************************				
			
        }

end


prog def myprint


        local level `e(level)'
        di
		di as text "Two-sample Instrumental variables (TS2SLS) regression"
  
        di
        di as text "First-stage F Results" _col(56) "Number of obs ="     /*
                */ as result %8.0f `e(n1)'
        di as text "{hline 23}" _col(56) "Wald chi2("`e(df_m)' ") =" /*
                */ as result %8.2f `e(chi2)'
        di as text "F(" %3.0f `e(numinst)' "," %6.0f         /*
                */ `e(df_r_first)' ") =" as result %8.2f        /*
                */ `e(F_first)' _col(56) as text "Prob > chi2   =" /*
                */ as result %8.4f chi2tail(`e(df_m)', `e(chi2)')     
        di as text "Prob > F      =" as result %8.4f            /*
                */ Ftail(`e(numinst)', `e(df_r_first)', `e(F_first)') /*
                */ _col(56) as text "R-squared     =" as result %8.4f `e(r2)'
        di as text       /*
                */ _col(56) as text "Adj R-squared =" as result %8.4f `e(r2_a)'
        di as text    /*
                */ _col(56) as text "Root MSE      =" as result         /*
                */ %8.4f `e(rmse)' 
        di
        eret display, level(`level')
        di as text "Instrumented:  " _c
        Disp `e(instd)'
        di as text "Instruments:   " _c
        if "`e(cons)'" == "yes" {
                Disp `e(insts)'
       }
       else {
               Disp `e(insts)' (No constant included)
       }
        di as text "Confidence set and p-value for "                    /*
                */ abbrev("`e(instd)'", 13) " "                         /*
                */ "are based on normal approximation,"
	    di as text "thus not robust to weak instruments."
        di as text "{hline 78}"
        
        
end

prog def myprint_robust


        local level `e(level)'
        di
		di as text "Two-sample Instrumental variables (TS2SLS) regression"
  
        di
        di as text _col(56) "Number of obs ="     /*
                */ as result %8.0f `e(n1)'
        di as text _col(56) "Wald chi2("`e(df_m)' ") ="  /*
                */ as result %8.2f `e(chi2)'
        di as text _col(56) as text "Prob > chi2   =" /*
                */ as result %8.4f chi2tail(`e(df_m)', `e(chi2)')     
        di as text _col(56) as text "R-squared     =" as result %8.4f `e(r2)'
 *       di as text       /*
 *              */ _col(56) as text "Adj R-squared =" as result %8.4f `e(r2_a)'
        di as text    /*
                */ _col(56) as text "Root MSE      =" as result         /*
                */ %8.4f `e(rmse)' 
        di
        eret display, level(`level')
        di as text "Instrumented:  " _c
        Disp `e(instd)'
        di as text "Instruments:   " _c
        if "`e(cons)'" == "yes" {
                Disp `e(insts)'
       }
       else {
               Disp `e(insts)' (No constant included)
       }
        di as text "Confidence set and p-value for "                    /*
                */ abbrev("`e(instd)'", 13) " "                         /*
                */ "are based on normal approximation,"
	    di as text "thus not robust to weak instruments."
        di as text "{hline 78}"
        
        
end

prog def myprintregions_robust

		args tsar_cset tsk_cset tsclr_cset tsar_p tsk_p tsclr_p 
        local level `e(level)'
        di _n _n
        di as text "{hline 78}"
        local what "`level'% confidence set"
        local plural

        local title : di "Weak IV Robust `what' and p-value"
        local titcol `=int((78 - length("`title'"))/2 + 1)'
        di as text _col(`titcol') "`title'"
        
        local title : di "for H0: _b[`e(instd)'] = " %-9.0g e(H0_b)
        local titcol `=int((78 - length("`title'"))/2 + 1)'
        di as text _col(`titcol') "`title'"
		
		local testnamelen = 5
		local name_tsclr : di "{txt}{ralign `testnamelen':Robust TSCLR}"
		local name_tsar : di "{txt}{ralign `testnamelen':Robust TSAR}"
		local name_tsk : di "{txt}{ralign `testnamelen':Robust TSK}"
		
		* calculate length of result parts
			local testlist "tsar tsk tsclr"
			foreach testname in `testlist' {
			   	local pvaltxt_`testname' : di "{res}" %8.4f ``testname'_p'
			   	local testtxt_`testname' "`pvaltxt_`testname''"
			}
	

		* print
        local titcol `=int((78 - length("`title'"))/2 + 1)'
		di			
        di as text "{hline 78}"
		local title : di "`level'% Confidence Set"
        local titcol `=int((70 - 15 - length("`title'"))/2) + 15'        
		di as text _col(2) "Test" _col(`titcol') "`title'" _col(71) "p-value"
		di as text "{hline 78}"
		di _col(2) "`name_tsclr'" _col(`titcol') "`tsclr_cset'" _col(70) "`testtxt_tsclr'{res}"
        di _col(2) "`name_tsar'" _col(`titcol') "`tsar_cset'" _col(70) "`testtxt_tsar'{res}"
		di _col(2) "`name_tsk'" _col(`titcol') "`tsk_cset'" _col(70) "`testtxt_tsk'{res}"
		di as text "{hline 78}"


end

prog def myprintregions_robust_just

		args tsclr_cset tsclr_p
        local level `e(level)'
        di _n _n
        di as text "{hline 78}"
        local what "`level'% confidence set"
        local title : di "Weak IV Robust `what' and p-value"
        local titcol `=int((78 - length("`title'"))/2 + 1)'
        di as text _col(`titcol') "`title'"
        local title : di "for H0: _b[`e(instd)'] = " %-9.0g e(H0_b)
        local titcol `=int((78 - length("`title'"))/2 + 1)'
        di as text _col(`titcol') "`title'"
		
		* column specs
			local testnamelen = 5
			
		* testnames
			local name_tsclr : di "{txt}{ralign `testnamelen':Robust TSCLR}"
			local testname "tsclr"
			local pvaltxt_`testname' : di "{res}" %8.4f ``testname'_p'
			local testtxt_`testname' "`pvaltxt_`testname''"
			
		* print
        local titcol `=int((78 - length("`title'"))/2 + 1)'
        
		di			
        di as text "{hline 78}"
		local title : di "`level'% Confidence Set"
        local titcol `=int((70 - 15 - length("`title'"))/2) + 15'        
		di as text _col(2) "Test" _col(`titcol') "`title'" _col(71) "p-value"
		di as text "{hline 78}"
		di _col(2) "`name_tsclr'" _col(`titcol') "`tsclr_cset'" _col(70) "`testtxt_tsclr'{res}"
		di as text "{hline 78}"
		di as text "Note: In the just identifed case, TSCLR = TSAR = TSK." 
end

program define myprintregions

        local level `e(level)'
        di _n _n
        di as text "{hline 78}"
        local what "`level'% confidence set"
        local title : di "Weak IV Robust `what's and p-values"
        local titcol `=int((78 - length("`title'"))/2 + 1)'
        di as text _col(`titcol') "`title'"
        
        local title : di "for H0: _b[`e(instd)'] = " %-9.0g e(H0_b)
        local titcol `=int((78 - length("`title'"))/2 + 1)'
        di as text _col(`titcol') "`title'"
 		        
        di
        di as text "{hline 78}"
        local title : di proper("`what'")
        local titcol `=int((70 - 15 - length("`title'"))/2) + 15'
        di as text _col(2) "Test" _col(`titcol') "`title'" _col(71) "p-value"
        di as text "{hline 78}" 
 
        local tsclrpval  "as result _col(72) %6.4f e(p_TSCLR)"
        if `e(TSCLR_type)' == 1 {
                local tsclrint "empty"
        }
        else if `e(TSCLR_type)' == 2 {
                local tsclrint : di "["                                    /*
                        */       %9.0g e(TSCLR_x1)                         /*
                        */       ", "                                   /*
                        */       %9.0g e(TSCLR_x2)                         /*
                        */       "]"
        }
        else if `e(TSCLR_type)' == 3 {
                local tsclrint "(-inf, +inf)"
        }
        else {
                        local tsclrint : di "(-inf, "                      /*
                                */       %9.0g e(TSCLR_x1)                 /*
                                */       "] U ["                        /*
                                */       %9.0g e(TSCLR_x2)                 /*
                                */       ", +inf)"  
        }
        local len `=length(`"`tsclrint'"')'
        local len `=int((70 - 15 - `len')/2)'
        local len `=`len'+15'
        di _col(2) "Benchmark TSCLR" _col(`len') as result "`tsclrint'" `tsclrpval'

                local tsarpval "as result _col(72) %6.4f e(p_TSAR)"
                if `e(TSAR_type)' == 1 {
                        local tsarint "empty"
                }
                else if `e(TSAR_type)' == 2 {
                        local tsarint : di "["                            /*
                                */       %9.0g e(TSAR_x1)                 /*
                                */       ", "                           /*
                                */       %9.0g e(TSAR_x2)                 /*
                                */       "]"
                }
                else if `e(TSAR_type)' == 3 {
                        local tsarint "(-inf, +inf)"
                }
                else {
                                local tsarint : di "(-inf, "              /*
                                        */       %9.0g e(TSAR_x1)         /*
                                        */       "] U ["                /*
                                        */       %9.0g e(TSAR_x2)         /*
                                        */       ", +inf)" 
                        
                       
                }
                local len `=length(`"`tsarint'"')'
                local len `=int((70 - 15 - `len')/2)'
                local len `=`len'+15'
              di _col(2) as text "Benchmark TSAR"                     /*
                        */      _col(`len') as result "`tsarint'" `tsarpval' 
		

                local tskpval "as result _col(72) %6.4f e(p_TSK)"
                if `e(TSK_type)' == 2 {
                        local tskint : di "["                            /*
                                */       %9.0g e(TSK_x1)                 /*
                                */       ", "                           /*
                                */       %9.0g e(TSK_x2)                 /*
                                */      "]"
                }
                else if `e(TSK_type)' == 3 {
                        local tskint "(-inf, +inf)"
                }
                else if `e(TSK_type)' == 4 {
                                local tskint : di "(-inf, "              /*
                                        */       %9.0g e(TSK_x1)         /*
                                        */       "] U ["                /*
                                        */       %9.0g e(TSK_x2)         /*
                                        */       ", +inf)" 
                        
                       
                }
                else if `e(TSK_type)' == 5 {
                                local tskint : di "(-inf, "              /*
                                        */       %6.0g e(TSK_x1)         /*
                                        */       "] U ["                /*
                                        */       %6.0g e(TSK_x2)         /*
                                        */       ", "                   /*
                                        */       %6.0g e(TSK_x3)         /*
                                        */       "] U ["                /*
                                        */       %6.0g e(TSK_x4)         /*
                                        */       ", +inf)"      
                        
                }
                else if `e(TSK_type)' == 6 & e(TSK_x1) < e(TSK_x3) {
                                local tskint : di "["                    /*
                                        */       %9.0g e(TSK_x1)         /*
                                        */       ", "                   /*
                                        */       %9.0g e(TSK_x2)         /*
                                        */       "] U ["                /*
                                        */       %9.0g e(TSK_x3)         /*
                                        */       ", "                   /*
                                        */       %9.0g e(TSK_x4)         /*
                                        */       "]"
                        
                        
                }
                else {
                                local tskint : di "["                    /*
                                        */       %9.0g e(TSK_x3)         /*
                                        */       ", "                   /*
                                        */       %9.0g e(TSK_x4)         /*
                                        */       "] U ["                /*
                                        */       %9.0g e(TSK_x1)         /*
                                        */       ", "                   /*
                                        */       %9.0g e(TSK_x2)         /*
                                        */       "]"                    
                        
                        
                }
                local len `=length(`"`tskint'"')'
                local len `=int((70 - 15 - `len')/2)'
                local len `=`len'+15'
            di _col(2) as text "Benchmark TSK"                         /*
                        */ _col(`len') as result "`tskint'" `tskpval' 
     
        di as text "{hline 78}"
end		

program define myprintregions_just

        local level `e(level)'
        di _n _n
        di as text "{hline 78}"
        local what "`level'% confidence set"
        local title : di "Weak IV Robust `what' and p-value"
        local titcol `=int((78 - length("`title'"))/2 + 1)'
        di as text _col(`titcol') "`title'"
        
        local title : di "for H0: _b[`e(instd)'] = " %-9.0g e(H0_b)
        local titcol `=int((78 - length("`title'"))/2 + 1)'
        di as text _col(`titcol') "`title'"
 		        
        di
        di as text "{hline 78}"
        local title : di proper("`what'")
        local titcol `=int((70 - 15 - length("`title'"))/2) + 15'
        di as text _col(2) "Test" _col(`titcol') "`title'" _col(71) "p-value"
        di as text "{hline 78}" 
 
        local tsclrpval  "as result _col(72) %6.4f e(p_TSCLR)"
        if `e(TSCLR_type)' == 1 {
                local tsclrint "empty"
        }
        else if `e(TSCLR_type)' == 2 {
                local tsclrint : di "["                                    /*
                        */       %9.0g e(TSCLR_x1)                         /*
                        */       ", "                                   /*
                        */       %9.0g e(TSCLR_x2)                         /*
                        */       "]"
        }
        else if `e(TSCLR_type)' == 3 {
                local tsclrint "(-inf, +inf)"
        }
        else {
                        local tsclrint : di "(-inf, "                      /*
                                */       %9.0g e(TSCLR_x1)                 /*
                                */       "] U ["                        /*
                                */       %9.0g e(TSCLR_x2)                 /*
                                */       ", +inf)"                 
        }
        local len `=length(`"`tsclrint'"')'
        local len `=int((70 - 15 - `len')/2)'
        local len `=`len'+15'
        di _col(2) "Benchmark TSCLR" _col(`len') as result "`tsclrint'" `tsclrpval'     
        di as text "{hline 78}"
        di as text "Note: In the just identifed case, TSCLR = TSAR = TSK." 
end		

program define myprintregions_footnote
        local points `e(points)'
		local grid `e(grid)'
        di as text "Confidence sets estimated for `points' points in `grid'."
end

	
/* Lifted straight out of ivreg.ado.    */
program define IsStop, sclass

        /* sic, must do tests one-at-a-time, 0 may be very large */
        if `"`0'"' == "[" {
                sret local stop 1
                exit
        }
        if `"`0'"' == "," {
                sret local stop 1
                exit
        }
        if `"`0'"' == "if" {
                sret local stop 1
                exit
        }
        if `"`0'"' == "in" {
                sret local stop 1
                exit
        }
        if `"`0'"' == "" {
                 sret local stop 1
                 exit
        }
        else     sret local stop 0

end

/* More codes borrowed from ivreg.ado.   */
program define Disp
        local first ""
        local piece : piece 1 64 of `"`0'"'
        local i 1
        while "`piece'" != "" {
                di in gr "`first'`piece'"
                local first "              "
                local i = `i' + 1
                local piece : piece `i' 64 of `"`0'"'
        }
        if `i'==1 { 
                di 
        }
end

***************************************	
*** Compute Test Statstics (robust) ***
***************************************
program computeivtests_robust, rclass
	syntax [, n1(name) n2(name) hateta(name) hatpi(name) sigmaf(name) sigmas(name) null(string) *]
    tempname SIGMA sqrtSIGMA TSARstat hatD ZZ PK TSKstat hatQt ZZZ ZZZZ sqrtZZZZ /*
	       */TSCLRstat rk tsclr_stat tsar_chi2 tsk_chi2 
	  mat `SIGMA' = `sigmaf' + `n1' / `n2' * (`null')^2 * `sigmas'
	  mat_inv_sqrt `SIGMA' `sqrtSIGMA'						
	  mat `TSARstat' = `n1' * (`hateta' - `hatpi' * `null') * invsym(`SIGMA') * (`hateta' - `hatpi' * `null')' 
	  mat `hatD' = -(`hatpi'' + (`n1'/`n2') * `null' * `sigmas' * invsym(`SIGMA') * (`hateta' - `hatpi' * `null')')		
	  mat `ZZ' = (`sqrtSIGMA' * `hatD')' * (`sqrtSIGMA' * `hatD')
	  mat `PK' = (`sqrtSIGMA' * `hatD') * invsym(`ZZ') *  (`sqrtSIGMA' * `hatD')' 
	  mat `TSKstat' = `n1' * (`sqrtSIGMA' * (`hateta' - `hatpi' * `null')')'  * (`sqrtSIGMA' * (`hateta' - `hatpi' * `null')')
	  mat `ZZZ' = (`n1'/`n2') * `sigmas' - (`n1'/`n2')^2 * `null'^2 * `sigmas' * invsym(`SIGMA') * `sigmas' 
	  mat `hatQt' = `n1' * `hatD'' * invsym(`ZZZ') * `hatD'
	  mat `ZZZZ' = (`TSARstat' + `hatQt') * (`TSARstat' + `hatQt') - 4 * `TSARstat' * `hatQt' + 4 * `TSKstat' * `hatQt'				
	  mat_inv_sqrt `ZZZZ' `sqrtZZZZ'
	  mat `TSCLRstat' = (`TSARstat' - `hatQt' + invsym(`sqrtZZZZ'))/2				
	  sca `rk' = `hatQt'[1,1]
	  sca `tsclr_stat' = `TSCLRstat'[1,1]
	  sca `tsar_chi2' = `TSARstat'[1,1]
	  sca `tsk_chi2' = `TSKstat'[1,1]
	* return tests in r()
	  return scalar tsar_chi2 = `tsar_chi2'
	  return scalar tsk_chi2 = `tsk_chi2'
	  return scalar tsclr_stat = `tsclr_stat'
	  return scalar rk = `rk'
end						

***************************************	
*** Compute P-values       (robust) ***
***************************************
program compute_pvals
	syntax [, null(string) rk(name) numk(string) level(string) ///
		tsar_p(name) tsar_chi2(name) tsar_df(name) ///
		tsk_p(name) tsk_chi2(name) tsk_df(string) ///
		tsclr_p(name) tsclr_stat(name) tsclr_df(name) ///
		tsar_r(name) tsk_r(name) tsclr_r(name)  *]
	scalar `tsk_df' = 1
	scalar `tsar_df' = `numk'
	scalar `tsar_p'= chi2tail(`tsar_df',`tsar_chi2')
	scalar `tsk_p'= chi2tail(1,`tsk_chi2')
	scalar `tsclr_df' = .
	* below is a subprogram for TSCLR p-value taken from Mikusheva and Poi's codes
			if `tsclr_stat'==.		scalar `tsclr_p'=.
			else {
				tsclr_pcal `numk' `rk' `tsclr_stat' `tsclr_p'
				* fix negative p-value approximations that occur because of rounding near zero 
					if `tsclr_p'<=-0.00001 & `tsclr_p'>-9999999		n di in red "error when approximating TSCLR p-value for null = `null'"	
					if `tsclr_p'<0								scalar `tsclr_p'=0.0000
			}
	* compute reject/dn reject binary
		scalar `tsar_r' = cond(`tsar_p'<=1-`level'/100,1,0)
		scalar `tsk_r' = cond(`tsk_p'<=1-`level'/100,1,0)
			if `tsclr_stat'==.		scalar `tsclr_r' = .
			else					scalar `tsclr_r' = cond(`tsclr_p'<=1-`level'/100,1,0)
end
 
exit
