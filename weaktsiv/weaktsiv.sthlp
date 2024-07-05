{smcl}
{* *! version 1.0.0  26feb2018}{...}
{cmd:help weaktsiv}{right: ({browse "https://doi.org/10.1177/1536867X19874235":SJ19-3: st0568})}
{hline}

{title:Title}

{p2colset 5 17 19 2}{...}
{p2col :{cmd:weaktsiv} {hline 2}}Two-sample instrumental-variables regression
with potentially weak instruments{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}
{cmd:weaktsiv} {it:depvar} {it:varlist_exog} {cmd:(}{it:varlist_endog} {cmd:=}
{it:varlist_iv}{cmd:)} {ifin} [{cmd:,} {cmd:nocons} {cmd:robust}
{opt level(#)} {opt test(#)} {opt points(#)}
{cmd:grid(}{it:#}{cmd:(}{it:#}{cmd:)}{it:#}{cmd:)}]

{phang2}
{it:depvar} is the outcome variable.

{phang2}
{it:varlist_exog} is the list of exogenous variables.

{phang2}
{it:varlist_endog} is the endogenous regressor of the model.

{phang2}
{it:varlist_iv} is the list of exogenous variables used together with
{it:varlist_exog} as instruments for {it:varlist_endog}.


{title:Description}

{p 4 4 2}
{cmd:weaktsiv} implements two-sample instrumental-variables (IV) regression
models with one endogenous regressor and potentially weak instruments of Choi,
Gu, and Shen (2018).

{p 4 4 2}
The default calculates the benchmark two-stage Anderson-Rubin (TSAR),
two-stage Kleibergen (TSK), and two-stage conditional likelihood-ratio (TSCLR)
tests and confidence sets.  If the {cmd:robust} option is specified,
{cmd:weaktsiv} calculates the fully robust versions of the tests and
confidence sets.

{p 4 4 2}
The command also reports results from two-sample two-stage least-squares
(TS2SLS) estimation, which is valid only with strong instrumental variables.


{title:Options}

{phang}
{cmd:nocons} suppresses the constant term in the regression model.

{phang}
{cmd:robust} provides versions of two-sample weak IV robust tests that are
also robust to heteroskedasticity and unequal moments of excluded instruments
and exogenous regressors across the two data samples.  {cmd:robust} also
reports Pacini and Windmeijer's (2016) robust standard error following the
TS2SLS estimation.

{phang}
{cmd:level(}{it:#}{cmd:)} sets the confidence level.  The default is
{cmd:level(95)}.

{phang}
{cmd:test(}{it:#}{cmd:)} sets the hypothesized value of the endogenous
variable's coefficient.  The default is {cmd:test(0)}.

{phang}
{opt points(#)} specifies the number of points used to create the grid for
confidence region calculation.  {cmd:points()} may be used together only with
the {cmd:robust} option and cannot be used together with the {cmd:grid()}
option.  The default is {cmd:points(100)}.

{phang}
{cmd:grid(}{it:#}{cmd:(}{it:#}{cmd:)}{it:#}{cmd:)} specifies the grid used for
confidence region calculation.  {cmd:grid()} may be used together only with
the {cmd:robust} option because the benchmark confidence regions are
calculated analytically.  The default uses the TS2SLS estimator plus or minus
two times the Pacini and Windmeijer (2016) standard error and 100 grid points.


{marker s_examples}{...}
{title:Examples}

{pstd}Load Currie and Yelowitz (2000) data{p_end}
{phang2}{bf:. {stata "use sample1.dta, clear":use sample1.dta}}{p_end}

{pstd}Implement the {cmd:weaktsiv} command in the just-identified case{p_end}
{phang2}{bf:. {stata "weaktsiv ry1 h* p* b* (ry2=z)":weaktsiv ry1 h* p* b* (ry2=z)}}{p_end}

{pstd}Load Olivetti and Paserman (2015) data{p_end}
{phang2}{bf:. {stata "use sample2.dta, clear":use sample2.dta}}{p_end}

{pstd}Implement the {cmd:weaktsiv} command in the overidentified case{p_end}
{phang2}{bf:. {stata "weaktsiv ry1 (ry2=z*), level(90)":weaktsiv ry1 (ry2=z*), level(90)}}{p_end}


{title:Stored results}

{pstd}
{cmd:weaktsiv} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(p_TSAR)}}TSAR test p-value{p_end}
{synopt:{cmd:e(p_TSK)}}TSK test p-value{p_end}
{synopt:{cmd:e(p_TSCLR)}}TSCLR test p-value{p_end}
{synopt:{cmd:e(TSAR_xi)}}endpoints of benchmark TSAR confidence sets{p_end}
{synopt:{cmd:e(TSK_xi)}}endpoints of benchmark TSK confidence sets{p_end}
{synopt:{cmd:e(TSCLR_xi)}}endpoints of benchmark TSCLR confidence sets{p_end}
{synopt:{cmd:e(level)}}confidence level for weak IV robust inference{p_end}
{synopt:{cmd:e(H0_b)}}value of beta under null for weak IV robust
inference{p_end}
{synopt:{cmd:e(n1)}}{it:#} of observations in sample 1 (the outcome sample){p_end}
{synopt:{cmd:e(n2)}}{it:#} of observations in sample 2 (the endogenous variable sample){p_end}
{synopt:{cmd:e(chi2)}}TS2SLS Wald statistic{p_end}
{synopt:{cmd:e(F_first)}}TS2SLS first-stage F{p_end}
{synopt:{cmd:e(numinst)}}number of instruments{p_end}
{synopt:{cmd:e(df_r_first)}}TS2SLS first-stage residual degrees of
freedom{p_end}
{synopt:{cmd:e(df_m)}}TS2SLS model degrees of freedom{p_end}
{synopt:{cmd:e(df_r)}}TS2SLS residual degrees of freedom{p_end}
{synopt:{cmd:e(r2)}}R^2{p_end}
{synopt:{cmd:e(r2_a)}}adjusted R^2{p_end}
{synopt:{cmd:e(mss)}}TS2SLS model  sum of squares{p_end}
{synopt:{cmd:e(rss)}}TS2SLS residual sum of squares{p_end}
{synopt:{cmd:e(rmse)}}TS2SLS root  mean squared errors{p_end}
{synopt:{cmd:e(points)}}{it:#} of grid points for robust confidence sets{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:weaktsiv}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(robust)}}whether the {cmd:robust} option is specified{p_end}
{synopt:{cmd:e(TSAR_type)}}type of benchmark  TSAR confidence set{p_end}
{synopt:{cmd:e(TSK_type)}}type of benchmark TSK confidence set{p_end}
{synopt:{cmd:e(TSCLR_type)}}type of benchmark TSCLR confidence set{p_end}
{synopt:{cmd:e(TSAR_cset)}}robust TSAR confidence set{p_end}
{synopt:{cmd:e(TSK_cset)}}robust TSK confidence set{p_end}
{synopt:{cmd:e(TSCLR_cset)}}robust TSCLR confidence set{p_end}
{synopt:{cmd:e(grid)}}grid range for robust  confidence set{p_end}
{synopt:{cmd:e(cons)}}whether constants are used{p_end}
{synopt:{cmd:e(instd)}}instrumented variable{p_end}
{synopt:{cmd:e(insts)}}instruments{p_end}
{synopt:{cmd:e(exog)}}exogenous variables{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}


{marker references}{...}
{title:References}

{marker CJS2018}{...}
{phang}
Choi, J., J. Gu, and S. Shen. 2018. Weak-instrument robust inference for
two-sample instrumental variables regression. {it:Journal of Applied Economics}
33: 109-125.

{phang}
Currie, J., and A. Yelowitz. 2000. Are public housing projects good for kids? 
{it:Journal of Public Economics} 75: 99-124.

{marker Inoue2010}{...}
{phang}
Inoue, A., and G. Solon. 2010. Two-sample instrumental variables estimators. 
{it:Review of Economics and Statistics} 92: 557-561.

{phang}
Olivetti, C., and M. D. Paserman. 2015. In the name of the son (and the
daughter): Intergenerational mobility in the United States, 1850-1940.
{it:American Economic Review} 105: 2695-2724.

{marker Pacini2016}{...}
{phang}
Pacini, D., and F. Windmeijer. 2016. Robust inference for the two-sample 2SLS
estimator. {it:Economics Letters} 146: 50-54.


{marker authors}{...}
{title:Authors}

{pstd}
Jaerim Choi{break}
University of Hawaii at Manoa{break}
Honolulu, HI{break}
choijm@hawaii.edu

{pstd}
Shu Shen{break}
University of California, Davis{break}
Davis, CA{break}
shushen@ucdavis.edu


{marker alsosee}{...}
{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume 19, number 3: {browse "https://doi.org/10.1177/1536867X19874235":st0568}{p_end}
