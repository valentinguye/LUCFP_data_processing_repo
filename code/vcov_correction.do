use "temp_data/temp_sample.dta", clear

ivregress 2sls lucpfip_pixelcount_total n_reachable_uml_lag1 (wa_ffb_price_imp1_lag1 = wa_prex_cpo_imp1_lag1)
estat vce
di e(rmse)
mat v2sls = e(V)
  
// First stage reg
reg wa_ffb_price_imp1_lag1 wa_prex_cpo_imp1_lag1 n_reachable_uml_lag1
predict double fit_endo, xb

// Second stage reg
reg lucpfip_pixelcount_total fit_endo n_reachable_uml_lag1
scalar rmsebyhand = e(rmse)

// the 'wrong' VCE, calculated from the instruments
mat vbyhand = e(V)
scalar dfk = e(df_r)

// the correct resids: orig regressors * second stage coeffs 
gen double eps2 = (lucpfip_pixelcount_total - _b[fit_endo]*wa_ffb_price_imp1_lag1 - _b[n_reachable_uml_lag1]*n_reachable_uml_lag1 - _b[_cons])^2
qui su eps2

// corrected RMSE, based on the correct resids
scalar rmsecorr = sqrt(r(sum) / dfk)

// corrected VCE, using the right s^2
mat vcorr = (rmsecorr / rmsebyhand)^2 * vbyhand
mat li vcorr

// check to see that it equals the real 2SLS VCE
mat diff = v2sls - vcorr
mat li diff
/*
* original version 
sysuse auto, clear
ivregress 2sls price headroom (weight = turn foreign)
estat vce
di e(rmse)
mat v2sls = e(V)
  
// First stage reg
qui reg weight turn foreign headroom
predict double what, xb

// Second stage reg
qui reg price what headroom
scalar rmsebyhand = e(rmse)

// the 'wrong' VCE, calculated from the instruments
mat vbyhand = e(V)
scalar dfk = e(df_r)

// the correct resids: orig regressors * second stage coeffs 
gen double eps2 = (price - _b[what]*weight - _b[headroom]*headroom - _b[_cons])^2
qui su eps2

// corrected RMSE, based on the correct resids
scalar rmsecorr = sqrt(r(sum) / dfk)

// corrected VCE, using the right s^2
mat vcorr = (rmsecorr / rmsebyhand)^2 * vbyhand
mat li vcorr

// check to see that it equals the real 2SLS VCE
mat diff = v2sls - vcorr
mat li diff
*/

*******************************************************
* NOW TEST ppmlhdfe 
cap ado uninstall ftools
cap ado uninstall reghdfe
cap ado uninstall ppmlhdfe

ssc install ftools
ssc install reghdfe
ssc install ppmlhdfe

clear all
ftools, compile
reghdfe, compile

* Test program
sysuse auto, clear
reghdfe price weight, a(turn)
ppmlhdfe price weight, a(turn)

* Install ftools (remove program if it existed previously)
cap ado uninstall ftools
net install ftools, from("https://raw.githubusercontent.com/sergiocorreia/ftools/master/src/")

* Install reghdfe
cap ado uninstall reghdfe
net install reghdfe, from("https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/src/")

* Install boottest (Stata 11 and 12)
if (c(version)<13) cap ado uninstall boottest
if (c(version)<13) ssc install boottest

* Install moremata (sometimes used by ftools but not needed for reghdfe)
cap ssc install moremata

* Install ivreg2, the core package
cap ado uninstall ivreg2
ssc install ivreg2

* Finally, install this package
cap ado uninstall ivreghdfe
net install ivreghdfe, from(https://raw.githubusercontent.com/sergiocorreia/ivreghdfe/master/src/)

ssc install ranktest


******************************************



use "temp_data/temp_sample.dta", clear

*ppmlhdfe 
ivreghdfe lucpfip_pixelcount_total n_reachable_uml_lag1 (wa_ffb_price_imp1_lag1 = wa_prex_cpo_imp1_lag1), absorb(parcel_id) small

ivreghdfe lucpfip_pixelcount_total n_reachable_uml_lag1 (wa_ffb_price_imp1_lag1 = wa_prex_cpo_imp1_lag1), absorb(parcel_id) cluster(parcel_id) 

scalar rmse = e(rmse)




* Now CONTROL FUNCTION APPROACH 
*https://www.statalist.org/forums/forum/general-stata-discussion/general/1381373-ivpoisson-with-panel-data-fixed-effects
use "temp_data/temp_sample.dta", clear
set seed 1234
capture drop v2h_fe_xt
xtset parcel_id year
xtreg wa_ffb_price_imp1_lag1 wa_prex_cpo_imp1_lag1 n_reachable_uml_lag1, fe
predict double v2h_fe_xt, e
xtpoisson lucpfip_pixelcount_total wa_ffb_price_imp1_lag1 v2h_fe_xt n_reachable_uml_lag1, fe vce(bootstrap)
mat li e(ci_normal)
mat li e(ci_bc)   

* with reghdfe, no vce(bootstrap), so use bootstrap function with a program 
* (inspired from D. Cameron's code  http://cameron.econ.ucdavis.edu/nhh2017/norway01_count_part2.pdf)
use "temp_data/temp_sample.dta", clear
capture program drop endogtwostep
program endogtwostep, eclass 
 version 14.0
 tempname b 
 capture drop v2h_fe_xt
 capture drop v2h_fe
 reghdfe wa_ffb_price_imp1_lag1 wa_prex_cpo_imp1_lag1 n_reachable_uml_lag1, absorb(parcel_id) residuals(v2h_fe)
 ppmlhdfe lucpfip_pixelcount_total wa_ffb_price_imp1_lag1 v2h_fe n_reachable_uml_lag1, absorb(parcel_id) vce(robust)
 matrix `b' = e(b)
 ereturn post `b'
 end

 bootstrap _b, reps(50) seed(1234) nowarn : endogtwostep

mat li e(ci_normal)
mat li e(ci_bc)    

/* NOTES: 
Clustering (both stages, and especially the 2nd) does not change the final SE bottstrap estimate
Both methods (xt* or *hdfe) return non bias-corrected CI. 
They don't yield the same results although the seed is the same... 
*/


* so after reghdfe, the e option does not work in predict 
use "temp_data/temp_sample.dta", clear
 capture drop v2h_fe
 reghdfe wa_ffb_price_imp1_lag1 wa_prex_cpo_imp1_lag1 n_reachable_uml_lag1, absorb(parcel_id) residuals(v2h_fe)
 ppmlhdfe lucpfip_pixelcount_total wa_ffb_price_imp1_lag1 v2h_fe n_reachable_uml_lag1, absorb(parcel_id) vce(robust)





sysuse auto, clear
ivreghdfe price weight (length=gear), absorb(turn trunk, tol(1e-6) accel(sd))
