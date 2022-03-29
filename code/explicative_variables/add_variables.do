

***** MACROECONOMIC VARIABLES *****
import excel "input_data/macro/fob_dom_prices_cpo.xlsx", firstrow clear
*sheet("Sheet1")

replace export_tax_effective = export_tax_effective/100

gen double spread_monthly = fob_blwn_cpo - dom_blwn_cpo  - export_tax_effective*fob_blwn_cpo
label variable spread_monthly "spread_cpo_fobblwn_domblwn"

* make annual average of monthhly data 
sort year 
bys year: egen spread = mean(spread_monthly)

/*
bys year: egen fob_blwn_cpo_y = mean(fob_blwn_cpo)
bys year: egen dom_blwn_cpo_y = mean(dom_blwn_cpo)
bys year: egen export_tax_effective_y = mean(export_tax_effective)
*/

duplicates drop year, force

keep year spread

drop if year > 2015

*drop Domesticprice FOBPrice
capture mkdir "temp_data/processed_macro"
save "temp_data/processed_macro/spread_cpo_fobblwn_domblwn.dta", replace

* read IBS panel 
use "temp_data/processed_mill_geolocalization/IBS_UML_panel.dta", clear

merge m:1 year using "temp_data/processed_macro/spread_cpo_fobblwn_domblwn.dta", nogenerate 
sort firm_id year 
*need to change the global for all new price variables each time (and don't forget to NOT DELFATE taxeffectiverate). 
*global prices ref_int_cpo_price cif_rtdm_cpo dom_blwn_cpo fob_blwn_cpo spread_int_dom_paspi rho dom_blwn_pko cif_rtdm_pko
global prices spread

*label variable spread_int_dom_paspi "Paspi calculation of fob_blwn_cpo minus dom_blwn_cpo"			

*DEFLATE 
foreach var of varlist $prices {
	quietly replace `var'=`var'/wpi_2010
}

/*
graph hbox cpo_price_imp2 in_tot_cpo_price_imp2 if in_tot_cpo_price_imp2<.& cpo_price_imp2<.
graph hbox cpo_price_imp1 in_tot_cpo_price_imp1 if in_tot_cpo_price_imp1<.& cpo_price_imp1<.
*celles qui ont à la fois un prix d'achat et un prix de vente de cpo ont le premier globalement plus bas. 
*/

/*
replace taxeffectiverate = taxeffectiverate/100

gen double spread1 = cif_rtdm_cpo - dom_blwn_cpo - taxeffectiverate*cif_rtdm_cpo
label variable spread1 "spread_cpo_rtdm_blwn"
gen double spread2 = fob_blwn_cpo - dom_blwn_cpo 
gen double spread3 = spread_int_dom_paspi - rho- taxeffectiverate*fob_blwn_cpo
label variable spread2 "spread_cpo_fobblwn_blwn"
label variable spread3 "spread"
gen double spread4 = ref_int_cpo_price - dom_blwn_cpo - taxeffectiverate*ref_int_cpo_price
label variable spread4 "spread_cpo_refpr_blwn"
gen double spread5 = spread_int_dom_paspi
label variable spread5 "spread_notax"
gen double spread6 = cif_rtdm_cpo - fob_blwn_cpo - taxeffectiverate*cif_rtdm_cpo
label variable spread6 "spread_cpo_rtdm_fobblwn"

tabstat spread1 spread2 spread3 spread4, statistics( mean ) by(year)
* all missing until 2007. 

*/

/*
forvalues s = 1/4{
gen double iv`s'_imp1 = prex_cpo_imp1*spread`s'
gen double iv`s'_imp2 = prex_cpo_imp2*spread`s'
}
*/

**** AVERAGE EXPORT SHARE OVER TIME ****
sort firm_id year
bys firm_id : gen double lag1_prex_cpo_imp1 = prex_cpo_imp1[_n-1] 
bys firm_id : gen double lag1_prex_cpo_imp2 = prex_cpo_imp2[_n-1] 

bys firm_id: egen double avg_prex_cpo_imp1 = mean(prex_cpo_imp1)
bys firm_id: egen double avg_prex_cpo_imp2 = mean(prex_cpo_imp2)



***** LOG VARIABLES *****
global to_log ffb_price_imp1 ffb_price_imp2 cpo_price_imp1 cpo_price_imp2 pko_price_imp1 pko_price_imp2 ///
out_val_cpo_imp1 out_val_cpo_imp2 out_val_pko_imp1 out_val_pko_imp2 revenue_total ///

*ref_int_cpo_price cif_rtdm_cpo dom_blwn_cpo fob_blwn_cpo spread_int_dom_paspi rho dom_blwn_pko cif_rtdm_pko 
*iv1_imp1 iv1_imp2 iv2_imp1 iv2_imp2 iv3_imp1 iv3_imp2 iv4_imp1 iv4_imp2

foreach var of varlist $to_log {
	gen `var'_ln = ln(`var')
}

sort firm_id year 



***** UML ESTABLISHMENT YEAR AND CONCENTRATION *****

merge m:1 trase_code year using "temp_data/processed_UML/UML_panel_valentin.dta", generate(merge_uml_panel) keepusing(est_year_imp concentration_10 concentration_30 concentration_50)
drop if merge_uml_panel == 2
drop merge_uml_panel
order est_year_imp, after(est_year)



***** SAVE *****
sort firm_id year 



* as for FFB, in_tot_ variable names are shortened. 
rename in_tot_* in_* 


save "temp_data/IBS_UML_panel_final.dta", replace 


export excel firm_id year uml_matched_sample geo_sample analysis_sample trase_code uml_id mill_name parent_co lat lon island_factor island_name district_name kec_name village_name ///
min_year est_year est_year_imp startYear max_year active industry_code ///
ffb_price_imp1 ffb_price_imp2 in_ton_ffb in_ton_ffb_imp1 in_ton_ffb_imp2 in_val_ffb in_val_ffb_imp1 in_val_ffb_imp2 ///
in_cpo_price_imp1 in_cpo_price_imp2 in_ton_cpo in_ton_cpo_imp1 in_ton_cpo_imp2 in_val_cpo in_val_cpo_imp1 in_val_cpo_imp2 ///
cpo_price_imp1 cpo_price_imp2 out_ton_cpo out_ton_cpo_imp1 out_ton_cpo_imp2 out_val_cpo out_val_cpo_imp1 out_val_cpo_imp2 prex_cpo prex_cpo_imp1 prex_cpo_imp2 out_cpo  ///
pko_price_imp1 pko_price_imp2 out_ton_pko out_ton_pko_imp1 out_ton_pko_imp2 out_val_pko out_val_pko_imp1 out_val_pko_imp2 prex_pko prex_pko_imp1 prex_pko_imp2 out_pko ///
out_ton_rpo out_ton_rpo_imp1 out_ton_rpo_imp2 out_val_rpo out_val_rpo_imp1 out_val_rpo_imp2 prex_rpo prex_rpo_imp1 prex_rpo_imp2 out_rpo ///
out_ton_rpko out_ton_rpko_imp1 out_ton_rpko_imp2 out_val_rpko out_val_rpko_imp1 out_val_rpko_imp2 prex_rpko prex_rpko_imp1 prex_rpko_imp2 out_rpko ///
using "temp_data/IBS_UML_panel_final_IOvar.xlsx", firstrow(variables) replace 




***** MERGE TFP VARIABLES FROM SEBI *****





