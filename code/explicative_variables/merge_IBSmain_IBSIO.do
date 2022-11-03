*****MERGING IBS_OUTPUTS WITH IBS_INPUTS*****

use "temp_data/processed_IBS/prepared_IO/IBS_outputs_prep.dta", clear

merge 1:1 firm_id year using "temp_data/processed_IBS/prepared_IO/IBS_inputs_prep.dta", generate(merge_codes_io) update
*option update for industry_code as only overlapping variable, firm_id_year being removed. Replaces missing values of same-named variables in master with values from using


save "temp_data/processed_IBS/prepared_IO/IBSIO.dta", replace 
*IBSIO is the dataset of IBS establishments that have at least one *commodity* codes of interest. Filtering was not based on *industry* code


*****MERGING IBS WITH IBSIO***** 
/* 
new way to merge: take rothenberg variables (kabu_code, startYear and EKSPOR PRPREX), add variables from final panel (i.e. actual order from Sebastian to BPS),
and keep only palm related establishments (either selected based on commodity code - all observations in IBSIO.dta - or on industry code from either Rothenberg or IBS_final_panel). 
and add the variables needed in other projects, that are in IBS shipment raw data but were not prepared by Sebastian in IBS_final_panel.dta
*/


use "input_data/IBS_rothenberg/si_panel.dta", clear
*this is Rothenberg's data, that go only until 2012. 
rename PSID firm_id
sort firm_id year
keep if year > 1997 
gen industry_code = prod5_rev3
replace industry_code = prod5_rev4 if mi(prod5_rev3)

/*
rename EKSPOR export_dummy
rename PRPREX export_pct 
*/
* variables for RSPO economics project 
rename ZPIVCU bonus_prod
rename ZNIVCU bonus_oth
rename ZPPVCU pension_prod
rename ZNPVCU pension_oth
rename IRRVCU royalties
rename IOTVCU oth_expenses
rename YISVCU mfc_services

keep firm_id year kabu_code industry_code startYear EKSPOR PRPREX ///
		bonus_prod bonus_oth pension_prod pension_oth royalties oth_expenses mfc_services

merge 1:1 firm_id year using "input_data/IBS_kraus/IBS_final_panel.dta", generate(merge_fp) keepusing(export_dummy export_dummy_imp export_pct export_pct_imp /// 
	 pct_own_cent_gov_imp pct_own_loc_gov_imp pct_own_nat_priv_imp pct_own_for_imp  /// 
	 revenue_total revenue_total_imp1 revenue_total_imp2 revenue_total_imp3 ///
	 value_added_self value_added_self_imp1 value_added_self_imp2 ///
	 inv_tot inv_tot_imp fc_add fc_add_imp ///
	 materials_tot materials_tot_imp1 materials_tot_imp2 materials_tot_imp3 ///
	 elec_qty elec_qty_imp1 elec_qty_imp2 elec_qty_imp3  ///
	 workers_total workers_total_imp1 workers_total_imp2 workers_total_imp3 ///
	 workers_prod workers_other workers_total_imp1 workers_prod_imp1 workers_other_imp1 workers_total_imp2 workers_prod_imp2 workers_other_imp2 ///
	 wage_prod wage_oth wage_prod_imp wage_oth_imp wage_prod_imp1 wage_prod_imp2 wage_oth_imp1 wage_oth_imp2 wage_prod_imp3 wage_oth_imp3 ///
	 kbli1 kbli2 ///
	 fc_est_tot_imp6co )

merge 1:1 firm_id year using "temp_data/processed_IBS/prepared_IO/IBSIO.dta", generate(merge_io)

gen palm_oil = (industry_code == 10431 | industry_code==15141 | industry_code==31151 | industry_code==10432 | industry_code==15144 | industry_code==31154)

keep if merge_io == 2 | merge_io == 3 | (merge_io == 1 & palm_oil == 1)

drop palm_oil 

* add the selection of variables available in IBS shipments but not available in IBS_final_panel.dta (Sebastian's panel)
* currently it's only ICOVCU: the gifts and donations
merge 1:1 firm_id year using "temp_data/processed_IBS/raw_IBS_shipment_panel_98_15.dta", generate(merge_raw)

keep if merge_raw == 1 | merge_raw == 3 

save "temp_data/processed_IBS/IBS_PO_98_15.dta", replace 
*save "$base_path_wd/build/input/IBS_1998.dta", replace
*saveold C:/Users/guyv/desktop/IBS_1998old, version(12) replace


/* 

previous way to merge IO and main, until 30th of October 2019
use "$base_path_wd/download/input/IBS_panel_pre_tfp.dta", clear

drop workers_total workers_prod workers_other export_dummy inv_tot fc_add elec_qty inv_tot_imp fc_add_imp kbli1 kbli2 elec_qty_imp2 workers_total_imp2 workers_prod_imp2 workers_other_imp2 workers_total_imp1 workers_prod_imp1 workers_other_imp1 elec_qty_imp1 elec_qty_imp3 materials_tot_imp3 fc_land_est_imp1 fc_est_tot_imp1co fc_est_tot_imp2co fc_est_tot_imp3co fc_est_tot_imp4co fc_est_tot_imp1cd fc_est_tot_imp2cd fc_est_tot_imp3cd fc_est_tot_imp4cd fc_est_tot_imp5co fc_est_tot_imp5cd fc_est_tot_imp6co fc_est_tot_imp6cd fc_est_tot_imp7co fc_est_tot_imp7cd fc_est_tot_imp8co fc_est_tot_imp8cd fc_est_tot_imp9co fc_est_tot_imp10co fc_est_tot_imp11co elec_qty_ln workers_total_ln workers_prod_ln workers_other_ln elec_qty_imp1_ln workers_total_imp1_ln workers_prod_imp1_ln workers_other_imp1_ln elec_qty_imp2_ln workers_total_imp2_ln workers_prod_imp2_ln workers_other_imp2_ln workers_total_imp3_ln workers_prod_imp3_ln workers_other_imp3_ln elec_qty_imp3_ln fc_est_tot_imp1co_ln fc_est_tot_imp1cd_ln fc_est_tot_imp2co_ln fc_est_tot_imp2cd_ln fc_est_tot_imp3co_ln fc_est_tot_imp3cd_ln fc_est_tot_imp4co_ln fc_est_tot_imp4cd_ln fc_est_tot_imp5co_ln fc_est_tot_imp5cd_ln fc_est_tot_imp6co_ln fc_est_tot_imp6cd_ln fc_est_tot_imp7co_ln fc_est_tot_imp7cd_ln fc_est_tot_imp8co_ln fc_est_tot_imp8cd_ln fc_est_tot_imp9co_ln fc_est_tot_imp10co_ln fc_est_tot_imp11co_ln fc_add_ln inv_tot_ln fc_add_imp_ln inv_tot_imp_ln fc_land_est_imp1_ln
gen palm_oil = 1 if (industry_code == 10431 | industry_code==15141 | industry_code==31151 | industry_code==10432 | industry_code==15144 | industry_code==31154) & year >1997
*We add refining industry_code, just to see what we get. 
* not filtering to these industry codes before merging, and filtering after like "keep if palm_oil == 1 or if commo_code was of interest (i.e. all those from using) 
* does not make a difference in the number of obs. in the end in this case."
keep if palm_oil ==1
drop palm_oil 
sort firm_id year

merge 1:1 firm_id year using $base_path_wd/download/output/IBSIO.dta

save C:/Users/guyv/ownCloud/opalval/download/output/IBS_1998.dta, replace 
save C:/Users/guyv/ownCloud/opalval/build/input/IBS_1998.dta, replace
saveold C:/Users/guyv/desktop/IBS_1998old, version(12) replace
*/

