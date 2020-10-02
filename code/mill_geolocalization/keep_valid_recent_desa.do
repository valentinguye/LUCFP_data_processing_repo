use "temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta", clear

*keep valid and most recent mills
g valid_desa_id = (!mi(desa_id))
codebook firm_id if valid_desa_id == 1

by firm_id: egen most_recent_vld = max(year) if valid_desa_id == 1 

g polygon_4match = (valid_desa_id == 1 & year == most_recent_vld)

keep if polygon_4match == 1 

destring desa_id, replace

keep firm_id year min_year max_year industry_code desa_id prov_name district_name kec_name village_name workers_total_imp3 /// 
in_ton_ffb_imp1 in_ton_ffb_imp2 out_ton_cpo_imp1 out_ton_cpo_imp2 out_ton_pko_imp1 out_ton_pko_imp2 out_ton_rpo_imp1 out_ton_rpo_imp2 out_ton_rpko_imp1 out_ton_rpko_imp2  ///
avg_out_ton_cpo_imp2 last_out_ton_cpo_imp2 avg_cpo_price_imp2 ///
avg_in_tot_ton_cpo_imp2 last_in_tot_ton_cpo_imp2 avg_in_ton_ffb_imp2 last_in_ton_ffb_imp2 avg_ffb_price_imp2

/*recast float firm_id year workers_total_imp3 /// 
in_ton_ffb_imp1 in_ton_ffb_imp2 out_ton_cpo_imp1 out_ton_cpo_imp2 out_ton_pko_imp1 out_ton_pko_imp2 out_ton_rpo_imp1 out_ton_rpo_imp2 out_ton_rpko_imp1 out_ton_rpko_imp2, force 
 */
sort firm_id year

capture mkdir "temp_data/processed_mill_geolocalization"

save "temp_data/processed_mill_geolocalization/IBSmills_valid_desa.dta", replace

