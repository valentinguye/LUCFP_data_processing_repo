/*
This script makes several subsampling for different purposes. 
*/
use "temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta", clear
******************************************************************************************************************************
***** Prepare manual matching sheet ******************************************************************************************************

/* Keep only those mills that we want to match manually. 
These are the 59 mills (out of 972 distinct ones before 2011) that had not been sent to spatial join in R because 
*all* their annual observations were originally mi(desa_id) or flagged for misreporting desa_id in a split. 
PLUS the *76* mills that could actually not be matched to a desa polygon with their desa_id.   
*/
sort firm_id year 

*this will be useful to filter only obs. with no invalid desa_id that are appearing before 2011 (and excluding those that appeared before 1998 and come back wierdly after 2010)
bys firm_id (year): egen min_year2 = min(year)
** those with no annual valid desa_id observation (59 mills)
g valid_desa_id = (!mi(desa_id))
bys firm_id (year) : egen any_valid_desa_id_2 = total(valid_desa_id) 

*codebook firm_id if any_valid_desa_id_2 == 0 & min_year2 < 2011 
/*
browse firm_id year desa_id valid_desa_id any_valid_desa_id if any_valid_desa_id == 0 
sum max_year if any_valid_desa_id == 0
browse firm_id year any_valid_desa_id
*/
** and those with at least one valid desa_id obs. but that were matched to no desa polygon. 
merge m:1 firm_id using "temp_data/processed_mill_geolocalization/pre2011_bad_desa_id.dta", generate(merge_baddesa) keepusing(firm_id)
*codebook year if merge_baddesa == 3

keep if (any_valid_desa_id_2 == 0 & min_year2 < 2011 ) | merge_baddesa == 3 
* there remains 135 mills with 671 obs. 
*codebook firm_id 

rename kec_name subdistrict_name 
keep firm_id year min_year max_year workers_total_imp3 district_name subdistrict_name village_name
order district_name, before(subdistrict_name)
gen mill_name = ""
gen double latitude = . 
gen double longitude = .
sort firm_id year
export excel using "temp_data/processed_mill_geolocalization/mills_to_georeference_pre2011.xls", firstrow(variables) replace
******************************************************************************************************************************
******************************************************************************************************************************



******************************************************************************************************************************
***** Prepare noto sheet ******************************************************************************************************
use "temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta", clear

merge m:m firm_id using "temp_data/processed_mill_geolocalization/noto.dta", generate(merge_noto)

*bys firm_id (year): egen min_year = min(year)
*drop if workers_total_imp3 >=. 
keep if merge_noto == 3 

rename kec_name subdistrict_name 
keep firm_id year district_name subdistrict_name village_name desa_id workers_total_imp3 min_year est_year parent_co mill_name lat lon grp ///
in_ton_ffb_imp1 in_ton_ffb_imp2 out_ton_cpo_imp1 out_ton_cpo_imp2 out_ton_pko_imp1 out_ton_pko_imp2 out_ton_rpo_imp1 out_ton_rpo_imp2 out_ton_rpko_imp1 out_ton_rpko_imp2 ///
pct_own_cent_gov_imp pct_own_loc_gov_imp pct_own_nat_priv_imp pct_own_for_imp

order in_ton_ffb_imp1 in_ton_ffb_imp2 out_ton_cpo_imp1 out_ton_cpo_imp2 out_ton_pko_imp1 out_ton_pko_imp2 out_ton_rpo_imp1 out_ton_rpo_imp2 out_ton_rpko_imp1 out_ton_rpko_imp2 /// 
pct_own_cent_gov_imp pct_own_loc_gov_imp pct_own_nat_priv_imp pct_own_for_imp, after(lon)
order district_name, before(subdistrict_name)
order grp, before(firm_id)
sort grp firm_id year 

export excel using "temp_data/processed_mill_geolocalization/noto.xls", firstrow(variables) replace
******************************************************************************************************************************
******************************************************************************************************************************










