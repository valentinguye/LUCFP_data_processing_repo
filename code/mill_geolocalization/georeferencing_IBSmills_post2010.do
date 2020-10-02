/* This was how mills_to_georeference was produced:

import excel "$base_path_wd\build\input\mill_geolocalization\province_district_code_names_2011_15.xlsx", sheet("Sheet1") firstrow clear
gen id = _n
reshape long bps_ name_ , i(id) j(year)
egen year_prov_distr = concat(year bps_)
* there are some duplicates, I don't know why, just remove them.
duplicates drop year_prov_distr, force
save "$base_path_wd\build\input\mill_geolocalization\province_district_code_names_2011_15.dta", replace 

use "$base_path_wd\build\output\IBS_1998_cleaned.dta", clear
sort firm_id year 
*browse firm_id year workers_total_imp3 out_ton_cpo out_ton_pko out_ton_rpo out_ton_rpko if firm_id == 5283
*72100 | firm_id == 72984 | firm_id == 71678

gen post_2010 = (min_year > 2010)
drop if post_2010 == 0 
drop if workers_total_imp3 >=. 
keep firm_id year province district min_year workers_total_imp3 

//bys firm_id: egen max_year = max(year)

tostring province, generate(province_str)
tostring district, generate(district_str)
replace district_str = "0" + district_str if district < 10 
egen prov_distr_str = concat(province_str district_str)
destring prov_distr_str, generate(prov_distr) force 
egen year_prov_distr = concat(year prov_distr)if prov_distr <.
drop province province_str district district_str prov_distr_str 

merge m:1 year_prov_distr using "$base_path_wd\build\input\mill_geolocalization\province_district_code_names_2011_15.dta"
drop if _merge == 2
drop _merge id bps_ prov_distr
gen mill_name = ""
gen desa_name = ""
gen double latitude = . 
gen double longitude = .
sort firm_id year 
order min_year, after(year)
export excel using "$base_path_wd\build\input\mill_geolocalization\mills_to_georeference.xls", firstrow(variables) replace

**************************
 Now we would rather make it like : 
*/
use "temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta", clear
sort firm_id year 
gen post_2010 = (min_year > 2010)
drop if post_2010 == 0 
drop if workers_total_imp3 >=. 
keep firm_id year min_year district_name workers_total_imp3 
export excel using "temp_data/processed_mill_geolocalization/mills_to_georeference_post2010.xls", firstrow(variables) replace

*There is no district_name change between the two versions. 

