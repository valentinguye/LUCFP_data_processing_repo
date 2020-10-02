use "input_data/IBS_IO/IBS_outputs.dta", clear

sort firm_id year

destring out_commodity_code_kki , generate(commo_code)


* Select records of palm biomass residues produced and sold (strings have been identifier with strpos(out_commodity_bah and _en, "") with different entries like sawit, waste, shell, kernel, etc..
* Codes are their associated codes)
keep if out_commodity_bah == "Bungkil kelapa sawit" | ///
out_commodity_bah == "Ampas palma lainnya"  | /// 
out_commodity_bah == "Bungkil biji kelapa sawit"  | /// 
out_commodity_bah == "Hasil ikutan/sisa industri minyak goreng dari minyak kelapa sawit"  | /// 
out_commodity_bah == "Hasil ikutan/sisa lainnya minyak kasar dari nabati dan hewani"  | /// 
out_commodity_bah == "Hasil ikutan/sisa lainnya minyak kasar dari nabati dan hewan" | /// 
out_commodity_bah == "Margarine padat lainnya"  | /// 
out_commodity_bah == "Cangkang sawit" | /// 
out_commodity_bah == "Ampas buah kelapa sawit" | ///
out_commodity_code_kki == "011349801"| ///
out_commodity_code_kki == "151419802" | ///
out_commodity_code_kki == "151419803" | ///
out_commodity_code_kki == "151419890" | ///
out_commodity_code_kki == "151420102" | ///
out_commodity_code_kki == "151449800" | ///
out_commodity_code_kki == "153249804" 

* There are 537 such records, 
egen firm_id_year_str = concat(firm_id year)
destring firm_id_year_str, generate(firm_id_year)
drop firm_id_year_str
* distributed over 485 firm_id x year records, and 214 firm_id

merge m:1 firm_id year using "temp_data/IBS_UML_panel_final.dta", gen(merge_ibs_final)
codebook firm_id if merge_ibs_final == 3 & uml_matched_sample == 1
* 103 IBS-UML mills have at least one record of biomass residue sold. Summing to 252 records. 

codebook firm_id if merge_ibs_final == 3 & analysis_sample == 1
* 147 geolocalized IBS mills have at least one record of biomass residue sold. Summing to 392 records. 

*** Check electricity data 
* add electricity sold from Rothenberg (not ordered by Sebastian)
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

keep firm_id year OELKHU YELVCU

merge 1:1 firm_id year using "input_data/IBS_kraus/IBS_panel_pre_clean_year.dta", generate(merge_preclean) keepusing(generator_capacity electricity_PLN_qty electricity_PLN_value electricity_non_PLN_qty electricity_non_PLN_value electricity_self_generated expenses_building expenses_land expenses_indirect_taxes expenses_loans)

merge 1:1 firm_id year using "temp_data/IBS_UML_panel_final.dta", generate(merge_mills)

drop if merge_mills == 1

rename OELKHU electricity_sold_qty
rename YELVCU electricity_sold_value

* is_mill
codebook electricity_sold_qty electricity_sold_value generator_capacity electricity_PLN_qty electricity_PLN_value electricity_non_PLN_qty electricity_non_PLN_value electricity_self_generated if analysis_sample == 1

