/* This script filters IBS output dataset by commodity (outputted) and keeps those relating to crude palm oil processing.  */

use "input_data/IBS_IO/IBS_outputs.dta", clear

sort firm_id year

destring out_commodity_code_kki , generate(commo_code)
replace commo_code = commo_code-151400000
egen firm_id_year_str = concat(firm_id year)
destring firm_id_year_str, generate(firm_id_year)
drop firm_id_year_str


***CLEANING THE COMMODITY VARIABLES*** 

	**This manually inspects the cases with incoherence between the three different variables for the commodity. 
	gen flag_incoherence = ((commo_code == 10103 & out_commodity_en != "Crude palm kernel oin")|(commo_code == 10102 & out_commodity_en != "Crude palm oil"))
	gen flag_incoherence2 = ((out_commodity_en == "Crude palm kernel oin" & commo_code != 10103 )|(out_commodity_en == "Crude palm oil" & commo_code != 10102))
	gen flag_incoherence3 = ((commo_code == 10103 & out_commodity_bah != "Minyak kasar/mentah biji kelapa sawit")|(commo_code == 10102 & out_commodity_bah != "Minyak kasar/mentah kelapa sawit"))
	gen flag_incoherence4 = ((out_commodity_bah == "Minyak kasar/mentah biji kelapa sawit" & commo_code != 10103 )|(out_commodity_bah == "Minyak kasar/mentah kelapa sawit" & commo_code != 10102))
	gen flag_incoherence5 = ((commo_code == 40101 & out_commodity_en != "Palm cooking oil")|(commo_code == 40102 & out_commodity_en != "Cookin oil from palm kernel"))
	gen flag_incoherence6 = ((out_commodity_en == "Palm cooking oil" & commo_code != 40101 )|(out_commodity_en == "Cookin oil from palm kernel" & commo_code != 40102))
	gen flag_incoherence7 = ((commo_code == 40101 & out_commodity_bah != "Minyak goreng kelapa sawit")|(commo_code == 40102 & out_commodity_bah != "Minyak goreng inti kelapa sawit"))
	gen flag_incoherence8 = ((out_commodity_bah == "Minyak goreng kelapa sawit" & commo_code != 40101 )|(out_commodity_bah == "Minyak goreng inti kelapa sawit" & commo_code != 40102))

	replace commo_code = . if commo_code ==10102 & out_commodity_bah =="COPEX"
	replace commo_code = . if firm_id==44984 & year==2003 & commo_code == 40102
	replace commo_code = . if flag_incoherence5==1
	replace commo_code = 40101 if flag_incoherence6==1
	replace commo_code = 40101 if flag_incoherence8==1

	drop flag_incoherence flag_incoherence2 flag_incoherence3 flag_incoherence4 flag_incoherence5 flag_incoherence6 flag_incoherence7 flag_incoherence8

	**Now, we add the obs. that have neither of the three variables filled with a standard value but are still likely to be one of our commodities. 
	g low_out_commodity_bah = lower(out_commodity_bah)
	replace commo_code = 10102 if (low_out_commodity_bah=="minyak sawit (cpo)" | low_out_commodity_bah=="cpo")
	replace commo_code = 10103 if low_out_commodity_bah=="c-pko"
	drop low_out_commodity_bah

	g low_out_commodity_en = lower(out_commodity_en)
	replace commo_code = 10102 if (low_out_commodity_en=="crude palm oil" | low_out_commodity_en=="cpo")
	drop low_out_commodity_en



***EXTRACTING MILLS***
keep if (commo_code==10102|commo_code== 10103|commo_code==40101|commo_code==40102)

*Just to make it clean
replace out_commodity_en = "Crude palm oil" if commo_code==10102  
replace out_commodity_en = "Crude palm kernel oin" if commo_code==10103 
replace out_commodity_en = "Palm cooking oil" if commo_code==40101  
replace out_commodity_en = "Cookin oil from palm kernel" if commo_code==40102 
	
*create industry_code variable*

/* This was the code running until 30th october 2019. 
destring out_industry_code_rev3, replace 
destring out_industry_code_rev4, replace 
gen rev3 = out_industry_code_rev3
replace rev3 = 0 if rev3 >=. 
gen rev4 = out_industry_code_rev4
replace rev4 = 0 if rev4 >=. 
gen accountrev3 = (year==2009)
gen industry_code = rev3*accountrev3 + rev4 
drop out_industry_code_rev3 out_industry_code_rev4 rev3 rev4 accountrev3
*/

destring out_industry_code_rev3, replace 
destring out_industry_code_rev4, replace 
gen industry_code = out_industry_code_rev3 
replace industry_code = out_industry_code_rev4 if year > 2009
drop out_industry_code_rev3 out_industry_code_rev4 

*browse firm_id year industry_code if mi(industry_code)


***CLEANING THE MEASUREMENT UNIT***

	gen out_measurement_ton = lower(out_measurement_unit)
	gen double out_qty_ton = out_output_qty
	

	* "*" is 159 instances, replace by "." and infer kg if in the same basket of all values with KG? Or delete ?
	replace out_measurement_ton = " " if out_measurement_unit == "*"
	replace out_measurement_ton="ton" if out_measurement_ton ==" "& (out_qty_ton==0 | out_qty_ton>=.) 
	*gen measure_missing = (out_measurement_ton==" ")
	*egen any_measure_missing = max(measure_missing), by(firm_id)
	
	*Case by case, missing values are either replaced when obvious (big change in output quantity, or unit always used in other years) or removed. 
	*out_measurement_unit_code is of no help. 
	replace out_measurement_ton = "kg" if firm_id_year == 21172011
	replace out_measurement_ton = "kg" if firm_id_year==104512001
	replace out_measurement_ton = "ton" if firm_id_year==563132010
	replace out_measurement_ton = "kg" if firm_id_year==682452011
	replace out_measurement_ton = "ton" if firm_id_year==703792011
	replace out_qty_ton = . if out_measurement_ton==" "
		
	* BATANG is 1 instance and means "tige" --> delete. 
	replace out_qty_ton = . if out_measurement_ton=="batang"

	* BUAH is 7 instances and means fruit apparently --> delete
	replace out_qty_ton = . if out_measurement_ton=="buah"

	* kg is 11,919 instances --> convert in TON
	replace out_qty_ton = out_qty_ton/1000 if out_measurement_ton=="kg"
	replace out_measurement_ton = "ton" if out_measurement_ton == "kg"

	*kw is 4 instances: all from this big refinery which is 47919. 
	replace out_qty_ton = . if out_measurement_ton =="kw"
	
	* LITER is 471 instances. Can we convert in KG ? 
		* 1 litre of Palm Oil is not up to 1 kg. It's about 912.28 Grams. from http://www.webconversiononline.com/mobile/weightof.aspx?quantity=1&measure=liter&ingredient=palmoil
		* Weight of 1 Liter Palm kernel is nearly  921.43 Gram. from http://www.webconversiononline.com/weightof.aspx?quantity=1&measure=liter&ingredient=palmkernel
		* We use rather https://de.scribd.com/doc/296169844/Density-Table-for-Palm-Oil-Products for its featuring densities for all 4 commo. 
	replace out_qty_ton = out_qty_ton*0.9026/1000 if out_measurement_ton =="liter"&commo_code==10102
	replace out_qty_ton = out_qty_ton*0.9099/1000 if out_measurement_ton =="liter"&commo_code==10103
	replace out_qty_ton = out_qty_ton*0.9000/1000 if out_measurement_ton =="liter"&commo_code==40101
	replace out_qty_ton = out_qty_ton*0.9099/1000 if out_measurement_ton =="liter"&commo_code==40102
	replace out_measurement_ton = "ton" if out_measurement_ton =="liter"

	* M3 is 13 instances. It may be converted as m3 is likely meaning 1000L 
	replace out_qty_ton = out_qty_ton*0.9026 if out_measurement_ton =="m3"&commo_code==10102
	replace out_qty_ton = out_qty_ton*0.9099 if out_measurement_ton =="m3"&commo_code==10103
	replace out_qty_ton = out_qty_ton*0.9000 if out_measurement_ton =="m3"&commo_code==40101
	replace out_qty_ton = out_qty_ton*0.9099 if out_measurement_ton =="m3"&commo_code==40102
	replace out_measurement_ton = "ton" if out_measurement_ton =="m3"

	* METER is 5 instances --> delete
	replace out_qty_ton = . if out_measurement_ton=="meter"

	* POND is 1 instance, it is RPO and can be converted
	replace out_qty_ton = out_qty_ton*0.45359237/1000 if out_measurement_ton=="pond"
	replace out_measurement_ton = "ton" if out_measurement_ton == "pond"




***LONG TO WIDE***
	**Multiple lines per commodity. 
		*DUPLICATES*
		/*First remove actual duplicates within duplicates of commo_code (and hence within the same year; those across years will be treated further in the cleaning step).
		because these duplicates, once aggregated, wouldn't be visible anymore. 
		The question is, within firm_id year commo_code, do we accept duplicates of qty that differ in value? I would say yes, like : they roughly calculated that half of 
		their production of this commo was sold to buyer A at a certain price, and the other half (equal then) sold to B at a different price. 
		What about duplicates in terms of value that differ in qty? What is the story behind this? They are only three, let's keep them. 
		Should we even remove duplicates in only qty and value? may be a total sale was broken up and all parts are equal but real? 
		Keeping them in the aggregation risks to produce overestimated values, while removing them risks to produce underestimates. But this should not influence the price, 
		or variable of interest, as a ratio of these cells multiplied or divided by the same factor. 
		
		duplicates tag firm_id year commo_code out_qty_ton out_output_value, gen(tag_du_firm_commo_qty_val)
		gen flag_du_firm_commo_qty_val= (tag_du_firm_commo_qty_val>0)
		egen any_du_firm_commo_qty_val= max(flag_du_firm_commo_qty_val), by(firm_id year commo_code)
		replace out_qty_ton=. if any_du_firm_commo_qty_val==1

		duplicates tag firm_id year commo_code out_qty_ton , gen(tag_du_firm_commo_qty)
		gen flag_du_firm_commo_qty= (tag_du_firm_commo_qty>0)

		duplicates tag firm_id year commo_code out_output_value , gen(tag_du_firm_commo_val)
		gen flag_du_firm_commo_val= (tag_du_firm_commo_val>0)

		count if flag_du_firm_commo_qty ==1
		count if flag_du_firm_commo_val==1
		count if flag_du_firm_commo_qty_val ==1
		count if flag_du_firm_commo_qty_val != flag_du_firm_commo_val

		replace out_qty_ton = . if flag_du_firm_commo_qty_val==1
		*/

		*Flag duplicates of commo_code to keep track of them for later. 
		duplicates tag firm_id year commo_code, generate(tag_multioutput)
		
		*aggregate by commo_code 
		by firm_id year commo_code, sort: egen double ag_out_qty_ton = total(out_qty_ton) 
		by firm_id year commo_code, sort: egen double ag_out_qty = total(out_output_qty) 
		by firm_id year commo_code, sort: egen double ag_out_value = total(out_output_value) 

		*average out_prex across duplicates. 
		gen weight_prex = out_qty_ton/ag_out_qty_ton
		gen normalized_out_prex = weight_prex*out_prex if out_prex <.
		by firm_id year commo_code, sort: egen double avg_out_prex = total(normalized_out_prex) if out_prex<. 
		drop weight_prex normalized_out_prex
		/*this prevents  out_prex missings from entering as zeros in the avg. However, weights are computed with account for outputs that have a missing out_prex.
		Finally, missing or zero out_qty_ton should give zero weight (the last condition on out_prex being nonmissing lets the sum operate on normalized_out_prex that 
		are zero solely because of out_qty_ton being zero or missing)
		sometimes ag_out_qty_ton is 0 and then weight_prex is missing, and that is why you have more mssing weighted_out_prex than missing out_prex even though there is no case of 
		out_qty_ton missing & out_prex<. 
		If you restrict to out_prex<. then you account as zeros obs. that have valid out_prex but missing weight (either because out_qty_ton is zero or missing) and that would not be accounted for 
		if you were restricted on normalized_out_prex<. */

		*drop these duplicates - and not replace them by missing of course otherwise reshape is not possible. 
		duplicates drop firm_id year commo_code, force 
		/*this droped arbitrary values (depending on how we would sort the data, as this command keeps the first occurence) in the sense that there is no reason to keep one 
		duplicate of out_qty_ton, out_output_qty or out_output_value´rather than an other. What we want is the aggregation only, hence these variables bring no additional information to the ag_ vars.*/
		drop out_qty_ton out_output_qty out_output_value out_prex
		

	**Reshape and clean a bit**
	drop out_industry_code_15 out_measurement_unit_code out_commodty_general_new out_commodity_desc_kbki out_commodity_code_hs2012 out_commodity_code_kbki out_commodity_bah out_commodity_en 
	drop exp_country_all exp_country_1_code exp_country_1_name exp_country_2_code exp_country_2_name exp_country_3_code exp_country_3_name
	drop out_measurement_ton out_commodity_code_kki out_commodity_general
	reshape wide out_measurement_unit ag_out_qty ag_out_qty_ton ag_out_value avg_out_prex tag_multioutput, i(firm_id year) j(commo_code)

	 rename *10102 *_cpo
	 rename *10103 *_pko
	 rename *40101 *_rpo
	 rename *40102 *_rpko

	 order out_measurement_unit_cpo, after(avg_out_prex_cpo)
	 order out_measurement_unit_pko, after(avg_out_prex_pko)
	 order out_measurement_unit_rpo, after(avg_out_prex_rpo)
	 order out_measurement_unit_rpko, after(avg_out_prex_rpko)

	 order ag_out_qty_cpo, before(out_measurement_unit_cpo)
	 order ag_out_qty_pko, before(out_measurement_unit_pko)
	 order ag_out_qty_rpo, before(out_measurement_unit_rpo)
	 order ag_out_qty_rpko, before(out_measurement_unit_rpko)

	 order industry_code, after(year)

drop firm_id_year

capture mkdir "temp_data/processed_IBS"
capture mkdir "temp_data/processed_IBS/prepared_IO"

save "temp_data/processed_IBS/prepared_IO/IBS_outputs_prep.dta", replace 
