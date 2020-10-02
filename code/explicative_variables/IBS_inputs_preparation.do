/* This script filters IBS input dataset by commodity (inputed) and keeps those relating to crude palm oil processing.  */

use "input_data/IBS_IO/IBS_inputs.dta", clear

sort firm_id year
replace in_commodity_code_kki = "." if in_commodity_code_kki == ""
gen in_commodity_bah_low = lower(in_commodity_bah)

egen firm_id_year_str = concat(firm_id year)
destring firm_id_year_str, generate(firm_id_year)
drop firm_id_year_str


***CLEANING THE COMMODITY VARIABLES*** 

	**This manually inspects the cases with incoherence between the three different variables for the commodity. 
	gen flag_incoherence = ((in_commodity_code_kki == "011340101" & in_commodity_bah_low != "tandan buah segar")|(in_commodity_code_kki == "011340501" & in_commodity_bah_low != "kelapa sawit segar"))
	gen flag_incoherence2 = ((in_commodity_bah_low == "tandan buah segar" & in_commodity_code_kki != "011340101")|(in_commodity_bah_low == "kelapa sawit segar" & in_commodity_code_kki != "011340501"))
	
	*Let's add cpo as an input, just to track them.
	/*
	gen flag_incoherence3 = ((in_commodity_bah_low == "minyak kasar/mentah kelapa sawit" & in_commodity_code_kki != "151410102"))
	gen flag_incoherence4 = ((in_commodity_code_kki == "151410102" & in_commodity_bah_low != "minyak kasar/mentah kelapa sawit"))
	There is no such incoherence between the code and the name. 
	*/
	*total flag_incoherence
	*total flag_incoherence2
	replace in_commodity_code_kki="011340101" if flag_incoherence2==1 & in_commodity_bah_low == "tandan buah segar"
	replace in_commodity_code_kki="011340501" if flag_incoherence2==1 & in_commodity_bah_low == "kelapa sawit segar"
	replace in_commodity_code_kki="011340101" if in_commodity_bah_low == "tbs"

	drop flag_incoherence flag_incoherence2 

*** EXTRACTING MILLS***
*the third code is for cpo. 
keep if in_commodity_code_kki=="011340101"| in_commodity_code_kki=="011340501" | in_commodity_code_kki=="151410102" 

		*create industry_code variable (only available until 2011)*
		* many cases of missing on rev3 but not on rev4 and vice-versa. Therefore the followging procedure rather than assigning rev3 to 1998-2009 and rev4 to 2009-2015.
		* there is no case of both being non-missing. 
		destring in_industry_code_rev3, replace 
		destring in_industry_code_rev4 , replace 
		replace in_industry_code_rev3 = 0 if in_industry_code_rev3 >=. 
		replace in_industry_code_rev4 = 0 if in_industry_code_rev4 >=. 
		gen industry_code = in_industry_code_rev3 + in_industry_code_rev4 
		drop in_industry_code_rev3 in_industry_code_rev4 


***CLEANING THE MEASUREMENT UNIT***

	gen in_measurement_ton = lower(in_measurement_unit)
	gen double in_dom_qty_ton = in_dom_qty
	gen double in_imp_qty_ton = in_imp_qty
	gen double in_tot_qty_ton = in_tot_qty
	

	* "*" is 232 instances, replace by "." and infer kg if in the same basket of all values with KG? Or delete ?
	replace in_measurement_ton = " " if in_measurement_unit == "*"
	replace in_measurement_ton="ton" if in_measurement_ton==" " & (in_tot_qty_ton==0 | in_tot_qty_ton>=. )
	*(There is no case of missing or zero in_tot_qty while non missing or non zero in imp or dom.)
	
	*gen measure_missing = (in_measurement_ton==" ")
	*egen any_measure_missing = max(measure_missing), by(firm_id)
	
	*Case by case, missing values are either replaced when obvious (big change in output quantity, or unit always used in other years) or droped. 
	*in_measurement_unit_code is of no help. 
	replace in_measurement_ton = "kg" if firm_id_year == 20612014
	replace in_measurement_ton = "kg" if firm_id_year == 20672014
	replace in_dom_qty_ton = . if in_measurement_ton == " "
	replace in_imp_qty_ton = . if in_measurement_ton == " "
	replace in_tot_qty_ton = . if in_measurement_ton == " "
		
	* bal is 1 instance and I don't know what it can be --> delete. 
	replace in_dom_qty_ton = . if in_measurement_ton == "bal"
	replace in_imp_qty_ton = . if in_measurement_ton == "bal"
	replace in_tot_qty_ton = . if in_measurement_ton == "bal"
	replace in_measurement_ton = " " if in_measurement_ton == "bal"

	* buah is 10 instances and means fruit. There is no case where it can be infered that an actual weigh measure was meant --> delete
	replace in_dom_qty_ton = . if in_measurement_ton == "buah"
	replace in_imp_qty_ton = . if in_measurement_ton == "buah"
	replace in_tot_qty_ton = . if in_measurement_ton == "buah"
	replace in_measurement_ton = " " if in_measurement_ton == "buah"

	* feet is one instance and cannot be converted --> delete 
	replace in_dom_qty_ton = . if in_measurement_ton == "feet"
	replace in_imp_qty_ton = . if in_measurement_ton == "feet"
	replace in_tot_qty_ton = . if in_measurement_ton == "feet"
	replace in_measurement_ton = " " if in_measurement_ton == "feet"

	* iu is one instance and cannot be converted --> delete 
	replace in_dom_qty_ton = . if in_measurement_ton == "iu"
	replace in_imp_qty_ton = . if in_measurement_ton == "iu"
	replace in_tot_qty_ton = . if in_measurement_ton == "iu"
	replace in_measurement_ton = " " if in_measurement_ton == "iu"

	* kg is 4993 instances --> convert in TON
	replace in_dom_qty_ton = in_dom_qty_ton/1000 if in_measurement_ton=="kg"
	replace in_imp_qty_ton = in_imp_qty_ton/1000 if in_measurement_ton=="kg"
	replace in_tot_qty_ton = in_tot_qty_ton/1000 if in_measurement_ton=="kg"
	replace in_measurement_ton = "ton" if in_measurement_ton == "kg"
	
	* liter is 3 instances, we cannot convert. --> delete.
	replace in_dom_qty_ton = . if in_measurement_ton == "liter" & in_commodity_code_kki != "151410102"
	replace in_imp_qty_ton = . if in_measurement_ton == "liter" & in_commodity_code_kki != "151410102"
	replace in_tot_qty_ton = . if in_measurement_ton == "liter" & in_commodity_code_kki != "151410102"
	replace in_measurement_ton = " " if in_measurement_ton == "liter" & in_commodity_code_kki != "151410102"
	* for cpo it is convertible: 
	replace in_dom_qty_ton = in_dom_qty_ton*0.9026/1000 if in_measurement_ton == "liter" & in_commodity_code_kki == "151410102"
	replace in_imp_qty_ton = in_imp_qty_ton*0.9026/1000 if in_measurement_ton == "liter" & in_commodity_code_kki == "151410102"
	replace in_tot_qty_ton = in_tot_qty_ton*0.9026/1000 if in_measurement_ton == "liter" & in_commodity_code_kki == "151410102"
	replace in_measurement_ton = "ton" if in_measurement_ton == "liter" & in_commodity_code_kki == "151410102"

	* M3 is 14 instances. cannot be converted --> delete. 
	replace in_dom_qty_ton = . if in_measurement_ton == "m3" & in_commodity_code_kki != "151410102"
	replace in_imp_qty_ton = . if in_measurement_ton == "m3" & in_commodity_code_kki != "151410102"
	replace in_tot_qty_ton = . if in_measurement_ton == "m3" & in_commodity_code_kki != "151410102"
	replace in_measurement_ton = " " if in_measurement_ton == "m3" & in_commodity_code_kki != "151410102"
	* for cpo it is convertible: 
	replace in_dom_qty_ton = in_dom_qty_ton*0.9026 if in_measurement_ton == "m3" & in_commodity_code_kki == "151410102"
	replace in_imp_qty_ton = in_imp_qty_ton*0.9026 if in_measurement_ton == "m3" & in_commodity_code_kki == "151410102"
	replace in_tot_qty_ton = in_tot_qty_ton*0.9026 if in_measurement_ton == "m3" & in_commodity_code_kki == "151410102"
	replace in_measurement_ton = "ton" if in_measurement_ton == "m3" & in_commodity_code_kki == "151410102"

	* meter is 6 instances --> delete
	replace in_dom_qty_ton = . if in_measurement_ton == "meter"
	replace in_imp_qty_ton = . if in_measurement_ton == "meter"
	replace in_tot_qty_ton = . if in_measurement_ton == "meter"
	replace in_measurement_ton = " " if in_measurement_ton == "meter"

	* pasang is 1 instance (and the mill has only one obs. --> delete
	replace in_dom_qty_ton = . if in_measurement_ton == "pasang"
	replace in_imp_qty_ton = . if in_measurement_ton == "pasang"
	replace in_tot_qty_ton = . if in_measurement_ton == "pasang"
	replace in_measurement_ton = " " if in_measurement_ton == "pasang"

	* roll is 1 instance --> delete
	replace in_dom_qty_ton = . if in_measurement_ton == "roll"
	replace in_imp_qty_ton = . if in_measurement_ton == "roll"
	replace in_tot_qty_ton = . if in_measurement_ton == "roll"
	replace in_measurement_ton = " " if in_measurement_ton == "roll"

	* yard is 6 instances --> delete
	replace in_dom_qty_ton = . if in_measurement_ton == "yard"
	replace in_imp_qty_ton = . if in_measurement_ton == "yard"
	replace in_tot_qty_ton = . if in_measurement_ton == "yard"
	replace in_measurement_ton = " " if in_measurement_ton == "yard"

	replace in_dom_qty_ton = . if in_measurement_ton == "galon"
	replace in_imp_qty_ton = . if in_measurement_ton == "galon"
	replace in_tot_qty_ton = . if in_measurement_ton == "galon"
	replace in_measurement_ton = " " if in_measurement_ton == "galon"

	replace in_dom_qty_ton = . if in_measurement_ton == "batang"
	replace in_imp_qty_ton = . if in_measurement_ton == "batang"
	replace in_tot_qty_ton = . if in_measurement_ton == "batang"
	replace in_measurement_ton = " " if in_measurement_ton == "batang"

	replace in_dom_qty_ton = . if in_measurement_ton == "lembar"
	replace in_imp_qty_ton = . if in_measurement_ton == "lembar"
	replace in_tot_qty_ton = . if in_measurement_ton == "lembar"
	replace in_measurement_ton = " " if in_measurement_ton == "lembar"

	replace in_dom_qty_ton = . if in_measurement_ton == "m2"
	replace in_imp_qty_ton = . if in_measurement_ton == "m2"
	replace in_tot_qty_ton = . if in_measurement_ton == "m2"
	replace in_measurement_ton = " " if in_measurement_ton == "m2"

	
*We wont use in_tot_qty any more so let's drop them 
drop in_imp_qty in_dom_qty in_tot_qty 

***LONG TO WIDE***

/*Reshape is either after agregating multi-inputs per commo_code per obs. or after removing those.*/

	**Agregate multiinput**
	destring in_commodity_code_kki, generate(commo_code)

		*ACTUAL DUPLICATES*
		/*There is no such duplicates within year. There are some duplicates across years, but, well, let's us keep them for now, they can still provide some information for further tests of output. 
		This would do the job however: 
		duplicates tag firm_id commo_code in_dom_qty_ton in_imp_qty_ton in_tot_qty_ton in_dom_val in_imp_val in_tot_val, gen(du_firm_commo_in)
		browse if du_firm_commo_in>0
		drop du_firm_commo_in
		duplicates drop firm_id commo_code in_dom_qty_ton in_imp_qty_ton in_tot_qty_ton in_dom_val in_imp_val in_tot_val, force
		*/

		*MULTIINPUT DUPLICATES*
		*keep track of what was multiinput per commo_code
		duplicates tag firm_id year commo_code, generate(tag_multiinput)
		gen flag_multiinput =(tag_multiinput>0)
		drop tag_multiinput

		*browse in_tot_qty_ton if flag_multiinput==1 & commo_code == 151410102
		*browse in_tot_qty_ton if ((in_imp_qty_ton>0&in_imp_qty_ton<.) | (in_imp_val>0&in_imp_val<.)) & commo_code == 151410102

		*Keeping only total quantity and imported on for keeping track of this. 205 importers of cpo (some are multiinput)

		
		*Aggregate and drop now useless information
		by firm_id year commo_code, sort: egen double ag_in_dom_qty_ton = total(in_dom_qty_ton) 
		by firm_id year commo_code, sort: egen double ag_in_dom_val = total(in_dom_val) 

		by firm_id year commo_code, sort: egen double ag_in_imp_qty_ton = total(in_imp_qty_ton) 
		by firm_id year commo_code, sort: egen double ag_in_imp_val = total(in_imp_val) 

		by firm_id year commo_code, sort: egen double ag_in_tot_qty_ton = total(in_tot_qty_ton) 
		by firm_id year commo_code, sort: egen double ag_in_tot_val = total(in_tot_val) 

		drop in_dom_qty_ton in_dom_val in_imp_qty_ton in_imp_val in_tot_qty_ton in_tot_val

		*drop these duplicates - and not replace them by missing of course otherwise reshape is not possible. 
		duplicates drop firm_id year commo_code, force
		

********************************************************************************************
	**Remove multiple lines per commodity. 
	*destring in_commodity_code_kki, generate(commo_code)
	*duplicates tag firm_id year commo_code, generate(tag_multiinput)
	*gen flag_multiinput_per_commo_code =(tag_multiinput>0)
	*total flag_multiinput_per_commo_code

	*egen any_multiinput = max(flag_multiinput_per_commo_code), by(firm_id)
	*total any_multiinput
	*drop if flag_multiinput_per_commo_code>0
*********************************************************************************************

	**Reshape and clean a bit**
	drop in_commodity_bah in_commodity_code_kki in_commodity_en  
	drop in_commodity_bah_misc in_industry_code_12 ISIC5D in_measurement_unit_code in_commodity_bah_14 in_commodity_code_kbki in_commodity_desc_kbki in_commodity_bah_low 
	drop in_measurement_ton imp_country_1_code imp_country_1_name imp_country_2_code imp_country_2_name imp_country_3_code imp_country_3_name
	reshape wide in_measurement_unit ag_in_dom_qty_ton ag_in_dom_val ag_in_imp_qty_ton ag_in_imp_val ag_in_tot_val ag_in_tot_qty_ton flag_multiinput, i(firm_id year) j(commo_code)


	rename *11340101 *_ffb
	rename *11340501 *_kss
	rename *151410102 *_cpo

*we will only keep the information on distinction between dom/imp/tot for cpo, not for ffb and kss
drop ag_in_dom_qty_ton_ffb ag_in_dom_val_ffb ag_in_imp_qty_ton_ffb ag_in_imp_val_ffb ag_in_dom_qty_ton_kss ag_in_dom_val_kss ag_in_imp_qty_ton_kss ag_in_imp_val_kss

	**Use only one input commodity***
	/*Comme on ne se sert du input en l'occurrence que pour avoir un ordre de grandeur pour évaluer les ouputs, on va utiliser l'info des kelapa sawit segar en les convertissant: 
	https://www.researchgate.net/publication/319298126_Characteristics_of_Fresh_Fruit_Bunch_Yield_and_the_Physicochemical_Qualities_of_Palm_Oil_during_Storage_in_North_Sumatra_Indonesia
	According to Wikipedia, Tenera is the main species used. So use the fruit/bunch ratio of 66.5% 
	First consider that it is credible to have both commo as inputs, then let's add up*/
	gen flag_2_inputs =((ag_in_tot_qty_ton_ffb!=0&ag_in_tot_qty_ton_ffb<.) & (ag_in_tot_qty_ton_kss!=0&ag_in_tot_qty_ton_kss<.)) 
	replace ag_in_tot_qty_ton_ffb = ag_in_tot_qty_ton_ffb + ag_in_tot_qty_ton_kss*(1/0.665) if flag_2_inputs==1
	/*And for values*/
	gen flag_2_inputsval =((ag_in_tot_val_ffb!=0&ag_in_tot_val_ffb<.) & (ag_in_tot_val_kss!=0&ag_in_tot_val_kss<.)) 
	replace ag_in_tot_val_ffb = ag_in_tot_val_ffb + ag_in_tot_val_kss if flag_2_inputsval==1
	
	/*Then use the information on segar when tandan is missing*/
	gen flag_fruittobunch = ((ag_in_tot_qty_ton_ffb==0|ag_in_tot_qty_ton_ffb>=.) & (ag_in_tot_qty_ton_kss!=0&ag_in_tot_qty_ton_kss<.))
	replace ag_in_tot_qty_ton_ffb = ag_in_tot_qty_ton_kss*(1/0.665) if flag_fruittobunch==1
	/*And for values (there are different cases than for qty NO WE SHOULD NOT CONVERT FOR VALUES; VALUES ARE STILL VALUES. 
	gen flag_fruittobunch_val = ((ag_in_tot_val_ffb==0|ag_in_tot_val_ffb>=.) & (ag_in_tot_val_kss!=0&ag_in_tot_val_kss<.))
	replace ag_in_tot_val_ffb = ag_in_tot_val_kss if flag_fruittobunch_val==1
*/
	/*So we add 361 obs. by doing so */

	replace flag_multiinput_ffb = 1 if (flag_multiinput_ffb==1|flag_multiinput_kss==1)&flag_2_inputs==1
	replace flag_multiinput_ffb = flag_multiinput_kss if flag_fruittobunch==1
	drop in_measurement_unit_kss flag_multiinput_kss ag_in_tot_qty_ton_kss ag_in_tot_val_kss


rename ag_in_tot_qty_ton_ffb in_ton_ffb 
rename ag_in_tot_val_ffb in_val_ffb 
rename ag_in* in*
order flag_multiinput_ffb flag_multiinput_cpo, before(flag_2_inputs)
drop in_measurement_unit_ffb
drop in_measurement_unit_cpo

drop firm_id_year

capture mkdir "temp_data/processed_IBS"
capture mkdir "temp_data/processed_IBS/prepared_IO"

save "temp_data/processed_IBS/prepared_IO/IBS_inputs_prep.dta", replace 
