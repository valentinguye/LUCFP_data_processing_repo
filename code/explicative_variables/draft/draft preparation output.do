*rev3 is until 2010 (included) and rev4 is from 2010 (included):*
inspect industry_code_rev4 if year <= 2009
inspect industry_code_rev3 if year >= 2010
*commo_code remains with the former code base (15141) in all years (manually checked)*

*for unique values : codebook varname* 
*outputs09-15 has 865 unique values and IBS_panel_pre-tfp has 856 for the period 2009-2015. 
*Hence 9 mills have been already removed. 
*Also, unique obs. (some years) may have been removed within mills.


***CLEANING THE COMMIDITY VARIABLES***

	**What commodities?**
		*The commodities we want to keep are CPO, PKO, and cookin oil from palm kernel and palm cooking oil. 

	**what should be the commodity variable?**
		*we don't have exaclty the same number of obs. if we extract with the string or the code variable: 
		browse firm_id year commo_code out_commodity_en out_measurement_unit out_output_qty out_output_value out_prex 
		gen flag_commo = (out_commodity_en == "Crude palm kernel oin" | out_commodity_en== "Crude palm oil")
		gen flag_commo2 =(commo_code == 151410102 | commo_code == 151410103)
		gen flag_commobah = (out_commodity_bah == "Minyak kasar/mentah biji kelapa sawit" | out_commodity_bah == "Minyak kasar/mentah kelapa sawit")
		total flag_commo
		total flag_commo2
		total flag_commobah
		*There are more obs. with commo_code

		*More precisely: 
		gen flag_incoherence = ((commo_code == 151410103 & out_commodity_en != "Crude palm kernel oin")|(commo_code == 151410102 & out_commodity_en != "Crude palm oil"))
		gen flag_incoherence2 = ((out_commodity_en == "Crude palm kernel oin" & commo_code != 151410103 )|(out_commodity_en == "Crude palm oil" & commo_code != 151410102))
		total flag_incoherence
		total flag_incoherence2
		browse out_commodity_en out_commodity_bah commo_code if flag_incoherence2==1
		*There are 8 cases of correct code but uncorrect names, of which one is surely not CPO/PKO, it is "COPEX". 
			*replace commo_code = . if commo_code ==151410102 & out_commodity_bah =="COPEX"
		*There is one case of correct English name but incoherent code, this is firm 44984 in 2003.  
		browse out_commodity_en out_commodity_bah commo_code if firm_id==44984
		*We can't infer much from its production over time because it is a changing pattern. Moreover its qty is a duplicate three times.
		*this year the questionaire was such that either the name or the code could be informed by the respondent. 
		*out_commodity_code_kbki is not usefull here (always missing)
		*Let's put it aside. 
			*replace commo_code = . if firm_id==44984 & year==2003 & commo_code == 151440102

		*Checking with bahasa names there are two more incoherent cases (i.e. correct English name but uncorrect bahasa names)  
		gen flag_incoherence3 = ((commo_code == 151410103 & out_commodity_bah != "Minyak kasar/mentah biji kelapa sawit")|(commo_code == 151410102 & out_commodity_bah != "Minyak kasar/mentah kelapa sawit"))
		gen flag_incoherence4 = ((out_commodity_bah == "Minyak kasar/mentah biji kelapa sawit" & commo_code != 151410103 )|(out_commodity_bah == "Minyak kasar/mentah kelapa sawit" & commo_code != 151410102))
		total flag_incoherence3
		total flag_incoherence4
		browse out_commodity_en out_commodity_bah commo_code if flag_incoherence4 ==1
		*The two additional cases are alright. 
		*The incoherence4 is the same case as incoherence2: firm 44984 in 2003. 
		 

	**Check for the other commodities**
		*Ce serait bien d'extraire aussi les productions de cookin oil from palm kernel et Palm cooking oil, resp. 151440102 et 151440101
		gen flag_incoherence5 = ((commo_code == 151440101 & out_commodity_en != "Palm cooking oil")|(commo_code == 151440102 & out_commodity_en != "Cookin oil from palm kernel"))
		gen flag_incoherence6 = ((out_commodity_en == "Palm cooking oil" & commo_code != 151440101 )|(out_commodity_en == "Cookin oil from palm kernel" & commo_code != 151440102))
		total flag_incoherence5
		total flag_incoherence6
		browse year out_commodity_en out_commodity_bah commo_code if flag_incoherence6 ==1 
		*So for incoherence5 it is the 44984_2003 case. 
		*The second instance is palm stearin accounted as palm cooking oil, we don't account for this as it is a slightly different commodity. 
			*replace commo_code = . if flag_incoherence5==1
		*For incoherence6, there are 43 cases of palm cooking oil industry with too short codes, but starting all with 1514401 or 15144 so consider they are 151440101
			*replace commo_code = 151440101 if flag_incoherence6==1
		*out_commodity_code_kbki is not usefull here (always missing)	

		*checking with bahasa names changes nothing. 
		gen flag_incoherence7 = ((commo_code == 151440101 & out_commodity_bah != "Minyak goreng kelapa sawit")|(commo_code == 151440102 & out_commodity_bah != "Minyak goreng inti kelapa sawit"))
		gen flag_incoherence8 = ((out_commodity_bah == "Minyak goreng kelapa sawit" & commo_code != 151440101 )|(out_commodity_bah == "Minyak goreng inti kelapa sawit" & commo_code != 151440102))
		total flag_incoherence7
		total flag_incoherence8
		browse out_commodity_en out_commodity_bah commo_code if flag_incoherence8 ==1 
		*incoherence7 is not different than incoherence5 (i.e. with English names)
		*incoherence8 features one additional case of correct name but no commo_code at all. --> add code. 
			*replace commo_code=151440101 if flag_incoherence8==1


	**Now, we add the obs. that have neither of the three variables filled with a standard value but are still likely to be one of our commodities. 

g var1_lower = lower(out_commodity_bah)
g var2_lower = lower(out_commodity_en)

gen addobs = (commo_code>=. & (strpos(var1_lower, "pko") | strpos(var1_lower,"sawit") | strpos(var1_lower,"cpo") | strpos(var2_lower,"cpo") | strpos(var1_lower,"inti") |strpos(var1_lower, "mentah") | strpos(var2_lower, "palm")|strpos(var1_lower, "minyak")|strpos(var1_lower, "kelapa")|strpos(var2_lower, "kernel")))

*Be more precise as there is no mentah or kelapa interesting
gen addobs2 = (commo_code>=. &(strpos(var1_lower, "pko") | strpos(var1_lower,"sawit") | strpos(var1_lower,"cpo") | strpos(var2_lower,"cpo") | strpos(var2_lower, "palm")|strpos(var2_lower, "kernel")))
browse firm_id year out_measurement_unit out_output_qty out_output_value out_commodity_bah commo_code out_commodity_en if addobs2==1
*Only "minyak sawit (CPO)", "cpo" are trustfully 1514101032. Others can be refined or not, we don't know. 
	*replace commo_code = 151410102 if (var1_lower=="minyak sawit (cpo)" | var1_lower=="cpo")
*Plus we have one "C-PKO" case that we can reasonably assume is crude PKO. 
	*replace commo_code = 151410103 if var1_lower=="c-pko"
*out_commodity_code_kbki is not usefull here (always missing)































































	**Extracting With industry_code**
***********************************************************************************************************************************************************************
**Producing variable for industry codes**
destring out_industry_code_rev3, replace 
destring out_industry_code_rev4 , replace 

*merge rev3 and rev4*
gen rev3 = out_industry_code_rev3
replace rev3 = 0 if rev3 >=. 
gen rev4 = out_industry_code_rev4
replace rev4 = 0 if rev4 >=. 
gen accountrev3 = (year==2009)
gen industry_code = rev3*accountrev3 + rev4 
drop out_industry_code_rev3 out_industry_code_rev4 rev3 rev4 accountrev3 out_commodity_code_kki

gen CPO_mill = 1 if (industry_code == 10431 | industry_code==15141 | industry_code==31151) & year >1989
keep if CPO_mill ==1
*and then cleaning the commodity variable*
gen flag_incoherence3 = ((commo_code == 151440101 & out_commodity_en != "Palm cooking oil")|(commo_code == 151440102 & out_commodity_en != "Cookin oil from palm kernel"))
replace commo_code = . if flag_incoherence3==1
*add some more obs.* 
replace commo_code = 151410103 if (out_commodity_bah=="PKO")
replace commo_code = 151410103 if (out_commodity_bah=="INTI SAWIT")
*all other cases are not clearly enough CPO or PKO, or have missing out_output_qty value.*

*And finally keeping only commo of interest
gen leftbehind = (commo_code!=151410102&commo_code!= 151410103&commo_code!=151440101&commo_code!=151440102)
by firm_id_year, sort : egen float other_output = total(leftbehind)
keep if leftbehind==0
********************************************************************************************************************************************************************************
