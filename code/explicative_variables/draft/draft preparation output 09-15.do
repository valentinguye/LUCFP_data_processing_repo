******************************************************************************************Cleaning the output dataset (2009-2015)****************************************************
use "C:\Users\guyv\ownCloud\opalval\download\input\IBS_IO\IBS_outputs09-15.dta", clear
sort firm_id year
**Producing variable for industry codes**
destring out_industry_code_rev3, generate(industry_code_rev3)
destring out_industry_code_rev4 , generate(industry_code_rev4)
destring out_commodity_code_kki , generate(commo_code)
drop out_industry_code_rev3
drop out_industry_code_rev4
drop out_commodity_code_kki
*merge rev3 and rev4*
gen rev3 = industry_code_rev3
replace rev3 = 0 if rev3 >=. 
gen rev4 = industry_code_rev4
replace rev4 = 0 if rev4 >=. 
gen accountrev3 = (year==2009)
gen industry_code = rev3*accountrev3 + rev4 
drop industry_code_rev3 industry_code_rev4 rev3 rev4 accountrev3 

egen firm_id_year_str = concat(firm_id year)
destring firm_id_year_str, generate(firm_id_year)
** Producing variable for commodity codes** 

		***EXTRACTING MILLS***
*Either on industry code or on commodity code.*

	**With commodity code**
**************************************************************************************************
gen flag_incoherence3 = ((commo_code == 151440101 & out_commodity_en != "Palm cooking oil")|(commo_code == 151440102 & out_commodity_en != "Cookin oil from palm kernel"))
gen flag_incoherence4 = ((out_commodity_en == "Palm cooking oil" & commo_code != 151440101 )|(out_commodity_en == "Cookin oil from palm kernel" & commo_code != 151440102))
replace commo_code = . if flag_incoherence3==1
replace commo_code = . if commo_code ==151410102 & out_commodity_bah =="COPEX"
replace commo_code = 151440101 if flag_incoherence4==1
replace commo_code = 151410103 if (out_commodity_bah=="PKO")
replace commo_code = 151410103 if (out_commodity_bah=="INTI SAWIT")
keep if (commo_code==151410102|commo_code== 151410103|commo_code==151440101|commo_code==151440102)

replace out_commodity_en = "Crude palm oil" if commo_code==151410102  
replace out_commodity_en = "Crude palm kernel oin" if commo_code==151410103 
replace out_commodity_en = "Palm cooking oil" if commo_code==151440101  
replace out_commodity_en = "Cookin oil from palm kernel" if commo_code==151440102 
**************************************************************************************************

		***CLEANING THE MEASUREMENT UNIT***

	gen out_measurement_ton = out_measurement_unit
	gen out_qty_ton = out_output_qty

	codebook out_measurement_ton
	*9 unique values
	* "*" is 156 instances, replace by "." and infer kg if in the same basket of all values with KG? Or delete ?
	replace out_measurement_ton = "." if out_measurement_unit == "*"
	* BATANG is 1 instance and means tige --> delete. 
	drop if out_measurement_unit=="BATANG"
	* BUAH is 6 instances and means fruit apparently --> delete
	drop if out_measurement_unit=="BUAH"
	* KG is 5779 instances --> convert in TON
	replace out_qty_ton = out_qty_ton/1000 if out_measurement_unit=="KG"
	replace out_measurement_ton = "TON" if out_measurement_unit == "KG"
	* LITER is 469 instances. Can we convert in KG ? 
		* 1 litre of Palm Oil is not up to 1 kg. It's about 912.28 Grams. from http://www.webconversiononline.com/mobile/weightof.aspx?quantity=1&measure=liter&ingredient=palmoil
		* Weight of 1 Liter Palm kernel is nearly  921.43 Gram. from http://www.webconversiononline.com/weightof.aspx?quantity=1&measure=liter&ingredient=palmkernel
		* We use rather https://de.scribd.com/doc/296169844/Density-Table-for-Palm-Oil-Products for its featuring densities for all 4 commo. 
	replace out_qty_ton = out_qty_ton*0.9026/1000 if out_measurement_unit =="LITER"&commo_code==151410102
	replace out_qty_ton = out_qty_ton*0.9099/1000 if out_measurement_unit =="LITER"&commo_code==151410103
	replace out_qty_ton = out_qty_ton*0.9000/1000 if out_measurement_unit =="LITER"&commo_code==151440101
	replace out_qty_ton = out_qty_ton*0.9099/1000 if out_measurement_unit =="LITER"&commo_code==151440102
	replace out_measurement_ton = "TON" if out_measurement_unit =="LITER"
	* M3 is 13 instances. It may be converted as m3 is likely meaning 1000L 
	replace out_qty_ton = out_qty_ton*0.9026 if out_measurement_unit =="M3"&commo_code==151410102
	replace out_qty_ton = out_qty_ton*0.9099 if out_measurement_unit =="M3"&commo_code==151410103
	replace out_qty_ton = out_qty_ton*0.9000 if out_measurement_unit =="M3"&commo_code==151440101
	replace out_qty_ton = out_qty_ton*0.9099 if out_measurement_unit =="M3"&commo_code==151440102
	replace out_measurement_ton = "TON" if out_measurement_unit =="M3"
	* METER is 5 instances --> delete
	drop if out_measurement_unit=="METER"

	* POND is 1 instance, it is RPO and can be converted
	replace out_qty_ton = out_qty_ton*0.45359237/1000 if out_measurement_unit=="POND"
	replace out_measurement_ton = "TON" if out_measurement_unit == "POND"
	* TON is 928 instances. 


		***LONG TO WIDE***
	**Remove multiple lines per commodity. 
	duplicates tag firm_id year commo_code, generate(tag_multioutput_per_commo_code)
	gen flag_multioutput_per_commo_code =(tag_multioutput_per_commo_code>0)
	
	drop if flag_multioutput_per_commo_code>0

	**Reshape and clean a bit**
	drop  out_industry_code_15 out_industry_code_15 out_commodity_bah out_commodity_en   KBKI out_commodity_code_hs2012 out_commodity_code_kbki out_commodity_desc_kbki out_commodty_general_new firm_id_year_str flag_incoherence4 flag_incoherence3 tag_multioutput_per_commo_code flag_multioutput_per_commo_code
	rename out_measurement_unit_code out_measureunit_code 
	 *sinon ca fait un nom trop grand apparemment vu qu'il concatène les noms. 

	 reshape wide out_measurement_unit out_measurement_ton out_output_qty out_qty_ton out_output_value out_prex exp_country_all exp_country_1_code exp_country_1_name exp_country_2_code exp_country_2_name exp_country_3_code exp_country_3_name out_measureunit_code, i(firm_id year) j(commo_code)

	 rename *151410102 *_cpo
	 rename *151410103 *_pko
	 rename *151440101 *_rpo
	 rename *151440102 *_rpko
	 order out_measurement_ton_cpo, before(out_output_value_cpo)
	 order out_measurement_ton_pko, before(out_output_value_pko)
	 order out_measurement_ton_rpo, before(out_output_value_rpo)
	 order out_measurement_ton_rpko, before(out_output_value_rpko)

	 order out_qty_ton_cpo, before(out_output_value_cpo)
	 order out_qty_ton_pko, before(out_output_value_pko)
	 order out_qty_ton_rpo, before(out_output_value_rpo)
	 order out_qty_ton_rpko, before(out_output_value_rpko)

	 order industry_code, after(year)









	**With industry code** 
*******************************************************************************************************
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
******************************************************************************************************* 





*28.02.2019

**export countries**
*Pas très utile car très peu d'observations sur la grande base, et pas propre. 
*(les variables code ont des valeurs non numériques, donc le destring ne generate rien)

*On ferait bien de bien cleaner les quantités avant de merger, 
*comme ca on peut se servire de ce qu'on ne va pas garder. Non en fait on cleanera après le merge car il y a davantage d'info dans la grande base. 
*Et on ne perd pas trop d'info si on extrait les quatres codes. ? 

*NOTRE GROSSE QUESTION MAINTENANT C'EST COMMENT TRAITER LES OBS QUI ONT PLUSIEURS LIGNES POUR LA MEME COMMO ?* 
*checker si l'aggrégat est cohérent avec le input ou avec d'autres trcs ? 
*aller voir le questionnaire pour voir ce qui change --> rien 
*a-t-on le problème systématiquement pour les années autres que 2014 et 2015 comme pour 2055? 
*ou un pb de firm_id ? 

*comment extraire*
*comment merger*















*************************************************************************************************************************************************************************
DRAFT

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
*we don'thave exaclty the same number of obs. if we extract with the string or the code variable. 
*((7 more if we use the code rather than the name (on the raw data, 16 after the adding of PKO and INTI SAWIT)): 
browse firm_id year commo_code out_commodity_en out_measurement_unit out_output_qty out_output_value out_prex 
gen flag_commo = (out_commodity_en == "Crude palm kernel oin" | out_commodity_en== "Crude palm oil")
gen flag_commo2 =(commo_code == 151410102 | commo_code == 151410103)
gen flag_commobah = (out_commodity_bah == "Minyak kasar/mentah biji kelapa sawit" | out_commodity_bah == "Minyak kasar/mentah kelapa sawit")
total flag_commo
total flag_commo2
total flag_commobah
*So in 7 cases, the code is right but the name is not the standard one. To see these instances: 
gen flag_incoherence = ((commo_code == 151410103 & out_commodity_en != "Crude palm kernel oin")|(commo_code == 151410102 & out_commodity_en != "Crude palm oil"))
gen flag_incoherence2 = ((out_commodity_en == "Crude palm kernel oin" & commo_code != 151410103 )|(out_commodity_en == "Crude palm oil" & commo_code != 151410102))
total flag_incoherence
total flag_incoherence2
browse out_commodity_en out_commodity_bah commo_code if flag_incoherence ==1
*Checking with bahasa names changes nothing. 
drop flag_incoherence flag_incoherence2
gen flag_incoherence = ((commo_code == 151410103 & out_commodity_bah != "Minyak kasar/mentah biji kelapa sawit")|(commo_code == 151410102 & out_commodity_bah != "Minyak kasar/mentah kelapa sawit"))
gen flag_incoherence2 = ((out_commodity_bah == "Minyak kasar/mentah biji kelapa sawit" & commo_code != 151410103 )|(out_commodity_bah == "Minyak kasar/mentah kelapa sawit" & commo_code != 151410102))
total flag_incoherence
total flag_incoherence2
browse out_commodity_en out_commodity_bah commo_code if flag_incoherence ==1 
*En bahasi le nom est bien rempli et correspond au code. 
*So all obs. with one of the two commo value are alright in this respect: the code says what the name says. 
*DONC ON VA EXTRAIRE AVEC LE CODE PLUTOT*  


	**Check for the other commodities**
*Ce serait bien d'extraire aussi les productions de cookin oil from palm kernel et Palm cooking oil, resp. 151440102 et 151440101
gen flag_incoherence3 = ((commo_code == 151440101 & out_commodity_en != "Palm cooking oil")|(commo_code == 151440102 & out_commodity_en != "Cookin oil from palm kernel"))
gen flag_incoherence4 = ((out_commodity_en == "Palm cooking oil" & commo_code != 151440101 )|(out_commodity_en == "Cookin oil from palm kernel" & commo_code != 151440102))
total flag_incoherence3
total flag_incoherence4

*checking with bahasa names changes nothing. 
drop flag_incoherence3 flag_incoherence4
gen flag_incoherence3 = ((commo_code == 151440101 & out_commodity_bah != "Minyak goreng kelapa sawit")|(commo_code == 151440102 & out_commodity_bah != "Minyak goreng inti kelapa sawit"))
gen flag_incoherence4 = ((out_commodity_bah == "Minyak goreng kelapa sawit" & commo_code != 151440101 )|(out_commodity_bah == "Minyak goreng inti kelapa sawit" & commo_code != 151440102))
total flag_incoherence3
total flag_incoherence4
browse out_commodity_en out_commodity_bah commo_code if flag_incoherence3 ==1 | flag_incoherence4 ==1 

*Cette obs. incohérente est en crude palm stearin, commodidté qu'on ne veut pas garder car prix pas comparable, donc on supprime sa valeur en commo_code: 
*replace commo_code = . if flag_incoherence3==1

*dans le cas d'une extraction par les commodités, on veut en revanche extraire les 9 informations trouvées par flag_incoherence4: c'est de la palm cooking oil mais le code est mal entré
*replace commo_code = 151440101 if flag_incoherence4==1


	**Adding more obs.?**
*Now, are there other obs. we could add that don't have the commo_code? 
browse firm_id year out_measurement_unit out_output_qty out_output_value out_commodity_bah commo_code out_commodity_en if commo_code >=. & (out_commodity_bah=="PKO" | out_commodity_bah=="INTI SAWIT")
*only PKO and INTI SAWIT are consider clear enough to be converted into 151410103. 
* on ne peut pas garder les autres outputs car même si on a la quantité et la valeur. Cela ne correspondra pas à un seul prix, celui de la CPO. 
*replace commo_code = 151410103 if (out_commodity_bah=="PKO")
*replace commo_code = 151410103 if (out_commodity_bah=="INTI SAWIT")

*si on extrait avec le commo_code, on ne peut pas utiliser cette méthode de visualisation directe car le commo_code manque pour de nombreuses observations qui n'ont rien à voir avec l'huile de palme. 
*Donc on ne rajoute rien, seulement deux trois obs. avec PKO et INTI SAWIT se rajoutent avec la commande au-dessus. 
 browse firm_id year out_measurement_unit out_output_qty out_output_value out_commodity_bah commo_code out_commodity_en if commo_code >=. & (out_commodity_bah=="PKO" | out_commodity_bah=="INTI SAWIT")



 *What do we not take with us? THESE LINES WILL HAVE TO BE INCLUDED BEFORE THE KEEP IF.
 gen leftbehind = (commo_code!=151410102&commo_code!= 151410103&commo_code!=151440101&commo_code!=151440102)
 total leftbehind
 inspect leftbehind
*anyway we should flag obs that sell other stuff than these 4: 
by firm_id_year, sort : egen float other_output = total(leftbehind)
*we would also like to have the income these guys have from the sales of other outputs, and not just know that they have some additional incomes. 
*once this is done, let's get rid of these additional lines. - keep if leftbehind==0
*and how many of the remainings have other outputs? 
codebook firm_id_year if other_output>0
*-->351



		***CLEANING THE MEASUREMENT UNIT***

	gen out_measure_unit = out_measurement_unit
	gen out_qty_ton = out_output_qty

	codebook out_measure_unit
	*9 unique values
	* "*" is 156 instances, replace by "." and infer kg if in the same basket of all values with KG? Or delete ?
	replace out_measure_unit = "." if out_measurement_unit == "*"
	* BATANG is 1 instance and means tige --> delete. 
	drop if out_measurement_unit=="BATANG"
	* BUAH is 6 instances and means fruit apparently --> delete
	drop if out_measurement_unit=="BUAH"
	* KG is 5779 instances --> convert in TON
	replace out_qty_ton = out_qty_ton/1000 if out_measurement_unit=="KG"
	replace out_measure_unit = "TON" if out_measurement_unit == "KG"
	* LITER is 469 instances. Can we convert in KG ? 
		* 1 litre of Palm Oil is not up to 1 kg. It's about 912.28 Grams. from http://www.webconversiononline.com/mobile/weightof.aspx?quantity=1&measure=liter&ingredient=palmoil
		* Weight of 1 Liter Palm kernel is nearly  921.43 Gram. from http://www.webconversiononline.com/weightof.aspx?quantity=1&measure=liter&ingredient=palmkernel
		* We use rather https://de.scribd.com/doc/296169844/Density-Table-for-Palm-Oil-Products for its featuring densities for all 4 commo. 
	replace out_qty_ton = out_qty_ton*0.9026/1000 if out_measurement_unit =="LITER"&commo_code==151410102
	replace out_qty_ton = out_qty_ton*0.9099/1000 if out_measurement_unit =="LITER"&commo_code==151410103
	replace out_qty_ton = out_qty_ton*0.9000/1000 if out_measurement_unit =="LITER"&commo_code==151440101
	replace out_qty_ton = out_qty_ton*0.9099/1000 if out_measurement_unit =="LITER"&commo_code==151440102
	replace out_measure_unit = "TON" if out_measurement_unit =="LITER"
	* M3 is 13 instances. It may be converted as m3 is likely meaning 1000L 
	replace out_qty_ton = out_qty_ton*0.9026 if out_measurement_unit =="M3"&commo_code==151410102
	replace out_qty_ton = out_qty_ton*0.9099 if out_measurement_unit =="M3"&commo_code==151410103
	replace out_qty_ton = out_qty_ton*0.9000 if out_measurement_unit =="M3"&commo_code==151440101
	replace out_qty_ton = out_qty_ton*0.9099 if out_measurement_unit =="M3"&commo_code==151440102
	replace out_measure_unit = "TON" if out_measurement_unit =="M3"
	* METER is 5 instances --> delete
	drop if out_measurement_unit=="METER"

	* POND is 1 instance, it is RPO and can be converted
	replace out_qty_ton = out_qty_ton*0.45359237/1000 if out_measurement_unit=="POND"
	replace out_measure_unit = "TON" if out_measurement_unit == "POND"
	* TON is 928 instances. 

	codebook out_commodity_en out_output_qty if out_measurement_unit == "POND"

	*donc pour cleaner output on va avoir besoin de output_ton et output original, car les outliers peuvent venir d'un pb de conversion (i.e. la valeur était normale avant de passer de M3 à TON)
	*mais le pb peut venir d'ailleurs bien sûr, auquel cas il est bien de tout pouvoir comparer et donc d'avoir tout mis en tonne. 

	



	*
*before spreading we need to fix the pb of several lines/commo_code/obs. (obs=firm_id*year)



 		***NOTRE GROSSE QUESTION MAINTENANT C'EST COMMENT TRAITER LES OBS QUI ONT PLUSIEURS LIGNES POUR LA MEME COMMO ?** 
 	**are they many?
duplicates tag firm_id year commo_code, generate(tag_multioutput_per_commo_code)
gen flag_multioutput_per_commo_code =(tag_multioutput_per_commo_code>0)
*codebook firm_id_year if flag_multioutput_per_commo_code>0 & commo_code==151410102
*--> c'est donc 197 observations perdues si on supprime celles qui ont plusieurs output sous le même commo_code. 
*codebook firm_id_year if flag_multioutput_per_commo_code>0 
*--> 349 sur les 4 commo. 

	**checker si l'aggrégat est cohérent avec le input ou avec d'autres trucs ? --> pour ca il faut rajouter des variables agregated_qty et agregated_value. 
*Pour ca il faudrait d'abord cleaner output, car s'il y a des valeurs extrêmes comptées parmi plusieurs dans les sommes, elles vont les influencer définitivement, 
*qu'on les retire plus tard ou pas, les sommes les auront capturées. 
*Mais pour cleaner il faut avoir spreadé, et donc avoir une seul obs. par firm_id_year_commo. Donc on va procéder comme suit: 
*delete multioutput, spread, clean output, then compare distribution of cleaned output and outputs of deleted multioutput, clean multioutput, agregate. 

	** delete multioutput**
	drop if flag_multioutput_per_commo_code>0



	**SPREAD**
	 drop  out_industry_code_15 out_industry_code_15 out_commodity_bah out_commodity_en   KBKI out_commodity_code_hs2012 out_commodity_code_kbki out_commodity_desc_kbki out_commodty_general_new firm_id_year_str flag_incoherence4 flag_incoherence3 tag_multioutput_per_commo_code flag_multioutput_per_commo_code
	 rename out_measurement_unit_code out_measureunit_code 
	 *sinon ca fait un nom trop grand apparemment vu qu'il concatène les noms. 

	 reshape wide out_measurement_unit out_measure_unit out_output_qty out_qty_ton out_output_value out_prex exp_country_all exp_country_1_code exp_country_1_name exp_country_2_code exp_country_2_name exp_country_3_code exp_country_3_name out_measureunit_code, i(firm_id year) j(commo_code)

	 

	 **OUT_QTY CLEANING**
	 *We should first clean the measurement unit and convert the qty values corresponding. This can be done more easily before spreading.
	 *so, now that this is done, let's clean by commo. 
*CLEANING OUT_QTY_CPO*
graph hbox out_output_qty_cpo out_qty_ton_cpo
*--> biggest outliers go away with expressing all obs. in the same unit obviously. 
*See a bit closer 
graph hbox out_qty_ton_cpo if out_qty_ton_cpo <. 

*It seems there are still big outliers. Are they also those who produce refined palm oil ?
codebook firm_id_year if out_qty_ton_cpo > 15000000 & out_qty_ton_cpo<.
gen outlier_qty_ton_cpo =(out_qty_ton_cpo > 15000000 & out_qty_ton_cpo<.)
browse if outlier_qty_ton_cpo ==1
*--> NOPE 

sum out_qty_ton_cpo
total out_qty_ton_cpo if year == 2015
*--> total production of CPO in 2015 is 832 Mt, whereas supposed to be around 62Mt in 2015 and 55Mt in 2013 (cifor report 2017). C'est un indice de présence de gros outliers. 

total out_qty_ton_cpo if year == 2015 & outlier_qty_ton_cpo==0
*--> on approche de l'ordre de grandeur, mais on est encore au-dessus alors qu'on n'est pas supposé avoir toute la population des producteurs. 
*Donc il doit encore y avoir de gros outliers. 
total out_qty_ton_cpo if year == 2015 & out_qty_ton_cpo<200000
codebook firm_id if year == 2015 & out_qty_ton_cpo<200000 & out_qty_ton_cpo<. & out_qty_ton_cpo>0
codebook firm_id if year == 2015
*--> en 2015 sur les 643 producteurs qu'on observe, les 446 plus petits ont produit 20Mt de CPO, ce qui est assez crédible.     

hist out_qty_ton_cpo if out_qty_ton_cpo <. & out_qty_ton_cpo<150000 
*Il y a un pb près de zero
codebook firm_id_year if out_qty_ton_cpo==0
*mais ce n'est pas les zero. C'est beaucoup de petites productions. 

histogram out_qty_ton_cpo if out_qty_ton_cpo <.& out_qty_ton_cpo <150 & out_qty_ton_cpo >0, bin(10) frequency


codebook firm_id_year if out_qty_ton_rpo >0 & out_qty_ton_rpo<. & year == 2015
codebook firm_id_year if out_qty_ton_cpo ==0

*Est-ce que c'est des obs. qui ont des valeurs particulières pour les autres commo? 



*Est-ce que des obs. qui ont une valeur particulièrement petite ? Un prix ? 



*******************************Bon voilà on à peu près la distribution importante à retenir. On peut enquêter un peu plus sur ces petites valeurs qui forment le plus gros décile.
*Ensuite l'idée ce sera de comparer la distribution des valeurs multioutput avec celle-ci. Si c'est systématiquement plus petit par exemple ? 
*On prendra les décisions finales de retirer des obs. quand on aura mergé avec la grosse base. pour l'instant produisons juste plein de flag informant les outliers, les zeros, les missings. 




**Ces obs. très proches de zero mais positive**
	*Une série d'histogrammes avec un plafond de plus en plus petit montre qu'il y a anormalement beaucoup d'obs. juste après zero, jusqu'à 150 tonnes/an.
	*histogram out_qty_ton_cpo if out_qty_ton_cpo <.& out_qty_ton_cpo <150 & out_qty_ton_cpo >0, bin(10) frequency
	gen small_out_qty_ton_cpo = (out_qty_ton_cpo < 150 & out_qty_ton_cpo>0 & out_qty_ton_cpo<.)
	*codebook firm_id_year if small_out_qty_ton_cpo ==1
	*browse if small_out_qty_ton_cpo==1
	*--> c'est beaucoup de LITER
	*Mais pas du tout tous les LITER puisqu'il y en a 4 fois plus. Donc on ne peut pas juste supprimer les LITER. 
	*histogram out_qty_ton_cpo if out_measurement_unit_cpo=="LITER"&out_qty_ton_cpo>0, bin(10) frequency

**Est-ce que c'est des obs. qui ont des valeurs particulières pour les autres commo? 
	*Il y a des valeurs tès élevées de out_qty_ton_pko qui sont pourtant des small_outliers de cpo --> pb car c'est censé être proportionnel (~x10). 
	*Donc ces obs. ont soit une erreur de saisie pour les pko soit pour cpo. 
	*histogram out_qty_ton_pko if out_qty_ton_pko >0 & small_out_qty_ton_cpo==1, bin(10) frequency
	*browse if out_qty_ton_pko >20000 & out_qty_ton_pko<. & small_out_qty_ton_cpo==1
	gen incoh_qty_cpo_pko =(out_qty_ton_pko >20000 & out_qty_ton_pko<. & small_out_qty_ton_cpo==1)
	*on peut regarder plus généralement à cette incohérence : 
	browse if 10*out_qty_ton_pko > out_qty_ton_cpo & out_qty_ton_pko<. 
	*--> il y en a 186 obs. qui déclarent produire plus de pko que de cpo ce qui n'est pas très normal. 
	*--> 275 pour x3
	*--> il y en a 1431 qui produisent seulement 5 fois plus de cpo que de pko. 

	***************************************************** IL FAUT SE RENSEIGNER MIEUX SUR CE RATIO
	browse if out_qty_ton_pko > out_qty_ton_cpo & out_qty_ton_pko<. 

	*Le seuil de normalité pour out_qty_ton_pko est 50
	*histogram out_qty_ton_pko if  out_qty_ton_pko <150000 & out_qty_ton_pko >0 & small_out_qty_ton_cpo==1, bin(10) frequency
	*browse if small_out_qty_ton_cpo ==1 & out_qty_ton_pko<50
	gen small_out_qty_ton_pko = (out_qty_ton_pko<50 &out_qty_ton_pko <. & out_qty_ton_pko >0)
	gen small_out_qty_ton = (small_out_qty_ton_cpo==1 & small_out_qty_ton_pko==1)
	*browse if  small_out_qty_ton==1
	*--> 85 obs. 

**Est-ce que des obs. qui ont une valeur particulièrement petite ? Un prix ? 

**Surtout il va falloir checker dans l'historique de chaque mill, pck là ces valeurs extrêmes sont souvent pas des répétitions pour une même mill. Donc on peut s'attendre à voir une chute une année. 

** et enfin, comparer aux inputs devrait aider. 
************************************************************************************ DONC IL FAUT MERGER O & I ! **************************************************************************************
* DONC VOILA LE PROGRAMME: ON FAIT UN GRAND DATA SET AVEC I & O SUR 1989-2015. 
* PUIS ON FAIT PASSER UNE BATTERIE DE TESTS Á CES VALEURS petites, SI ELLES EN RATENT UN SEUL, ELLES SONT SUPPRIMÉES.
   
* 1. Est-ce que la production de PKO est plus grande que celle de CPO. (reste à définir le ratio...), le calculer pour chaque obs et voir sa distribution. 
* 2. Est-ce que pour celles qui ont un bon ratio, la production des DEUX est encore anormalement faible ? 
* 3. Est-ce que le le prix est anormal ? 
* 4. Est-ce que lorsque le prix est normal, la valeur dans l'absolue est anormale ? 
* 5. Est-ce que les valeurs du même individu sont très différentes les autres années ? 
* 6. Est-ce que le ratio avec l'input est anormal ? 
* 7. Est-ce que lorsque le ratio output/input est normal, l'input ET l'output dans l'absolu sont anormaux ? 



	 gen price_cpo = out_value_cpo/out_qty_cpo 
	 gen price_pko = out_value_pko/out_qty_pko 
	 gen price_rpo = out_value_rpo/out_qty_rpo 
	 gen price_rpko = out_value_rpko/out_qty_rpko 
	 
	 browse out_qty_cpo out_value_cpo price_cpo out_prexcpo
	 drop price_cpo
	 drop price_pko
	 drop price_rpo
	 drop price_rpko


	 





	**agregate** 
	tostring firm_id, generate(firm_id_str)
	tostring year, generate(year_str)
	gen short_commo_code = commo_code-151400000
	tostring short_commo_code, generate(short_commo_code_str) usedisplayformat
	egen firm_id_year_commo_str = concat(firm_id_str year_str short_commo_code_str)
	destring firm_id_year_commo_str, generate(firm_id_year_commo)
	sort firm_id_year_commo
	by firm_id_year_commo, sort : egen double agr_output_qty = total(out_output_qty) 
	browse firm_id_year_commo year firm_id out_output_qty agr_output_qty 
	*--> ca ne mache pas, il y a plein de trucs bizarres. 




by firm_id_year_commo: egen agr_output_qty = total(out_output_qty)
gen flag_bizarre = ( out_output_qty!= agr_output_qty & flag_multioutput_per_commo_code==0)
inspect flag_bizarre 
drop agr_output_qty flag_bizarre

sort firm_id year 
sort firm_id_year
by firm_id_year, sort: egen ag_in_tot_qty_ton = total(in_tot_qty_ton) if flag_multiinput_per_commo_code==1 



	**aller voir le questionnaire pour voir ce qui change --> rien 

	**a-t-on le problème systématiquement pour les années autres que 2014 et 2015 comme pour 2055? 
inspect flag_multioutput_per_commo_code if year >2013
*--> 358 positive donc encore un pb ces années là. 


	**ou un pb de firm_id ? 








		***CLEANING OUT_OUTPOUT_QTY***
	**Missing out_output_qty**
		browse firm_id year out_measurement_unit out_output_qty out_output_value out_commodity_bah commo_code out_commodity_en if (out_output_qty >=. | out_output_qty ==0)&(commo_code==151410102|commo_code== 151410103|commo_code==151440101|commo_code==151440102)
		gen flag_missing_in_commo = ((out_output_qty >=. | out_output_qty ==0)&(commo_code==151410102|commo_code== 151410103|commo_code==151440101|commo_code==151440102))
		gen flag_missing_in_CPO = ((out_output_qty >=. | out_output_qty ==0)&(commo_code==151410102))
		total flag_missing_in_commo
		*231 obs.
		total flag_missing_in_CPO
		*160 obs.
		* Ne pas les supprimer tout de suite, il y aura peut-être de l'info dans l'autre base (si sales est différent de out_output_value... mais que dire alors? 
		*Un autre check pour ces obs.: est ce que l'output est renseigné pour leurs autres commo? Est ce que la value est différente ? 

		***CLEANING OUT_MEASUREMENT_UNIT***
* on peut extrapoler (c'est tout le temps KG) si mm ordre de grandeur. 
 
