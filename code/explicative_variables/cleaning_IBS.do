/*STRUCTURE: 
** Some pre-processing
***DUPLICATES***
***DEFLATE***
***OUTLIERS***
	0. statistical outliers
	1. output/input ratio
	2. cpo/pko ratio
	3. variation rate

***PRICE***
***AGGREGATED OBS.***
***PERCENTAGE EXPORTED***
***REFINING***

*** /// ADD ALL AVAILABLE GEOGRAPHIC INFORMATION from IBS_base_prelu.dta and IBS_desa2000.dta \\\ *** 


Over all the process, imp1 variables are associated with stronger imputations than imp2. 
It is to say that imp1 sample will be more modified than imp2, in an attempt to reduce noise.
  
*/

***** FIRST AND LAST YEAR ***** 
* We want to have the first year for all mills, even those that appeared before we have economic information on them (i.e. before 1998). 
* keep the oldest year for all IBS firms that existed before 1998. 
use  "input_data/IBS_kraus/IBS_final_panel.dta", clear
keep firm_id year 
keep if year < 1998
bys firm_id (year): egen min_year = min(year)
keep if year == min_year 
* among these old firms, keep only palm oil mills. 
merge 1:m firm_id using "temp_data/processed_IBS/IBS_PO_98_15.dta", generate(merge_mill_id) keepusing(firm_id year)
keep if merge_mill_id == 3 
drop merge_mill_id min_year
sort firm_id year
duplicates drop firm_id year, force 
* reintroduce the mill dataset, and now compute the min year. 
append using "temp_data/processed_IBS/IBS_PO_98_15.dta"
sort firm_id year
order year, after(firm_id)
bys firm_id (year): egen min_year = min(year)
bys firm_id (year): egen max_year = max(year)
order max_year, after(year)
order min_year, after(year)

drop if year < 1998 

sort firm_id year 



order  ag_out_qty_ton_cpo ag_out_value_cpo avg_out_prex_cpo ag_out_qty_cpo tag_multioutput_cpo out_measurement_unit_cpo  ///
	   ag_out_qty_ton_pko ag_out_value_pko avg_out_prex_pko ag_out_qty_pko tag_multioutput_pko out_measurement_unit_pko  ///
	   ag_out_qty_ton_rpo ag_out_value_rpo avg_out_prex_rpo ag_out_qty_rpo tag_multioutput_rpo out_measurement_unit_rpo  ///
	   ag_out_qty_ton_rpko ag_out_value_rpko avg_out_prex_rpko ag_out_qty_rpko tag_multioutput_rpko out_measurement_unit_rpko, after(industry_code)

order in_ton_ffb in_val_ffb flag_multiinput_ffb ///
	  in_dom_qty_ton_cpo in_dom_val_cpo in_imp_qty_ton_cpo in_imp_val_cpo in_tot_qty_ton_cpo in_tot_val_cpo flag_multiinput_cpo, after(industry_code)

order EKSPOR PRPREX, before(export_pct)
order export_pct_imp, after(export_pct)

rename ag_out_qty_* out_* 
rename ag_out_value* out_val* 
rename *qty_ton_cpo *ton_cpo 

global output_vars out_ton_cpo out_val_cpo out_ton_pko out_val_pko out_ton_rpo out_val_rpo out_ton_rpko out_val_rpko in_ton_ffb in_val_ffb
*NOTE: we add input_vars in output_vars because the rationale of their cleaning will be the same. 
*this is a bit messy but now we have also added cpo as an input, distinctively from domestic, imported or total. 
global in_cpo_vars in_dom_ton_cpo in_dom_val_cpo in_imp_ton_cpo in_imp_val_cpo in_tot_ton_cpo in_tot_val_cpo

* actually not necessary, because duplicates are plausible  
* global rspoproj_vars gifts bonus_prod bonus_oth pension_prod pension_oth royalties oth_expenses mfc_services


*Creating indicators of validity for output variables: 1 for obs. that are neither 0 nor missing. 
foreach var of varlist $output_vars {
	gen valid_`var' = (`var' <. & `var'>0)
}

gen valid_cpo = (valid_out_ton_cpo==1 & valid_out_val_cpo==1)
gen valid_pko = (valid_out_ton_pko==1 & valid_out_val_pko==1)
gen valid_rpo = (valid_out_ton_rpo==1 & valid_out_val_rpo==1)
gen valid_rpko = (valid_out_ton_rpko==1 & valid_out_val_rpko==1)

gen valid_ffb = (in_ton_ffb<. & in_ton_ffb>0 & in_val_ffb<. & in_val_ffb>0)

foreach var of varlist $in_cpo_vars {
	gen valid_`var' = (`var' <. & `var'>0)
}
global breakdown in_dom in_imp in_tot 
foreach br of global breakdown{
	gen valid_`br'_cpo = (valid_`br'_ton_cpo==1 & valid_`br'_val_cpo==1)
}

/*
foreach var of varlist $rspoproj_vars {
	gen valid_`var' = (`var' <. & `var'>0)
}
*/


*Creating lags for all output variables. 
sort firm_id year
foreach var of varlist $output_vars {
	bys firm_id : gen double lag_`var' = `var'[_n-1] 
}
global lag_output_vars lag_out_ton_cpo lag_out_val_cpo lag_out_ton_pko lag_out_val_pko lag_out_ton_rpo lag_out_val_rpo lag_out_ton_rpko lag_out_val_rpko ///
lag_in_ton_ffb lag_in_val_ffb /// 

sort firm_id year
foreach var of varlist $in_cpo_vars {
	bys firm_id : gen double lag_`var' = `var'[_n-1] 
}
global lag_in_cpo_vars lag_in_dom_ton_cpo lag_in_dom_val_cpo lag_in_imp_ton_cpo lag_in_imp_val_cpo lag_in_tot_ton_cpo lag_in_tot_val_cpo

/*
sort firm_id year
foreach var of varlist $rspoproj_vars {
	bys firm_id : gen double lag_`var' = `var'[_n-1] 
}
global lag_rspoproj_vars lag_gifts lag_bonus_prod lag_bonus_oth lag_pension_prod lag_pension_oth lag_royalties lag_oth_expenses lag_mfc_services
*/





*********************************************************************************************************************************************************************
***DUPLICATES***
*********************************************************************************************************************************************************************


**WITHIN FIRM_ID duplicates in terms of out_ton_ and out_val_ . 

	/* Why are there such duplicates?
	Either BPS, the Indonesian stat agency, has filled some missing cells with the information from the same firm one other year before. 
	Or the mill itself did that because it did not know the exact amount at the time of the survey. 
	Or it is a data manipulation error. 
	In the two first cases, which are most likely, the copies still contain some information, if they are attempts of estimation of the real amount. 
	So we will have two extreme rationales: 
	_imp1 will consider that these estimations are too noisy and should be rather missing; while 
	_imp2 will keep some of them, considering they might still contain some actual info. 
	*/

	*_imp1: any copy is too noisy* 
		// So for out_ton_cpo_imp1 it is 832 real changes to missing, and for out_val_cpo_imp1 it is only 29. 
		foreach var of varlist $output_vars {
			gen double `var'_imp1 = `var'
			order `var'_imp1, after(`var')
			gen flag_du_firm_`var' = (`var'== lag_`var' & valid_`var'==1)
			replace `var'_imp1 = . if flag_du_firm_`var'==1 
			drop flag_du_firm_`var'
		}

		foreach var of varlist $in_cpo_vars {
			gen double `var'_imp1 = `var'
			order `var'_imp1, after(`var')
			gen flag_du_firm_`var' = (`var'== lag_`var' & valid_`var'==1)
			replace `var'_imp1 = . if flag_du_firm_`var'==1 
			drop flag_du_firm_`var'
		}

		/*
		foreach var of varlist $rspoproj_vars {
			gen double `var'_imp1 = `var'
			order `var'_imp1, after(`var')
			gen flag_du_firm_`var' = (`var'== lag_`var' & valid_`var'==1)
			replace `var'_imp1 = . if flag_du_firm_`var'==1 
			drop flag_du_firm_`var'
		}	
		*/

	*_imp2: only duplicates in terms of both qty and value are too noisy, but when either one varies, there will be some information in the price that will be worth the noise. 
		// So for out_ton_cpo_imp2 it is only 16 real changes to missing. 	
		foreach var of varlist $output_vars {
			gen double `var'_imp2 = `var'
			order `var'_imp2, after(`var'_imp1)
		}

		foreach var of varlist $in_cpo_vars {
			gen double `var'_imp2 = `var'
			order `var'_imp2, after(`var'_imp1)
		}

		global commodity cpo pko rpo rpko 
		foreach commo of global commodity {
			gen flag_du_firm_qty_val_`commo' = (lag_out_ton_`commo' == out_ton_`commo' & lag_out_val_`commo'==out_val_`commo' & valid_`commo'==1)
		}
		gen flag_du_firm_qty_val_ffb = (lag_in_ton_ffb == in_ton_ffb & lag_in_val_ffb==in_val_ffb & valid_ffb==1)

		global breakdown in_dom in_imp in_tot
		foreach br of global breakdown{
			gen flag_du_firm_qty_val_`br'_cpo = (lag_`br'_ton_cpo == `br'_ton_cpo & lag_`br'_val_cpo == `br'_val_cpo & valid_`br'_cpo==1)
		}
		
		global out_cpo_vars out_ton_cpo out_val_cpo
		foreach var of varlist $out_cpo_vars {
			replace `var'_imp2 = . if flag_du_firm_qty_val_cpo==1
		}
		
		global out_pko_vars out_ton_pko out_val_pko
		foreach var of varlist $out_pko_vars {
			replace `var'_imp2 = . if flag_du_firm_qty_val_pko==1
		}		

		global out_rpo_vars out_ton_rpo out_val_rpo
		foreach var of varlist $out_rpo_vars {
			replace `var'_imp2 = . if flag_du_firm_qty_val_rpo==1
		}		

		global out_rpko_vars out_ton_rpko out_val_rpko
		foreach var of varlist $out_rpko_vars {
			replace `var'_imp2 = . if flag_du_firm_qty_val_rpko==1
		}


		global in_ffb_vars in_ton_ffb in_val_ffb
		foreach var of varlist $in_ffb_vars{
				replace `var'_imp2 = . if flag_du_firm_qty_val_ffb==1
			}

		foreach br of global breakdown{
			replace `br'_ton_cpo_imp2 = . if flag_du_firm_qty_val_`br'_cpo
			replace `br'_val_cpo_imp2 = . if flag_du_firm_qty_val_`br'_cpo
		}

		drop flag_du_firm_qty_val_cpo flag_du_firm_qty_val_pko flag_du_firm_qty_val_rpo flag_du_firm_qty_val_rpko flag_du_firm_qty_val_ffb ///
			flag_du_firm_qty_val_in_dom_cpo flag_du_firm_qty_val_in_imp_cpo flag_du_firm_qty_val_in_tot_cpo
		
		foreach var of varlist $lag_output_vars{
			drop `var'
		}
		foreach var of varlist $lag_in_cpo_vars{
			drop `var'
		}

**WITHIN YEAR duplicates in terms of out_ton_ and out_val_  

	/* Why are there such duplicates? Probably because BPS wanted to fill some missings when it was their first year (so they could not impute from previous years within same firm).
	Then, they may have used the information from an already existing mill, or a mill appearing this year too but without missing. 
	Here, we would like to keep the originals. If they were startup firms too, then we cannot identify them, but if among a group of duplicates one and only one obs. was existing before 
	then it is likely to be the original. 
	In the more imputational _imp1 variables, we remove also such original obs. 

	Also, we keep duplicates in terms of either only qty or value, considering that within a whole cross-section this is possible. 
	*/
		bysort firm_id : gen lag_year = year[_n-1]
		gen nostartup = (lag_year<.)

	*For CPO
		local a year out_ton_cpo out_val_cpo
		egen du_id_year_cpo = group(`a')
		duplicates tag `a' if valid_cpo==1, gen(tag_du_year_cpo)

		bysort du_id_year_cpo : egen n_nostartup_cpo = total(nostartup)
		gen original_cpo = (n_nostartup_cpo ==1 & nostartup==1)
		
		global imp1_cpo_vars out_ton_cpo_imp1 out_val_cpo_imp1 
		foreach var of varlist $imp1_cpo_vars {
			replace `var' = . if tag_du_year_cpo>0 & tag_du_year_cpo<. 
		}
		global imp2_cpo_vars out_ton_cpo_imp2 out_val_cpo_imp2 
		foreach var of varlist $imp2_cpo_vars {
			replace `var' = . if tag_du_year_cpo>0 & tag_du_year_cpo<. & original_cpo==0 
		}

	*And for PKO
		local b year out_ton_pko out_val_pko
		egen du_id_year_pko = group(`b')
		duplicates tag `b' if valid_pko==1, gen(tag_du_year_pko)

		bysort du_id_year_pko : egen n_nostartup_pko = total(nostartup)
		gen original_pko = (n_nostartup_pko ==1 & nostartup==1)
		
		global imp1_pko_vars out_ton_pko_imp1 out_val_pko_imp1 
		foreach var of varlist $imp1_pko_vars {
			replace `var' = . if tag_du_year_pko>0 & tag_du_year_pko<. 
		}
		global imp2_pko_vars out_ton_pko_imp2 out_val_pko_imp2 
		foreach var of varlist $imp2_pko_vars {
			replace `var' = . if tag_du_year_pko>0 & tag_du_year_pko<. & original_pko==0 
		}

	*And for RPO 
		local c year out_ton_rpo out_val_rpo
		egen du_id_year_rpo = group(`c')
		duplicates tag `c' if valid_rpo==1, gen(tag_du_year_rpo)

		bysort du_id_year_rpo : egen n_nostartup_rpo = total(nostartup)
		gen original_rpo = (n_nostartup_rpo ==1 & nostartup==1)
		
		global imp1_rpo_vars out_ton_rpo_imp1 out_val_rpo_imp1 
		foreach var of varlist $imp1_rpo_vars {
			replace `var' = . if tag_du_year_rpo>0 & tag_du_year_rpo<. 
		}
		global imp2_rpo_vars out_ton_rpo_imp2 out_val_rpo_imp2 
		foreach var of varlist $imp2_rpo_vars {
			replace `var' = . if tag_du_year_rpo>0 & tag_du_year_rpo<. & original_rpo==0 
		}

	*And for RPKO
		local d year out_ton_rpko out_val_rpko
		egen du_id_year_rpko = group(`d')
		duplicates tag `d' if valid_rpko==1, gen(tag_du_year_rpko)

		bysort du_id_year_rpko : egen n_nostartup_rpko = total(nostartup)
		gen original_rpko = (n_nostartup_rpko ==1 & nostartup==1)
		
		global imp1_rpko_vars out_ton_rpko_imp1 out_val_rpko_imp1 
		foreach var of varlist $imp1_rpko_vars {
			replace `var' = . if tag_du_year_rpko>0 & tag_du_year_rpko<. 
		}
		global imp2_rpko_vars out_ton_rpko_imp2 out_val_rpko_imp2 
		foreach var of varlist $imp2_rpko_vars {
			replace `var' = . if tag_du_year_rpko>0 & tag_du_year_rpko<. & original_rpko==0 
		}

	*And for FFB
		local e year in_ton_ffb in_val_ffb
		egen du_id_year_ffb = group(`e')
		duplicates tag `e' if valid_ffb==1, gen(tag_du_year_ffb)

		bysort du_id_year_ffb : egen n_nostartup_ffb = total(nostartup)
		gen original_ffb = (n_nostartup_ffb ==1 & nostartup==1)
		
		global imp1_ffb_vars in_ton_ffb_imp1 in_val_ffb_imp1 
		foreach var of varlist $imp1_ffb_vars {
			replace `var' = . if tag_du_year_ffb>0 & tag_du_year_ffb<. 
		}
		global imp2_ffb_vars in_ton_ffb_imp2 in_val_ffb_imp2 
		foreach var of varlist $imp2_ffb_vars {
			replace `var' = . if tag_du_year_ffb>0 & tag_du_year_ffb<. & original_ffb==0 
		}

	*And for in_cpo we do it at the level of total input. 
		local f year in_tot_ton_cpo in_tot_val_cpo
		egen du_id_year_in_cpo = group(`f')
		duplicates tag `f' if valid_in_tot_cpo==1, gen(tag_du_year_in_tot_cpo)

		bysort du_id_year_in_cpo : egen n_nostartup_in_tot_cpo = total(nostartup)
		gen original_in_tot_cpo = (n_nostartup_in_tot_cpo ==1 & nostartup==1)
		
		global imp1_in_tot_cpo_vars in_tot_ton_cpo_imp1 in_tot_val_cpo_imp1 
		foreach var of varlist $imp1_in_tot_cpo_vars {
			replace `var' = . if tag_du_year_in_tot_cpo>0 & tag_du_year_in_tot_cpo<. 
		}
		global imp2_in_tot_cpo_vars in_tot_ton_cpo_imp2 in_tot_val_cpo_imp2 
		foreach var of varlist $imp2_in_tot_cpo_vars {
			replace `var' = . if tag_du_year_in_tot_cpo>0 & tag_du_year_in_tot_cpo<. & original_in_tot_cpo==0 
		}

drop lag_year nostartup du_id_year_cpo n_nostartup_cpo tag_du_year_cpo /// 
du_id_year_pko tag_du_year_pko n_nostartup_pko original_cpo original_pko  ///
du_id_year_rpo tag_du_year_rpo n_nostartup_rpo original_rpo /// 
du_id_year_rpko tag_du_year_rpko n_nostartup_rpko original_rpko ///
du_id_year_ffb tag_du_year_ffb n_nostartup_ffb original_ffb ///
du_id_year_in_cpo tag_du_year_in_tot_cpo n_nostartup_in_tot_cpo original_in_tot_cpo

/* 
browse firm_id year du_id_year_cpo tag_du_year_cpo if tag_du_year_cpo>0 & tag_du_year_cpo<.
One can see that in most cases these duplicates are duplicates within firm_id also. This means that our explanation "they filled missing startupers with existing firms" is not so true. 
So there may be another story behind these duplicates. Anyway, they are already removed from the duplicates treatment above then. Except for those that encompass startupers, then 
the program here saved 
*/









**********************************************************************************************************************************************************
***DEFLATE***
**********************************************************************************************************************************************************


*merge m:1 year using /Users/valentinguye/Desktop/wpi_palmoil_1998_2015.dta, nogenerate 
merge m:1 year using "input_data/macro/wpi_palmoil_1998_2015.dta", nogenerate 
*Change reference year from 2000 to 2010: 
replace wpi = wpi/290.133117675781
rename wpi wpi_2010
*codebook wpi_2010 if year == 2010

global monetary_vars in_val_ffb ///
					 out_val_cpo ///
					 out_val_pko ///
					 out_val_rpo ///
					 out_val_rpko ///
					 out_val_cpo_imp1 ///
					 out_val_pko_imp1 ///
					 out_val_rpo_imp1 ///
					 out_val_rpko_imp1 ///
					 in_val_ffb_imp1 ///
					 in_dom_val_cpo_imp1 ///
					 in_imp_val_cpo_imp1 ///
					 in_tot_val_cpo_imp1 ///
					 out_val_cpo_imp2 ///
					 out_val_pko_imp2 ///
					 out_val_rpo_imp2 ///
					 out_val_rpko_imp2 ///
					 in_val_ffb_imp2 ///
					 in_dom_val_cpo_imp2 ///
					 in_imp_val_cpo_imp2 ///
					 in_tot_val_cpo_imp2 /// 
					 	gifts ///
					 	bonus_prod ///
						bonus_oth ///
						pension_prod ///
						pension_oth ///
						royalties ///
						oth_expenses ///
						mfc_services 
/* 
Convert to full IDR (monetary values are measured in 1000 IDR in survey)
And deflate to 2010 price level. 
*/
foreach mon_var of varlist $monetary_vars {
	quietly replace `mon_var' = `mon_var'*1000
	quietly replace `mon_var' = `mon_var'/wpi_2010
}

// Convert to 2010 Dollars

foreach mon_var of varlist $monetary_vars {
	quietly replace `mon_var' = `mon_var'/9090.433333
}








**********************************************************************************************************************************************************
***OUTLIERS***
**********************************************************************************************************************************************************



/* 

The main idea behind removing qty and val extreme values and not merely extreme prices is that, for some reason, there are impossibly high qty and val, and we don't want to risk 
their ratio is deemed a normal price and these obs. used in the analyses. 

First, let us define otls. 
Statistically, they are those obs for which var > (p75+1.5iqr). And to be more precise we want to apply this definition within yearly cross-sections. 
Let us restrict this definition: otls are the obs that are statistical otls AND fail at least one test of likelihood. 
These tests use output/input ratio and CPO/PKO ratio and the variation_rate over time.
In other words, we will keep statistical otls that are however coherent with ALL the information we have. 
We don't extand the definition of otls (i.e. deem as otls obs. that are not statistical otls but are however not coherent with some of the information we have)
because this would require us to have high quality information on these other variables used (inputs and pko), which is not the case. This means 
that we don't see cases of extremely high say output/input ratio with normal output quantity as a hint of flawed measurement of output. We assume this is purely due to
input mismeasurement, and never because of a combination of mismeasurement of both. 

Regarding low otls, i.e. <(p25 - 1.5iqr). In theory, they should also be a concern. 
	...but there is none (for out_ton_cpo and out_val_cpo). For the ratios in the tests, we are not interested in spotting low otls because, since there are 
	only high otls for the variables of interest, they would only spot abnormally high secondary information (input, pko or previous). 

Regarding the test on variation_rate.  
 	Why is it ok to add it for cpo and pko qty and not for value nor for rpo and rpko qty? 
 	Because as such, this test is not very precise (and designing more sophisticated one would be time consuming and not necessarily worth it).
 	For out_ton_, including it makes an additional condition to the loosing of the otl definition, i.e. adding this test makes always less obs. "saved" from being statistical otls. 
 	On the other hand, adding it when it is the only test, in the out_val_ cases, means that adding this test makes always more obs. "saved", since otherwise you only use statistical 
 	otl definition. 

Regarding obs. that are value otl but not qty otl. 
	They have been deemed credible obs. with big quantities produced. So we keep them even if values are otl. If it appears that values are really too high, then the cleaning on 
	prices outliers will remove them. 


INFINE, 4 cases: 
	- Neither out_ton_cpo nor out_val_cpo are otls --> the obs. is kept. 
	- out_ton_cpo is an otl (it is a statistical otl and it failed one of tests 1,2,3) and out_val_cpo is not an otl --> obs is droped
	- out_ton_cpo is not an otl but out_val_cpo is --> obs is KEPT. 
	- out_ton_cpo and out_val_cpo are otls --> obs is droped. 


Some notes on the codes below. 
	- tabstat is not conditioned because we remove zeros (see just below) and missings don't impact the distribution calculation. 

	- Taking the `var'<. in the if and not in the () enables to interpret otl_`var'=0 as "the obs. is not an otl", and not as
	"the obs. is not an otl OR it is missing on out_ton_cpo". Recall, this "if" ensures that missings are not counted as otls.

Beside, we will only modify _imp1 and _imp2, and the procedures need to be done separately for each, because they don't have the same distributions
and hence not the same otl thresholds etc. 


In order for the 0 cells not to be accounted as such in the tabstats. This is justified by the fact that we won't be able to interpret and hence to use a price coming from a firm 
with 0 in one of these cells. 
*/

global imp_out_vars out_ton_cpo_imp1 out_val_cpo_imp1 /// 
					out_ton_pko_imp1 out_val_pko_imp1 /// 
					out_ton_rpo_imp1 out_val_rpo_imp1 /// 
					out_ton_rpko_imp1 out_val_rpko_imp1 /// 
					in_ton_ffb_imp1 in_val_ffb_imp1 ///
					in_dom_ton_cpo_imp1 in_dom_val_cpo_imp1 ///
					in_imp_ton_cpo_imp1 in_imp_val_cpo_imp1 ///
					in_tot_ton_cpo_imp1 in_tot_val_cpo_imp1 ///
					out_ton_cpo_imp2 out_val_cpo_imp2 /// 
					out_ton_pko_imp2 out_val_pko_imp2 /// 
					out_ton_rpo_imp2 out_val_rpo_imp2 /// 
					out_ton_rpko_imp2 out_val_rpko_imp2 /// 
					in_ton_ffb_imp2 in_val_ffb_imp2 ///
					in_dom_ton_cpo_imp2 in_dom_val_cpo_imp2 ///
					in_imp_ton_cpo_imp2 in_imp_val_cpo_imp2 ///
					in_tot_ton_cpo_imp2 in_tot_val_cpo_imp2 

foreach var of varlist $imp_out_vars {
	quietly replace `var'= . if `var' == 0 
}

/*

First let's compute the flags for outliers regarding additional information. 
This does not need to be done in a yearly  fashion. 

1. OUTPUT/INPUT RATIO
	
*/

	capture program drop output_input_ratio
	program define output_input_ratio
	gen oir_`1' = `1'/in_ton_ffb if `1'<. & in_ton_ffb <. & in_ton_ffb >0
	quietly tabstat oir_`1', statistics(iqr, p75) columns(variables) save
	matrix stats = r(StatTotal)
	local iqr_oir_`1' = stats[1,1]
	local p75_oir_`1' = stats[2,1]
	local otl_th_oir_`1' = `p75_oir_`1'' + 1.5*`iqr_oir_`1''
	gen otl_oir_`1' = (oir_`1'> `otl_th_oir_`1'') if oir_`1'<.
	drop oir_`1'
	*display `otl_th_oir_`1''
	end
	output_input_ratio out_ton_cpo_imp1
	output_input_ratio out_ton_cpo_imp2
	output_input_ratio out_ton_pko_imp1
	output_input_ratio out_ton_pko_imp2

	*browse out_ton_cpo oir if otl_oir==1 & otl_outtoncpo==0

	/*
	So among the 423 outtoncpo otls, 140 are also otls in terms of oir. Meaning that for these, both output is abnormal alone and it is abnormal with respect to input. 
	the remaining 283 outtoncpo otls are still to be evaluated with creterion below. 
	Regarding the 722 outtoncpo-non-otl that however are OIR otls, the only thing we can do is assuming this is because of an error on the sole input cell. 

	Note also that the median is .22 which is in line with the oil/bunch ratio of 21-23% in the literature. The mean is highly falsed by otls. 


2. CPO/PKO RATIO

 	*/
 	forvalues a = 1/2 {
		gen cpopko_ratio_imp`a' = out_ton_cpo_imp`a'/out_ton_pko_imp`a' if out_ton_cpo_imp`a'<. & out_ton_pko_imp`a'<.
		quietly tabstat cpopko_ratio_imp`a' if cpopko_ratio_imp`a'<., statistics(iqr, p75) columns(variables) save
		matrix stats = r(StatTotal)
		local iqr_cpopko_ratio = stats[1,1]
		local p75_cpopko_ratio = stats[2,1]
		local otl_th_cpopko_ratio = `p75_cpopko_ratio' + 1.5*`iqr_cpopko_ratio'
		gen otl_cpopko_ratio_imp`a' = (cpopko_ratio_imp`a' > `otl_th_cpopko_ratio') if cpopko_ratio_imp`a'<.
		drop cpopko_ratio_imp`a'
		*display `otl_th_cpopko_ratio'
	}

	/*	
	Here again, the median of cpopko_ratio of 4.5 is pretty broadly consistent with the literature giving an oil/bunch ratio between 3 and 4 times higher than the kernel/bunch ratio. 
	In Byerlee et al. (book) Table B2.1 the oil and kernel extraction rates are respectively 22.8 and 4.70 in Sime Darby plantations in 2012/2013 in Indonesia. 

	This is while still accounting for otls, the actual descriptive statistics will come later of course. 


3. VARIATION RATE
*/	
	sort firm_id year
	global imp_cpo_pko_ffb_ton_vars out_ton_cpo_imp1 out_ton_pko_imp1 out_ton_cpo_imp2 out_ton_pko_imp2 in_ton_ffb_imp1 in_ton_ffb_imp2
	foreach var of varlist $imp_cpo_pko_ffb_ton_vars {	
		bys firm_id : gen double lag_`var' = `var'[_n-1]	
		gen dt_`var' = (`var'-lag_`var')/(lag_`var') if `var'<. & lag_`var'<.
		quietly tabstat dt_`var', statistics(iqr p75) columns(variables) save
		matrix stats = r(StatTotal)
		local iqr_dt_`var' = stats[1,1]
		local p75_dt_`var' = stats[2,1]
		local otl_th_dt_`var'  = `p75_dt_`var'' + 1.5*`iqr_dt_`var''
		quietly gen otl_dt_`var' = (dt_`var' > `otl_th_dt_`var'') if dt_`var'<.
		drop lag_`var' dt_`var' 
		*display `otl_th_dt_`var''
	} 

	
/*


4. "STATISTICAL" OUTLIERS **

	This is to define statistical outliers ACROSS years	
		foreach var of varlist $imp_out_vars {
			quietly tabstat `var', statistics(iqr p75) columns(variables) save
			matrix stats = r(StatTotal)
			local iqr_`var' = stats[1,1]
			local p75_`var' = stats[2,1]
			local otl_th_`var' = `p75_`var'' + 1.5*`iqr_`var''
			quietly gen otl_`var' = (`var' > `otl_th_`var'') if `var'<.
			*display `otl_th_`var''
		}

	This is to define statistical outliers WITHIN years 	
	*/
	foreach var of varlist $imp_out_vars {
		quietly gen otl_`var' = 0 if `var'<.
		forvalues y = 1998/2015 {
			quietly tabstat `var' if year == `y', statistics(iqr p75) columns(variables) save
			matrix stats = r(StatTotal)
			local iqr_`var' = stats[1,1]
			local p75_`var' = stats[2,1]
			local otl_th_`var' = `p75_`var'' + 1.5*`iqr_`var''
			quietly replace otl_`var' = 1 if year == `y' & `var' > `otl_th_`var'' & `var'<.
			*display `otl_th_`var''
		}
	}
/*
browse if otl_out_ton_cpo_imp1==1 & out_measurement_unit_cpo!="Ton"
it does not seem that the mesurement unit conversions made in IBS_output_preparation.do have triggered unproportionnaly more otls. 
browse firm_id if firm_id == 2117 | firm_id == 10451 | firm_id == 56313 | firm_id == 68245 | firm_id == 70379



**DEFINE WHAT IS NOT AN OUTLIER IN THE END**

We should be be careful that when additional information is missing (there is no input data, or pko data, or it is the first observed year for this mill) 
the obs. is neither automatically discarded or kept: the statistical definition of otl should decide alone. 
generate otl_cpo = (otl_out_ton_cpo==1 & (otl_oir==1 | otl_cpopko_ratio==1 | otl_dt_out_ton_cpo==1))
for instance, this synthax is not good because statistical otls that happen to have no information on input, pko and t-1 out_ton_cpo fill non of the conditions in the second parentheses
thus they would be kept while we want the statistical otl definition to rule them out in this case. 
*/

* For CPO and PKO quantities
	forvalues a = 1/2 {
		global imp`a'_cpo_pko_ton_vars out_ton_cpo_imp`a' out_ton_pko_imp`a'
		foreach var of global imp`a'_cpo_pko_ton_vars {
			generate not_otl_`var' = (otl_`var'==0 | ///
				(otl_`var'==1 & otl_oir_`var'==0 & otl_cpopko_ratio_imp`a'==0 & otl_dt_`var'==0)) ///
				if otl_`var' <. 
			replace `var' = . if not_otl_`var'==0 
		}
	}
			/* This was just to create a flag for some checks. 
			forvalues a = 1/2 {
				global imp`a'_cpo_pko_ton_vars out_ton_cpo_imp`a' out_ton_pko_imp`a'
				foreach var of global imp`a'_cpo_pko_ton_vars {
					gen saved_`var' = (otl_`var'==1 & otl_oir_`var'==0 & otl_cpopko_ratio_imp`a'==0 & otl_dt_`var'==0) if otl_`var'<. 
					inspect saved_`var'
				}
			}
			*/
	/*
	With this code:
		Suppose you are missing on qty_cpo, then your otl_out_ton_cpo is neither 1 nor 0, and get a 0 here. 
		Suppose you are missing on val_cpo, then you should be discarded, but later, not on the behalf of being an out_ton_cpo otl.   
		Suppose you are missing on additional information (input, pko or lag), then you cannot have related flags equal to zero and what you get depends only on the first line. 
	So this synthax enables to interpret not_otl_outtoncpo==0 as otls among out_ton_cpo-valid information.

		*Refined outliers removal*
		We impute on a case by case basis among the remaining extreme values. 
			*For CPO
			There is still one extreme value for out_ton_cpo_imp2: It would not be used for final analyses because its value is extreme and removed. 
			But let us still remove it so it does not flaw aggregation checks implying qty for instance. 
			*/
			replace out_ton_cpo_imp2 = . if out_ton_cpo_imp2> 1000000 & out_ton_cpo_imp2<.
			*browse year if out_ton_cpo_imp2> 200000 & out_ton_cpo_imp2<.
			*the other extreme values "saved" with out_ton_cpo_imp2 > 150,000 let's just keep them in imp2, they are not in imp1. 

			*For PKO, there remain important outliers (5 in imp2 out of which 2 in imp1). They are all outliers on val_imp1 and 2 anyway. Let's just remove those from imp1. 
			replace out_ton_pko_imp1 = . if out_ton_pko_imp1 > 200000 & out_ton_pko_imp1<. 


* For CPO and PKO VALUES 
/*
	Outliers are merely statistical outliers only among obs. whose quantities have been deemed non-credible (i.e. outliers). 
*/	
	forvalues a = 1/2 {
		global commodity cpo pko 
		foreach commo of global commodity {
			replace out_val_`commo'_imp`a' = . if otl_out_val_`commo'_imp`a'==1 & not_otl_out_ton_`commo'_imp`a'==0
			*gen flag_otl_val_`commo'_imp`a' = (otl_out_val_`commo'_imp`a'==1 & not_otl_out_ton_`commo'_imp`a'==1)
			*quietly drop not_otl_out_ton_`commo'_imp`a'
		}
	}
/*
browse firm_id year out_val_cpo_imp1 out_ton_cpo out_ton_cpo_imp1 out_ton_cpo_imp2 avg_out_prex_cpo out_val_rpo out_rpo out_ton_rpko revenue_total revenue_total_imp1 revenue_total_imp2 otl_out_ton_cpo_imp1 ///
if  out_val_cpo_imp1>600000000 & out_val_cpo_imp1<.
*/
replace out_val_cpo_imp1 = . if out_val_cpo_imp1>600000000
replace out_val_cpo_imp2 = . if out_val_cpo_imp2>600000000

*graph hbox out_val_cpo_imp1
/*

To summarize: 
# removed obs. respectively for: 	out_ton_cpo_imp1 ; out_ton_pko_imp1 ; out_ton_cpo_imp2 ; out_ton_pko_imp2 ; out_val_cpo_imp1 ; out_val_pko_imp1 ; out_val_cpo_imp2 ; out_val_pko_imp2
Within year otl & with deflation	337					369					388					452				406					534					406					534
Within year otl & w/o deflation 	all the same
Across year otl & with deflation	329					392					384					478				426					565					426					565
Across year otl & w/o deflation 	329					392					384					478				357					517					357					517

What's to take away: 
Many qty outliers are duplicates on qty (likely within firm_id, since this is where most duplicates where) that kept only in imp2. 
Second line being the same as first line is logic: deflation does not change qty, nor does it change the identification of outliers for values if this is made within year anyway. 
If one defines outliers across the whole panel, more values are identified as being outliers after being deflated. 
Defining outliers across the whole panel rather than within each year does not lead to a very different number of obs. being removed. 

With the new code that keeps values outliers if they are not qty outliers, the first row (within year & with deflation) is 186 279 217 344. 
Hence many obs. are kept this way. 

*For RPO AND RPKO 
For now the only purpose is saying how much cpo weighs in one unit's production. 
We areconcerned that variance can be very strong for these variables because we can have both big refineries and many small mills doing only a bit of rpo. 
Hence we remove first very strong outliers. They produce more than 15 million ton in one year, which is not so credible. 
*/
replace out_ton_rpo_imp1 = . if out_ton_rpo_imp1 > 15000000


	/*
	forvalues a = 1/2 {
		forvalues y = 1998/2015{
				quietly tabstat out_ton_rpo_imp`a' if year == `y' , format statistics(sum) columns(variables) save
				matrix stats = r(StatTotal)
				global sample`a'00_out_ton_rpo_`y' = stats[1,1]
		}	
	}


	forvalues year = 1998/2015{ 
		display ${sample100_out_ton_rpo_`year'}, ${sample200_out_ton_rpo_`year'}
	}


browse firm_id year otl_out_ton_rpo_imp1  otl_out_ton_rpo_imp2 otl_out_ton_rpko_imp1 otl_out_ton_rpko_imp2 ///
out_ton_rpo_imp1  out_ton_rpo_imp2 out_ton_rpko_imp1 out_ton_rpko_imp2 ///
 if otl_out_ton_rpo_imp1==1 | otl_out_ton_rpo_imp2==1 |otl_out_ton_rpko_imp1==1|otl_out_ton_rpko_imp2==1

forvalues a=1/2 {
	global imp`a'_rpo_rpko_ton_vars out_ton_rpo_imp`a' out_ton_rpko_imp`a'
	foreach var of global imp`a'_rpo_rpko_ton_vars {
		replace `var' = . if otl_`var'==1 
	}
}

*/

/*
For FFB 
	for qty, saving some statistical outliers when they have a normal output/input ratio is not a good idea since this input is in the denumerator. 
	The following procedure using the input/output ratio saves 24 and 38 obs. for imp1 and 2 resp. compared to merely removing statistical outliers. 
	Let's do this because we seemingly have some big mills in the sample. 
	*/
	forvalues a = 1/2 {
		gen double ior_imp`a' = in_ton_ffb_imp`a'/out_ton_cpo_imp`a' if out_ton_cpo_imp`a'<. & in_ton_ffb_imp`a' <.
		quietly tabstat ior_imp`a', statistics(iqr p75) columns(variables) save
		matrix stats = r(StatTotal)
		local iqr_ior_imp`a' = stats[1,1]
		local p75_ior_imp`a' = stats[2,1]
		local otl_th_ior_imp`a' = `p75_ior_imp`a'' + 1.5*`iqr_ior_imp`a''
		gen otl_ior_imp`a' = (ior_imp`a'>`otl_th_ior_imp`a'') if ior_imp`a'<.
		drop ior_imp`a'
	}

	forvalues a = 1/2 {
		generate not_otl_ton_ffb_imp`a' = (otl_in_ton_ffb_imp`a'==0 | ///
			(otl_in_ton_ffb_imp`a'==1 & otl_ior_imp`a'==0 & ///
			otl_dt_in_ton_ffb_imp`a'==0)) ///
			if otl_in_ton_ffb_imp`a' <. 
		replace  in_ton_ffb_imp`a' = . if not_otl_ton_ffb_imp`a'==0 
	}


	*values 
	forvalues a = 1/2 {
		replace in_val_ffb_imp`a' = . if otl_in_val_ffb_imp`a'==1 & not_otl_ton_ffb_imp`a'==0
		drop not_otl_ton_ffb_imp`a'
	}

/*
*For Input CPO 
	We don't use any input/output ratio to confirm high values because cpo can be processed into many different goods.
	*/
	forvalues a = 1/2 {
		local imp_in_cpo in_dom_ton_cpo_imp`a' in_dom_val_cpo_imp`a' in_imp_ton_cpo_imp`a' in_imp_val_cpo_imp`a' in_tot_ton_cpo_imp`a' in_tot_val_cpo_imp`a'
		foreach var of varlist `imp_in_cpo'{
			replace `var' = . if otl_`var'==1
		}
	}

*graph hbox out_ton_cpo_imp1 if otl_out_ton_cpo_imp1 ==0
*codebook out_ton_cpo_imp1 out_ton_cpo_imp2

*sum out_ton_cpo out_ton_cpo_imp1 out_ton_cpo_imp2, detail


*graph hbox out_val_cpo 
*sum out_val_cpo, detail
*/

drop otl_oir_out_ton_cpo_imp1 ///
otl_oir_out_ton_cpo_imp2 ///
otl_oir_out_ton_pko_imp1 ///
otl_oir_out_ton_pko_imp2 ///
otl_cpopko_ratio_imp1 ///
otl_cpopko_ratio_imp2 ///
otl_dt_out_ton_cpo_imp1 ///
otl_dt_out_ton_pko_imp1 ///
otl_dt_out_ton_cpo_imp2 ///
otl_dt_out_ton_pko_imp2 ///
otl_dt_in_ton_ffb_imp1 /// 
otl_dt_in_ton_ffb_imp2 /// 
not_otl_out_ton_cpo_imp1 /// 
not_otl_out_ton_pko_imp1 /// 
not_otl_out_ton_cpo_imp2 /// 
not_otl_out_ton_pko_imp2 /// 
otl_ior_imp1 otl_ior_imp2

foreach var of varlist $imp_out_vars{
	drop otl_`var'
}




*******************************************************************************************************************************************************************************************************
***PRICE***
*******************************************************************************************************************************************************************************************************




*CPO price
forvalues a = 1/2 {
	gen cpo_price_imp`a' = out_val_cpo_imp`a'/out_ton_cpo_imp`a'
	order cpo_price_imp`a', before(out_ton_cpo)
}

*PKO price 
forvalues a = 1/2 {
	gen pko_price_imp`a' = out_val_pko_imp`a'/out_ton_pko_imp`a'
	order pko_price_imp`a', before(out_ton_pko)
}

*FFB price 
forvalues a = 1/2 {
	gen ffb_price_imp`a' = in_val_ffb_imp`a'/in_ton_ffb_imp`a'
	order ffb_price_imp`a', before(in_ton_ffb)
}

*in_cpo prices 
forvalues a = 1/2 {
	foreach br of global breakdown{
		gen `br'_cpo_price_imp`a' = `br'_val_cpo_imp`a'/`br'_ton_cpo_imp`a'
		order `br'_cpo_price_imp`a', before(`br'_ton_cpo)
	}	
}


/* 
There are many otls in the price. We attribute them to extremely low quantities which we have not treated as outliers supra. 
Indeed it is difficult to say until which minimum production of a commodity a mill can go, it can credibly produce very little sometimes. 
For imputation, we couldn't use input or byproducts because the mismeasurement can as well be rather there. 
So by removing price upper outliers we remove obs. that have either mismeasurement of quantity (too low) relative to value, or mismeasurement of value (too high though not outlier) relative to 
a true small quantity. 
And by removing price lower outlier we remove obs. that have either mismeasurement of value (too low) relative to quantity, or mismeasurement of quantity (too high though not outlier) relative to 
a true small value. 

For these two processes, is there a reason to spot outliers within years? To me no, because extreme values are only deemed as such with respect to a contemporaneous information.

*/

global imp_prices cpo_price_imp1 pko_price_imp1 ffb_price_imp1 cpo_price_imp2 pko_price_imp2 ffb_price_imp2 ///
	 in_dom_cpo_price_imp1 in_dom_cpo_price_imp2 in_imp_cpo_price_imp1 in_imp_cpo_price_imp2 in_tot_cpo_price_imp1 in_tot_cpo_price_imp2
foreach var of varlist $imp_prices{
	quietly tabstat `var', statistics(iqr p75 p25) columns(variables) save
	matrix stats = r(StatTotal)
	local iqr_`var' = stats[1,1]
	local p75_`var' = stats[2,1]
	local p25_`var' = stats[3,1]
	local otl_uth_`var' = `p75_`var'' + 1.5*`iqr_`var''
	local otl_lth_`var' = `p25_`var'' - 1.5*`iqr_`var''
	quietly gen u_otl_`var' = (`var'>=`otl_uth_`var'') if `var'<.
	quietly gen l_otl_`var' = (`var'<=`otl_lth_`var'') if `var'<.
	*display `otl_lth_`var''
	*count if u_otl_`var'==1 
	*count if l_otl_`var'==1 
	quietly replace `var' = . if u_otl_`var' ==1 
	quietly replace `var' = . if l_otl_`var' ==1 
	drop u_otl_`var'
	drop l_otl_`var'
}

/*Among those obs. that are val outliers but qty not outliers (resp. 72 and 167 for cpo) there are resp. 31 and 115 price upper outliers, i.e. obs. with really too high values. 
And among this same group the obs. "saved" are 15 (resp 20) among price not upper outlier and 1 (resp 2) among price outliers. 
This shows that our "saving" obs. saved val outliers that are yet credible price wise. 
*/

*************************************************************************************************************************************************************************************
***AGGREGATED OBS.***
*************************************************************************************************************************************************************************************







/*
/*At the end of IBS_output preparation, FOR CPO, there are 317 obs. that come from an aggregation. 
The remaining are either zero (not coming from an aggregation) or missing (the obs. is only producing other commodities, 1210 instances). 
Now (once merged), the number of missing has raised of course (no matches) and the number of obs. coming from an aggregation is still 317. 
Among which 65 are otls. Hence an important proportion, raising awareness that other aggregated obs. that are not deemed as otls are problematic as they are likely to be inflated. 
If one takes the statistical definition of otls, there are still 65 instances coming from aggregation. 
*/
inspect tag_multioutput_cpo
count if tag_multioutput_cpo >0 & tag_multioutput_cpo<. & otl_cpo == 1 
count if tag_multioutput_cpo >0 & tag_multioutput_cpo<. & (otl_out_ton_cpo == 1 | otl_out_val_cpo==1)
count if tag_multioutput_cpo >0 & tag_multioutput_cpo<. & (otl_out_ton_cpo == 1)
count if tag_multioutput_cpo >0 & tag_multioutput_cpo<. & (otl_out_val_cpo == 1)
count if tag_multioutput_cpo >0 & tag_multioutput_cpo<. & (otl_cpo_price == 1)
browse if tag_multioutput_cpo >0 & tag_multioutput_cpo<. & otl_cpo == 1 
browse if tag_multioutput_cpo >0 & tag_multioutput_cpo<.
/*
For remaining price otls, 7 obs come from aggregation and are price otls, of the 190 price otls at this stage. And these are not otls in terms of qty or value. 
In fact, they are small qtys with normal values.   
*/
browse if tag_multioutput_cpo >0 & tag_multioutput_cpo<. & (otl_cpo_price == 1)

/*
Voilà, sachant cela et avec la réponse de BPS on décidera de garder ces 247 (317-70) obs. ou non.  
On peut aussi checker comment ca se passe pour pko. */
inspect tag_multioutput_pko
*There are 326 obs. coming from aggregation of pko. 
count if tag_multioutput_pko >0 & tag_multioutput_pko<. & otl_pko == 1 
*Parmi lesquelles 179 ont été considérées comme otls ! 
count if tag_multioutput_pko >0 & tag_multioutput_pko<. & (otl_out_ton_pko == 1 | otl_out_val_pko==1)
*184 si on considère les définitions statistiques d'otl (qui produisent davantage d'otls)
count if tag_multioutput_pko >0 & tag_multioutput_pko<. & (otl_out_ton_pko == 1)
*160 obs. venant de l'aggrégation de pko sont des otls en terme de quantité 
count if tag_multioutput_pko >0 & tag_multioutput_pko<. & (otl_out_val_pko == 1)
*et 155 en termes de value. 
count if tag_multioutput_pko >0 & tag_multioutput_pko<. & (otl_pko_price == 1)
*seulement 6 en termes de prix: lorsque l'aggrégation a produit des valeurs très élevées pour quantité et value, le ratio lui n'est pas affecté. 
browse if tag_multioutput_pko >0 & tag_multioutput_pko<. & otl_pko == 1 
browse if tag_multioutput_pko >0 & tag_multioutput_pko<.
*Pour être vraiment rigoureux il faudrait faire un test de comparaison des distributions de ces 247 obs. et de l'ensemble des obs nettes d'otls. 
gen flag_multioutput = (tag_multioutput_cpo>0) if tag_multioutput_cpo<. 
ksmirnov out_ton_cpo, by(flag_multioutput) exact
ksmirnov out_val_cpo, by(flag_multioutput) exact
*/

*** Construct some cross-sectional indicators of mills cpo-related activities
** inputs
*cpo
bys firm_id: egen avg_in_tot_ton_cpo_imp2 = mean(in_tot_ton_cpo_imp2)

bys firm_id: gen templast_in_tot_ton_cpo_imp2 = in_tot_ton_cpo_imp2 if year == max_year 
bys firm_id: egen last_in_tot_ton_cpo_imp2 = total(templast_in_tot_ton_cpo_imp2)
replace last_in_tot_ton_cpo_imp2 = . if last_in_tot_ton_cpo_imp2 == 0 
drop templast_in_tot_ton_cpo_imp2

*ffb 
bys firm_id: egen avg_in_ton_ffb_imp2 = mean(in_ton_ffb_imp2)

bys firm_id: gen templast_in_ton_ffb_imp2 = in_ton_ffb_imp2 if year == max_year 
bys firm_id: egen last_in_ton_ffb_imp2 = total(templast_in_ton_ffb_imp2)
replace last_in_ton_ffb_imp2 = . if last_in_ton_ffb_imp2 == 0 
drop templast_in_ton_ffb_imp2

bys firm_id: egen avg_ffb_price_imp2 = mean(ffb_price_imp2)

** outputs
bys firm_id: egen avg_out_ton_cpo_imp2 = mean(out_ton_cpo_imp2)

bys firm_id: gen templast_out_ton_cpo_imp2 = out_ton_cpo_imp2 if year == max_year 
bys firm_id: egen last_out_ton_cpo_imp2 = total(templast_out_ton_cpo_imp2)
replace last_out_ton_cpo_imp2 = . if last_out_ton_cpo_imp2 == 0 
drop templast_out_ton_cpo_imp2

** price
bys firm_id: egen avg_cpo_price_imp1 = mean(cpo_price_imp1)
bys firm_id: egen avg_cpo_price_imp2 = mean(cpo_price_imp2)





*save "$base_path_wd\build\temp\IBS_1998_cleaned_nogeo.dta", replace


*use  "$base_path_wd\build\temp\IBS_1998_cleaned_nogeo.dta", clear

**************************************************************************************************************************************************************
***PERCENTAGE EXPORTED***
**************************************************************************************************************************************************************


***** STRUCTURE OF THE SECTION *****
	**** VARIABLE DESIGNATIONS
	**** ROUTINE CLEANING
	**** IMPUTATIONS 


***** VARIABLE DESIGNATIONS *****

/*
There are several variables at our disposition that relate with export. They contain different types of information, and their availability varies across years.
The support excel file summarizes their annual availability.
Based on that, we keep only one dummy variable and one total share variable. 

In what follows, export_dummy and export_dummy_imp refer to whether a firm exports or not, EKSPOR is the equivalent in Rothenberg.
export_pct, export_pct_imp or total/average/mean export share refer to the firm-level averaged percentage of all its products exported. 
_imp variables differ in that full duplicates in the main IBS dataset are removed from the former. Let us keep these duplicates for now (i.e. don't use _imp)
PRPREX is the same but from Rothenberg dataset.
prex refers to commodity specific export shares. 


EXPORT DUMMY
There is no information in export_dummy that there is not in EKSPOR. 
Moreover, we would rather keep all possible information and hence not remove within firm duplicates of export_dummy (i.e. we don't need _imp either)	
*/	
drop export_dummy export_dummy_imp
* EKSPOR is either 1 if there was some export, or 2 if there was NO export. 
replace EKSPOR = . if EKSPOR == 0 
replace EKSPOR = 0 if EKSPOR == 2

/*
EXPORT SHARE
The export_pct is the same as the PRPREX from 1998 to 2009 (with one exception of 1 missing in 2006 for PRPREX, that is a 0 in export_pct)
Then, in 2010-2012, all missings in PRPREX are counted as zeros in export_pct. 
Then in 2013-2015, all PRPREX are missing, and not export_pct. 
Hence, we can keep using only export_pct, with zeros being backed by EKSPOR.
Code to see it:
gen pos_PRPREX =(PRPREX >0 & !mi(PRPREX)) if !mi(PRPREX)
gen pos_export_pct =(export_pct >0 & !mi(export_pct)) if !mi(export_pct)
gen pos_prex_cpo =(prex_cpo >0 & !mi(prex_cpo)) if !mi(prex_cpo)
gen pos_prex_pko =(prex_pko >0 & !mi(prex_pko)) if !mi(prex_pko)
tabulate year pos_PRPREX , missing
tabulate year pos_export_pct , missing 
tabulate year pos_prex_cpo, missing 
tabulate year pos_prex_pko, missing 
*/
drop PRPREX export_pct_imp

gen export_pct_imp = export_pct
order export_pct_imp, after(export_pct)


/*
COMMODITY EXPORT SHARES 
imp1 and imp2 are not related to qty and value imp1 and imp2 above. But the same logic applies: imp1 is subject to stronger imputations than imp2, i.e. imp1_prex_vars take value on 
a more modified sample. 
The rationale is: imp1 makes strong imputations in the aim to get as many informative observations as possible. 
imp2 makes imputations in the aim to get straightforward available information. 
no imp is the least imputed version, necessary to see how raw data look like.
*/
rename avg_out_prex* prex*
global prex_vars prex_cpo prex_pko prex_rpo prex_rpko

foreach var of varlist $prex_vars{
	gen double `var'_imp1 = `var'
	gen double `var'_imp2 = `var'
	order `var'_imp2, after(`var')
	order `var'_imp1, after(`var')
}



**** ROUTINE CLEANING ****

*Some cleaning that turn some obs. down to missing or rescale (but does not change zeros nor missings to anythin)

	*For some reason some export_pct are higher than 100. It is not clear whether they should be replaced with 100 or divided by 10 (or 100), so remove them. 
	replace export_pct_imp = . if export_pct >100 | export_pct <0 

	*Also, because of the avg over multioutput, some prex are slightly higher than 100.
	foreach var of varlist $prex_vars{
		replace `var'_imp1 = 100 if `var'> 100 & `var' <. 
		replace `var'_imp2 = 100 if `var'> 100 & `var' <. 
	}
	foreach var of varlist $prex_vars{
		replace `var'_imp1 = . if  `var' <0 
		replace `var'_imp2 = . if  `var' <0 
	}

	*Also, this is not so likely that firms reported export shares of less than 1%, it is thus likely that these answers were expressed in fractional form. 
	replace export_pct_imp = export_pct_imp*100 if export_pct < 1 & export_pct > 0 
	foreach var of varlist $prex_vars{
		replace `var'_imp1 = `var'*100 if `var' <1 & `var' >0 
		replace `var'_imp2 = `var'*100 if `var' <1 & `var' >0 
	}
/*
	For export shares exactly equal to one, this is a bit less straightforward, because this could mean 100%, or 1 in the EKSPOR answer which means some export, or just 1% is exported.
	export shares equal to 2 could also be answers reported in the EKSPOR coding, meaning "no export". 
	So let us not change impute anything.
					foreach var of varlist $prex_vars {
						count if `var' == 1
						gen `var'_eq1 = (`var'==1)
						sort firm_id year
						by firm_id : egen any_`var'_eq1 = max(`var'_eq1)
					}
					browse prex_cpo if any_prex_cpo_eq1==1
*/
		

*browse firm_id year industry_code EKSPOR export_pct export_pct_imp prex_cpo prex_cpo_imp1 prex_cpo_imp2 if EKSPOR == 0 & export_pct == 0  & prex_cpo > 0 & !mi(prex_cpo)


**** IMPUTATIONS ****

/*
MAIN PROBLEMS:
	
	ZEROS: The problem is raised empirically and theoretically. 
			Theoretically, there are two cases one can think of: 
			* It is likely that some actual exporters report null export shares for they really believe their outputs are not exported 
			because they sell them to domestic actors and have no downstream visibility. (Eventhough the questionaire asks "by own or others") 
			* or that those same firms, because of this uncertainty, leave a blank which might be in turn assigned a zero throughout the data collection and processing. 
				those latter zeros are called * SPURIOUS ZEROS * in what follows. 
			Empirically: 
			The number of zeros in all variables in all (available) years is striking at a first glance. 
			The annual effective (i.e. weighted by output) export shares are much below official national aggregates. 

	MISSINGS: Two distinct causes. 
			First, some variables have just not been delivered at all some years:
				EKSPOR: all obs. missing in 1998, 2013-2015
		        export_pct: all obs. missing in 2001-2003, 2005, 2007
		        prex: all obs. missing in 2006 and many missing in 2013-2015 
		    Second, for the remaining years, each variable still has missings. The reason for those may be numerous. One possible explanation is that some firms 
		    just did not answer those export questions because they did not know whether or to what extent their outputs were exported. 

TIME VARYING MISMEASUREMENT. 
		/*\ These causes being distinct is important, because it means that missing observations should not be looked at identically, depending on the variable and the year. \*/
			Therefore, the imputation process should account for the nature of the missing observation. 
			I.e. the imputations need to be conditionnal on available variables, and not homogeneously on both the dummy and the total share.   
			In particular, the years 2013-2015 heavily rely on imputations because they have very few non missing prex observations. They also have no available dummy. 
			Thus, conditionning on the dummy would result in heavily misleading imputations for these years. 
			 


MAIN RATIONALES HELPING TO HANDLE THESE PROBLEMS.
The availability of three different types of variables (dummy, total share and commodity share) makes it possible to identify some misreportings and impute some other values. 
The main rationales are (broadly here, more in details within the code below): 
Positive export shares are more credible information than the dummy, as they require an active reporting effort while the latter may be subject to 
ticking box mistakes. 
The dummy may adress the problem of missings, because it can clame for no export activity at all. It is easy to imagine that a firm that knows it does not 
export any of its output may report that in the dummy (which comes first in the questionaires) and then just leave blanks in the following questions on shares. 
The dummy may also help confirming that null export shares are actual believes of the respondents that their outputs weren't exported, rather than spurious zeros.


SUMMARIZED IMPUTATION PROCESS:

*** Turn missing export shares to zeros.
	- export_pct_imp: when the dummy is null.
	- imp1: if at least one of the dummy or the total export share is null, and the other one does not contradict.
	- imp2: idem, but in years when both variables are available, both must report null export.  

*** Turn null export shares to missing
	- export_pct_imp: when the dummy is positive.
	- imp1: when the dummy discredits it, or the total export share discredits it, or both are missing. 
	- imp2: never. 

*** Turn missing or null commodity export shares to positive. 

  ** = 100% for all commodities' imp1 and imp2: when the total export share is 100%.

  ** For CPO
	- imp1 and imp2 = export_pct_imp if CPO is the only commodity produced.
	- imp1 and imp2 = (total export - other commo export)/CPO output when only one other commodity is produced and it has a known export share. 
	- imp1 only = export_pct_imp in remaining cases, i.e. when imp1 is still missing because other commodities are produced but their export shares are unknown. 

  ** For PKO
  	idem. 

  ** Regarding RPO and RPKO, these are firms that are more likely to export other products that are not under our radar, because they are bigger factories and because there are many different
	 processed palm oil products. 
	 So the imputation is likely to be even less robust, so we don't do it at all. 

*** FILL GAPS 
	- imp1 = average of previous and next year when missing.
	- imp2 is not changed. 

****************************************************

One of the reasons behind the remaining spread between our averages and national ones (collected at traders in harbours directly) is respondents underestimating their 
export shares because they just don't know for some part of their productions. 

*/


*** TURN MISSINGS TO ZEROS
/*
It is pretty straightforward to impute that if a firm reports to have no export activity at all, then its total and commodity export shares are null. Also, it is flexible regarding which 
variable is available and hence minimizes the dependency to variables'availability.
However, imp2 is more conservative and requires more evidence of null export,  This limits the risk of a spurious zero at the level of total export share to change a missing into a 
new spurious zero at the commodity level.
Note that this implies that imp2 makes no such imputations in years where either the dummy or the total share is missing. Therefore the imputations are more dependent on 
variables' availability. 
Knowing that in 2013-2015 in particular there is almost no null commodity export shares, this choice will result in an underestimated number of null prex in these years. 
*/

replace export_pct_imp = 0 if mi(export_pct_imp) & EKSPOR == 0 
* (condition on export_pct_imp rather export_pct because we prefer to trust EKSPOR where export_pct was > 100 or negative, and EKSPOR is null)


foreach var of varlist $prex_vars{
	replace `var'_imp1 = 0 if mi(`var') & ((EKSPOR == 0 & (export_pct == 0 | mi(export_pct))) | (mi(EKSPOR) & export_pct == 0))
	replace `var'_imp2 = 0 if mi(`var') & (EKSPOR == 0 & export_pct == 0) | ///
										  (mi(EKSPOR) & export_pct == 0 & (year == 1998 | year > 2012)) | ///
										  (EKSPOR == 0 & mi(export_pct) & (year == 2001 | year == 2002 | year == 2003 | year == 2005 | year == 2007))

}
* (This excludes the case where both variables are missing, in which case we have no reason to think that prex is zero rather than missing.) 

* for CPO this turns 3353 missing obs. to ZERO in imp1 and 3013 in imp2, THIS IS IMPORTANT. 


*** TURN ZEROS TO MISSINGS
/* 
When the dummy indicates export activity, but the total export share is reported a zero, there might be two scenarii: 
- either the dummy reflects reality, and the export share is a spurious zero
- or the export share reflects reality and the dummy has been mistakenly ticked. 
We choose to trust more in the former scenario, i.e. that mistakenly ticking a box is less likely than spurious zero cases.
We also remove null total shares that are not credible in light of positive commodity shares. 
The consequence of that is that those 120 likely spurious zeros won't be available for later imputation on commodity export shares based on total.   
*/
replace export_pct_imp = . if export_pct == 0 & ///
								(EKSPOR == 1 | ///
								(prex_cpo > 0 & !mi(prex_cpo)) | ///
								(prex_pko > 0 & !mi(prex_pko)) | /// 
								(prex_rpo > 0 & !mi(prex_rpo)) | /// 
								(prex_rpko > 0 & !mi(prex_rpko)))

* we could impute total export shares from commodity export shares in different cases, but as the total export share is not planned to be used in the analysis, we don't. 

/* 
When the commodity export share is null but there are signs of firm level export because of the dummy and/or the total share, it can be either that: 
- the truly positive share of commodity exported was not reported by the respondent and this ended in a spurious zero.
- null prex reflects reality: there was some firm level export but this commodity was indeed not exported at all.  
 
We trust the former more in imp1 and the latter more in imp2, i.e. we discredit null prex_ variables only in imp1 and not in imp2. 
BUT, some imp2 that remain 0 here will actually be turned to positive values below, when they are in obvious cases of 
being equal to the total share (either because the total share is 100% or because there is only one commo exported, but not when there is another commo exported 
because it is precisely less clear whether the null commodity export share is true or not then).

Also, in imp1, null prex are deemed spurious zeros if both the dummy and total export share are missing. 
*/
foreach var of varlist $prex_vars{
	replace `var'_imp1 = . if `var' == 0 & ((EKSPOR == 1) | (export_pct > 0 & !mi(export_pct)) | (mi(EKSPOR) & mi(export_pct)))
}
* For CPO this is 803 zeros that are turned to missing. 





*** TURN MISSING OR NULL TO POSITIVE ***

/*
we use the already imputed version of export_pct, i.e. export_pct_imp so that imputations on prex_vars account for the prior routine cleaning.
The imputations made on export_pct_imp, from missing to zeros and from zeros to missing, won't affect the following imputations, as they only use positive 
total export shares that were affected by neither of them.  

For CPO, 
Among obs. that have prex_cpo missing and export_pct non missing, we can make increasingly strong imputations:

	1. If the total export share is 100%, all the commodity shares are necessarily 100% too (for those commodities that are produced). 

	2. Suppose you produce only CPO. Then the total export share is equal to the commodity export share. replace prex_cpo_imp1 and prex_cpo_imp2 with export_pct if prex_cpo <. 

	3. If you produce CPO and other commodities - FOR WHICH PREX_ IS KNOWN - we can impute prex_cpo rather precisely (not exactly because other products than those 4 may be exported). 
		There will be a different formula depending on whether the firm produces only one other commo, 2, or 3 other commodities.  
	
	4. If you produce CPO and other commodities - FOR WHICH PREX_ OR NECESSARY INFO IS MISSING (so that we cannot compute prex_cpo) - there are two extreme options: 
		a) we use the firm level export share as a proxy for the commodity level export share, which is the same as assuming that all commodities produced 
			were exported in the same proportions (once weighted with the amounts produced), this is _imp1. 
		b) we are conservative and don't impute more. 
		 
	These options don't depend on whether the firm produces 1, 2, or 3 other commodities. 

1. and 2. make the same modifications to imp1 and imp2. 
3. modifies a subgroup of missing prex when export_pct is not missing. It makes the same imputations on imp1 and on imp2, but this does not affect both equally. 
	This is because it also affects the initially null prex that have first been turned to missing in imp1, but not in imp2 because these may actually mean that the respondents
	knew one thing for sure: that this commodity was not exported, while still exporting an other product.  
4. modifies an other subgroup of missing prex when export_pct is not missing (those that missed some information to enter the imputation in 3.)
Finally: 
THE DIFFERENCE BETWEEN _IMP1 AND _IMP2 IS THAT _IMP2 IS NOT IMPUTED TO BE EQUAL TO THE TOTAL SHARE UNLESS 
IT IS 100; OR THERE IS ONLY ONE COMMODITY; OR IT IS MISSING (AND NOT NULL) AND THERE IS ENOUGH INFORMATION ON THE OTHER EXPORTED COMMODITY.

Then repeat this for PKO. 

*/

** For all commodities
		* 1. 
		global commo cpo pko rpo rpko 

		foreach commo of global commo {
			forvalues a = 1/2{
				replace prex_`commo'_imp`a' = 100 if (mi(prex_`commo') | prex_`commo' == 0) & /// replaces those missing but also those with a null commodity export share. 
				(valid_out_ton_`commo'==1 | valid_out_val_`commo'==1) & /// 
				export_pct_imp == 100
			}	
		} 
		*This yields 130 real changes. 

	**For CPO**

		global imp_prex_cpo prex_cpo_imp1 prex_cpo_imp2
		
		* 2.

		foreach var of varlist $imp_prex_cpo {
			replace `var' = export_pct_imp if ( /// 
			(mi(`var') | `var' == 0)  /// replaces those missing but also those with a null commodity export share.  
			& !mi(export_pct_imp) & export_pct_imp > 0 /// AND export_pct_imp > 0, otherwise missing prex that have not been replaced with 0s in imp2 supra would get replaced here. 
			& (valid_out_ton_cpo==1 | valid_out_val_cpo==1) ///
			& valid_out_ton_pko == 0 & valid_out_val_pko == 0 ///
			& valid_out_ton_rpo == 0 & valid_out_val_rpo == 0 ///
			& valid_out_ton_rpko == 0 & valid_out_val_rpko == 0 ///
				)
		}

		/* 
		This yields 127 real changes. 
		For cpo vars we don't use imp variables just for the sake of code writing simplicity. But in the end the duplicates and outliers lines of cpo vars will not be taken into account anyway. 
		We condition on validity of qty OR value so that the aggregates on export shares can be calculated on as many information as possible. 
		We condition on invalidity of both qty and value at non-imp state because we want to be sure that there is no sign of the production of something else than CPO. 
		For instance, suppose you have a positive cpo qty (and you want your prex_cpo to be replaced by export_pct only if you produced only cpo) and an outlying pko qty. 
		Then it means you probably also produced pko. But if one looks at you out_ton_pko_imp vars, you look like you don't produce pko. 
		Or, your raw pko info is duplicated from your previous for some reason. It is likely that you are a producer of pko, but this might not appear in your imp vars. 
		
		* 3. 
		The choice of val and imp2 is because this is where there are the least missings, and also not the outliers. Repeating the command with ton and imp2 makes no new change. 
		
		*/	
		foreach var of varlist $imp_prex_cpo {
			global second_commo pko rpo rpko 
			foreach commo of global second_commo {
				local all_unproduced_commo : list global(second_commo) - local(commo)
				local unproduced_commo1 : word 1 of `all_unproduced_commo'
				local unproduced_commo2 : word 2 of `all_unproduced_commo'
				replace `var' = (export_pct_imp*(out_val_cpo_imp2 + out_val_`commo'_imp2) - (prex_`commo'*out_val_`commo'_imp2))/out_val_cpo_imp2 if /// 
				mi(`var')  /// prex_cpo_imp is missing, but not null because this would make these imputations also on null imp2 that have been deemed as true zero in step *** turning zeros to missings ***
				& export_pct_imp <. & export_pct_imp > 0 /// and it is available in export_pct_imp and positive - i.e. it maks sense to make imputation from it. 
				& prex_`commo' <. 								/// and the information on the prex_ of the other commo is available. 
				& out_val_cpo_imp2 <. 							 /// and there is a sign of production in cpo 
				& out_val_`commo'_imp2<.						 /// and there is a sign of production in the second commo 
				& out_val_`unproduced_commo1'_imp2 >=. & out_val_`unproduced_commo2'_imp2>=. 
				* and there is no sign of production for the two other commos. 	
				*di "`unproduced_commo2'"
			}
		}

				* As a sign that this is not so accurate (may be there were other products exported?) this gives several prex_cpo_imp > 100. 
				replace prex_cpo_imp1 = 100 if prex_cpo_imp1> 100 & prex_cpo_imp1 <. 
				replace prex_cpo_imp2 = 100 if prex_cpo_imp2> 100 & prex_cpo_imp2 <. 	
				*Just for generality: 
				replace prex_cpo_imp1 = 0 if prex_cpo_imp1<0  
				replace prex_cpo_imp2 = 0 if prex_cpo_imp2<0  	
		/*
		Among these obs. that have a sign of cpo and pko and non missing prex_pko, there is no obs that has also sign of a production for other commo, so there is no need to go for thoses cases. 


		* 4. 
		Proxy the cpo export share with the total export share in remaining cases (only for imp1)	
		The conditioning on missing is on prex_cpo_imp1 and not just prex_cpo in order to make the imputation also for obs. that were switched from 0 to missing because there was some sign of export. 
		*/
			replace prex_cpo_imp1 = export_pct_imp if mi(prex_cpo_imp1) & !mi(export_pct_imp) & /// 
				(industry_code==15141 | industry_code== 10431 | industry_code== 31151) & ///
				(valid_out_ton_cpo==1 | valid_out_val_cpo==1) 	 & ///
				(valid_out_ton_pko == 1 | valid_out_val_pko == 1 | ///
				valid_out_ton_rpo == 1 | valid_out_val_rpo == 1  | ///
				valid_out_ton_rpko == 1 | valid_out_val_rpko == 1)	& /// 
				(prex_pko>=. | out_val_pko_imp2 >=.) & (prex_rpo>=. | out_val_rpo_imp2>=.) & (prex_rpko>=. | out_val_rpko_imp2>=.)	 
		/*
		80 real changes made (non zeros). 
		So here with this, we replace prex_cpo_imp1 with the export_pct in the cases that are not treated above, it is to say: 
		this is not case 1. because there is some sign of another commodity being produced 
		this is not case 2. because the prex_vars are missing. (last line)
		Also, we restrict this to firms that have a crude vegetable oil industry_code in order to limit the risk of multiple commodities and hence commo export share different than 
		total export share. 
		*/


	**For PKO**

		* 2.

		global imp_prex_pko prex_pko_imp1 prex_pko_imp2

		foreach var of varlist $imp_prex_pko {
			replace `var' = export_pct_imp if ( /// 
			(mi(`var') | `var' == 0)  /// replaces those missing but also those with a null commodity export share.  
			& !mi(export_pct_imp) & export_pct_imp > 0 ///
			& (valid_out_ton_pko==1 | valid_out_val_pko==1) ///
			& valid_out_ton_cpo == 0 & valid_out_val_cpo == 0 ///
			& valid_out_ton_rpo == 0 & valid_out_val_rpo == 0 ///
			& valid_out_ton_rpko == 0 & valid_out_val_rpko == 0 ///
				)
		}
		/* 
		

		* 3. 

		*/		
		foreach var of varlist $imp_prex_pko {
			global second_commo cpo rpo rpko 
			foreach commo of global second_commo {
				local all_unproduced_commo : list global(second_commo) - local(commo)
				local unproduced_commo1 : word 1 of `all_unproduced_commo'
				local unproduced_commo2 : word 2 of `all_unproduced_commo'
				replace `var' = (export_pct_imp*(out_val_pko_imp2 + out_val_`commo'_imp2) - (prex_`commo'*out_val_`commo'_imp2))/out_val_pko_imp2 if /// 
				mi(`var') /// prex_pko is missing
				& export_pct_imp <. & export_pct_imp > 0 ///  and it is available in export_pct_imp 
				& prex_`commo' <. 								/// and the information on the prex_ of the other commo is available. 
				& out_val_pko_imp2 <. 							 /// and there is a sign of production in pko 
				& out_val_`commo'_imp2<.						 /// and there is a sign of production in the second commo 
				& out_val_`unproduced_commo1'_imp2 >=. & out_val_`unproduced_commo2'_imp2>=. 
				* and there is no sign of production for the two other commos. 	
				*di "`commo'"
			}
		}

		/*  425 real changes made when CPO is the other commo.
			Several of these obs. are negative on prex_pko_imp, showing that these cases are not so exact imputation. There were other products weighing in the firm level export_pct
			or this is the sign of some other measurement error. 
		*/
				replace prex_pko_imp1 = 100 if prex_pko_imp1> 100 & prex_pko_imp1 <. 
				replace prex_pko_imp2 = 100 if prex_pko_imp2> 100 & prex_pko_imp2 <. 
				replace prex_pko_imp1 = 0 if prex_pko_imp1<0 
				replace prex_pko_imp2 = 0 if prex_pko_imp2<0 
		

		/*Among these obs. that have a sign of pko and cpo and non missing prex_cpo or prex_rpo, there is no obs that has also sign of a production for other commo, 
		so there is no need to go for thoses cases. 


		* 4. 

		Old code
		// First transform any remaining prex_pko missing with 0 if total export_pct is null or dummy is 2 (no export at all)
		forvalues a = 1/2 {
			replace prex_pko_imp`a' = 0 if prex_pko_imp`a' >=. & export_pct_imp == 0 | export_dummy_imp==2
		}
		*/

			replace prex_pko_imp1 = export_pct_imp if mi(prex_pko_imp1) & !mi(export_pct_imp) & ///
				(industry_code==15141 | industry_code== 10431 | industry_code== 31151) & ///
				prex_cpo>=. & prex_rpo>=. & prex_rpko>=. 		 & /// 
				(valid_out_ton_pko==1 | valid_out_val_pko==1) 	 & ///
				(valid_out_ton_cpo == 1 | valid_out_val_cpo == 1 | ///
				valid_out_ton_rpo == 1 | valid_out_val_rpo == 1  | ///
				valid_out_ton_rpko == 1 | valid_out_val_rpko == 1)	& /// 
				(prex_cpo>=. | out_val_cpo_imp2 >=.) & (prex_rpo>=. | out_val_rpo_imp2>=.) & (prex_rpko>=. | out_val_rpko_imp2>=.)	



*** FILL GAPS 
* Only for imp1 and only when the gap is a missing, not when it's a zero (which is only 3 and 4 cases for cpo and pko resp. anyways).
global commodity cpo pko
foreach commo of global commodity{
	*forvalues a = 1/2{
		local a 1 
		sort firm_id year
		by firm_id: gen double lag_prex_`commo'_imp`a' = prex_`commo'_imp`a'[_n-1]
		by firm_id: gen double lead_prex_`commo'_imp`a' = prex_`commo'_imp`a'[_n+1]

		*gen flag_0prex_`commo'_imp`a' = (prex_`commo'_imp`a'==0 & lag_prex_`commo'_imp`a' >0 & lag_prex_`commo'_imp`a'<. & lead_prex_`commo'_imp`a' >0 & lead_prex_`commo'_imp`a'<.)
		*replace prex_`commo'_imp`a' = lead_prex_`commo'_imp`a' if flag_0prex_`commo'_imp`a'==1 & lead_prex_`commo'_imp`a'==lag_prex_`commo'_imp`a'

		gen flag_missprex_`commo'_imp`a' = (prex_`commo'_imp`a'>=.& lag_prex_`commo'_imp`a'<. & lead_prex_`commo'_imp`a'<.)
		replace prex_`commo'_imp`a' = (lag_prex_`commo'_imp`a' + lead_prex_`commo'_imp`a')/2 if flag_missprex_`commo'_imp`a'==1 
		*count if flag_missprex_`commo'_imp`a'==1 
		*count if flag_missprex_`commo'_imp`a'==1 & lag_prex_`commo'_imp`a' == lead_prex_`commo'_imp`a'
		drop lag_prex_`commo'_imp`a' lead_prex_`commo'_imp`a' 
	*}
}	
*there are no cases of gaps of 2 years. flag_missprex_`commo'_imp`a'







/* DRAFT NOTES ON THIS SECTION: 
			However, regarding the diversity of cases, this would mean applying a variety of different rationales to annual imputations, making the tractability of 
			those choices too intricate. 
			We rather choose to apply homogeneous data processing over years acknowledge that some structural flaws in available data induce time varying noise.

		Let's just summarize what this section did regarding null prex and null export_pct. 
		- It replaced any missing prex with export_pct when the latter was null, following the argument that a firm may have not informed the prex for each commodity when it is 
		exporting nothing anyway, and just responded 0 for the total export share. 
		- It never replaced prex with export_pct when the former was positive and the latter null, assuming there is no reason for a firm to respond a positive value on a commo line if the export is null, 
		so the error has to come from export_pct... 
		- It never replaced prex with export_pct when the former was null and the latter positive bc it could be positive for other commo but null for one in particular and hence positive on average. 		

		// PRPREX
			g PR = (PRPREX > 0) if !mi(PRPREX)
			tabulate year PR, missing
			/*What one can see is that there are years when there are no missings (and #zeros might be more than the actual instances of null export.) 
			That's 1998-2009 bar the 5 years when full variable is missing. 
			In 2010, 2011, 2012, that's the other way round - no zero and only positive or missing values. 
			If we use these, then it is equivalent in terms of imputation strength: we wuold have to assume that all missings are null exports. 
			*/

		// export_pct

			g ep = (export_pct > 0) if !mi(export_pct)
			tabulate year ep, missing
			/*
			Same: only zeros and no missings (beside positive values). With a small exception in 2015 with 60 missing values. 
			And the 5 years of missing still. 
			*/

	*combien sont les proportions positives alors que la dummy dit qu'il n'y avait pas d'export
	count if EKSPOR == 0 & export_pct >0 & !mi(export_pct)  
	browse firm_id year export_pct if EKSPOR == 0 & export_pct >0 & !mi(export_pct
		*pour ceux là il y a beaucoup de cas où prex est aussi nul ou missing. Il

		browse firm_id year export_pct prex_cpo if cpo == 1 & EKSPOR == 0

			gen export_pct_conf = export_pct
			replace export_pct_conf = . if export_pct == 0 & EKSPOR != 

			g  cpo = (prex_cpo > 0) if !mi(prex_cpo)
			tabulate year cpo, missing

			g  pko = (prex_pko > 0) if !mi(prex_pko)
			tabulate year pko, missing
*/






************************************************************************************************************
*** ADD AND CLEAN GEOGRAPHIC VARIABLES ***
************************************************************************************************************

**** Attribute as much geographic information as possible.  
sort firm_id year 

merge 1:1 firm_id year using "input_data/IBS_kraus/IBS_base_prelu.dta", generate(merge_base_prelu) keepusing(island_number island_name island_factor province district /// 
	desa_code desa_name kab_code kab_name kec_code kec_name prov_code prov_name desa_id district_id desa_code_2000 district_code_1993) update
* all master's are matched because IBS_base_prelu is on 1998-2015 although desa_id, like other variables, is available only until 2010. 
* 1220 obs. that were missing on either province or district have been updated (0 nonmissing conflict). 

*remove obs. from using that have not matched (i.e. IBS firms that are no mills)
drop if merge_base_prelu == 2
drop merge_base_prelu

* We don't add desa_id and desa_code_2000 for 368 and 367 obs. resp. from IBS_desa2000 (other variables are not updated by IBS_desa2000) 
* because these are obs already removed purposedly when desa_fl (string_fl) ==1. 

*** ISLANDS 

replace island_name = "Sumatra" if island_name == "Riau"
replace island_factor = 8 if island_factor == 6

* maxmode handles missings (if there is the same amount of missing as of one island, then the mode is the island string.)
bys firm_id (year): egen island_mode = mode(island_name), missing maxmode

replace island_name = island_mode 

drop island_mode

*** DISTRICTS 
** For 2000-2010 use district_id already available. For other years, build it. 
tostring province, generate(province_str)
tostring district, generate(district_str)
*apparently some dots are not converted to missings, (problem does not occur for kab_code) so: 
replace district_str = "" if district_str == "." 
replace province_str = "" if province_str == "." 
replace province_str = "0" + province_str if province < 10 
replace district_str = "0" + district_str if district < 10 
replace province_str = prov_code if !mi(prov_code) & mi(province_str)
replace district_str = kab_code if !mi(kab_code) & mi(district_str)

egen prov_distr = concat(province_str district_str)
replace prov_distr = "" if mi(province) | mi(district_str)
replace prov_distr = district_id if year > 1999 & year < 2011 
replace prov_distr = "" if length(prov_distr) < 4
/* About the difference between district_id - concatenation of prov_code and kec_code - and prov_distr, concatenation of district and province.
91 cases (+7 where it's just that district_id has not the right length), most of the times it is the first obs of a mill, and the following codes are 
the same as the district, not the same as the kab_code. This suggests that the prov_distr is more accurate, and that we should not remove it in these cases, which will 
be done though, because the district mode_fl is based on district_id crosswalked. But not so important because in these cases we necessarily have information on the following 
years (and we need only one correct). 
gen prov_distr_di =((district_str!=kab_code &!mi(kab_code)) | (province_str != prov_code &!mi(prov_code)))
bys firm_id (year): gen any_di = 1 if sum(prov_distr_di) >0
*/
** Add names
egen year_prov_distr = concat(year prov_distr) if !mi(prov_distr)
merge m:1 year_prov_distr using "temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta", generate( /// 
	merge_district_names) keepusing(name_)
drop if merge_district_names == 2 
drop year_prov_distr
rename name_ district_name  
replace district_name = "*" + kab_name if mi(district_name) & !mi(kab_name)

** Clean 
* remove likely wrong district ids and names 
bys firm_id (year): egen district_code_1993_mode = mode(district_code_1993)
g district_mode_fl = (district_code_1993 != district_code_1993_mode & !mi(district_code_1993) & !mi(district_code_1993_mode))
replace prov_distr = "" if district_mode_fl==1
replace district_name = "" if district_mode_fl==1
drop district_code_1993_mode 
/* Don't use firm_district_1993_fl_md because apparently it is a flag on firms that have more than 2 instances of collapsed 2000 district code and _mode 
(from district_mode.do from Sebi). And we want to allow for long periods or repeated years with un-actualised district codes after a split. 
*/

* some don't match even though they have a deemed valid prov_distr. Let's flag them... 
g prov_distr_suspect_fl = 1 if merge_district_names == 1 & !mi(prov_distr)
drop merge_district_names
/* most of them have inherited a name from kab_name 
browse prov_distr district_id if prov_distr_fl_suspect == 1
*/

*** VILLAGES
** Use desa_id imported from IBS_base_prelu
egen year_desa_id = concat(year desa_id) if !mi(desa_id)

** Add names
merge m:1 year_desa_id using  "temp_data/processed_indonesia_spatial/desa_code_names_98_2014.dta", generate( /// 
	merge_village_names) keepusing(nm)
drop if merge_village_names == 2 
rename nm village_name  
replace village_name = "*" + desa_name if mi(village_name) & !mi(desa_name)
count if merge_village_names ==1 & year > 1999 & year < 2011


** Clean 
bys firm_id (year): egen desa_code_2000_mode = mode(desa_code_2000)
g village_mode_fl = 1 if desa_code_2000 != desa_code_2000_mode & !mi(desa_code_2000) & !mi(desa_code_2000_mode)

replace desa_id = "" if village_mode_fl==1
replace village_name = "" if village_mode_fl==1
drop desa_code_2000_mode  
/* the string mismatch flag has already been taken into account by using desa_id from IBS_base_prelu from the beginning rather than desa_id from IBS_desa2000
*/

* some don't match even though they have a deemed valid desa_id. Let's flag them... 
g desa_id_suspect_fl = 1 if merge_village_names == 1 & !mi(desa_id)
drop merge_village_names
 
*774 obs (in 2000-2010) don't match a village code from crosswalk eventhough they have a desa_id
* among them 557 are replaced by desa_name 



order desa_id, after(industry_code)
sort firm_id year 



save "temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta", replace




/*


********* EXCEL EXPORT FOR ROBERT ***********************
export excel firm_id year industry_code desa_id ///
 ffb_price_imp1 ffb_price_imp2 in_ton_ffb_imp1 in_ton_ffb_imp2 in_val_ffb_imp1 in_val_ffb_imp2 ///
 in_tot_cpo_price_imp1 in_tot_cpo_price_imp2 in_tot_ton_cpo_imp1 in_tot_ton_cpo_imp2 in_tot_val_cpo_imp1 in_tot_val_cpo_imp2  ///
 cpo_price_imp1 cpo_price_imp2 out_ton_cpo_imp1 out_ton_cpo_imp2 out_val_cpo_imp1 out_val_cpo_imp2 prex_cpo_imp1 prex_cpo_imp2 ///
 pko_price_imp1 pko_price_imp2 out_ton_pko_imp1 out_ton_pko_imp2 out_val_pko_imp1 out_val_pko_imp2 prex_pko_imp1 prex_pko_imp2 ///
 out_ton_rpo_imp1 out_ton_rpo_imp2 out_val_rpo_imp1 out_val_rpo_imp2 prex_rpo_imp1 prex_rpo_imp2 ///
 out_ton_rpko_imp1 out_ton_rpko_imp2 out_val_rpko_imp1 out_val_rpko_imp2 prex_rpko_imp1 prex_rpko_imp2 ///
 province district ///
 pct_own_cent_gov_imp pct_own_loc_gov_imp pct_own_nat_priv_imp pct_own_for_imp  ///
 export_pct_imp ///
 revenue_total_imp1 revenue_total_imp2 materials_tot_imp1 materials_tot_imp2 revenue_total_imp3 ///
 workers_total_imp3 workers_prod_imp3 workers_other_imp3 ///
 value_added_self_imp1 value_added_self_imp2 value_added_self ///
 using "$base_path_wd\build\input\IBS_1998_cleaned.xls", firstrow(variables) replace

saveold C:/Users/guyv/desktop/IBS_1998_cleaned.dta, version(12) replace 

*/

***CLEANING IN_TOT_QTY_TON***

*Bunch weight
*23-27 kg
*Fruit/bunch
*60-65 %
*Oil/bunch
*21-23 %
*Kernel/bunch
*5-7 %
*Mesocarp/bunch
*44-46 %
*Mesocarp/fruit
*71-76 %
*Kernel/fruit
*21-22
*Shell/fruit
*10-11
*However, such high yields are rarely achieved in practice because climatic conditions are usually less than ideal. Rainfall is erratic in Central and West Africa and hence the tree suffer water-related stresses. 


*Malaysia: 33.32% the oil in bunch content in ripe bunch

*from http://jopr.mpob.gov.my/wp-content/uploads/2013/07/joprv11n1-p21.pdf 
*--> Oil/bunch (of mean weight of 13.8kg) btw 13% and 43% with mean 25%
*for Kernel it is rather 7%

*https://www.researchgate.net/publication/319298126_Characteristics_of_Fresh_Fruit_Bunch_Yield_and_the_Physicochemical_Qualities_of_Palm_Oil_during_Storage_in_North_Sumatra_Indonesia
 *the  ratio  of  fruit  to  bunch(67.7%), bref voir le screenshot. 
