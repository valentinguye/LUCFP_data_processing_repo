*¨Prepare most recent version of UML. This version has establishment year variable, est_year, taking values earlier than 2004. 
* Such values come from different sources: "We have found some establishment dates going back to the 80s and up to 2019 based on satellite imagery, 
* journal articles, company reports, mill installation contractor websites, and government websites"

import excel  "input_data/uml/mills_estyear_clean.xlsx", firstrow clear 
rename latitude lat
rename longitude lon

rename mill_id trase_code

* correct the one duplicate case on coordinates with information from mills_20200129.xlsx
replace lat = 1.078	 if uml_id == "PO1000004435"
replace lon = 100.261 if uml_id == "PO1000004435"
duplicates list lat lon 

* add the guys that are not here but were in traseMills_capEstyear.xlsx though.  
merge 1:1 trase_code using "temp_data/processed_UML/traseMills_capEstyear_modified.dta", nogenerate
* 13 mills were not in mills_estyear_clean.xlsx but are in traseMills_capEstyear_selected.dta

* round coordinates
replace lat = round(lat, 0.001) 
replace lon = round(lon, 0.001) 
duplicates list lat lon 


save "temp_data/processed_UML/mills_estyear_clean_modified.dta", replace 

// use "$base_path_wd\build\input\mills_estyear_clean.dta", clear

/* use "$base_path_wd\build\input\mill_geolocalization\traseMills_capEstyear_selected.dta", clear
browse if trase_code == "M-01200"

parent_co	mill_name	lat	lon
PT Pelita Agung Agrindustri	PT. Pelita Agung Agro Industri	1.427	101.189
 */


* read most recent version of UML, that does not have est_year variable.
use "temp_data/processed_UML/mills_20200129_modified.dta", clear 

* add the est_year variable, and the 10 mills that are in mills_estyear_clean.dta (either from mills_estyear_clean.xlsx or from traseMills_capEstyear_selected.dta)
merge 1:1 trase_code using "temp_data/processed_UML/mills_estyear_clean_modified.dta", nogenerate

* add the min_year variable from IBS 
merge 1:1 trase_code using "temp_data/processed_mill_geolocalization/IBS_UML_cs.dta", nogenerate keepusing(min_year)

gen est_year_imp = est_year

*** MILLS WITH DOUPTFUL ESTABLISHMENT YEAR

** mills that appeared in IBS prior to UML establishment year 
count if !mi(est_year) & min_year < est_year
* 183 

** est_year = 2004 
count if est_year == 2004
* 234
count if est_year == 2004 & !mi(min_year)
* 199
count if est_year == 2004 & min_year == 2004
* 13
count if est_year == 2004 & min_year < 2004
* 99

/*  
234 UML mills have an establishment year in 2004; Out of which, 199 are matched with IBS; out of which 13 appear "indeed" in IBS in 2004 for the first time, 
99 appear before 2004, and 86 after 2004. 

183 mills have an establishment year that is questioned by their first year of appearance in IBS (min_year). 99 of these are mills with 
establishment year in 2004. 
*/

* we deem that the earliest year of est_year and min_year is the best approximation of the true establishment year.  
replace est_year_imp = min_year if !mi(est_year_imp) & min_year < est_year_imp

* this replaced the 99 est_year = 2004 cases, plus another 84 cases of mills appearing in IBS prior to UML est_year. 


*** MILLS WITHOUT ESTABLISHMENT YEAR 
count if mi(est_year) 
* 337  
count if mi(est_year) & !mi(min_year)
* 82
count if mi(est_year) & !mi(min_year) & min_year ==2004
* 4
count if mi(est_year) & !mi(min_year) & min_year > 2004
* 47

/* 
So 337 UML mills don't have an est_year. We approximate their establishment years with the first years of appearance in IBS. 
Out of 337 mills with missing establishment years, 82 can be approximated this way, of which 47 with a first year in IBS later than 2004, 4 in 2004
and 31 before 2004. 
*/

replace est_year_imp = min_year if mi(est_year_imp) & !mi(min_year)

save "temp_data/processed_UML/UML_valentin_imputed_est_year.dta", replace 
