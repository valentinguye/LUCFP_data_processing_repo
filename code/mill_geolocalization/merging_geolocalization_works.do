/*
This script merges together the spreadsheets used for manual geolocalization 
(through triangular matching between IBS, Manufacturing directories (MD) and UML.) 
Then, it makes and documents additional manual modifications to the data, in order to resolve conflicts (mill name or coordinates duplicates).

It needs: 
* input: temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta
*		 input_data/manually_matched_ibs_uml/mills_to_georeference_post2010_done.xls 	 sheet("mills to georef")
*  	     input_data/manually_matched_ibs_uml/mills_to_georeference_pre2011_done.xlsx	 sheet("Mills to georef")
*  	     input_data/manually_matched_ibs_uml/noto_done.xls 	sheet("Sheet1")
*  		 input_data/uml/traseMills_capEstyear.xlsx 	sheet("traseMills_capEstyear")
*		 input_data/manually_matched_ibs_uml/matching_unref/md_millMatching_modified.xlsx
* 		 input_data/manually_matched_ibs_uml/overall_btw_conflicts_done.xlsx


For this purpose, the first step is to prepare (in a .dta harmonized format) all 
the Excel spreadsheets containing the outputs of the manual works. 

	- oto is already in .dta from R file, it was not processed in an excel format. 
	- for noto, pre2011 and post2010 we want to check for each subset that there is only one mill_name and one coordinate value within firm_id
		(Duplicated mill_name and coordinate value across firm_id will be checked later, once they are all together.)
	- unref is in Excel cross-sectional format, with matched trase_code (named mill_id originally) and without coordinates, so this is to be fixed. 
*/



***** PREPARE NOT-ONE-TO-ONE (NOTO) MANUALLY MATCHED GROUP *****
import excel "input_data/manually_matched_ibs_uml/noto_done.xls", sheet("Sheet1") firstrow clear
destring lat, replace
destring lon, replace
* also, in case the excel sheet produces new empty lines: 
drop if mi(firm_id)
sort firm_id year 

*drop those that could not be geolocalized (remember that some obs. - in green in excel basically - have coordinates but no names
* and others have names from directories but no coordinates could be found. We want to keep both kinds.)
drop if mi(lat) & mi(mill_name)

*check that there is only one mill_name and one coordinate value within firm_id 
by firm_id mill_name, sort: gen nvals = (_n ==1)
by firm_id: replace nvals = sum(nvals)
count if nvals != 1
drop nvals
by firm_id lat, sort: gen nvals = (_n ==1)
by firm_id: replace nvals = sum(nvals)
count if nvals != 1
drop nvals

sort firm_id year 
duplicates drop firm_id, force

/* 
duplicates list mill_name
duplicates list lat lon
check duplicated mill_name and coordinate value across firm_id (4 instances here, which are the 4 cases 
of firm_id that need to be changed because they are found to be the mill as another firm_id)
*/

save "temp_data/processed_mill_geolocalization/noto_done.dta", replace 
*that's 102 mills with coordinates + 1 (AGROINDO INDAH PERKASA) with name but no coordinates. 


***** PREPARE PRE 2011 MANUALLY MATCHED GROUP *****
import excel "input_data/manually_matched_ibs_uml/mills_to_georeference_pre2011_done.xlsx", sheet("Mills to georef") firstrow clear
rename latitude lat
rename longitude lon
destring lat, replace
destring lon, replace
drop if mi(firm_id)

*drop those that could not be geolocalized (remember that some obs. - in green in excel basically - have coordinates but no names
* and others have names from directories but no coordinates could be found. We want to keep both kinds.)
drop if mi(lat) & mi(mill_name)
* it is 123 mills (485 obs.) that could not be geolocalized. 

*check that there is only one mill_name and one coordinate value within firm_id 
by firm_id mill_name, sort: gen nvals = (_n ==1)
by firm_id: replace nvals = sum(nvals)
by firm_id: replace nvals = nvals[_N]
count if nvals != 1
drop nvals
by firm_id lat, sort: gen nvals = (_n ==1)
by firm_id: replace nvals = sum(nvals)
count if nvals != 1
drop nvals

sort firm_id year 
duplicates drop firm_id, force

/* 
duplicates list mill_name
duplicates list lat lon
across firm_id there is one duplicated mill_name (BONANZA MEGAH LTD, PT) and the only duplicated coordinates are missings. 
*/

save "temp_data/processed_mill_geolocalization/mills_to_georeference_pre2011_done.dta", replace 
/*that's 38 mills with coordinates + 14 with names only ( 
SEGAR KEMBANG SEJATI, PT 
SEI SEKALA 
LEMBAH KRYA, PT 
PMK SUMBER WARAS 
PRISCOLIN, PT 
INDOSCO UTAMA PT 
BINA KARYA PRIMA, PT 
BONANZA MEGAH LTD, PT 
SUKUR, PT 
SAMPURNA SPUTNIK PT 
KETUPAT MAS,PT 
DUTA SUMBER NABATI, PT 
KURNIA TUNGGAL NUGRAHA, PT 
BONANZA MEGAH LTD, PT 
)
*/

***** PREPARE POST 2010 MANUALLY MATCHED GROUP ***** 
import excel "input_data/manually_matched_ibs_uml/mills_to_georeference_post2010_done.xlsx", sheet("mills to georef") firstrow clear
sort firm_id year
rename latitude lat
rename longitude lon
destring lat, replace
destring lon, replace
drop if mi(firm_id)

*drop those that could not be geolocalized or obs. of geolocalized mills for which the info was not written
* (remember that some obs. - in green in excel basically - have coordinates but no names
* and others have names from directories but no coordinates could be found. We want to keep both kinds.)
drop if mi(lat) & mi(mill_name)

*check that there is only one mill_name and one coordinate value within firm_id 
by firm_id mill_name, sort: gen nvals = (_n ==1)
by firm_id: replace nvals = sum(nvals)
by firm_id: replace nvals = nvals[_N]
sort firm_id year 
count if nvals != 1
drop nvals
by firm_id lat, sort: gen nvals = (_n ==1)
by firm_id: replace nvals = sum(nvals)
count if nvals != 1
drop nvals

sort firm_id year 
duplicates drop firm_id, force

save "temp_data/processed_mill_geolocalization/mills_to_georeference_post2010_done.dta", replace 
/*That is 98 mills with name and coordinates, and 4 mills with only names
TUNAS LESTARI SEJATI, 
PTSUMATERA JAYA AGRO INDUSTRI, 
PTANAM KOTO, 
PTANJ AGRI SIASIS
*/


capture mkdir "temp_data/processed_UML"


***** PREPARE UML FIRST VERSION (2019)
import excel "input_data/uml/traseMills_capEstyear.xlsx", sheet("traseMills_capEstyear") firstrow clear
keep trase_code parent_co mill_name latitude longitude est_year
rename latitude lat
rename longitude lon
destring lat, replace
destring lon, replace
* they are replaced as double. Let us round these coordinates, so that we can match them more easily in the merges that follow.
replace lat = round(lat, 0.001) 
replace lon = round(lon, 0.001) 

codebook trase_code
save "temp_data/processed_UML/traseMills_capEstyear_modified.dta", replace 

***** PREPARE UML SECOND VERSION (2020) *****
import excel "input_data/uml/mills_20200129.xlsx", firstrow clear
rename latitude lat
rename longitude lon
destring lat, replace
destring lon, replace

/* lat lon do not identify uniquely observations in this UML. 
M-00299 and M-00722 have the same coordinates. 
We could either drop one. 
duplicates tag lat lon, generate(latlon_du)
browse if latlon_du == 1
duplicates drop lat lon, force 

or replace with the coordinates in the former UML version, where they have different ones. 
NO duplicates in the version from 12/02/2020
*/

* round coordinates
replace lat = round(lat, 0.001) 
replace lon = round(lon, 0.001) 

codebook trase_code
save "temp_data/processed_UML/mills_20200129_modified.dta", replace 




***** PREPARE UNREF PROGRAMATICALLY MATCHED GROUP *****
* just rename uml's mill_id var and add the coordinates to this matching work between unref ibs and UML. 
import excel "input_data/manually_matched_ibs_uml/matching_unref/md_millMatching_modified.xlsx", firstrow clear
keep firm_id mill_id 
rename mill_id trase_code
*duplicates tag trase_code , generate(trase_code_dup)
*browse if trase_code_dup > 0
count if !mi(trase_code)
* 151 (149 - 1 (67252) +1 (55845) +1 (66743), the other cases commented are either changes in matched mill_id, or spurious matches that were already spotted by Jason. 
drop if mi(trase_code)

* add coordinates to all manually matched ibs firm_id, and parent_co and mill names in particular to the 4 mills that are not in the 
* earlier version of UML
merge 1:1 trase_code using "temp_data/processed_UML/mills_20200129_modified.dta", keepusing(trase_code parent_co mill_name lat lon)

keep if _merge == 3
drop _merge


* For those which were already found in the earlier version of UML, we prefer their parent_co and mill names because that's those that were used in the manual works
* So this update replace replaces 147 mills' names, and the 4 unmatched from master are those mills that are only in latest version of UML
merge 1:1 trase_code using "temp_data/processed_UML/traseMills_capEstyear_modified.dta", keepusing(trase_code parent_co mill_name) update replace 
drop if _merge == 2 
drop _merge

save "temp_data/processed_mill_geolocalization/md_millMatching_modified.dta", replace 




***** MERGE FIRST WAVES OF MANUAL WORKS TOGETHER *****
* (md_millMatching_modified (i.e. unref) is merged later) 

use "temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta", clear
sort firm_id year

merge m:1 firm_id using "temp_data/processed_mill_geolocalization/oto.dta", generate(merge_oto) keepusing(parent_co mill_name lat lon) update

merge m:1 firm_id using "temp_data/processed_mill_geolocalization/noto_done.dta", generate(merge_noto) keepusing(parent_co mill_name lat lon) update

merge m:1 firm_id using "temp_data/processed_mill_geolocalization/mills_to_georeference_pre2011_done.dta", generate(merge_pre2011) keepusing(mill_name lat lon) update

merge m:1 firm_id using "temp_data/processed_mill_geolocalization/mills_to_georeference_post2010_done.dta", generate(merge_post2010) keepusing(mill_name lat lon) update


***** MANUAL MODIFICATIONS / CONFLICT RESOLUTION ***** 

**** CORRECTIONS AFTER DOUBLE CHECKING OF AUTOMATIC SPATIAL MATCHES (OTO) ****
/* the double check is performed in double_checking_manual_matches.R
Case by case innvestigation is performed for firms which spatially matched UML mill_name and parent_co have no word in common with at least one of 
manufacturing directory company names matched by number of workers and district. 

When possible, we replace in priority with UML values in the format of mills_20200129.xlsx

2042 was spatially matched with M-01027. In MD it matches convincingly 	SOCFINDO AEK LOBA PT which is M-00803 in UML, which location is right near 
the village of 2042, AEK LOBA AFDELING I. M-00803 is matched with no other ibs firm. 
Therefore, 2042 is more likely M-00803 than M-01027. 
*/
replace mill_name = "AEK LOBA" if firm_id == 2042
replace parent_co = "SOCFIN INDONESIA" if firm_id == 2042
replace lat = 2.651401827 if firm_id == 2042
replace lon = 99.61735134 if firm_id == 2042


/* 
2057 was spatially matched with M-00207. In MD it matches convincingly TOLAN TIGA INDONESIA PT which is M-00796 in UML, which location is right near 
the village of 2057, PERKEBUNAN PERLABIAN. M-00796 is matched with no other ibs firm. 
Therefore, 2057 is more likely M-00796 than M-00207.  
*/
replace mill_name = "PERLABIAN" if firm_id == 2057
replace parent_co = "TOLAN TIGA INDONESIA" if firm_id == 2057
replace lat = 2.060466692	 if firm_id == 2057
replace lon = 100.0823895 if firm_id == 2057

/*
2113 was spatially matched with M-00539, named TELUK DALAM in UML. In MD, it matches with a PADASA ENAM UTAMA. 
In Google Maps, the dormitories just near M-00539 coordinates are called PADASA ENAM UTAMA, and other UML mills with this name 
are not in the same area. Therefore, we believe that M-00539 is the correct match with 2113, and the MD match too, but MD and UML 
do not correspond for some reason (PADASA ENAM UTAMA is a group name). 
*/

* 50453 matches nothing in MD; and its spatial match is removed in a manual excel work. 

/*
52424 was spatially matched to M-00455, which parent_co and mill_name are SAWIT ANUGERAH SEJAHTERA. In MD, the convincing match is 
ADEI PLANTATION AND INDUSTRY PT, which is a group name, and none of other UML references known under this name are located in the area of 
52424's village. 
*/

* For 52429, the uml parent_co actually matches the MD company name, it was just not detected in the R algorithm

* 52429 matches with nothing in MD. 
* 52436 idem 
* 52471 idem (and also this one's spatial match is actually attributed to another firm, as for 50453)
* 55630 idem 
* 66780 idem

* 54281: its match in MD is in 2012-2014 while this firm has only one IBS observation in 2005, therefore the spatial match is more trustworthy than the MD match. 

* 55779 was spatially matched to a mill called SISIRAU in UML, and its only match in MD is called SISIRAN PT, we believe there is no conflict. 

/* 55842 was spatially matched with M-00288. In MD, it matches with EKA KARYA PRIMA PT, which adress is indeed the village where 55842 is. 
In this village, there is also M-00288 and no other mill, and on Google Maps there is something called Eka Karya Prima, not on a mill. 
So the correct mill matching is M-00288, and MD reference does not designate another mill, but something related to M-00288.
*/

/* 67172 the MD match is named BAKRIE SUMATERA PLANTATIONS TBK PT, which address is in the village of 67172 (and hence its uml spatial match)
 Therefore we could believe it is the plantation related to the spatial match M-01004. 
 BUT BAKRIE SUMATERA PLANTATIONS is the name of a parent_co in UML, which mill, M-00096, is right outside the border of 67172 village PERKEBUNAN SEI BALAI
 And it matching between 67172 and MD is convincing in terms of workers and address. Moreover, M-00096 is not matched with another IBS firm. 
 Therefore, 67172 is rather M-00096 than M-01004. 
*/
replace mill_name = "UNIT KISARAN" if firm_id == 67172
replace parent_co = "BAKRIE SUMATERA PLANTATIONS" if firm_id == 67172
replace lat = 3.040709837 if firm_id == 67172
replace lon = 99.58323885 if firm_id == 67172

/*
69488 was spatially matched with M-00403. But in MD, it matches with MAZUMA AGRO INDONESIA PT (only on one number of workers, repeatedly) which is M-00723
which is located right outside the border of 69488 village SUNGAI KORANG. M-00723 is matched with no other ibs firm. 
Therefore, 69488 is more likely M-00723 than M-00403. 
*/
replace mill_name = "MAZUMA AGRO INDONESIA" if firm_id == 69488
replace parent_co = "MAZUMA AGRO INDONESIA" if firm_id == 69488
replace lat = 1.128939824 if firm_id == 69488
replace lon = 100.1552193 if firm_id == 69488



**** MANUAL WORK FOR CONFLICTING CASES WHERE DIFFERENT firm_id HAVE THE SAME MILL_NAME OR SAME COORDINATES ****

/*firm_id that were found to be the same mill as another firm_id while performing the manual geolocalization.
The oldest firm_id is kept
*/
replace firm_id = 10464 if firm_id == 68838 
replace firm_id = 36434 if firm_id == 69486
replace firm_id = 52333 if firm_id == 55834
replace firm_id = 66746 if firm_id == 70267  

/*
For BONANZA MEGAH LTD, there are two different refineries in the directories indeed, with different number of workers which each match with one of 18244 and 60758
but they have exactly the same adress. Just had a prime in front of the string. 
*/
replace mill_name = "'BONANZA MEGAH LTD, PT" if firm_id == 60758

/* 
For Baras case, there are two mills in directories, and according to this document and to google image, these mills have different capacities. 
According to number of workers and all cpo and pko variables available, we assume that 34164 is the biggest and 67078 is the smallest one. This also 
makes sense in terms of apparent age difference. 
http://ispo-org.or.id/images/notifikasi/notif_152_PT%20Unggul%20Widya%20Teknologi%20Lestari.pdf
*/
replace mill_name = "Agribaras" if firm_id == 67078
replace parent_co = "PT Unggul Widya Teknologi Lestari" if firm_id == 67078
replace lat = -1.521139 if firm_id == 67078
replace lon = 119.433556 if firm_id == 67078

/*
From directories, PT. Agri Andalas is both 36728 and 69026. 
*/
replace firm_id = 36728 if firm_id == 69026

/*
The case of 55627, 66728 and PT. Bumi Pratama Khatulistiwa: it is seemingly all the same mill but id are following but there is a duplicate in 2008
(because number of workers are the same in 2008 and many info are missing in 2008 for 55627.)
*/
drop if firm_id == 55627 & year == 2008
replace firm_id = 55627 if firm_id == 66728

/*
5283 and 69025 are also likely the same mill, namely PT. Dario Dharma Pratama - Ipuh 
because number of workers match for both in contemporaneous directories. 
Moreover, 70449 and 70450 are less likely Ipuh because they match only on one number of workers while 69025 matches on 133 and 136. 
*/
replace firm_id = 5283 if firm_id == 69025 

replace mill_name = "" if firm_id == 70449 
replace parent_co = "" if firm_id == 70449 
replace lat = . if firm_id == 70449 
replace lon = . if firm_id == 70449 

replace mill_name = "" if firm_id == 70450 
replace parent_co = "" if firm_id == 70450 
replace lat = . if firm_id == 70450 
replace lon = . if firm_id == 70450 
/*
2125 and 2127 in Torgamba. 
2125 matches in directories with TASIK RAJA, PT, which is in heilmayr mills but outside of BUKIT TUJUH village. 2127 remains PT. Torganda - PKS Sibisa Mangantur, the only
mill from Torganda in this area (2127 matches with a TORGANDA, PT in directories, the only one in this region.)
*/
replace mill_name = "PT. Tasik Raja" if firm_id == 2125
replace parent_co = "PT Tasik Raja" if firm_id == 2125
replace lat = 1.669279923 if firm_id == 2125
replace lon = 100.1615075 if firm_id == 2125


/*
in some years there are two PASANGKAYU, PT in the same desa in the directories (and at the same time 48970 and 67074 overlap)
There is only one big mill in this area. Hence, it is likely all the same mill. 
In the four years (2006-2009) when the two firm_id overlap in IBS, quantities are all duplicated. As if they had just replicated the quantities from 2005. 
This is one more hint that it is the same mill, with valid obs. from 1998 to 2005 for 48970 and then from 2006 for 67074. 
*/
drop if firm_id == 48970 & (year == 2006 | year == 2007 | year == 2008 | year == 2009)
replace firm_id = 48970 if firm_id == 67074 


/*
PT Sawita Leidong Jaya matches with 52434 on several numbers of workers, while 76315 only on 105 (because only one observation in 2015)
*/

replace mill_name = "" if firm_id == 76315 
replace parent_co = "" if firm_id == 76315 
replace lat = . if firm_id == 76315 
replace lon = . if firm_id == 76315 

/*
70396 has 2 matches with directories while 511439 has much more. 
*/
replace mill_name = "" if firm_id == 70396 
replace parent_co = "" if firm_id == 70396 
replace lat = . if firm_id == 70396 
replace lon = . if firm_id == 70396 

/*
68248 (from oto) clearly matches Citra Riau Sarana in Muara Langsat, which is CRS - Teso II in ref_mills. 
72997 has observations only since 2013 and matches CRS - Pks Jake in 2015. But Pks Jake has directories since 2006 so it is necessarily 
CRS - Teso I (because Teso I is estimated to appear in 2004 and Teso III in 2013 - and google earth 2007 images confirm Teso III was not there but Teso I was). 
As we don't take appearance year in IBS as an trustworthy indicator of actual birthdate, we can suppose that 72997 is CRS - Teso 1. 
*/
replace mill_name = "PT. Citra Riau Sarana - Teso I" if firm_id == 72997 
replace parent_co = "PT Citra Riau Sarana" if firm_id == 72997 
replace lat = -0.212423211 if firm_id == 72997 
replace lon = 101.4734808 if firm_id == 72997 



/*
68244 matches several times and more precisely on geography with Hasil than 76312
*/
replace mill_name = "" if firm_id == 76312 
replace parent_co = "" if firm_id == 76312 
replace lat = . if firm_id == 76312 
replace lon = . if firm_id == 76312 

/*
55727 matched spatially with PT. Kalimantan Sawit Abadi - Natai Baru, but this mill (found in directory) never matches with in terms of workers with 55727. 
Rather, 55727 is found to match sawit sumber mas sarana in directories, in Sulung desa. . 
*/
replace mill_name = "PT. Sawit Subermas Sarana - Sulung Mill" if firm_id == 55727 
replace parent_co = "PT Sawit Sumbermas Sarana" if firm_id == 55727 
replace lat = -2.304451535 if firm_id == 55727 
replace lon = 111.6132597 if firm_id == 55727 

/*
77946 seems to be the same mill as 75041 because their years are continuous, and they have the same number of workers in the same district. 
*/
replace firm_id = 75041 if firm_id == 77946

/*
54271 cannot be assumed to be the same mill as 72995 because there are at least two mills in Pelintung according to imagery.
72995 matching on number of workers is more trustworthy than the mere spatial match of 54271. 
*/
replace mill_name = "" if firm_id == 54271 
replace parent_co = "" if firm_id == 54271 
replace lat = . if firm_id == 54271 
replace lon = . if firm_id == 54271 

/*
Clearly, the automatic spatial match between 52344 and Begerpang mill is right. 
The match in 2015 on 99 for 72956 likely links to a head quarter office rather because we cannot find two Lonsum mills in deli serdang.  
*/
replace mill_name = "" if firm_id == 72956 
replace parent_co = "" if firm_id == 72956 
replace lat = . if firm_id == 72956 
replace lon = . if firm_id == 72956

/*
51434 matches with masterindo on several numbers of workers while 72994 only on 132.
72994 matches with Swastisiddhi Amagra matches 2013 ibs with 2014 directory. 
*/ 		
replace mill_name = "PT. Swastisidi Amagra" if firm_id == 72994 
replace parent_co = "PT. Swastisiddhi Amagra" if firm_id == 72994 
replace lat = 0.210324795 if firm_id == 72994 
replace lon = 101.2753159 if firm_id == 72994

/*
55635 and 72781 are likely both PT. Poliplant Sejahtera, because they both match it in their contemporaneous directories 
and don't have overlapping years. 
*/
replace firm_id = 55635 if firm_id == 72781

/*
68247 is a better match than 76319
*/
replace mill_name = "" if firm_id == 76319 
replace parent_co = "" if firm_id == 76319 
replace lat = . if firm_id == 76319 
replace lon = . if firm_id == 76319



sort firm_id year 
*we check for duplicates firm_id coordinates (or very close ones, relating to the same mill) on R. These cases are: 

/*
4780 (found with spatial matching) is clearly PTP Nusantara Sungai Lengi based on number of workers. 
The other one, 56417, is matching with another PTP Nusantara VII in directories, but it cannot be found in ref, nor on 
google maps/earth, nor on internet. Although it produces CPO in IBS, i tmight be the plantations associated to Sungai Lengi,
as is stated here: http://www.ptpn7.com/displaycontent.aspx?topic=Distrik%20Muara%20Enim 
*/
replace mill_name = "" if firm_id == 56417 
replace parent_co = "" if firm_id == 56417 
replace lat = . if firm_id == 56417 
replace lon = . if firm_id == 56417

/*
The case of 36342 and 67075 with Letawa, PT. 
36342 has observations from 1998 to 2009, and all duplicates from 2006 on. It matches directories in 2003 and 2006. 
In these years, there are only one Letawa in directories. 
Then, in 2009 and 2010, there are 2 Letawa in directories, in the same desa, with different numbers of workers that match 
36342 and 67075 respectively. 
Then, there is only one Letawa again, that matches 67075 now. 
Therefore, it seems likely that it has been the same mill from the beginning, but with some change in firm_id. 
The overlapping years of firm_id are 2006-2009. 
For 36342, these are duplicates of 2006. 
For 67075, the quantity is misreported (>100k tons of CPO.)
Therefore, we believe the most credible observation of Letawa is 36342 until 2006 included, and 67075 afterwards
*/
drop if firm_id == 36342 & year > 2006 
drop if firm_id == 67075 & year == 2006 
replace firm_id = 36342 if firm_id == 67075 

/*
see https://www.industryabout.com/country-territories-3/1937-indonesia/vegetable-oil-industry/30072-bsp-agro-mitra-madani-palm-oil-mill
https://www.tuv.com/media/indonesia/brochure_2/system_clients_2016/per_2017/update_public_summary_rspo_part_1/RSPO_ASA4_Public_Summary_PT_Agrowiyana_BSP_Jambi.pdf
67238 matches Agro MITRA MADANI, which cannot be located where we thought it was according to these two sources, because this 
is alreadz where PT Trimitra Lestari is. But these sources are trustworthy (TÜV especially). Other documents online listing
palm oil mills give the same coordinates to both mills, pointing at the same factory.
What happens is that there is another mill very closeby, a bit more west. According to google names, this is Mitra Madani 
(they also give the name of Agro Wiyana, which relates to Mitra Madani according to above tüv source.) and the one at which 
all documents point at is Trimitra Lestari. 
*/
replace lat = -1.078053 if firm_id == 67238 
replace lon = 103.110272 if firm_id == 67238


/*
69531 is PT. Agro Sarimas Indonesia.
54274 might also be PT. Agro Sarimas Indonesia because it is observed in same desa years continuously before,
but that is likely rather a coco mill than a palm mill
https://www.gmdu.net/corp-414564.html and we have no palm data in IBS. 
May be it was a conversion of activity. Let us remove mill name and coordinates to 54274 and assume it indeed was not a 
palm oil mill at this time.  
*/
replace mill_name = "" if firm_id == 54274 
replace parent_co = "" if firm_id == 54274 
replace lat = . if firm_id == 54274 
replace lon = . if firm_id == 54274

/*
56234 is a more precise match with Andalas than 76289. 
*/
replace mill_name = "" if firm_id == 76289 
replace parent_co = "" if firm_id == 76289 
replace lat = . if firm_id == 76289 
replace lon = . if firm_id == 76289

/*
70390 and 76329 have both only one match with directories when looking only to n-1, n and n+1 year directories, 
which RA did in the first place. But actually, 70390 has two matches when looking at n-2 year directories. 
Therefore we assume 70390 is PT. Bina Sawit Nusantara 
*/
replace mill_name = "PT. Bina Sawit Nusantara" if firm_id == 70390 
replace parent_co = "PT. Bina Sawit Nusantara" if firm_id == 70390 
replace lat = 0.165799445 if firm_id == 70390 
replace lon = 101.3260084 if firm_id == 70390

/*
55775 has better matches with PT. Socfin - Lae Butar (Desa rimo in directories).
71593 might rather be PT. Delima Makmur.   
*/
replace mill_name = "PT. Delima Makmur" if firm_id == 71593 
replace parent_co = "PT. Delima Makmur" if firm_id == 71593 
replace lat = 2.245276912 if firm_id == 71593 
replace lon = 98.03032975 if firm_id == 71593

/*
2062 is clearly PTP Nusantara Sei Silau. 
In 2015 directory there is Sei Silau and SEI S´, each matching with one of these 2 IBS mills. 
2062 has no observation in 2015, and many missings in 2014. 
72957 has obs. in 2013-2015 but mainly missing in 2013, 2014. 
They have the same number of workers in IBS (185 in 2014 for 2062 and in 2015 for 72957)
That looks like a firm_id change that resulted in some information lost. 
*/
drop if firm_id == 72957 & (year == 2013 | year == 2014)
replace firm_id = 2062 if firm_id == 72957

/*
49003 better matches Sei Mangkei 
*/
replace mill_name = "" if firm_id == 72959 
replace parent_co = "" if firm_id == 72959 
replace lat = . if firm_id == 72959 
replace lon = . if firm_id == 72959

/*
For 2119, obs. on years 2004 and 2005 are duplicates of year 2003 (for quantities but not for values, as is often the case)
72959 goes from 2004, with same number of workers at that time. 
*/
drop if firm_id == 2119 & (year == 2004 | year == 2005)
replace firm_id = 2119 if firm_id == 52344

/*
browse firm_id year workers_total_imp3 district_name kec_name village_name  mill_name parent_co lat lon ///
merge_oto merge_noto merge_pre2011 merge_post2010 in_ton_ffb ffb_price_imp1 out_ton_cpo out_ton_cpo_imp2 out_val_cpo out_val_cpo_imp2 ///
cpo_price_imp2 out_ton_pko out_ton_pko_imp2 out_val_pko out_val_pko_imp2 pko_price_imp2     ///
if firm_id == 2119 | firm_id == 52344
*/
sort firm_id year 
/*
*check these conflicting cases after the manual modifications

bys mill_name: egen n_dif_id_per_mill_name = nvals(firm_id)
count if n_dif_id_per_mill_name > 1 & !mi(mill_name) 
*drop n_dif_id_per_mill_name

bys lon lat: egen n_dif_id_per_coord = nvals(firm_id)
count if n_dif_id_per_coord > 1 & !mi(lat) 
*drop n_dif_id_per_coord

OR : 
duplicates tag mill_name if !mi(mill_name), generate(du_mill_name) 
duplicates tag firm_id mill_name if !mi(mill_name), generate(du_firm_id_mill_name) 
sort mill_name firm_id year 

browse firm_id year workers_total_imp3 district_name kec_name village_name  mill_name parent_co lat lon ///
merge_oto merge_noto merge_pre2011 merge_post2010 du_mill_name du_firm_id_mill_name if du_mill_name != du_firm_id_mill_name
*/

**** ADD COORDINATES TO UNREF firm_id ****

merge m:1 firm_id using "temp_data/processed_mill_geolocalization/md_millMatching_modified.dta", ///
 generate(merge_unref_geo) keepusing(trase_code parent_co mill_name lat lon) update

sort firm_id year 

* we indeed have 151 new firm_id being updated on their names, lat lon 
*codebook firm_id if merge_unref_geo == 4
save "temp_data/processed_mill_geolocalization/merge_geoloc_works_temp.dta", replace 



**** THE OVERALL CONFLICT RESOLUTION WAS MADE MANUALLY IN AN EXCEL SPREADSHEET. IT IS MERGED BACK WITH FULL PANEL HERE ****

** check the conflicts brought by this new wave of geolocalized mills
/*
bys mill_name: egen n_dif_id_per_mill_name = nvals(firm_id)
gen name_du = n_dif_id_per_mill_name > 1 & !mi(mill_name) 
*codebook firm_id if name_du == 1

bys lon lat: egen n_dif_id_per_coord = nvals(firm_id)
gen coord_du = n_dif_id_per_coord > 1 & !mi(lat) 
*codebook firm_id if coord_du == 1

keep if name_du == 1 | coord_du == 1
 sort mill_name firm_id year 
export excel mill_name parent_co trase_code lat lon firm_id year min_year workers_total_imp3 district_name kec_name village_name     ///
merge_oto merge_noto merge_pre2011 merge_post2010 merge_unref_geo n_dif_id_per_mill_name name_du n_dif_id_per_coord coord_du ///
using "C:\Users\GUYE\Desktop\opalval\build\input\mill_geolocalization\overall_btw_conflicts.xlsx", firstrow(variables) replace 
*/

** introduce the spreadsheet used to resolve overall between conflicts into a .dta 
import excel "input_data/manually_matched_ibs_uml/overall_btw_conflicts_done.xlsx", firstrow clear  
destring lat lon, replace dpcomma

*keep only the merging key variables and the variables that have been modified manually and hence that we don't want to update
keep firm_id year mill_name parent_co lat lon trase_code 

* this merge keeps unchanged the above variables for this block of obs., and adds all the rest of the panel. 
merge 1:1 firm_id year using "temp_data/processed_mill_geolocalization/merge_geoloc_works_temp.dta", generate(merge_overall_cfl) 
order year, before(mill_name)
order firm_id, before(year)
sort firm_id year 



**** HANDLE CASES WHERE SEVERAL firm_id WERE FOUND TO ACTUALLY REFER TO THE SAME MILL (THE firm_id CHANGED OVER TIME) ****

/*
Code used to spot them
*drop n_dif_id_per_mill_name name_du n_dif_id_per_coord coord_du

bys mill_name: egen n_dif_id_per_mill_name = nvals(firm_id)
gen name_du = n_dif_id_per_mill_name > 1 & !mi(mill_name) 
codebook firm_id if name_du == 1

bys lon lat: egen n_dif_id_per_coord = nvals(firm_id)
gen coord_du = n_dif_id_per_coord > 1 & !mi(lat) 
codebook firm_id if coord_du == 1

sort mill_name firm_id year 
browse if name_du == 1 | coord_du == 1 
*/

/*
50450 56305 are likely one same mill (the last year of the former is the first year of the latter, both with 1828 workers). 
As usual: keep the line with the most information on quantities and values and keep the earlier firm_id. 
Here, the year with the most information is the first of 56305. The two firm_id have common values in 2006 on establishment level 
(like workers or export_pct) but not at commodity variable level. 
*/
drop if firm_id == 50450 & year == 2006
replace firm_id = 50450 if firm_id == 56305

* 68951 has only one record in 2008 that matches workers and kecamatan of SURYA RAYA LESTARI II PT in 2009 MD. 
replace firm_id = 68951 if firm_id == 69449



**** AND RESOLVE NEW CONFLICTS ARISING FROM NEW MATCHES MADE THROUGH THE OVERALL CONFLICT RESOLUTION ****

/* 76306 matches only one number of workers vs. 2 for 70386 */
replace mill_name = "" if firm_id == 76306
replace lat = . if firm_id == 76306
replace lon = . if firm_id == 76306


/*
56383 is initially a noto, that was matched to PT. Inti Indosawit Subur - Tungkal Ulu in the noto manual resolution, 
but then matched rather to Rudy Agung AgraLaksana in the overall conflict resolution - which remains the more convincing match. 
(comment in overall_btw_resolution_done.xlsx is:
"This one is more likely Rudy Agung AgraLaksana, as they match two different workers 
values - 104 in 2010 with 2014 directory and 75 in 2009 with 2010 directory. It was not prefered over PT. Inti Indosawit Subur - Tungkal Ulu 
because I thought I could not see it in 2006 on Google earth Pro timelaps, but actually using
 GEE  LANDSAT/LT05/C01/T1_SR 2004-2006 it is less clear that is was not there at that time. Finally, the MD address matches the IBS kecamatan."

69008 is a oto. Looking at workers and reported desa in IBS and MD, it rather seems that it is PALMA ABADI PT in MD, which is 
in turn M-00551 which exact location is indeed just on the other side of the village reported in IBS and MD by 69008, explaining
why it was spatially mismatched.
*/
replace mill_name = "PT. Palma Abadi" if firm_id == 69008
replace parent_co = "PT. Palma Abadi" if firm_id == 69008
replace lat = -1.321715223 if firm_id == 69008
replace lon = 103.2653722 if firm_id == 69008
replace trase_code = "M-00551" if firm_id == 69008


**** RESOLVE CONFLICTS FOUND FROM INTERSECTIONS OF 500m BUFFERS AROUND COORDINATES ****

* 40139 70397
/* 40139 matches better in workers and location with Sei Pagar than 70397
70397 was thought to be Sei Pagar because it had one n-1 worker match in MD. 
But actually it is more likely that it is GANDA BUANINDO PT now that we look at all possible MD for all ibs years (then it matches 
on 2 different numbers of workers, in three years.)
*/
replace mill_name = "PT. Ganda Buanindo" if firm_id == 70397
replace parent_co = "PT Ganda Buanindo" if firm_id == 70397
replace trase_code = "M-00237" if firm_id == 70397
replace lat = 0.007161207 if firm_id == 70397
replace lon = 101.2395325 if firm_id == 70397

* 44938 67076 those are Suryaraya Lestari 1
/* The decision taken here is to remove the four last years of 44938 because 3 of these 4 years are duplicates to previous year(s) for quantity variables.
ATTENTION: This assumes that those two firms are the same one while this is not so obvious: 
Suryaraya Lestari 1 is only one mill in the UML, while it is two different ones in MD 2008, 2009 and 2010 with exactly the same name but 
different number of workers, and slightly different addresses; and two IBS firms (44938 and 67076) match them respectively, and  
look like different facilities, that both produce CPO and PKO, but in different quantities and values. 
*/
drop if firm_id == 44938 & (year == 2006 | year == 2007 | year == 2008 | year == 2009)
replace firm_id = 44938 if firm_id == 67076


* 44993 55831 are actually two separate mills from imagery! 

* 51385 76247 
/* look like the same mill: 76247 is just the 2015 record, with the same number of workers as years before in 51385, and matching to only one
MD mill, the same 51385 matches with, and that has been matched in md_millMatching_modified to M-00422 */
replace mill_name = "PT. Langkat Nusantara Kepong" if firm_id == 76247
replace parent_co = "PT Langkat Nusantara Kepong" if firm_id == 76247
replace trase_code = "M-00422" if firm_id == 76247
replace lat = 3.560677077563 if firm_id == 76247
replace lon = 98.39074238997 if firm_id == 76247
replace firm_id = 51385 if firm_id == 76247

* 56377 67238 
/* it does not look like they refer to the same mill for 2007 record, because they each have different cpo info. 
and 56377 is a slightly better match with AGRO MITRA MADANIsince this is PT. Agrowiyana according to Jason (md_millMatching_modified)
and the desa (Brasau) reported in IBS by 56377 straight next to this mill's location. */
replace mill_name = "" if firm_id == 67238
replace parent_co = "" if firm_id == 67238
replace trase_code = "" if firm_id == 67238
replace lat = . if firm_id == 67238
replace lon = . if firm_id == 67238

* 69489 71626
/* The argument that was used to geolocalize 71626 was (from mills_to_georeference_post2010.xlsx) was: 
"Match found in directory but not in ref_mills. Mill name and coordinates from google maps (looking around desa rondaman lombang)."
Which makes it a weaker candidate than 69489 that matches on three numbers of workers and on village location
*/
replace mill_name = "" if firm_id == 71626
replace parent_co = "" if firm_id == 71626
replace trase_code = "" if firm_id == 71626
replace lat = . if firm_id == 71626
replace lon = . if firm_id == 71626



***** ADD MORE UML VARIABLES ***** 
/* 

/// ATTENTION \\\ mill_name does not identify uniquely firm_id! There is PT. Torganda for instance (and non UML don't always have a name anyways)

* donc là on en a 151 qui ont un trase_code depuis mill_mdMatching, de base. 
Puis 8 qui ont été ajoutés manuellement dans overall_btw_conflicts
Puis +3 avec les changements manuels dans la résolution des nouveaux conflits générés (ceux à qui on enlève n'en ont pas en fait)
Puis avec le merge_trase_code on ajoute le trase_code à 295 firm_id. 
*/
codebook firm_id if !mi(trase_code)
*gen nm_tc = !mi(trase_code)
codebook firm_id if !mi(lat)
* add to all these geolocalized ibs mills their trase_code and est_year variables from traseMills_capEstyear

* It is necessary to round coordinates to match them (some decimals have been lost for some reason (manual work or data wrangling) 
* and therefore reduce the number of matches on coordinates.  
replace lat = round(lat, 0.001) 
replace lon = round(lon, 0.001) 

merge m:1 lat lon using "temp_data/processed_UML/traseMills_capEstyear_modified.dta", generate(merge_trase_code) keepusing(trase_code est_year) update 
drop if merge_trase_code == 2
*tab merge_trase_code nm_tc
/*those that already have a trase_code but though don't match with traseMills_capEstyear are those with a recent trase_code, 
or Suryaraya Lestari 2 coordinates from traseMills_capEstyear which are wrong. 
*/

merge m:1 lat lon using "temp_data/processed_UML/mills_20200129_modified.dta", generate(merge_trase_code2) keepusing(trase_code active uml_id) update 
drop if merge_trase_code2 == 2

sort firm_id year

**** FIRMS THAT WERE MANUALLY GEOLOCALIZED AND HAVE NO UML MILL WITHIN 1km
* (this is spotted in spatial_routines.R)
* 1763 was indeed not in UML, Jason has added it (cf his email)
* 4782 is a refinery 
replace mill_name = "" if firm_id == 4782
replace parent_co = "" if firm_id == 4782
replace lat = . if firm_id == 4782
replace lon = . if firm_id == 4782

* 4783 is a refinery 
replace mill_name = "" if firm_id == 4783
replace parent_co = "" if firm_id == 4783
replace lat = . if firm_id == 4783
replace lon = . if firm_id == 4783

* 40024 is a refinery 
replace mill_name = "" if firm_id == 40024
replace parent_co = "" if firm_id == 40024
replace lat = . if firm_id == 40024
replace lon = . if firm_id == 40024

* 66897 is a refinery 
replace mill_name = "" if firm_id == 66897
replace parent_co = "" if firm_id == 66897
replace lat = . if firm_id == 66897
replace lon = . if firm_id == 66897

/* 
* 68193, MD name: PABRIK KELAPA SAWIT TEUPAH, found at 2.398 96.423
The IBS shows input and output data for this facility for 2008 and 2009 and MD has info for this in 2009. 
However, the mill you located looks to have only started to be constructed in 2010 according to: 
https://www.infosawit.com/news/4042/semenjak-2010-pks-milik-kabupaten-simeulue-belum-rampung
https://aceh.tribunnews.com/2016/01/27/polres-simeulue-akan-panggil-rekanan-pembangunan-pabrik-kelapa-sawit
Also, the mill was never fully operational according to those articles as well. I can't find any other information based on the name of the facility in the MD though.
*/
replace mill_name = "" if firm_id == 68193
replace parent_co = "" if firm_id == 68193
replace lat = . if firm_id == 68193
replace lon = . if firm_id == 68193

/* 
* 72948 is actually a mill, the comment from mills_to_georeference_post2010 is:
match in 2014 directory (and still present in 2015 with still 131 workers). 
coordinates found with google maps following a search for Perk. Simpang Gambir, the name of the desa given in the directory, leading to this photo: 
https://www.google.com/maps/place/Perk.+Simpang+Gambir,+Lingga+Bayu,+Kabupaten+de+Mandailing+Natal,+Sumatra+du+Nord,+Indon%C3%A9sie/@0.5917666,99.3055242,3a,75y,90t/data=!3m8!1e2!3m6!1sAF1QipOdlj0PGQH9vTLKR4IKMiBcNGR3LNzHI0CAzMk1!2e10!3e12!6shttps:%2F%2Flh5.googleusercontent.com%2Fp%2FAF1QipOdlj0PGQH9vTLKR4IKMiBcNGR3LNzHI0CAzMk1%3Dw160-h120-k-no!7i1280!8i960!4m13!1m7!3m6!1s0x302a2d3dddda127f:0x8e3d68fb54b397d4!2sPerk.+Simpang+Gambir,+Lingga+Bayu,+Kabupaten+de+Mandailing+Natal,+Sumatra+du+Nord,+Indon%C3%A9sie!3b1!8m2!3d0.5888083!4d99.3506604!3m4!1s0x302a2d3dddda127f:0x8e3d68fb54b397d4!8m2!3d0.5888083!4d99.3506604

UPDATE:
Jason has added it to UML and it will thus be matched now. 
 */


* Firms without trase_code that have been given coordinates pointing at UML mills have been corrected for in their respective manual spreadsheets.

*** Firms that have been given a name from MD but that could not be given coordinates. 
/*
gen no_coord = (mi(lat) & !mi(mill_name))
keep if no_coord == 1
keep firm_id year mill_name parent_co lat lon workers_total_imp3 ///
merge_oto merge_noto merge_pre2011 merge_post2010 merge_unref_geo merge_overall_cfl merge_trase_code merge_trase_code2 /// 
district_name kec_name village_name 

export excel using "C:\Users\GUYE\Desktop\opalval\build\input\mill_geolocalization\names_but_no_coord_raw.xlsx", firstrow(variables) replace 

browse firm_id year mill_name parent_co lat lon ///
merge_oto merge_noto merge_pre2011 merge_post2010 merge_unref_geo merge_overall_cfl merge_trase_code merge_trase_code2 /// 
district_name kec_name village_name if no_coord == 1

1760 I don't understand where this mill_name comes from; 
2035 idem
*/
/* 3783 - Matched with 'LEMBAH KRYA, PT' from the MD that shows the product as Minyak Kelapa (Coconut oil), also no FFB input or CPO output in 
the IBS_UML_panel.xlsx. 
http://scholar.unand.ac.id/15586/2/BAB%20I%20PENDAHULUAN.pdf
https://opencorpdata.com/id
https://m2indonesia.com/informasi/perusahaan/profil-perusahaan-lembah-krya-pt-unit-minyak-kelapa-padang-provinsi-sumatera-barat.htm
*/

/*
56374 - Matched with 'Agroindo Indah Perkasa' from MD. The inputs for this facility are both FFB and CPO so possibly a refinery/kernel crusher? 
Also, the main product based on 
https://companychambers.blogspot.com/2018/11/profil-dan-alamat-lengkap-pt-agroindo.html
https://www.smartbisnis.co.id/directory/agroindo-indah-perkasa 
states that it is Minyak Goreng / Curah (Bulk Cooking Oil). 
*/

/*
70406 - Matched with 'Tunas Lestari Sejati, PT' from MD. 
Could not find any specific information about this being a mill. The best I could find is this article 
https://jambi-independent.co.id/read/2018/10/09/30309/perusahaan-bermasalah-dengan-petani
saying it is a plantation company. 
*/



***** HOMOGENIZE MILL AND PARENT COMPANY NAME VARIABLES *****

* Remove names of refineries or names that are mistakes. 
replace mill_name = "" if firm_id == 1760
replace mill_name = "" if firm_id == 2035
replace mill_name = "" if firm_id == 3783
replace mill_name = "" if firm_id == 4556
replace mill_name = "" if firm_id == 10466
replace mill_name = "" if firm_id == 10470
replace mill_name = "" if firm_id == 10471
replace mill_name = "" if firm_id == 18244
replace mill_name = "" if firm_id == 33757
replace mill_name = "" if firm_id == 34160
replace mill_name = "" if firm_id == 43180
replace mill_name = "" if firm_id == 51460
replace mill_name = "" if firm_id == 56374
replace mill_name = "" if firm_id == 60758
replace mill_name = "" if firm_id == 71677
replace mill_name = "" if firm_id == 72952

* give the names from most recent UML to all obs. that have a trase_code 
merge m:1 trase_code using "temp_data/processed_UML/mills_20200129_modified.dta", ///
 generate(merge_uml_names) keepusing(mill_name parent_co lat lon) update replace 
 drop if merge_uml_names == 2
/* this also replaced all the coordinates, which are anyways the same in all cases but one, M-00455, where traseMills_capEstyear points to an older facility,
while the more recent UML points to a more recent facility. */

/*
codebook trase_code 
codebook mill_name 
codebook mill_name if !mi(trase_code)
codebook firm_id if !mi(trase_code) 
codebook firm_id if !mi(lat)
*/

/* Test that all firm_id have only one unique set of coordinates 
*drop n_lat_per_firm lat_error
bys firm_id: egen n_lat_per_firm = nvals(lat)
gen lat_error = n_lat_per_firm > 1 & !mi(lat)
codebook firm_id if lat_error == 1
browse firm_id year lat lon lat_error if lat_error == 1

* drop n_lon_per_firm lon_error
bys firm_id: egen n_lon_per_firm = nvals(lon)
gen lon_error = n_lon_per_firm > 1 & !mi(lon)
codebook firm_id if lon_error == 1
browse firm_id year lon lon lon_error if lon_error == 1
*/

/* Test that all coordinates have only one firm_id 
bys lon lat: egen n_dif_id_per_coord = nvals(firm_id)
gen coord_du = n_dif_id_per_coord > 1 & !mi(lat) 
codebook firm_id if coord_du == 1
browse firm_id year lon lat if coord_du == 1
* there is only SURYARAYA LESTARI 1 
drop n_dif_id_per_coord coord_du
*/

/* Test that all trase_code have only one unique set of coordinates 
bys trase_code: egen n_lat_per_trase_code = nvals(lat)
gen lat_error = n_lat_per_trase_code > 1 & !mi(trase_code)
codebook trase_code if lat_error == 1
browse firm_id year lat lon lat_error if lat_error == 1
drop n_lat_per_trase_code lat_error

* Test that all coordinates have only one trase_code 
bys lon lat: egen n_dif_tc_per_coord = nvals(trase_code)
gen coord_du = n_dif_tc_per_coord > 1 & !mi(lat) 
codebook trase_code if coord_du == 1
browse firm_id year trase_code mill_name lon lat if coord_du == 1 
drop n_dif_id_per_coord coord_du

The only cases so far are the three mills that are not in UML. 
*/

* One last cleaning action: a mill geolocalized in Papua has island_name in Maluku, it is 46510
replace island_name = "Papua" if firm_id == 46510
replace island_factor = 5 if island_factor == 4

* I further checked that no firm has missing island_name while it has more precise (province, district, village) info) 


*** Distinguish geolocalized IBS mills, mills, and the rest. 
* at this stage, only mills matched with UML have coordinates
gen uml_matched_sample = (!mi(lat))

sort firm_id year
global mill_commo out_ton_cpo out_val_cpo out_ton_pko out_val_pko in_ton_ffb in_val_ffb
foreach var of varlist $mill_commo{
	bys firm_id: egen sum_`var' = sum(`var')
}

* The criterion to be a mill is to have at least one record of producing CPO or PKO, either in quantity or in value AND one record of sourcing FFB. 
* The AND condition removes 184 establishments from the is_mill sample with respect to the OR condition (producing either cpo or pko OR sourcing ffb)
* In the subsample of 470 mills matched with UML, the difference between the AND and the OR conditions is only 11 mills 
* (indicating that in IBS, almost all establishments credibly recognized as mills declare some FFB sourcing) 
* However, we are not interested only in mills declaring buying FFB. We also want to look at potentially integrated plantation-mill systems
* in which FFB is not traded, and the only output is CPO. 
* Therefore, we use the OR condition. 

gen is_mill = (((sum_out_ton_cpo >0 & !mi(sum_out_ton_cpo)) | ///
			    (sum_out_val_cpo >0 & !mi(sum_out_val_cpo)) | ///
			    (sum_out_ton_pko >0 & !mi(sum_out_ton_pko)) | ///
			    (sum_out_val_pko >0 & !mi(sum_out_val_pko))) | ///
			   ((sum_in_ton_ffb >0 & !mi(sum_in_ton_ffb)) | ///
			    (sum_in_val_ffb >0 & !mi(sum_in_val_ffb)))) 
* also, do not count mills that satisfy the above conditions but are in Java or Bali; 
replace is_mill = 0 if island_name == "Java" | island_name == "Bali Nusa Tenggara" 

* Of the 470 geolocalized mills, 2 have only missing info, 7 have missing FFB info, and 2 are in Java. 

* Now, add the village centroid coordinates.

merge m:1 firm_id using "temp_data/processed_mill_geolocalization/IBSmills_desacentro.dta", generate(merge_desa_centro) update
* nonmissing conflicts are those that have true coordinates 
* missing updated are those that have not been geo-localized through UML matching, and are now approximately geo-localized with their village centroids
* not matched from using come from the fact that some firm_id have been removed above in this script. We do not want to reintroduce them
drop if merge_desa_centro == 2

* round these new coordinates too
replace lat = round(lat, 0.001) if merge_desa_centro == 4
replace lon = round(lon, 0.001) if merge_desa_centro == 4

gen geo_sample = (!mi(lat))

gen analysis_sample = (geo_sample == 1 & is_mill == 1)

*codebook firm_id if analysis_sample == 1
* Now 587 establishments are identified as palm oil mills and are geo-localized. 466 precisely through matching with UML, and 121 with their valid, most recent, village centroids  



** remake min_year variables (to update for firm_id changes)
bys firm_id: egen minmin_year = min(min_year)
drop min_year 
rename minmin_year min_year
order min_year, after(lon)
sort firm_id year 
* and max_year
bys firm_id: egen maxmax_year = max(max_year)
drop max_year 
rename maxmax_year max_year
order max_year, after(min_year)
sort firm_id year 

** some tidying
order trase_code, after(year)
order uml_id, after(trase_code)

order desa_id, after(desa_code)
order kabu_code, after(kab_code)

order est_year, after(min_year)
order startYear, after(est_year)
order active, after(max_year)
label var min_year "first ibs year"
label var est_year "UML estimated first year"
label var startYear "Rothenberg estimated first year (uncleaned)"

order island_factor, after(lon)
order island_name, after(island_factor)
order district_name, after(island_name)
order kec_name, after(district_name)
order village_name, after(kec_name)

* save a panel complete version 
sort firm_id year 
save "temp_data/processed_mill_geolocalization/IBS_geo_panel.dta", replace







codebook firm_id if !mi(lat)
codebook firm_id if !mi(uml_matched_sample)
codebook firm_id if !mi(trase_code)
codebook firm_id if !mi(mill_name)

export excel firm_id year uml_matched_sample is_mill geo_sample analysis_sample trase_code uml_id mill_name parent_co lat lon district_name kec_name village_name /// 
min_year est_year startYear max_year active industry_code ///
ffb_price_imp1 ffb_price_imp2 in_ton_ffb in_ton_ffb_imp1 in_ton_ffb_imp2 in_val_ffb in_val_ffb_imp1 in_val_ffb_imp2 flag_multiinput_ffb ///
in_dom_cpo_price_imp1 in_dom_cpo_price_imp2 in_dom_ton_cpo in_dom_ton_cpo_imp1 in_dom_ton_cpo_imp2 in_dom_val_cpo in_dom_val_cpo_imp1 in_dom_val_cpo_imp2 ///
in_imp_cpo_price_imp1 in_imp_cpo_price_imp2 in_imp_ton_cpo in_imp_ton_cpo_imp1 in_imp_ton_cpo_imp2 in_imp_val_cpo in_imp_val_cpo_imp1 in_imp_val_cpo_imp2 ///
in_tot_cpo_price_imp1 in_tot_cpo_price_imp2 in_tot_ton_cpo in_tot_ton_cpo_imp1 in_tot_ton_cpo_imp2 in_tot_val_cpo in_tot_val_cpo_imp1 in_tot_val_cpo_imp2 ///
cpo_price_imp1 cpo_price_imp2 out_ton_cpo out_ton_cpo_imp1 out_ton_cpo_imp2 out_val_cpo out_val_cpo_imp1 out_val_cpo_imp2 prex_cpo prex_cpo_imp1 prex_cpo_imp2 out_cpo ///
pko_price_imp1 pko_price_imp2 out_ton_pko out_ton_pko_imp1 out_ton_pko_imp2 out_val_pko out_val_pko_imp1 out_val_pko_imp2 prex_pko prex_pko_imp1 prex_pko_imp2 out_pko ///
out_ton_rpo out_ton_rpo_imp1 out_ton_rpo_imp2 out_val_rpo out_val_rpo_imp1 out_val_rpo_imp2 prex_rpo prex_rpo_imp1 prex_rpo_imp2 out_rpo tag_multioutput_rpo ///
out_ton_rpko out_ton_rpko_imp1 out_ton_rpko_imp2 out_val_rpko out_val_rpko_imp1 out_val_rpko_imp2 prex_rpko prex_rpko_imp1 prex_rpko_imp2 out_rpko ///
EKSPOR export_pct export_pct_imp revenue_total pct_own_cent_gov_imp pct_own_loc_gov_imp pct_own_nat_priv_imp pct_own_for_imp workers_total_imp3 ///
using "temp_data/processed_mill_geolocalization/IBS_UML_panel.xlsx", firstrow(variables) replace 

* save a cross-sectional selected version of only matched ibs-uml mills
keep if uml_matched_sample == 1
keep firm_id year is_mill trase_code uml_id mill_name parent_co lat lon island_factor island_name district_name kec_name village_name /// 
min_year est_year startYear max_year active industry_code ///
avg_in_tot_ton_cpo_imp2 last_in_tot_ton_cpo_imp2 avg_in_ton_ffb_imp2 last_in_ton_ffb_imp2 avg_ffb_price_imp2 ///
avg_out_ton_cpo_imp2 last_out_ton_cpo_imp2 avg_cpo_price_imp1 avg_cpo_price_imp2

duplicates drop firm_id, force
sort firm_id
save "temp_data/processed_mill_geolocalization/IBS_UML_cs.dta", replace


/*
browse firm_id year merge_overall_cfl mill_name parent_co trase_code lat lon ///
merge_oto merge_noto merge_pre2011 merge_post2010 merge_unref_geo workers_total_imp3 district_name kec_name village_name ///
min_year n_dif_id_per_mill_name name_du n_dif_id_per_coord coord_du ///
export_pct export_pct_imp revenue_total pct_own_cent_gov_imp pct_own_loc_gov_imp pct_own_nat_priv_imp pct_own_for_imp ///
 ffb_price_imp1 ffb_price_imp2 in_ton_ffb in_ton_ffb_imp1 in_ton_ffb_imp2 in_val_ffb in_val_ffb_imp1 in_val_ffb_imp2 ///
 cpo_price_imp1 cpo_price_imp2 out_ton_cpo out_ton_cpo_imp1 out_ton_cpo_imp2 out_val_cpo out_val_cpo_imp1 out_val_cpo_imp2 prex_cpo prex_cpo_imp1 prex_cpo_imp2 ///
 pko_price_imp1 pko_price_imp2 out_ton_pko out_ton_pko_imp1 out_ton_pko_imp2 out_val_pko out_val_pko_imp1 out_val_pko_imp2 prex_pko prex_pko_imp1 prex_pko_imp2 ///
 if !mi(mill_name) & mi(uml_matched_sample)
*/
