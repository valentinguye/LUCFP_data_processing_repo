*****************************************************************************************************************************************
* Master do file for the LUCFP project. R scripts are called from here, but this requires user to specify the local path to R executable
*****************************************************************************************************************************************

/////!!!!!!\\\\\/////!!!!!!\\\\\/////!!!!!!\\\\\
*** SPECIFY YOUR LOCAL PATH TO YOUR R PROGRAM 
global Rterm_path "C:\Program Files\R\R-3.6.1\bin\R.exe"
/////!!!!!!\\\\\/////!!!!!!\\\\\/////!!!!!!\\\\\

* Set permanent R options to be applied at each call of an R script from Stata. 
* vanilla is necessary for the R script call to work, but it prevents R from reading its .Rprofile file at start up 
* see https://stat.ethz.ch/R-manual/R-devel/library/base/html/Startup.html 
* and https://stackoverflow.com/questions/12540138/how-can-i-make-my-r-session-vanilla
global Rterm_options "--vanilla"


******* NOTHING IS USER SPECIFIC BELOW THIS LINE ********

*** THE PROJECT WORKING DIRECTORY DOES NOT NEED TO BE MANUALLY SPECIFIED, AS IT IS HOME OF THE PRESENT master.do, i.e. ~/LUCFP/data_processing/    
/* 
However, note that if you run the present .do file from a code editor like Sublime Text 3, 
or if you set a working directory prior to running it (in your profile.do for instance), you need to set the current working directory to 
your local absolute path to the present master.do 
*/
cd "$LUCFP"

* This is where all the processed data are stored through the following script execution. 
* Subsequent directories are generated within scripts.  
capture mkdir "temp_data"


*** PACKAGES 
** Stata Packages
ssc install rsource, replace

* Note on rsource: the working directory returned by getwd() at the start of R scripts called from rsource 
* is the current working directory in Stata. 

** R Packages
rsource using "install_R_project_packages.R"
* see the project README and/or comments within the script for more info


***** PREPARE DATA ***** 

**** explicative variables 

	*** extract palm oil mills from IBS main input/output datasets (IBS_IO) - and clean measurement unit and reshape. 

	do code/explicative_variables/IBS_inputs_preparation.do
	* input: input_data/IBS_IO/IBS_inputs    downloaded here from C:/Users/guyv/ownCloud/opal (2)/download/output/IBS_IO
	* output: temp_data/processed_IBS/prepared_IO/IBS_inputs_prep.dta

	do code/explicative_variables/IBS_outputs_preparation.do
	* input: input_data/IBS_IO/IBS_outputs   downloaded here from C:/Users/guyv/ownCloud/opal (2)/download/output/IBS_IO
	* output: temp_data/processed_IBS/prepared_IO/IBS_outputs_prep.dta

	** merge them together and with main IBS. 
	do code/explicative_variables/merge_IBSmain_IBSIO.do
	* input: temp_data/processed_IBS/prepared_IO/IBS_outputs_prep.dta
	*		 temp_data/processed_IBS/prepared_IO/IBS_inputs_prep.dta
	*		 input_data/si_panel.dta  		  (Rothenberg data) uploaded here from C:/Users/guyv/ownCloud/opal (2)/download/input/IBS_full/si_panel.dta
	* 		 input_data/IBS_final_panel.dta  (Sebastian cleaned data) uploaded here from C:/Users/guyv/ownCloud/opal (2)/build/output/IBS_final_panel.dta

	* output: temp_data/processed_IBS/IBS_PO_98_15.dta

	*** clean IBS mill dataset. 
	** some preparatory work with district and village crosswalks
	do code/explicative_variables/prepare_crosswalks.do
	* input 	input_data/indonesia_spatial/District-Proliferation-Crosswalk_complete.csv
	*           input_data/indonesia_spatial/desa_crosswalk_1998_2014.csv 

	* output:   temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta
	*		    temp_data/processed_indonesia_spatial/desa_code_names_98_2014.dta
	
	** clean IBS_PO_98_15 (including geographic variables)
	do code/explicative_variables/cleaning_IBS.do
	* input: temp_data/processed_IBS/IBS_PO_98_15.dta
	*		 input_data/IBS_kraus/IBS_final_panel.dta
	*        input_data/IBS_kraus/IBS_base_prelu.dta
	*        temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta
	*        temp_data/processed_indonesia_spatial/desa_code_names_98_2014.dta

	* output: temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta

**** mill geolocalization 

	*** automatic matching with Heilmayr mill list. 

		** keep mill obs. of most recent year with a valid desa_id
		do code/mill_geolocalization/keep_valid_recent_desa.do
		* input: temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta
		* output: temp_data/processed_mill_geolocalization/IBSmills_valid_desa.dta

		** give IBS mills a desa geometry and centroid
		rsource using "code/mill_geolocalization/make_IBSmills_desageom.R"
		* input: temp_data/processed_mill_geolocalization/IBSmills_valid_desa.dta
		* 		 input_data/indonesia_spatial/village_shapefiles/desa_1998_crosswalked.shp 
		*		 ... 
		*        input_data/indonesia_spatial/village_shapefiles/desa_2009_crosswalked.shp
		* 		 input_data/indonesia_spatial/village_shapefiles/desa_map_2010/indo_by_desa_2010.shp
		
		* output: temp_data/processed_mill_geolocalization/IBSmills_desageom.Rdata
		*		  temp_data/processed_mill_geolocalization/IBSmills_desacentro.dta

		** perform spatial matching and output distinct subsets that will need different processing
		rsource using "code/mill_geolocalization/IBS_UML_matching.R"
		* input: temp_data/processed_mill_geolocalization/IBSmills_desageom.Rdata
		*  		 input_data/uml/traseMills_capEstyear.xlsx
		
		* output: temp_data/processed_mill_geolocalization/pre2011_bad_desa_id.dta
		* 		  temp_data/processed_mill_geolocalization/unreferenced_mill_desa.shp (driver geojson)
		*		  temp_data/processed_mill_geolocalization/ibs_unref.dta
		* 		  temp_data/processed_mill_geolocalization/oto.dta
		* 		  temp_data/processed_mill_geolocalization/noto.dta


	*** manual matching / geo-localization

		** Take sub-dataset of mills with observations only since 2011 included and thus surely no desa_id.
		do code/mill_geolocalization/georeferencing_IBSmills_post2010.do
		* input: temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta	 
		
		* output: temp_data/processed_mill_geolocalization/mills_to_georeference_post2010.xls

		** Take sub-datasets of mills with at least one observation before 2011 but that cannot match automatically (no valid desa_id, or conflicting matches - "noto"). 
		do code/mill_geolocalization/georeferencing_IBSmills_pre2011.do
		* input: temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta
		* 		 temp_data/processed_mill_geolocalization/pre2011_bad_desa_id.dta
		* 		 temp_data/processed_mill_geolocalization/noto.dta

		* output: temp_data/processed_mill_geolocalization/mills_to_georeference_pre2011.xls
		* 		  temp_data/processed_mill_geolocalization/noto.xls  

		** MANUAL WORK; **NO CODE** **NOTE THAT OUTPUT IS IN READ-ONLY INPUT DATA**
		* input: temp_data/processed_mill_geolocalization/mills_to_georeference.xls
		* 		 temp_data/processed_mill_geolocalization/mills_to_georeference_pre2011.xls
		* 		 temp_data/processed_mill_geolocalization/noto.xls 

		* output:  input_data/manually_matched_ibs_uml/mills_to_georeference_post2010_done.xls
		* 		   input_data/manually_matched_ibs_uml/mills_to_georeference_pre2011_done.xlsx
		* 		   input_data/manually_matched_ibs_uml/noto_done.xls

	*** Automatic matching of remaining unreferenced mills. 
		rsource using "code/mill_geolocalization/match_unref.R"
		* input:	input_data/manufacturing_directories/direktori_industri_merged_cleaned.xlsx
		* 			temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta
		*			temp_data/processed_mill_geolocalization/ibs_unref.Rdata
		*			input_data/manually_matched_ibs_uml/matching_unref/wtn_cfl_done.xlsx
		*			input_data/manually_matched_ibs_uml/matching_unref/btw_cfl_done.xlsx

		* output:	temp_data/processed_mill_geolocalization/wtn_cfl.xlsx
		*			temp_data/processed_mill_geolocalization/btw_cfl.xlsx
		*			temp_data/processed_mill_geolocalization/named_unref.xlsx		

		* This named_unref.xlsx was sent to Jason Benedict who matched these MD names to UML names and outputed md_millMatching.xlsx
		* Some modifications were made manually to md_millMatching.xlsx, resulting in md_millMatching_modified.xlsx


	*** merge automatic and manual works - and resolve remaining conflicts.  
		do code/mill_geolocalization/merging_geolocalization_works.do
		* input: temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta
		*		 input_data/manually_matched_ibs_uml/mills_to_georeference_post2010_done.xls 	 sheet("mills to georef")
		*  	     input_data/manually_matched_ibs_uml/mills_to_georeference_pre2011_done.xlsx	 sheet("Mills to georef")
		*  	     input_data/manually_matched_ibs_uml/noto_done.xls 	sheet("Sheet1")
		*  		 input_data/uml/traseMills_capEstyear.xlsx 	sheet("traseMills_capEstyear")
		*		 input_data/uml/mills_20200129.xlsx
		*		 input_data/manually_matched_ibs_uml/matching_unref/md_millMatching_modified.xlsx
		* 		 input_data/manually_matched_ibs_uml/overall_btw_conflicts_done.xlsx
		*		 temp_data/processed_mill_geolocalization/IBSmills_desacentro.dta


		* output: temp_data/processed_mill_geolocalization/mills_to_georeference_post2010_done.dta 	 
		*  	      temp_data/processed_mill_geolocalization/mills_to_georeference_pre2011_done.dta	 
		*  	      temp_data/processed_mill_geolocalization/noto_done.dta 	
		* 		  temp_data/processed_UML/traseMills_capEstyear_modified.dta
		*		  temp_data/processed_UML/mills_20200129_modified.dta
		* 		  temp_data/processed_mill_geolocalization/matching_unref/md_millMatching_modified.dta
		* 		  temp_data/processed_mill_geolocalization/merge_geoloc_works_temp.dta

		*  	      temp_data/processed_mill_geolocalization/IBS_UML_panel.dta
		*  	      temp_data/processed_mill_geolocalization/IBS_UML_panel.xlsx
		*  		  temp_data/processed_mill_geolocalization/IBS_UML_cs.dta


**** Make additional modifications to IBS, for which IBS-UML matching was necessary. 

	** impute establishment year from est_year variable in UML and min_year variable in IBS (requires both have been matched already)
	do code/explicative_variables/impute_UML_est_year.do
	* input		input_data/uml/mills_estyear_clean.xlsx
	*			temp_data/processed_UML/traseMills_capEstyear_modified.dta
	*			temp_data/processed_UML/mills_20200129_modified.dta
	*			temp_data/processed_mill_geolocalization/IBS_UML_cs.dta
	
	* output 	temp_data/processed_UML/mills_estyear_clean_modified.dta
	*			temp_data/processed_UML/UML_valentin_imputed_est_year.dta


	** Compute concentration variable for all UML mills. 
	* /!\ 10:40 hours to be executed.
	rsource using "code/explicative_variables/mill_concentration.R"
	* input 	temp_data/processed_UML/UML_valentin_imputed_est_year.dta

	* output 	temp_data/processed_UML/UML_panel_valentin.dta

	** add international and domestic prices; compute log variables; add mill concentration and imputed establishment year variables 
	do code/explicative_variables/add_variables.do
	* input 	input_data/macro/prices_exp.xlsx		 
	*  	      	temp_data/processed_mill_geolocalization/IBS_UML_panel.dta
	*			temp_data/processed_UML/UML_panel_valentin.dta

	* output 	temp_data/processed_macro/prices_exp.dta
	*			temp_data/IBS_UML_panel_final.dta


**** Prepare Indonesian island shapes
		rsource using "code/prepare_island_sf.R"
		* input 	input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp
		*
		* output 	temp_data/processed_indonesia_spatial/island_sf


**** Build outcome variables at the parcel level

	*** prepare gfc data in a separate script to isolate use of gfcanalysis package. 
		rsource using "code/outcome_variables/prepare_gfc.R"
		* input: 	input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp
		*			INTERNET CONNEXION
		
		* output: 	temp_data/processed_lu/gfc_data_Indonesia_30th.tif 
		*			temp_data/processed_lu/gfc_data_Indonesia_60th.tif 
		*			temp_data/processed_lu/gfc_data_Indonesia_90th.tif 

	*** prepare parcel maps of annual LUC from primary forest to industrial plantations (LUCPFIP)
		rsource using "code/outcome_variables/prepare_lucpfip.R" 
		* input:	temp_data/processed_indonesia_spatial/island_sf
		* 			input_data/austin_plantation_maps/IIASA_indo_oilpalm_map/oilpalm_2015_WGS1984.tif
		*			input_data/margono_primary_forest/timeseq_change00_12.tif
		*			temp_data/IBS_UML_panel_final.dta

		* output: 
		*		Main raster outputs 
		*			temp_data/processed_lu/annual_maps/lucpfip_ISLAND_TYPE_YEAR.tif for ISLAND = ("Sumatra, Kalimantan, Papua), TYPE = (intact, degraded, total) and YEAR = (2001-2018)
		*			temp_data/processed_lu/annual_maps/parcel_lucpfip_ISLAND_PS_TYPE.tif for ISLAND = ("Sumatra, Kalimantan, Papua), PS = 3km, TYPE = (intact, degraded, total) and YEAR = (2001-2018)
		
		*		Dataframe outputs 
		*			temp_data/processed_parcels/lucpfip_panel_ISLAND_PS_CR_TYPE.rds for ISLAND = ("Sumatra, Kalimantan, Papua), PS = 3km, CR = (10km_IBS_CR, 30km_IBS_CR, 50km_IBS_CR, 10km_UML_CR, 30km_UML_CR, 50km_UML_CR) TYPE = (intact, degraded, total) and YEAR = (2001-2018)
		*			temp_data/processed_parcels/lucpfip_panel_PS_CR.rds for PS = 3km, CR = (10km_IBS_CR, 30km_IBS_CR, 50km_IBS_CR, 10km_UML_CR, 30km_UML_CR, 50km_UML_CR) : each one has rows from three islands and columns for three forest definitions.  

	*** prepare parcel maps of annual LUC from 30, 60 and 90% tree cover forest outside 2000 industrial oil palm plantations (LUCFIP)
		rsource using "code/outcome_variables/prepare_lucfip" 
		* input:	temp_data/processed_indonesia_spatial/island_sf
		*			temp_data/processed_lu/gfc_data_Indonesia_30th.tif
		*			temp_data/processed_lu/gfc_data_Indonesia_60th.tif
		*			temp_data/processed_lu/gfc_data_Indonesia_90th.tif
		* 			input_data/austin_plantation_maps/IIASA_indo_oilpalm_map/oilpalm_2000_WGS1984.tif
		*			temp_data/processed_lu/austin_ioppm_2015_Sumatra_aligned.tif
		*			temp_data/processed_lu/austin_ioppm_2015_Kalimantan_aligned.tif
		*			temp_data/processed_lu/austin_ioppm_2015_Papua_aligned.tif
		*			temp_data/IBS_UML_panel_final.dta

		* output: 
		*		Main raster outputs 
		*			temp_data/processed_lu/annual_maps/lucfip_ISLAND_TH_YEAR.tif for ISLAND = ("Sumatra", "Kalimantan", "Papua"), TH = (30th, 60th, 90th)and YEAR = (2001-2018)
		*			temp_data/processed_lu/annual_maps/parcel_lucfip_ISLAND_PS_TH.tif for ISLAND = ("Sumatra, Kalimantan, Papua), PS = 3km, TH = (30th, 60th, 90th) and YEAR = (2001-2018)
		
		*		Dataframe outputs 
		*			temp_data/processed_parcels/lucfip_panel_ISLAND_PS_CR_TH.rds for ISLAND = ("Sumatra, Kalimantan, Papua), PS = 3km, CR = (10km_IBS_CR, 30km_IBS_CR, 50km_IBS_CR, 10km_UML_CR, 30km_UML_CR, 50km_UML_CR) TH = (30, 60, 90) and YEAR = (2001-2018)
		*			temp_data/processed_parcels/lucfip_panel_PS_CR.rds for PS = 3km, CR = (10km_IBS_CR, 30km_IBS_CR, 50km_IBS_CR, 10km_UML_CR, 30km_UML_CR, 50km_UML_CR) : each one has rows from three islands and columns for three forest definitions.  

	*** prepare parcel maps of annual LUC from primary forest to small and medium size plantations (LUCPFSP) 
		rsource using "code/outcome_variables/prepare_lucpfsp"
		* input: 

		* output:	temp_data/processed_parcels/lucfip_panel_ISLAND_PS_CR_TH.rds for ISLAND = ("Sumatra, Kalimantan, Papua), PS = 3km, CR = (10km_IBS_CR, 30km_IBS_CR, 50km_IBS_CR, 10km_UML_CR, 30km_UML_CR, 50km_UML_CR) TH = (30, 60, 90) and YEAR = (2001-2018)
		*			temp_data/processed_parcels/lucfip_panel_PS_CR.rds for PS = 3km, CR = (10km_IBS_CR, 30km_IBS_CR, 50km_IBS_CR, 10km_UML_CR, 30km_UML_CR, 50km_UML_CR) : each one has rows from three islands and columns for three forest definitions.  

	*** prepare driving travel time between every pairs of parcels and mills. 
	capture "code/programs/make_osrm_durations"
	/*  
	It is not sourced here, because it requires a local instance of OSRM to be set up. See within the script for more details. 
	For this reason, the output being not easily reproducible, it is saved in input_data, and not temp_data. 
	
	This is run by island (SUMATRA AND KALIMANTAN ONLY BECAUSE NOT WORKING IN PAPUA),
	for any existing mill up to 2015 (geolocalized IBS or UML) and for annually existing mills (geolocalized IBS only so far). 
	*/


	*** Program that selects outcome variable raster grid cells into data frames, based on OSRM travel time constraints - for either IBS or UML mills. 
	/* This corresponds to the 3rd parts of prepare_* scripts, that selected grid cells based on a catchment radius (i.e. distance) constraint. 
	But unlike the code in these prepare_* scripts, here only total primary forest is kept (and 30th forest in preparation) and this script gathers 
	parcels from industrial, small and medium sized plantations.  	
	*/
	rsource using "code/programs/make_osrm_CA"
	* input: 		temp_data/processed_indonesia_spatial/island_sf
	*				temp_data/IBS_UML_panel_final.dta
	*				input_data/uml/mills_20200129.xlsx
	*				temp_data/processed_lu/parcel_lucpfip_ISLAND_PS_TYPE for ISLAND = c("Sumatra", "Kalimantan", "Papua"); PS = 3km, TYPE = c("intact", "degraded", "total", "30th" "60th", "90th")
	*				input_data/local_osrm_outputs/osrm_driving_durations_ISLAND_PS_SAMPLE for ISLAND = ("Sumatra", "Kalimantan") PS = 3km and SAMPLE = (IBS, UML)

	* output:		temp_data/processed_parcels/lucpfSIZEp_panel_ISLAND_PS_CA.rds for SIZE = ("s", "m", "i"); ISLAND = ("Sumatra", "Kalimantan"); PS = 3km, CA = (2h_IBS_CA, 4h_IBS_CA, 6h_IBS_CA, 2h_UML_CA, 4h_UML_CA, 6h_UML_CA)  : each one is for total primary forest only
	*				temp_data/processed_parcels/lucpfp_panel_PS_CA.rds for PS = 3km, CA = (2h_IBS_CA, 4h_IBS_CA, 6h_IBS_CA, 2h_UML_CA, 4h_UML_CA, 6h_UML_CA) : each one is for total primary forest only, and has rows from two islands and columns for three plantation sizes  
	*				temp_data/processed_parcels/lucpfip_panel_DYNA_ISLAND_PS_CA.rds for DYNA = ("replace", "rapid", "slow"); ISLAND = ("Sumatra", "Kalimantan"); PS = 3km, CA = (2h_IBS_CA, 4h_IBS_CA, 6h_IBS_CA, 2h_UML_CA, 4h_UML_CA, 6h_UML_CA)  : each one is for total primary forest only
	*				temp_data/processed_parcels/lucpfip_panel_dynamics_ISLAND_PS_CA.rds for ISLAND = ("Sumatra", "Kalimantan"); PS = 3km, CA = (2h_IBS_CA, 4h_IBS_CA, 6h_IBS_CA, 2h_UML_CA, 4h_UML_CA, 6h_UML_CA)  : each one is for total primary forest only



	*** prepare parcel maps and dataframes of baseline forest extents according to different definitions
	* This script outputs these baseline forest covers for parcels within CR and within CA (OSRM duration based). 
	rsource using "code/outcome_variables/prepare_2000_forest_extents.R"
	* input: 	temp_data/processed_indonesia_spatial/island_sf
	*			temp_data/processed_lu/gfc_data_Indonesia_TH.tif for TH = (30th, 60th, 90th)
	*			temp_data/processed_lu/austin_ioppm_2000_ISLAND_aligned.tif for ISLAND = ("Sumatra", "Kalimantan", "Papua")
	*			temp_data/processed_lu/margono_primary_forest_ISLAND_aligned.tif for ISLAND = ("Sumatra", "Kalimantan", "Papua")
	*			temp_data/processed_mill_geolocalization/IBS_UML_panel.dta
	*			input_data/uml/mills_20200129.xlsx
	*			input_data/local_osrm_outputs/osrm_driving_durations_ISLAND_PS_SAMPLE for ISLAND = ("Sumatra", "Kalimantan") PS = 3km and SAMPLE = (IBS, UML)

	* output: 	temp_data/processed_lu/gfc_fc2000_ISLAND_TH.tif for ISLAND = ("Sumatra", "Kalimantan", "Papua") and TH = (30th, 60th, 90th)
	*			temp_data/processed_lu/gfc_fc2000_ISLAND_TH_prj.tif for ISLAND = ("Sumatra", "Kalimantan", "Papua") and TH = (30th, 60th, 90th)
	*			temp_data/processed_lu/gfc_fc2000_outside_ip_ISLAND_TH.tif for ISLAND = ("Sumatra", "Kalimantan", "Papua") and TH = (30th, 60th, 90th)
	*			temp_data/processed_lu/margono_TYPE_primary_forest_ISLAND_aligned.tif for TYPE = c("intact", "degraded", "total") ISLAND = ("Sumatra", "Kalimantan", "Papua") 
	*			temp_data/processed_lu/parcel_fc2000_ISLAND_30th_PS.tif for ISLAND = ("Sumatra", "Kalimantan", "Papua") and PS = 3km
	*			temp_data/processed_lu/parcel_pfc2000_ISLAND_total_PS.tif for ISLAND = ("Sumatra", "Kalimantan", "Papua") and PS = 3km
	*			temp_data/processed_lu/parcel_fc2000_ISLAND_30th_PS_SAMPLE_masked.tif for ISLAND = ("Sumatra", "Kalimantan", "Papua"), PS = 3km and SAMPLE = (IBS, UML)
	*			temp_data/processed_lu/parcel_pfc2000_ISLAND_total_PS_SAMPLE_masked.tif for ISLAND = ("Sumatra", "Kalimantan", "Papua"), PS = 3km and SAMPLE = (IBS, UML)
	*			temp_data/processed_lu/fc2000_cs_ISLAND_30th_PS_CR.tif for ISLAND = ("Sumatra", "Kalimantan", "Papua"), PS = 3km and CR = (10km_IBS_CR, 30km_IBS_CR, 50km_IBS_CR, 10km_UML_CR, 30km_UML_CR, 50km_UML_CR)
	*			temp_data/processed_lu/pfc2000_cs_ISLAND_total_PS_CR.tif for ISLAND = ("Sumatra", "Kalimantan", "Papua"), PS = 3km and CR = (10km_IBS_CR, 30km_IBS_CR, 50km_IBS_CR, 10km_UML_CR, 30km_UML_CR, 50km_UML_CR)
	* 		main output: 
	*			temp_data/processed_parcels/baseline_fc_cs_PS_CR.rds for PS = 3km and CR = (10km_IBS_CR, 30km_IBS_CR, 50km_IBS_CR, 10km_UML_CR, 30km_UML_CR, 50km_UML_CR)
	*			temp_data/processed_parcels/baseline_fc_cs_PS_CA.rds for PS = 3km and CA = (2h_IBS_CA, 4h_IBS_CA, 6h_IBS_CA, 2h_UML_CA, 4h_UML_CA, 6h_UML_CA)

	*** prepare emissions from LUCPFIP 
		rsource using "code/outcome_variables/prepare_emissions.R"




**** build explicative variables at the parcel level

	/* Distribute geolocalized IBS mill variables to square parcels with distance-weighted averages.  
	 Parcels are grouped in dataframes according to their size (3x3km only currently - i.e. PS = 3000 (meters)) and 
	 how many they are (what is the catchment radius: CR = 10km, 30km, 50km).  
	
	Note that prepare_OUTCOME_VAR scripts outputted dataframes for catchment radii of both IBS and UML mill sets, and hence this is specified in file names. 
	However, in further RHS-related scripts, we always use parcels within IBS catchment radii, and hence do not specify it in file names.  
	 */	

	 * /!\ ~ 61 hours to execute 
	rsource using "code/explicative_variables/wa_at_parcels_distances.R"
		* input   	temp_data/IBS_UML_panel_final.dta
		*			temp_data/processed_parcels/lucpfip_panel_PS_CR.rds for PS = 3km, CR = (10km_IBS_CR, 30km_IBS_CR, 50km_IBS_CR)

		* output 	temp_data/processed_parcels/temp_cs_wa_explanatory/cs_wa_explanatory_PS_CR_YEAR.rds for PS = 3km, CR = (10CR, 30CR, 50CR) and YEAR = (1998-2015)
		*			temp_data/processed_parcels/wa_panel_parcels_PS_CR.rds for PS = 3km, CR = (10CR, 30CR, 50CR)



	/* Or ditribute to the same square parcels but within a travel time catchment AREA and with DURATION-weighted averages. 
	 It is different than wa_at_parcels_distances.R because it runs the analysis by ISLAND. 
	 And because it uses an additional input: the duration matrices annually computed by a program (make_osrm_durations.R) that requires a local OSRM server. 
	 Contrary to its "distance" counterparts, it requires the parcel template (lucpfip_panel typically) the ISLAND LEVEL. 
	 The end of the script stacks the two islands together. 
	 */
	rsource using "code/explicative_variables/wa_at_parcels_durations.R"
	* input 	temp_data/IBS_UML_panel_final.dta
	*			temp_data/processed_parcels/lucpfip_panel_ISLAND_PS_CA_total.rds for ISLAND = (Sumatra, Kalimantan); PS = 3km, CA = (2h_CA, 4h_CA, 6h_CA)
	*			input_data/local_osrm_outputs/osrm_driving_durations_ISLAND_PS_TT_YEAR for ISLAND = (Sumatra, Kalimantan); PS = 3km, TT = (2h_IBS, 4h_IBS, 6h_IBS); YEAR = 1998:2015

	* output 	temp_data/processed_parcels/temp_cs_wa_explanatory/cs_wa_explanatory_ISLAND_PS_CA_YEAR.rds for ISLAND = (Sumatra, Kalimantan); PS = 3km, CA = (2h_CA, 4h_CA, 6h_CA) and YEAR = (1998-2015)
	*			temp_data/processed_parcels/wa_panel_parcels_ISLAND_PS_CA.rds for PS = 3km, CA = (2h_CA, 4h_CA, 6h_CA)
	*			temp_data/processed_parcels/wa_panel_parcels_PS_CA.rds for PS = 3km, CA = (2h_CA, 4h_CA, 6h_CA) 





	/* add information on the number of UML mills that are reachable from each parcel in each year within 10, 30 and 50km. 
	  and on the share of IBS sample in this total number of reachable UML mills. 
	  AND add island and district variables
	  AND legal or certified land use attributes
	  AND add time series of international and domestic prices and export tax and spreads. 
	  /!\ ~1h */
	rsource using "code/explicative_variables/add_CR_parcel_variables.R"	
		* input:	temp_data/processed_UML/UML_valentin_imputed_est_year.dta
		*			temp_data/processed_parcels/wa_panel_parcels_PS_CR.rds  	 for PS = 3km, CR = (10CR, 30CR, 50CR)
		*			temp_data/processed_indonesia_spatial/island_sf
		*			input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp
		*			input_data/indonesia_spatial/district_shapefiles/district_2015_base2000.shp
		*			temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta
		*			temp_data/processed_parcels/baseline_fc_cs_PS_CR.rds for PS = 3km, CR = (10km_IBS_CR, 30km_IBS_CR, 50km_IBS_CR)
		*			input_data/RSPO_supply_bases/RSPO-certified_oil_palm_supply_bases_in_Indonesia.shp
		*			input_data/oil_palm_concessions
		*			input_data/kawasan_hutan/Greenorb_Blog/final/KH-INDON-Final.shp

		*output: 	temp_data/processed_parcels/parcels_panel_reachable_uml_PS_CR.rds for PS = 3km, CR = (10CR, 30CR, 50CR)
		*			temp_data/processed_parcels/parcels_panel_geovars_PS_CR.rds for PS = 3km, CR = (10CR, 30CR, 50CR)
		*			temp_data/processed_parcels/parcels_panel_w_dyn_PS_CR.rds   for PS = 3km, CR = (10CR, 30CR, 50CR)
		*			temp_data/processed_parcels/parcels_panel_final_PS_CR.rds   for PS = 3km, CR = (10CR, 30CR, 50CR)

	/* Same as above, but on the duration-based CA parcels. 
	In particular, this differs from the CR counterpart because this one does not compute n_reachable_uml nor sample_coverage variables.
	*/
	rsource using "code/explicative_variables/add_CA_parcel_variables.R"
		* input:	input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp
		*			input_data/indonesia_spatial/district_shapefiles/district_2015_base2000.shp
		*			temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta
		*			input_data/RSPO_supply_bases/RSPO-certified_oil_palm_supply_bases_in_Indonesia.shp
		*			input_data/oil_palm_concessions
		*			input_data/kawasan_hutan/Greenorb_Blog/final/KH-INDON-Final.shp
		*   	 	temp_data/processed_parcels/wa_panel_parcels_PS_CA.rds for PS = 3km, CA = (2h_CA, 4h_CA, 6h_CA) 

		* output: 	temp_data/processed_parcels/parcels_panel_geovars_PS_CA for PS = 3km, CA = (2h_CA, 4h_CA, 6h_CA) 
		*			temp_data/processed_parcels/parcels_panel_w_dyn_PS_CA for PS = 3km, CA = (2h_CA, 4h_CA, 6h_CA) 
		*			temp_data/processed_parcels/parcels_panel_final_PS_CA for PS = 3km, CA = (2h_CA, 4h_CA, 6h_CA) 



**** Merge LHS and RHS (!)
	* Add also baseline forest extents variables computed in prepare_2000_forest_extents.R
		rsource using "code/merge_lhs_rhs_CR_parcels.R"
		* input: 	temp_data/processed_parcels/OV_panel_PS_CR.rds 		for OV = (lucpfip, lucfip, ... ); PS = 3km, CR = (10km_IBS_CR, 30km_IBS_CR, 50km_IBS_CR)
		*			temp_data/processed_parcels/parcels_panel_final_PS_CR.rds for PS = 3km, CR = (10CR, 30CR, 50CR)
		*			temp_data/processed_parcels/baseline_fc_cs_PS_CR.rds for PS = 3km and CR = (10km_IBS_CR, 30km_IBS_CR, 50km_IBS_CR)


		* output: 	temp_data/panel_parcels_ip_final_PS_CR.rds  for PS = 3km, CR = (10CR, 30CR, 50CR)

	* And same as above, but for duration-based CA parcels. 
		rsource using "code/merge_lhs_rhs_CA_parcels.R"
		* input: 	input/processed_parcels/lucpfp_panel_PS_CA.rds for PS = 3km, CA = (2h_IBS_CA, 4h_IBS_CA, 6h_IBS_CA) : each one is for total primary forest only, and has rows from two islands and columns for three plantation sizes  
		*			temp_data/processed_parcels/parcels_panel_final_PS_CA for PS = 3km, CA = (2h_CA, 4h_CA, 6h_CA) 
		* 			temp_data/processed_parcels/baseline_fc_cs_PS_CA.rds for PS = 3km and CA = (2h_IBS_CA, 4h_IBS_CA, 6h_IBS_CA)

		* output: 	temp_data/panel_parcels_ip_final_PS_CA.rds for PS = 3km and CA = (2h_IBS_CA, 4h_IBS_CA, 6h_IBS_CA)



***** ANALYSIS *****
* note: we do the conversion from pixel counts to areas in here, 
* depending on what we need. This is not done in earlier stages, so it will have to be done before analysis in further scripts too. 
**** Descriptive statistics

	*** Outcome variables
	rsource using "code/descriptive_statistics/lucfp_des_stats.R"
		* input 	temp_data/processed_mill_geolocalization/IBS_UML_panel.dta
		*			input_data/uml/mills_20200129.xlsx
		*			temp_data/processed_lu/gfc_data_Indonesia_30th.tif, gfc_data_Indonesia_60th.tif, gfc_data_Indonesia_90th.tif 
		*			temp_data/processed_lu/austin_ioppm_2000_",island,"_aligned.tif
		*			temp_data/processed_lu/margono_primary_forest_island_aligned.tif

		*			temp_data/processed_lu/parcel_lucfip_ISLAND_PS_TH.tif for ISLAND = ("Sumatra", "Kalimantan", "Papua") PS = 3km, TH = (30th, 60th, 90th)
		*			temp_data/processed_lu/parcel_lucpfip_ISLAND_PS_TH.tif for ISLAND = ("Sumatra", "Kalimantan", "Papua") PS = 3km, TH = (30th, 60th, 90th)

		*			input_data/GEE_outputs 	/!\ these data come from the execution of code at 
		*			with inputs from this script's parts: ##### Prepare polygons of interest ##### ; 
		*												  #### Prepare 30, 60, 90% forest cover outside industrial plantations in 2000 ####
		*											      #### Prepare intact, degraded and total primary forest cover in 2000 ####


		* Main output 	temp_data/processed_indonesia_spatial/mill_cr_prj_ISLAND
		*				temp_data/processed_lu/gfc_fc2000_outside_ip_ISLAND_TH.tif
		*				temp_data/processed_lu/margono_TYPE_primary_forest_ISLAND_aligned.tif	

		*				LateX forest extent 2000 and LUCFIP descriptive statistics table codes. 


	*** Explicative variables
	rsource using "code/descriptive_statistics/rhs_tables.R"
	* input 	temp_data/IBS_UML_panel_final.dta
	*			temp_data/panel_parcels_ip_final_PS_CR.rds for PS = 3km and CR = (10CR, 30CR, 50CR)

	* output 	LateX IBS and grid cell descriptive statistics table codes 

	rsource using "code/descriptive_statistics/rhs_maps.R"
	* input 	temp_data/processed_parcels/wa_panel_parcels_reachable_uml_PS_CR.rds for PS = 3km and CR = (10CR, 30CR, 50CR)
	*			temp_data/IBS_UML_panel_final.dta

	* output 	outputs/figures/ 		(for exemple accu_lucpfip, tm_n_reachable_uml, tm_cpo_price_imp1 ... ) 









