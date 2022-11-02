library(foreign)
library(dplyr)
library(stringr)

# Give here the names of the IBS variables to be added to workstream 
names_raw_vars <- c("ICOVCU")

# stock pre-processed cross sections in:
cs_list <- list()
i <- 1
for(year in 1998:2015){
  
  ibs <- read.dbf(paste0("C:/Users/GUYE/Desktop/LUCFP/data_processing/input_data/IBS_shipment/IBS_",year,".dbf"), as.is = TRUE)
  
  year_key <- str_sub(year, start = -2, end = -1)
  
  names(ibs) <- gsub(pattern = year_key, replacement = "", x = names(ibs))
  
  ibs$year <- year 
  
  # make them all the same class 
  ibs <- mutate(ibs, across(.cols = everything(), 
                            ~as.numeric(.)))
  # this throws warnings 'NAs introduits lors de la conversion automatique" which normal: just saying it recognizes <NA> in character class and converts into NA numeric class
  
  ibs <- ibs[,c("PSID", "year", names_raw_vars)]
  
  cs_list[[i]] <- ibs
  
  i <- i + 1
}

panel <- bind_rows(cs_list)











# code to check which variables have been shiped by BPS. 
for(YRno in 0:15){
  year <- 2000+YRno
  ibs <- read.dbf(paste0("C:/Users/GUYE/Desktop/LUCFP/data_processing/input_data/IBS_shipment/IBS_",year,".dbf"), as.is = TRUE)
  
  
lookfor <- c("ZPIVCU", "ZNIVCU", # bonus in cash and in kind and others
             "ZPPVCU", "ZNPVCU",  # pension allowance and insurance
             #"IINVCU", # test
             # "OUTPUT", # test
             "ICOVCU",  # gifts 
             "IRRVCU", # "royalties"
             "IOVCU", # others, among other expenses
             "YISVCU") # manufacturing services

if(YRno<10){YRno <- paste0("0",YRno)}
lookfor <- paste0(lookfor, YRno)

# lookfor[lookfor %in% names(ibs)] %>% print()
print(year)
print(nrow(ibs))
print(summary(ibs$ICOVCU))

}


