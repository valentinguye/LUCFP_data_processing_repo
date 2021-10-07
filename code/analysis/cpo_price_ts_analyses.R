RHS_50 <-  readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_w_dyn_",
                                    parcel_size/1000,"km_",
                                    "50CR.rds")))


prices <- RHS_50[,c("lonlat", "year", "cpo_price_imp1")]

nrow(prices)/length(unique(prices$year))

sample_ids <- prices$lonlat %>% sample(size = 700)
prices <- prices[prices$lonlat %in% sample_ids,]

fp <- prices[,c("lonlat", "year")]

length_panel <- length(unique(prices$year))
fp <- mutate(fp, year = year + length_panel*3)
fp$cpo_price_imp1 <- NA
prices <- rbind(prices, fp)
prices <- dplyr::arrange(prices, lonlat, year)
prices <- prices[,c("cpo_price_imp1")]
prices_ts <- ts(prices)

pacf_14 <- pacf(prices_ts, na.action = na.pass, lag.max = 14)
acf_14 <- acf(prices_ts, na.action = na.exclude, lag.max = 14)



# there should be no differencing, as the trend over the time serie means nothing, because of how it's constructed. 
# yes, there could be differencing, because differencing is not removing the mean of the full serie 
# but removing the value in the previous period. 
arima_800 <- arima(prices_ts, order = c(8,0,0)) 
arima_400 <- arima(prices_ts, order = c(4,0,0)) 

arima_008 <- arima(prices_ts, order = c(0,0,8)) 
arima_004 <- arima(prices_ts, order = c(0,0,4)) 

arima_808 <- arima(prices_ts, order = c(8,0,8)) 
arima_404 <- arima(prices_ts, order = c(4,0,4)) 

# as in Sukiyono
arima_104 <- arima(prices_ts, order = c(1,0,4)) 


arima_test <- arima(prices_ts, order = c(6,0,6)) 


ar_test7 <- ar(prices_ts, order.max = 7, na.action = na.pass)
