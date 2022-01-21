library(raster)
library(sp)
library(tidyverse)
library(ncdf4)
library(lubridate)
library(geosphere)

#' USA rice irrigated
lat <- 34.5852
lon <- -91.7517
lat_m <- 3827752.872
lon_m <- 614484.527

#' Set flux footprint extent
US_HRA_xmin <- lon_m-150
US_HRA_xmax <- lon_m+150
US_HRA_ymin <- lat_m-150
US_HRA_ymax <- lat_m+150
US_HRA_FFP <- extent(c(US_HRA_xmin, US_HRA_xmax,
                       US_HRA_ymin, US_HRA_ymax))

#' Read NDVI files
NDVI_dir_list <- list.dirs("E:\\USArice\\PySEBAL_data\\SEBAL_out\\")
NDVI_dir_list <- grep(pattern="vegetation", x=NDVI_dir_list, value=TRUE)
NDVI_list <- list.files(NDVI_dir_list, pattern="_NDVI", all.files=TRUE, full.names = TRUE)
NDVI_stack <- stack(NDVI_list)

#' Read Tact files
Tact_dir_list <- list.dirs("E:\\USArice\\PySEBAL_data\\SEBAL_out\\")
Tact_dir_list <- grep(pattern="evapotranspiration", x=Tact_dir_list, value=TRUE)
Tact_list <- list.files(Tact_dir_list, pattern="_Tact", all.files=TRUE, full.names = TRUE)
Tact_stack <- stack(Tact_list)

#' Read Tpot files
Tpot_dir_list <- list.dirs("E:\\USArice\\PySEBAL_data\\SEBAL_out\\")
Tpot_dir_list <- grep(pattern="evapotranspiration", x=Tpot_dir_list, value=TRUE)
Tpot_list <- list.files(Tpot_dir_list, pattern="_Tpot", all.files=TRUE, full.names = TRUE)
Tpot_stack <- stack(Tpot_list)

#' Read LAI files
LAI_dir_list <- list.dirs("E:\\USArice\\PySEBAL_data\\SEBAL_out\\")
LAI_dir_list <- grep(pattern="vegetation", x=LAI_dir_list, value=TRUE)
LAI_list <- list.files(LAI_dir_list, pattern="_lai", all.files=TRUE, full.names = TRUE)
LAI_stack <- stack(LAI_list)

#' Read SAVI files
SAVI_dir_list <- list.dirs("E:\\USArice\\PySEBAL_data\\SEBAL_out\\")
SAVI_dir_list <- grep(pattern="vegetation", x=LAI_dir_list, value=TRUE)
SAVI_list <- list.files(SAVI_dir_list, pattern="_SAVI", all.files=TRUE, full.names = TRUE)
SAVI_stack <- stack(SAVI_list)

#' Read biomass files
biomass_dir_list <- list.dirs("E:\\USArice\\PySEBAL_data\\SEBAL_out\\")
biomass_dir_list <- grep(pattern="biomass", x=biomass_dir_list, value=TRUE)
biomass_list <- list.files(biomass_dir_list, pattern="_Biomass_production", all.files=TRUE, full.names = TRUE)
biomass_stack <- stack(biomass_list)

#' Extract NDVI vals from stack
NDVI_soy_FFP <- crop(NDVI_stack, US_HRA_FFP)
NDVI_vals <- cellStats(NDVI_soy_FFP, stat=mean, na.rm=TRUE)


#' Calculate fAPAR
fPAR_vals <- -0.161 + 1.257*NDVI_vals
fPAR_vals <- as_tibble_col(matrix(unlist(fPAR_vals)))

#' Calculate moisture stress biomass
Tact_soy_FFP <- crop(Tact_stack, US_HRA_FFP)
Tact_vals <- cellStats(Tact_soy_FFP, stat=mean, na.rm=TRUE)
Tpot_soy_FFP <- crop(Tpot_stack, US_HRA_FFP)
Tpot_vals <- cellStats(Tpot_soy_FFP, stat=mean, na.rm=TRUE)

Moisture_stress_vals <- Tact_vals/Tpot_vals
Moisture_stress_vals <- as_tibble_col(matrix(unlist(Moisture_stress_vals)))

#' Get tower data and supp. RH data
US_HRA_daily_data <- read.csv("data/US_HRA.csv")

Rh_dir_list <- list.dirs("E:\\USArice\\PySEBAL_data\\Meteo\\")
Rh_list <- list.files(Rh_dir_list, pattern="_Rh_avg", all.files=TRUE, full.names = TRUE)
Rh_stack_list <- lapply(Rh_list, raster)

Rh_stack <- stack(Rh_stack_list)
Rh_stack <- projectRaster(Rh_stack, crs=crs(NDVI_stack))
Rh_stack <- crop(Rh_stack, US_HRA_FFP)

Rh_vals <- cellStats(Rh_stack, stat=mean, na.rm=TRUE)
Rh_dates <- as.character(as.POSIXct(str_sub(names(Rh_stack),2,10), format=("%Y_%m%d")))
Rh_df <- bind_cols(Rh_dates, Rh_vals)
names(Rh_df) <- c("date", "Rh")

#' Get tower values
US_HRA_daily_data$Date <- as.character(as.Date(as.character(US_HRA_daily_data$TIMESTAMP), format="%Y%m%d"))
US_HRA_daily_df <- US_HRA_daily_data %>% 
  dplyr::select(Date, TA_F, SW_IN_F, LW_IN_F, GPP_DT) #Add RH if available

names(US_HRA_daily_df) <- c("date", "daily_Ta", "daily_Fsd", "daily_Fld",
                            "daily_GPP")

########################## Calculate thermal time ##############################
US_HRA_daily_df$dateD <- as.Date(US_HRA_daily_df$date)
US_HRA_daily_df$DOY <-  as.numeric(strftime(US_HRA_daily_df$dateD, format = "%j"))

US_HRA_daily_df <- US_HRA_daily_df %>% 
  mutate_at(vars(dateD), funs(year, month, day)) %>% 
  mutate(DAS = if_else(DOY >= 105 & DOY <= 305, DOY-105,0))

USA_daily_GDD <- US_HRA_daily_df %>%
  filter(DAS > 0) %>% 
  group_by(year) %>% 
  mutate(GDD = if_else(DAS >0, cumsum(daily_Ta-10), 0)) %>% 
  ungroup()

USA_daily_GDD <- USA_daily_GDD %>% 
  dplyr::select("dateD", "GDD")

US_HRA_daily_df <- left_join(US_HRA_daily_df, USA_daily_GDD, by="dateD")

#' Get times for Landsat dates
US_HRA_time <- str_sub(names(Tpot_stack),16,25)
US_HRA_time <- as_tibble_col(US_HRA_time)
US_HRA_time <- as.data.frame(US_HRA_time)
US_HRA_time$Day_times <- as.POSIXct(US_HRA_time$value, format=("%Y_%m_%d"))
US_HRA_time$date <- as.character(as.Date(US_HRA_time$Day_times+86400))
names(US_HRA_time) <- c("value", "time", "date")

#' Join tower dates with Landsat dates

US_HRA_towervals <- US_HRA_time  %>% 
  left_join(US_HRA_daily_df, by="date") %>% 
  left_join(Rh_df, by="date")

#' Define constants
Th <- 35
Kt <- 23
Tl <- 0
rl <- 130
Jarvis_coeff <- (Th-Kt)/(Kt-Tl)

#' Calculate PAR
US_HRA_towervals$PAR <- US_HRA_towervals$daily_Fsd*0.48*0.0864
US_HRA_towervals$fPAR <- fPAR_vals$value
US_HRA_towervals$moisture_stress <- Moisture_stress_vals$value

#' Calculate vapor stress biomass
US_HRA_towervals$esat <- 0.6108*exp(17.27*US_HRA_towervals$daily_Ta/
                                      (US_HRA_towervals$daily_Ta+237.3))
US_HRA_towervals$eact <- US_HRA_towervals$Rh*US_HRA_towervals$esat/100

US_HRA_towervals$vapor_stress <- if (0.88-0.183*log(US_HRA_towervals$esat-US_HRA_towervals$eact) >1){
  1
} else {
  0.88-0.183*log(US_HRA_towervals$esat-US_HRA_towervals$eact)
}

#' Calculate heat stress biomass
US_HRA_towervals$heat_stress <- ((US_HRA_towervals$daily_Ta-Tl)*
                                   (Th-US_HRA_towervals$daily_Ta)^Jarvis_coeff)/
  ((Kt-Tl)*((Th-Kt)^Jarvis_coeff))

#' Calculate LUE
US_HRA_towervals$APAR <- US_HRA_towervals$fPAR*US_HRA_towervals$PAR
US_HRA_towervals$LUE <- US_HRA_towervals$daily_GPP/
  (US_HRA_towervals$APAR)
US_HRA_towervals$LUEmax <- US_HRA_towervals$LUE/
  (US_HRA_towervals$heat_stress*US_HRA_towervals$vapor_stress*US_HRA_towervals$moisture_stress)

#' Get SEBAL GPP
biomass_FFP <- crop(biomass_stack, US_HRA_FFP)
biomass_vals <- cellStats(biomass_FFP, stat=mean, na.rm=TRUE)
biomass_vals <- as_tibble_col(biomass_vals)

LAI_FFP <- crop(LAI_stack, US_HRA_FFP)
LAI_vals <- cellStats(LAI_FFP, stat=mean, na.rm=TRUE)
LAI_vals <- as_tibble_col(LAI_vals)

SAVI_FPP <- crop(SAVI_stack, US_HRA_FFP)
SAVI_vals <- cellStats(SAVI_FPP, stat=mean, na.rm=TRUE)
SAVI_vals <- as_tibble_col(SAVI_vals)

#DOY <- as_tibble_col(str_sub(names(NDVI_stack),24,26))
Site <- as_tibble_col("US_HRA")

US_HRA_MLvals <- bind_cols(Site, US_HRA_towervals, biomass_vals,
                           NDVI_vals, LAI_vals, SAVI_vals)
names(US_HRA_MLvals)[1:2] <- c("Site", "origin_date")
names(US_HRA_MLvals)[27:30] <- c("biomass", "NDVI_vals", "LAI_vals",
                                 "SAVI_vals")
#US_HRA_MLvals$DOY <- as.numeric(US_HRA_MLvals$DOY)

#'Add latitude
US_HRA_MLvals <- US_HRA_MLvals %>% 
  mutate(Lat = lat) %>% 
  mutate(Dlength = daylength(Lat, DOY))

US_HRA_MLvals$GPP_pred <- US_HRA_MLvals$biomass/10

ggplot(data= US_HRA_MLvals, aes(x=daily_GPP, y=GPP_pred))+geom_point()+
  geom_abline(intercept = 0, slope=1)
US_HRA_daily_df$dateD <- as.Date(US_HRA_daily_df$date)
ggplot(data=US_HRA_daily_df, aes(x=dateD, y=daily_GPP))+geom_point()+
  scale_x_date(date_labels="%Y%m", date_breaks = "6 months")
ggplot(data=US_HRA_MLvals, aes(x=date, y=daily_GPP))+geom_point()
US_HRA_MLvals <- US_HRA_MLvals %>% 
  mutate_at(vars(date), funs(year, month, day))

ggplot(data= US_HRA_MLvals, aes(x=daily_GPP, y=GPP_pred, col=Crop))+geom_point()+
  geom_abline(intercept = 0, slope=1)

US_HRA_MLvals <- US_HRA_MLvals %>%   
  mutate(Crop = "Rice")

ggplot(data= US_HRA_MLvals, aes(x=daily_GPP, y=GPP_pred, col=Crop))+geom_point()+
  geom_abline(intercept = 0, slope=1)

############################## Landsat 7 ######################################

#' Calculate and add GCVI = (NIR/Green) - 1
#' LE07 = (B4 / B2 ) -1
#' LC08 = (B5 / B3) -1
US_HRA_LE07_dir_list <- list.dirs("E:\\USArice\\PySEBAL_data\\Satellite_data\\gapfilled")

#' Get metadata to convert to TOA reflectance
US_HRA_LE07_MTL_list <- grep(pattern="LE07", x=US_HRA_LE07_dir_list, value=TRUE)
US_HRA_LE07_MTL_list <- list.files(US_HRA_LE07_MTL_list, pattern="MTL", all.files=TRUE, full.names = TRUE)

#' Write a for loop to store the conversion values for each date (row) in a column for each band
US_HRA_LE07_conv_vals <- data.frame(
  B1_MULT_val=rep(NA,length(US_HRA_LE07_MTL_list)), B1_ADD_val=rep(NA,length(US_HRA_LE07_MTL_list)),
  B2_MULT_val=rep(NA,length(US_HRA_LE07_MTL_list)), B2_ADD_val=rep(NA,length(US_HRA_LE07_MTL_list)),
  B3_MULT_val=rep(NA,length(US_HRA_LE07_MTL_list)), B3_ADD_val=rep(NA,length(US_HRA_LE07_MTL_list)),
  B4_MULT_val=rep(NA,length(US_HRA_LE07_MTL_list)), B4_ADD_val=rep(NA,length(US_HRA_LE07_MTL_list)),
  B5_MULT_val=rep(NA,length(US_HRA_LE07_MTL_list)), B5_ADD_val=rep(NA,length(US_HRA_LE07_MTL_list)),
  B7_MULT_val=rep(NA,length(US_HRA_LE07_MTL_list)), B7_ADD_val=rep(NA,length(US_HRA_LE07_MTL_list)))

for (i in 1:length(US_HRA_LE07_MTL_list)) {
  f <- readLines(US_HRA_LE07_MTL_list[i])
  
  B1_MULT <- grep("REFLECTANCE_MULT_BAND_1", f, value=TRUE)
  US_HRA_LE07_conv_vals$B1_MULT_val[i] <- as.numeric(sub(".*=", "", B1_MULT))
  B1_ADD <- grep("REFLECTANCE_ADD_BAND_1", f, value=TRUE)
  US_HRA_LE07_conv_vals$B1_ADD_val[i] <- as.numeric(sub(".*=", "", B1_ADD))
  
  B2_MULT <- grep("REFLECTANCE_MULT_BAND_2", f, value=TRUE)
  US_HRA_LE07_conv_vals$B2_MULT_val[i] <- as.numeric(sub(".*=", "", B2_MULT))
  B2_ADD <- grep("REFLECTANCE_ADD_BAND_2", f, value=TRUE)
  US_HRA_LE07_conv_vals$B2_ADD_val[i] <- as.numeric(sub(".*=", "", B2_ADD))
  
  B3_MULT <- grep("REFLECTANCE_MULT_BAND_3", f, value=TRUE)
  US_HRA_LE07_conv_vals$B3_MULT_val[i] <- as.numeric(sub(".*=", "", B3_MULT))
  B3_ADD <- grep("REFLECTANCE_ADD_BAND_3", f, value=TRUE)
  US_HRA_LE07_conv_vals$B3_ADD_val[i] <- as.numeric(sub(".*=", "", B3_ADD))
  
  B4_MULT <- grep("REFLECTANCE_MULT_BAND_4", f, value=TRUE)
  US_HRA_LE07_conv_vals$B4_MULT_val[i] <- as.numeric(sub(".*=", "", B4_MULT))
  B4_ADD <- grep("REFLECTANCE_ADD_BAND_4", f, value=TRUE)
  US_HRA_LE07_conv_vals$B4_ADD_val[i] <- as.numeric(sub(".*=", "", B4_ADD))
  
  B5_MULT <- grep("REFLECTANCE_MULT_BAND_5", f, value=TRUE)
  US_HRA_LE07_conv_vals$B5_MULT_val[i] <- as.numeric(sub(".*=", "", B5_MULT))
  B5_ADD <- grep("REFLECTANCE_ADD_BAND_5", f, value=TRUE)
  US_HRA_LE07_conv_vals$B5_ADD_val[i] <- as.numeric(sub(".*=", "", B5_ADD))
  
  B7_MULT <- grep("REFLECTANCE_MULT_BAND_7", f, value=TRUE)
  US_HRA_LE07_conv_vals$B7_MULT_val[i] <- as.numeric(sub(".*=", "", B7_MULT))
  B7_ADD <- grep("REFLECTANCE_ADD_BAND_7", f, value=TRUE)
  US_HRA_LE07_conv_vals$B7_ADD_val[i] <- as.numeric(sub(".*=", "", B7_ADD))
  
}

US_HRA_LE07_B4_file_list <- list.files(US_HRA_LE07_dir_list, pattern="_B4", all.files=TRUE, full.names = TRUE)
US_HRA_LE07_B4_rast <- lapply(US_HRA_LE07_B4_file_list, raster) 
US_HRA_LE07_B4_crop <- lapply(X= US_HRA_LE07_B4_rast, FUN=crop, y=US_HRA_FFP)
US_HRA_LE07_B4_stack <- stack(US_HRA_LE07_B4_crop)

US_HRA_LE07_B2_file_list <- list.files(US_HRA_LE07_dir_list, pattern="_B2", all.files=TRUE, full.names = TRUE)
US_HRA_LE07_B2_rast <- lapply(US_HRA_LE07_B2_file_list, raster) 
US_HRA_LE07_B2_crop <- lapply(X= US_HRA_LE07_B2_rast, FUN=crop, y=US_HRA_FFP)
US_HRA_LE07_B2_stack <- stack(US_HRA_LE07_B2_crop)

US_HRA_LE07_GCVI <- (US_HRA_LE07_B4_stack/US_HRA_LE07_B2_stack) -1
US_HRA_LE07_GCVI <- as_tibble_col(cellStats(US_HRA_LE07_GCVI, stat=mean, na.rm=TRUE))

#' Get times for Landsat 7 dates
USA_time <- str_sub(US_HRA_LE07_B2_file_list,67,74)
USA_time <- as_tibble_col(USA_time)
USA_time <- as.data.frame(USA_time)
USA_time$Day_times <- as.POSIXct(USA_time$value, format=("%Y%m%d"))
USA_time$date <- as.character(as.Date(USA_time$Day_times+86400))
names(USA_time) <- c("value", "time", "date")
US_HRA_LE07_GCVI <- bind_cols(USA_time$date, US_HRA_LE07_GCVI$value)
names(US_HRA_LE07_GCVI) <- c("date", "GCVI")

#' Get blue band = B1
US_HRA_LE07_B1_file_list <- list.files(US_HRA_LE07_dir_list, pattern="_B1", all.files=TRUE, full.names = TRUE)
US_HRA_LE07_B1_rast <- lapply(US_HRA_LE07_B1_file_list, raster) 
US_HRA_LE07_B1_crop <- lapply(X= US_HRA_LE07_B1_rast, FUN=crop, y=US_HRA_FFP)
US_HRA_LE07_B1_stack <- stack(US_HRA_LE07_B1_crop)
US_HRA_LE07_blue <- as_tibble_col(cellStats(US_HRA_LE07_B1_stack, stat=mean, na.rm=TRUE))

#' Get green band = B2
US_HRA_LE07_B2_file_list <- list.files(US_HRA_LE07_dir_list, pattern="_B2", all.files=TRUE, full.names = TRUE)
US_HRA_LE07_B2_rast <- lapply(US_HRA_LE07_B2_file_list, raster) 
US_HRA_LE07_B2_crop <- lapply(X= US_HRA_LE07_B2_rast, FUN=crop, y=US_HRA_FFP)
US_HRA_LE07_B2_stack <- stack(US_HRA_LE07_B2_crop)
US_HRA_LE07_green <- as_tibble_col(cellStats(US_HRA_LE07_B2_stack, stat=mean, na.rm=TRUE))

#' Get red band = B3
US_HRA_LE07_B3_file_list <- list.files(US_HRA_LE07_dir_list, pattern="_B3", all.files=TRUE, full.names = TRUE)
US_HRA_LE07_B3_rast <- lapply(US_HRA_LE07_B3_file_list, raster) 
US_HRA_LE07_B3_crop <- lapply(X= US_HRA_LE07_B3_rast, FUN=crop, y=US_HRA_FFP)
US_HRA_LE07_B3_stack <- stack(US_HRA_LE07_B3_crop)
US_HRA_LE07_red <- as_tibble_col(cellStats(US_HRA_LE07_B3_stack, stat=mean, na.rm=TRUE))

#' Get nir band = B4
US_HRA_LE07_B4_file_list <- list.files(US_HRA_LE07_dir_list, pattern="_B4", all.files=TRUE, full.names = TRUE)
US_HRA_LE07_B4_rast <- lapply(US_HRA_LE07_B4_file_list, raster) 
US_HRA_LE07_B4_crop <- lapply(X= US_HRA_LE07_B4_rast, FUN=crop, y=US_HRA_FFP)
US_HRA_LE07_B4_stack <- stack(US_HRA_LE07_B4_crop)
US_HRA_LE07_nir <- as_tibble_col(cellStats(US_HRA_LE07_B4_stack, stat=mean, na.rm=TRUE))

#' Get swir1 band = B5
US_HRA_LE07_B5_file_list <- list.files(US_HRA_LE07_dir_list, pattern="_B5", all.files=TRUE, full.names = TRUE)
US_HRA_LE07_B5_rast <- lapply(US_HRA_LE07_B5_file_list, raster) 
US_HRA_LE07_B5_crop <- lapply(X= US_HRA_LE07_B5_rast, FUN=crop, y=US_HRA_FFP)
US_HRA_LE07_B5_stack <- stack(US_HRA_LE07_B5_crop)
US_HRA_LE07_swir1 <- as_tibble_col(cellStats(US_HRA_LE07_B5_stack, stat=mean, na.rm=TRUE))

#' Get swir2 band = B7
US_HRA_LE07_B7_file_list <- list.files(US_HRA_LE07_dir_list, pattern="_B7", all.files=TRUE, full.names = TRUE)
US_HRA_LE07_B7_rast <- lapply(US_HRA_LE07_B7_file_list, raster) 
US_HRA_LE07_B7_crop <- lapply(X= US_HRA_LE07_B7_rast, FUN=crop, y=US_HRA_FFP)
US_HRA_LE07_B7_stack <- stack(US_HRA_LE07_B7_crop)
US_HRA_LE07_swir2 <- as_tibble_col(cellStats(US_HRA_LE07_B7_stack, stat=mean, na.rm=TRUE))

#' Convert to reflectance
US_HRA_LE07_blue$refl <- US_HRA_LE07_blue$value* US_HRA_LE07_conv_vals$B1_MULT_val+
  US_HRA_LE07_conv_vals$B1_ADD_val

US_HRA_LE07_green$refl <- US_HRA_LE07_green$value* US_HRA_LE07_conv_vals$B2_MULT_val+
  US_HRA_LE07_conv_vals$B2_ADD_val

US_HRA_LE07_red$refl <- US_HRA_LE07_red$value* US_HRA_LE07_conv_vals$B3_MULT_val+
  US_HRA_LE07_conv_vals$B3_ADD_val

US_HRA_LE07_nir$refl <- US_HRA_LE07_nir$value* US_HRA_LE07_conv_vals$B4_MULT_val+
  US_HRA_LE07_conv_vals$B4_ADD_val

US_HRA_LE07_swir1$refl <- US_HRA_LE07_swir1$value* US_HRA_LE07_conv_vals$B5_MULT_val+
  US_HRA_LE07_conv_vals$B5_ADD_val

US_HRA_LE07_swir2$refl <- US_HRA_LE07_swir2$value* US_HRA_LE07_conv_vals$B7_MULT_val+
  US_HRA_LE07_conv_vals$B7_ADD_val

#' Combine band data with GCVI
US_HRA_LE07_dat <- bind_cols(US_HRA_LE07_GCVI, US_HRA_LE07_blue$refl, US_HRA_LE07_green$refl,
                              US_HRA_LE07_red$refl, US_HRA_LE07_nir$refl, US_HRA_LE07_swir1$refl,
                              US_HRA_LE07_swir2$refl)

names(US_HRA_LE07_dat) <- c("date", "GCVI", "blue", "green", "red", "nir", "swir1", "swir2")

#' Add sensor
US_HRA_LE07_dat <- US_HRA_LE07_dat %>% 
  mutate(sensor = "Ls7")


######################## Repeat for Landsat 8 #################################
US_HRA_LC08_dir_list <- list.dirs("E:\\USArice\\PySEBAL_data\\Satellite_data\\")

#' Get metadata to convert to TOA reflectance
US_HRA_LC08_MTL_list <- grep(pattern="LC08", x=US_HRA_LC08_dir_list, value=TRUE)
US_HRA_LC08_MTL_list <- list.files(US_HRA_LC08_MTL_list, pattern="MTL", all.files=TRUE, full.names = TRUE)

#' Write a for loop to store the conversion values for each date (row) in a column for each band
US_HRA_LC08_conv_vals <- data.frame(
  B2_MULT_val=rep(NA,length(US_HRA_LC08_MTL_list)), B2_ADD_val=rep(NA,length(US_HRA_LC08_MTL_list)),
  B3_MULT_val=rep(NA,length(US_HRA_LC08_MTL_list)), B3_ADD_val=rep(NA,length(US_HRA_LC08_MTL_list)),
  B4_MULT_val=rep(NA,length(US_HRA_LC08_MTL_list)), B4_ADD_val=rep(NA,length(US_HRA_LC08_MTL_list)),
  B5_MULT_val=rep(NA,length(US_HRA_LC08_MTL_list)), B5_ADD_val=rep(NA,length(US_HRA_LC08_MTL_list)),
  B6_MULT_val=rep(NA,length(US_HRA_LC08_MTL_list)), B6_ADD_val=rep(NA,length(US_HRA_LC08_MTL_list)),
  B7_MULT_val=rep(NA,length(US_HRA_LC08_MTL_list)), B7_ADD_val=rep(NA,length(US_HRA_LC08_MTL_list)))

for (i in 1:length(US_HRA_LC08_MTL_list)) {
  f <- readLines(US_HRA_LC08_MTL_list[i])
  B2_MULT <- grep("REFLECTANCE_MULT_BAND_2", f, value=TRUE)
  US_HRA_LC08_conv_vals$B2_MULT_val[i] <- as.numeric(sub(".*=", "", B2_MULT))
  B2_ADD <- grep("REFLECTANCE_ADD_BAND_2", f, value=TRUE)
  US_HRA_LC08_conv_vals$B2_ADD_val[i] <- as.numeric(sub(".*=", "", B2_ADD))
  
  B3_MULT <- grep("REFLECTANCE_MULT_BAND_3", f, value=TRUE)
  US_HRA_LC08_conv_vals$B3_MULT_val[i] <- as.numeric(sub(".*=", "", B3_MULT))
  B3_ADD <- grep("REFLECTANCE_ADD_BAND_3", f, value=TRUE)
  US_HRA_LC08_conv_vals$B3_ADD_val[i] <- as.numeric(sub(".*=", "", B3_ADD))
  
  B4_MULT <- grep("REFLECTANCE_MULT_BAND_4", f, value=TRUE)
  US_HRA_LC08_conv_vals$B4_MULT_val[i] <- as.numeric(sub(".*=", "", B4_MULT))
  B4_ADD <- grep("REFLECTANCE_ADD_BAND_4", f, value=TRUE)
  US_HRA_LC08_conv_vals$B4_ADD_val[i] <- as.numeric(sub(".*=", "", B4_ADD))
  
  B5_MULT <- grep("REFLECTANCE_MULT_BAND_5", f, value=TRUE)
  US_HRA_LC08_conv_vals$B5_MULT_val[i] <- as.numeric(sub(".*=", "", B5_MULT))
  B5_ADD <- grep("REFLECTANCE_ADD_BAND_5", f, value=TRUE)
  US_HRA_LC08_conv_vals$B5_ADD_val[i] <- as.numeric(sub(".*=", "", B5_ADD))
  
  B6_MULT <- grep("REFLECTANCE_MULT_BAND_6", f, value=TRUE)
  US_HRA_LC08_conv_vals$B6_MULT_val[i] <- as.numeric(sub(".*=", "", B6_MULT))
  B6_ADD <- grep("REFLECTANCE_ADD_BAND_6", f, value=TRUE)
  US_HRA_LC08_conv_vals$B6_ADD_val[i] <- as.numeric(sub(".*=", "", B6_ADD))
  
  B7_MULT <- grep("REFLECTANCE_MULT_BAND_7", f, value=TRUE)
  US_HRA_LC08_conv_vals$B7_MULT_val[i] <- as.numeric(sub(".*=", "", B7_MULT))
  B7_ADD <- grep("REFLECTANCE_ADD_BAND_7", f, value=TRUE)
  US_HRA_LC08_conv_vals$B7_ADD_val[i] <- as.numeric(sub(".*=", "", B7_ADD))
  
}

#' Get blue band = B2
US_HRA_LC08_B2_file_list <- grep(pattern="LC08", x=US_HRA_LC08_dir_list, value=TRUE)
US_HRA_LC08_B2_file_list <- list.files(US_HRA_LC08_B2_file_list, pattern="B2", all.files=TRUE, full.names = TRUE)
US_HRA_LC08_B2_rast <- lapply(US_HRA_LC08_B2_file_list, raster) 
US_HRA_LC08_B2_crop <- lapply(X= US_HRA_LC08_B2_rast, FUN=crop, y=US_HRA_FFP)
US_HRA_LC08_B2_stack <- stack(US_HRA_LC08_B2_crop)
US_HRA_LC08_blue <- as_tibble_col(cellStats(US_HRA_LC08_B2_stack, stat=mean, na.rm=TRUE))

#' Get green band = B3
US_HRA_LC08_B3_file_list <- grep(pattern="LC08", x=US_HRA_LC08_dir_list, value=TRUE)
US_HRA_LC08_B3_file_list <- list.files(US_HRA_LC08_B3_file_list, pattern="B3", all.files=TRUE, full.names = TRUE)
US_HRA_LC08_B3_rast <- lapply(US_HRA_LC08_B3_file_list, raster) 
US_HRA_LC08_B3_crop <- lapply(X= US_HRA_LC08_B3_rast, FUN=crop, y=US_HRA_FFP)
US_HRA_LC08_B3_stack <- stack(US_HRA_LC08_B3_crop)
US_HRA_LC08_green <- as_tibble_col(cellStats(US_HRA_LC08_B3_stack, stat=mean, na.rm=TRUE))

#' Get red band = B4
US_HRA_LC08_B4_file_list <- grep(pattern="LC08", x=US_HRA_LC08_dir_list, value=TRUE)
US_HRA_LC08_B4_file_list <- list.files(US_HRA_LC08_B4_file_list, pattern="B4", all.files=TRUE, full.names = TRUE)
US_HRA_LC08_B4_rast <- lapply(US_HRA_LC08_B4_file_list, raster) 
US_HRA_LC08_B4_crop <- lapply(X= US_HRA_LC08_B4_rast, FUN=crop, y=US_HRA_FFP)
US_HRA_LC08_B4_stack <- stack(US_HRA_LC08_B4_crop)
US_HRA_LC08_red <- as_tibble_col(cellStats(US_HRA_LC08_B4_stack, stat=mean, na.rm=TRUE))

#' Get nir band = B5
US_HRA_LC08_B5_file_list <- grep(pattern="LC08", x=US_HRA_LC08_dir_list, value=TRUE)
US_HRA_LC08_B5_file_list <- list.files(US_HRA_LC08_B5_file_list, pattern="B5", all.files=TRUE, full.names = TRUE)
US_HRA_LC08_B5_rast <- lapply(US_HRA_LC08_B5_file_list, raster) 
US_HRA_LC08_B5_crop <- lapply(X= US_HRA_LC08_B5_rast, FUN=crop, y=US_HRA_FFP)
US_HRA_LC08_B5_stack <- stack(US_HRA_LC08_B5_crop)
US_HRA_LC08_nir <- as_tibble_col(cellStats(US_HRA_LC08_B5_stack, stat=mean, na.rm=TRUE))

#' Get swir1 band = B6
US_HRA_LC08_B6_file_list <- grep(pattern="LC08", x=US_HRA_LC08_dir_list, value=TRUE)
US_HRA_LC08_B6_file_list <- list.files(US_HRA_LC08_B6_file_list, pattern="B6", all.files=TRUE, full.names = TRUE)
US_HRA_LC08_B6_rast <- lapply(US_HRA_LC08_B6_file_list, raster) 
US_HRA_LC08_B6_crop <- lapply(X= US_HRA_LC08_B6_rast, FUN=crop, y=US_HRA_FFP)
US_HRA_LC08_B6_stack <- stack(US_HRA_LC08_B6_crop)
US_HRA_LC08_swir1 <- as_tibble_col(cellStats(US_HRA_LC08_B6_stack, stat=mean, na.rm=TRUE))

#' Get swir2 band = B7
US_HRA_LC08_B7_file_list <- grep(pattern="LC08", x=US_HRA_LC08_dir_list, value=TRUE)
US_HRA_LC08_B7_file_list <- list.files(US_HRA_LC08_B7_file_list, pattern="B7", all.files=TRUE, full.names = TRUE)
US_HRA_LC08_B7_rast <- lapply(US_HRA_LC08_B7_file_list, raster) 
US_HRA_LC08_B7_crop <- lapply(X= US_HRA_LC08_B7_rast, FUN=crop, y=US_HRA_FFP)
US_HRA_LC08_B7_stack <- stack(US_HRA_LC08_B7_crop)
US_HRA_LC08_swir2 <- as_tibble_col(cellStats(US_HRA_LC08_B7_stack, stat=mean, na.rm=TRUE))

#' Convert to reflectance
US_HRA_LC08_blue$refl <- US_HRA_LC08_blue$value* US_HRA_LC08_conv_vals$B2_MULT_val+
  US_HRA_LC08_conv_vals$B2_ADD_val

US_HRA_LC08_green$refl <- US_HRA_LC08_green$value* US_HRA_LC08_conv_vals$B3_MULT_val+
  US_HRA_LC08_conv_vals$B3_ADD_val

US_HRA_LC08_red$refl <- US_HRA_LC08_red$value* US_HRA_LC08_conv_vals$B4_MULT_val+
  US_HRA_LC08_conv_vals$B4_ADD_val

US_HRA_LC08_nir$refl <- US_HRA_LC08_nir$value* US_HRA_LC08_conv_vals$B5_MULT_val+
  US_HRA_LC08_conv_vals$B5_ADD_val

US_HRA_LC08_swir1$refl <- US_HRA_LC08_swir1$value* US_HRA_LC08_conv_vals$B6_MULT_val+
  US_HRA_LC08_conv_vals$B6_ADD_val

US_HRA_LC08_swir2$refl <- US_HRA_LC08_swir2$value* US_HRA_LC08_conv_vals$B7_MULT_val+
  US_HRA_LC08_conv_vals$B7_ADD_val

US_HRA_LC08_GCVI <- (US_HRA_LC08_B4_stack/US_HRA_LC08_B2_stack) -1
US_HRA_LC08_GCVI <- as_tibble_col(cellStats(US_HRA_LC08_GCVI, stat=mean, na.rm=TRUE))


#' Get times for Landsat 8 dates
USA_time <- str_sub(US_HRA_LC08_B2_file_list,58,65)
USA_time <- as_tibble_col(USA_time)
USA_time <- as.data.frame(USA_time)
USA_time$Day_times <- as.POSIXct(USA_time$value, format=("%Y%m%d"))
USA_time$date <- as.character(as.Date(USA_time$Day_times+86400))
names(USA_time) <- c("value", "time", "date")
US_HRA_LC08_GCVI <- bind_cols(USA_time$date, US_HRA_LC08_GCVI$value)
names(US_HRA_LC08_GCVI) <- c("date", "GCVI")

#' Combine band data with GCVI
US_HRA_LC08_dat <- bind_cols(US_HRA_LC08_GCVI, US_HRA_LC08_blue$refl, US_HRA_LC08_green$refl,
                              US_HRA_LC08_red$refl, US_HRA_LC08_nir$refl, US_HRA_LC08_swir1$refl,
                              US_HRA_LC08_swir2$refl)

names(US_HRA_LC08_dat) <- c("date", "GCVI", "blue", "green", "red", "nir", "swir1", "swir2")

#' Add sensor
US_HRA_LC08_dat <- US_HRA_LC08_dat %>% 
  mutate(sensor = "Ls8")

########################## Combine two satellite datasets #######################
US_HRA_sat_dat <- bind_rows(US_HRA_LE07_dat, US_HRA_LC08_dat)

#' Join tower dates with Landsat dates

US_HRA_MLvals <- US_HRA_MLvals  %>% 
  left_join(US_HRA_sat_dat, by="date")


#' Export data
write.csv(US_HRA_MLvals, "MLdata/US_HRA_MLvals.csv")

