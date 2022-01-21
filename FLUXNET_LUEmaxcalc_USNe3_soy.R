library(raster)
library(sp)
library(tidyverse)
library(ncdf4)
library(lubridate)
library(geosphere)

#' USA maize_cont_irrigated
lat <- 41.1797
lon <- -96.4397
lat_m <- 4561866.086
lon_m <- 714748.908

#' Set flux footprint extent
USA_Ne3_xmin <- lon_m-150
USA_Ne3_xmax <- lon_m+150
USA_Ne3_ymin <- lat_m-150
USA_Ne3_ymax <- lat_m+150
USA_Ne3_FFP <- extent(c(USA_Ne3_xmin, USA_Ne3_xmax,
                          USA_Ne3_ymin, USA_Ne3_ymax))

#' Read NDVI files
NDVI_dir_list <- list.dirs("E:\\Omaha\\PySEBAL_data\\SEBAL_out_soy\\")
NDVI_dir_list <- grep(pattern="vegetation", x=NDVI_dir_list, value=TRUE)
NDVI_list <- list.files(NDVI_dir_list, pattern="_NDVI", all.files=TRUE, full.names = TRUE)
NDVI_stack <- stack(NDVI_list)

#' Read Tact files
Tact_dir_list <- list.dirs("E:\\Omaha\\PySEBAL_data\\SEBAL_out_soy\\")
Tact_dir_list <- grep(pattern="evapotranspiration", x=Tact_dir_list, value=TRUE)
Tact_list <- list.files(Tact_dir_list, pattern="_Tact", all.files=TRUE, full.names = TRUE)
Tact_stack <- stack(Tact_list)

#' Read Tpot files
Tpot_dir_list <- list.dirs("E:\\Omaha\\PySEBAL_data\\SEBAL_out_soy\\")
Tpot_dir_list <- grep(pattern="evapotranspiration", x=Tpot_dir_list, value=TRUE)
Tpot_list <- list.files(Tpot_dir_list, pattern="_Tpot", all.files=TRUE, full.names = TRUE)
Tpot_stack <- stack(Tpot_list)

#' Read LAI files
LAI_dir_list <- list.dirs("E:\\Omaha\\PySEBAL_data\\SEBAL_out_soy\\")
LAI_dir_list <- grep(pattern="vegetation", x=LAI_dir_list, value=TRUE)
LAI_list <- list.files(LAI_dir_list, pattern="_lai", all.files=TRUE, full.names = TRUE)
LAI_stack <- stack(LAI_list)

#' Read SAVI files
SAVI_dir_list <- list.dirs("E:\\Omaha\\PySEBAL_data\\SEBAL_out_soy\\")
SAVI_dir_list <- grep(pattern="vegetation", x=LAI_dir_list, value=TRUE)
SAVI_list <- list.files(SAVI_dir_list, pattern="_SAVI", all.files=TRUE, full.names = TRUE)
SAVI_stack <- stack(SAVI_list)

#' Read biomass files
biomass_dir_list <- list.dirs("E:\\Omaha\\PySEBAL_data\\SEBAL_out_soy\\")
biomass_dir_list <- grep(pattern="biomass", x=biomass_dir_list, value=TRUE)
biomass_list <- list.files(biomass_dir_list, pattern="_Biomass_production", all.files=TRUE, full.names = TRUE)
biomass_stack <- stack(biomass_list)

#' Extract NDVI vals from stack
NDVI_soy_FFP <- crop(NDVI_stack, USA_Ne3_FFP)
NDVI_vals <- cellStats(NDVI_soy_FFP, stat=mean, na.rm=TRUE)


#' Calculate fAPAR
fPAR_vals <- -0.161 + 1.257*NDVI_vals
fPAR_vals <- as_tibble_col(matrix(unlist(fPAR_vals)))

#' Calculate moisture stress biomass
Tact_soy_FFP <- crop(Tact_stack, USA_Ne3_FFP)
Tact_vals <- cellStats(Tact_soy_FFP, stat=mean, na.rm=TRUE)
Tpot_soy_FFP <- crop(Tpot_stack, USA_Ne3_FFP)
Tpot_vals <- cellStats(Tpot_soy_FFP, stat=mean, na.rm=TRUE)

Moisture_stress_vals <- Tact_vals/Tpot_vals
Moisture_stress_vals <- as_tibble_col(matrix(unlist(Moisture_stress_vals)))

#' Get tower data and supp. RH data
USA_daily_data <- read.csv("data/USA_Ne3.csv")

Rh_dir_list <- list.dirs("E:\\Omaha\\PySEBAL_data\\Meteo\\")
Rh_list <- list.files(Rh_dir_list, pattern="_Rh_avg", all.files=TRUE, full.names = TRUE)
Rh_stack <- lapply(Rh_list, raster)

Rh_stack <- stack(Rh_stack)
Rh_stack <- projectRaster(Rh_stack, crs=crs(NDVI_stack))
Rh_stack <- crop(Rh_stack, USA_Ne3_FFP)

Rh_vals <- cellStats(Rh_stack, stat=mean, na.rm=TRUE)
Rh_dates <- as.character(as.POSIXct(str_sub(names(Rh_stack),2,10), format=("%Y_%m%d")))
Rh_df <- bind_cols(Rh_dates, Rh_vals)
names(Rh_df) <- c("date", "Rh")

#' Get tower values
USA_daily_data$Date <- as.character(as.Date(as.character(USA_daily_data$TIMESTAMP), format="%Y%m%d"))
USA_daily_df <- USA_daily_data %>% 
  dplyr::select(Date, TA_F, SW_IN_F, LW_IN_F, GPP_DT_VUT_MEAN) #Add RH if available

names(USA_daily_df) <- c("date", "daily_Ta", "daily_Fsd", "daily_Fld",
                          "daily_GPP")

#' Get times for Landsat dates
USA_time <- str_sub(names(Tpot_stack),16,25)
USA_time <- as_tibble_col(USA_time)
USA_time <- as.data.frame(USA_time)
USA_time$Day_times <- as.POSIXct(USA_time$value, format=("%Y_%m_%d"))
USA_time$date <- as.character(as.Date(USA_time$Day_times+86400))
names(USA_time) <- c("value", "time", "date")

########################## Calculate thermal time ##############################
USA_daily_df$dateD <- as.Date(USA_daily_df$date)
USA_daily_df$DOY <-  as.numeric(strftime(USA_daily_df$dateD, format = "%j"))

USA_daily_df <- USA_daily_df %>% 
  mutate_at(vars(dateD), funs(year, month, day)) %>% 
  mutate(DAS = if_else(DOY >= 135 & DOY <= 305, DOY-135,0))

USA_daily_GDD <- USA_daily_df %>%
  filter(DAS > 0) %>% 
  group_by(year) %>% 
  mutate(GDD = if_else(DAS >0, cumsum(daily_Ta-10), 0)) %>% 
  ungroup()

USA_daily_GDD <- USA_daily_GDD %>% 
  dplyr::select("dateD", "GDD")

USA_daily_df <- left_join(USA_daily_df, USA_daily_GDD, by="dateD")

#' Join tower dates with Landsat dates
USA_towervals <- USA_time  %>% 
  left_join(USA_daily_df, by="date") %>% 
  left_join(Rh_df, by="date")

#' Define constants
Th <- 35
Kt <- 23
Tl <- 0
rl <- 130
Jarvis_coeff <- (Th-Kt)/(Kt-Tl)

#' Calculate PAR
USA_towervals$PAR <- USA_towervals$daily_Fsd*0.48*0.0864
USA_towervals$fPAR <- fPAR_vals$value
USA_towervals$moisture_stress <- Moisture_stress_vals$value

#' Calculate vapor stress biomass
USA_towervals$esat <- 0.6108*exp(17.27*USA_towervals$daily_Ta/
                               (USA_towervals$daily_Ta+237.3))
USA_towervals$eact <- USA_towervals$Rh*USA_towervals$esat/100

USA_towervals$vapor_stress <- if (0.88-0.183*log(USA_towervals$esat-USA_towervals$eact) >1){
  1
} else {
  0.88-0.183*log(USA_towervals$esat-USA_towervals$eact)
}

#' Calculate heat stress biomass
USA_towervals$heat_stress <- ((USA_towervals$daily_Ta-Tl)*
                            (Th-USA_towervals$daily_Ta)^Jarvis_coeff)/
  ((Kt-Tl)*((Th-Kt)^Jarvis_coeff))

#' Calculate LUE
USA_towervals$APAR <- USA_towervals$fPAR*USA_towervals$PAR
USA_towervals$LUE <- USA_towervals$daily_GPP/
  (USA_towervals$APAR)
USA_towervals$LUEmax <- USA_towervals$LUE/
  (USA_towervals$heat_stress*USA_towervals$vapor_stress*USA_towervals$moisture_stress)

#' Get SEBAL GPP
biomass_FFP <- crop(biomass_stack, USA_Ne3_FFP)
biomass_vals <- cellStats(biomass_FFP, stat=mean, na.rm=TRUE)
biomass_vals <- as_tibble_col(biomass_vals)

LAI_FFP <- crop(LAI_stack, USA_Ne3_FFP)
LAI_vals <- cellStats(LAI_FFP, stat=mean, na.rm=TRUE)
LAI_vals <- as_tibble_col(LAI_vals)

SAVI_FPP <- crop(SAVI_stack, USA_Ne3_FFP)
SAVI_vals <- cellStats(SAVI_FPP, stat=mean, na.rm=TRUE)
SAVI_vals <- as_tibble_col(SAVI_vals)

#DOY <- as_tibble_col(str_sub(names(NDVI_stack),24,26))
Site <- as_tibble_col("USA_Ne3")

USA_Ne3_MLvals <- bind_cols(Site, USA_towervals, biomass_vals,
                              NDVI_vals, LAI_vals, SAVI_vals)
names(USA_Ne3_MLvals)[1:2] <- c("Site", "origin_date")
names(USA_Ne3_MLvals)[27:30] <- c("biomass", "NDVI_vals", "LAI_vals",
                                    "SAVI_vals")
#USA_Ne3_MLvals$DOY <- as.numeric(USA_Ne3_MLvals$DOY)

#'Add latitude
USA_Ne3_MLvals <- USA_Ne3_MLvals %>% 
  mutate(Lat = lat) %>% 
  mutate(Dlength = daylength(Lat, DOY))

USA_Ne3_MLvals$GPP_pred <- USA_Ne3_MLvals$biomass/10

ggplot(data= USA_Ne3_MLvals, aes(x=daily_GPP, y=GPP_pred))+geom_point()+
  geom_abline(intercept = 0, slope=1)
USA_daily_df$dateD <- as.Date(USA_daily_df$date)
ggplot(data=USA_daily_df, aes(x=dateD, y=daily_GPP))+geom_point()+
  scale_x_date(date_labels="%Y%m", date_breaks = "6 months")
ggplot(data=USA_Ne3_MLvals, aes(x=date, y=daily_GPP))+geom_point()
USA_Ne3_MLvals <- USA_Ne3_MLvals %>% 
 mutate_at(vars(date), funs(year, month, day)) %>% 
 filter(year %in% c(2001,2002,2003,2004,2005,2006,2007,2008,2009, 2010, 2011, 2012))
USA_Ne3_MLvals <- USA_Ne3_MLvals %>%   
 mutate(Crop = if_else(year %in% c("2001","2003","2005","2007","2009", "2011"),
                        "Maize", "Soy"))

ggplot(data= USA_Ne3_MLvals, aes(x=daily_GPP, y=GPP_pred, col=Crop))+geom_point()+
  geom_abline(intercept = 0, slope=1)

USA_Ne3_MLvals <- USA_Ne3_MLvals %>%   
  filter(Crop == "Soy")
#USA_Ne3_MLvals <- USA_Ne3_MLvals %>%   
#  mutate(Crop = "Maize")

ggplot(data= USA_Ne3_MLvals, aes(x=daily_GPP, y=GPP_pred, col=Crop))+geom_point()+
  geom_abline(intercept = 0, slope=1)

############################## Landsat 7 ######################################

#' Calculate and add GCVI = (NIR/Green) - 1
#' LE07 = (B4 / B2 ) -1
#' LC08 = (B5 / B3) -1
USA_Ne3_LE07_dir_list <- list.dirs("E:\\Omaha\\PySEBAL_data\\Satellite_data\\gapfilled")

#' Get metadata to convert to TOA reflectance
USA_Ne3_LE07_MTL_list <- grep(pattern="LE07", x=USA_Ne3_LE07_dir_list, value=TRUE)
USA_Ne3_LE07_MTL_list <- list.files(USA_Ne3_LE07_MTL_list, pattern="MTL", all.files=TRUE, full.names = TRUE)

#' Write a for loop to store the conversion values for each date (row) in a column for each band
USA_Ne3_LE07_conv_vals <- data.frame(
  B1_MULT_val=rep(NA,length(USA_Ne3_LE07_MTL_list)), B1_ADD_val=rep(NA,length(USA_Ne3_LE07_MTL_list)),
  B2_MULT_val=rep(NA,length(USA_Ne3_LE07_MTL_list)), B2_ADD_val=rep(NA,length(USA_Ne3_LE07_MTL_list)),
  B3_MULT_val=rep(NA,length(USA_Ne3_LE07_MTL_list)), B3_ADD_val=rep(NA,length(USA_Ne3_LE07_MTL_list)),
  B4_MULT_val=rep(NA,length(USA_Ne3_LE07_MTL_list)), B4_ADD_val=rep(NA,length(USA_Ne3_LE07_MTL_list)),
  B5_MULT_val=rep(NA,length(USA_Ne3_LE07_MTL_list)), B5_ADD_val=rep(NA,length(USA_Ne3_LE07_MTL_list)),
  B7_MULT_val=rep(NA,length(USA_Ne3_LE07_MTL_list)), B7_ADD_val=rep(NA,length(USA_Ne3_LE07_MTL_list)))

for (i in 1:length(USA_Ne3_LE07_MTL_list)) {
  f <- readLines(USA_Ne3_LE07_MTL_list[i])
  
  B1_MULT <- grep("REFLECTANCE_MULT_BAND_1", f, value=TRUE)
  USA_Ne3_LE07_conv_vals$B1_MULT_val[i] <- as.numeric(sub(".*=", "", B1_MULT))
  B1_ADD <- grep("REFLECTANCE_ADD_BAND_1", f, value=TRUE)
  USA_Ne3_LE07_conv_vals$B1_ADD_val[i] <- as.numeric(sub(".*=", "", B1_ADD))
  
  B2_MULT <- grep("REFLECTANCE_MULT_BAND_2", f, value=TRUE)
  USA_Ne3_LE07_conv_vals$B2_MULT_val[i] <- as.numeric(sub(".*=", "", B2_MULT))
  B2_ADD <- grep("REFLECTANCE_ADD_BAND_2", f, value=TRUE)
  USA_Ne3_LE07_conv_vals$B2_ADD_val[i] <- as.numeric(sub(".*=", "", B2_ADD))
  
  B3_MULT <- grep("REFLECTANCE_MULT_BAND_3", f, value=TRUE)
  USA_Ne3_LE07_conv_vals$B3_MULT_val[i] <- as.numeric(sub(".*=", "", B3_MULT))
  B3_ADD <- grep("REFLECTANCE_ADD_BAND_3", f, value=TRUE)
  USA_Ne3_LE07_conv_vals$B3_ADD_val[i] <- as.numeric(sub(".*=", "", B3_ADD))
  
  B4_MULT <- grep("REFLECTANCE_MULT_BAND_4", f, value=TRUE)
  USA_Ne3_LE07_conv_vals$B4_MULT_val[i] <- as.numeric(sub(".*=", "", B4_MULT))
  B4_ADD <- grep("REFLECTANCE_ADD_BAND_4", f, value=TRUE)
  USA_Ne3_LE07_conv_vals$B4_ADD_val[i] <- as.numeric(sub(".*=", "", B4_ADD))
  
  B5_MULT <- grep("REFLECTANCE_MULT_BAND_5", f, value=TRUE)
  USA_Ne3_LE07_conv_vals$B5_MULT_val[i] <- as.numeric(sub(".*=", "", B5_MULT))
  B5_ADD <- grep("REFLECTANCE_ADD_BAND_5", f, value=TRUE)
  USA_Ne3_LE07_conv_vals$B5_ADD_val[i] <- as.numeric(sub(".*=", "", B5_ADD))
  
  B7_MULT <- grep("REFLECTANCE_MULT_BAND_7", f, value=TRUE)
  USA_Ne3_LE07_conv_vals$B7_MULT_val[i] <- as.numeric(sub(".*=", "", B7_MULT))
  B7_ADD <- grep("REFLECTANCE_ADD_BAND_7", f, value=TRUE)
  USA_Ne3_LE07_conv_vals$B7_ADD_val[i] <- as.numeric(sub(".*=", "", B7_ADD))
  
}

USA_Ne3_LE07_B4_file_list <- list.files(USA_Ne3_LE07_dir_list, pattern="_B4", all.files=TRUE, full.names = TRUE)
USA_Ne3_LE07_B4_rast <- lapply(USA_Ne3_LE07_B4_file_list, raster) 
USA_Ne3_LE07_B4_crop <- lapply(X= USA_Ne3_LE07_B4_rast, FUN=crop, y=USA_Ne3_FFP)
USA_Ne3_LE07_B4_stack <- stack(USA_Ne3_LE07_B4_crop)

USA_Ne3_LE07_B2_file_list <- list.files(USA_Ne3_LE07_dir_list, pattern="_B2", all.files=TRUE, full.names = TRUE)
USA_Ne3_LE07_B2_rast <- lapply(USA_Ne3_LE07_B2_file_list, raster) 
USA_Ne3_LE07_B2_crop <- lapply(X= USA_Ne3_LE07_B2_rast, FUN=crop, y=USA_Ne3_FFP)
USA_Ne3_LE07_B2_stack <- stack(USA_Ne3_LE07_B2_crop)

USA_Ne3_LE07_GCVI <- (USA_Ne3_LE07_B4_stack/USA_Ne3_LE07_B2_stack) -1
USA_Ne3_LE07_GCVI <- as_tibble_col(cellStats(USA_Ne3_LE07_GCVI, stat=mean, na.rm=TRUE))

#' Get times for Landsat 7 dates
USA_time <- str_sub(USA_Ne3_LE07_B2_file_list,65,72)
USA_time <- as_tibble_col(USA_time)
USA_time <- as.data.frame(USA_time)
USA_time$Day_times <- as.POSIXct(USA_time$value, format=("%Y%m%d"))
USA_time$date <- as.character(as.Date(USA_time$Day_times+86400))
names(USA_time) <- c("value", "time", "date")
USA_Ne3_LE07_GCVI <- bind_cols(USA_time$date, USA_Ne3_LE07_GCVI$value)
names(USA_Ne3_LE07_GCVI) <- c("date", "GCVI")

#' Get blue band = B1
USA_Ne3_LE07_B1_file_list <- list.files(USA_Ne3_LE07_dir_list, pattern="_B1", all.files=TRUE, full.names = TRUE)
USA_Ne3_LE07_B1_rast <- lapply(USA_Ne3_LE07_B1_file_list, raster) 
USA_Ne3_LE07_B1_crop <- lapply(X= USA_Ne3_LE07_B1_rast, FUN=crop, y=USA_Ne3_FFP)
USA_Ne3_LE07_B1_stack <- stack(USA_Ne3_LE07_B1_crop)
USA_Ne3_LE07_blue <- as_tibble_col(cellStats(USA_Ne3_LE07_B1_stack, stat=mean, na.rm=TRUE))

#' Get green band = B2
USA_Ne3_LE07_B2_file_list <- list.files(USA_Ne3_LE07_dir_list, pattern="_B2", all.files=TRUE, full.names = TRUE)
USA_Ne3_LE07_B2_rast <- lapply(USA_Ne3_LE07_B2_file_list, raster) 
USA_Ne3_LE07_B2_crop <- lapply(X= USA_Ne3_LE07_B2_rast, FUN=crop, y=USA_Ne3_FFP)
USA_Ne3_LE07_B2_stack <- stack(USA_Ne3_LE07_B2_crop)
USA_Ne3_LE07_green <- as_tibble_col(cellStats(USA_Ne3_LE07_B2_stack, stat=mean, na.rm=TRUE))

#' Get red band = B3
USA_Ne3_LE07_B3_file_list <- list.files(USA_Ne3_LE07_dir_list, pattern="_B3", all.files=TRUE, full.names = TRUE)
USA_Ne3_LE07_B3_rast <- lapply(USA_Ne3_LE07_B3_file_list, raster) 
USA_Ne3_LE07_B3_crop <- lapply(X= USA_Ne3_LE07_B3_rast, FUN=crop, y=USA_Ne3_FFP)
USA_Ne3_LE07_B3_stack <- stack(USA_Ne3_LE07_B3_crop)
USA_Ne3_LE07_red <- as_tibble_col(cellStats(USA_Ne3_LE07_B3_stack, stat=mean, na.rm=TRUE))

#' Get nir band = B4
USA_Ne3_LE07_B4_file_list <- list.files(USA_Ne3_LE07_dir_list, pattern="_B4", all.files=TRUE, full.names = TRUE)
USA_Ne3_LE07_B4_rast <- lapply(USA_Ne3_LE07_B4_file_list, raster) 
USA_Ne3_LE07_B4_crop <- lapply(X= USA_Ne3_LE07_B4_rast, FUN=crop, y=USA_Ne3_FFP)
USA_Ne3_LE07_B4_stack <- stack(USA_Ne3_LE07_B4_crop)
USA_Ne3_LE07_nir <- as_tibble_col(cellStats(USA_Ne3_LE07_B4_stack, stat=mean, na.rm=TRUE))

#' Get swir1 band = B5
USA_Ne3_LE07_B5_file_list <- list.files(USA_Ne3_LE07_dir_list, pattern="_B5", all.files=TRUE, full.names = TRUE)
USA_Ne3_LE07_B5_rast <- lapply(USA_Ne3_LE07_B5_file_list, raster) 
USA_Ne3_LE07_B5_crop <- lapply(X= USA_Ne3_LE07_B5_rast, FUN=crop, y=USA_Ne3_FFP)
USA_Ne3_LE07_B5_stack <- stack(USA_Ne3_LE07_B5_crop)
USA_Ne3_LE07_swir1 <- as_tibble_col(cellStats(USA_Ne3_LE07_B5_stack, stat=mean, na.rm=TRUE))

#' Get swir2 band = B7
USA_Ne3_LE07_B7_file_list <- list.files(USA_Ne3_LE07_dir_list, pattern="_B7", all.files=TRUE, full.names = TRUE)
USA_Ne3_LE07_B7_rast <- lapply(USA_Ne3_LE07_B7_file_list, raster) 
USA_Ne3_LE07_B7_crop <- lapply(X= USA_Ne3_LE07_B7_rast, FUN=crop, y=USA_Ne3_FFP)
USA_Ne3_LE07_B7_stack <- stack(USA_Ne3_LE07_B7_crop)
USA_Ne3_LE07_swir2 <- as_tibble_col(cellStats(USA_Ne3_LE07_B7_stack, stat=mean, na.rm=TRUE))

#' Convert to reflectance
USA_Ne3_LE07_blue$refl <- USA_Ne3_LE07_blue$value* USA_Ne3_LE07_conv_vals$B1_MULT_val+
  USA_Ne3_LE07_conv_vals$B1_ADD_val

USA_Ne3_LE07_green$refl <- USA_Ne3_LE07_green$value* USA_Ne3_LE07_conv_vals$B2_MULT_val+
  USA_Ne3_LE07_conv_vals$B2_ADD_val

USA_Ne3_LE07_red$refl <- USA_Ne3_LE07_red$value* USA_Ne3_LE07_conv_vals$B3_MULT_val+
  USA_Ne3_LE07_conv_vals$B3_ADD_val

USA_Ne3_LE07_nir$refl <- USA_Ne3_LE07_nir$value* USA_Ne3_LE07_conv_vals$B4_MULT_val+
  USA_Ne3_LE07_conv_vals$B4_ADD_val

USA_Ne3_LE07_swir1$refl <- USA_Ne3_LE07_swir1$value* USA_Ne3_LE07_conv_vals$B5_MULT_val+
  USA_Ne3_LE07_conv_vals$B5_ADD_val

USA_Ne3_LE07_swir2$refl <- USA_Ne3_LE07_swir2$value* USA_Ne3_LE07_conv_vals$B7_MULT_val+
  USA_Ne3_LE07_conv_vals$B7_ADD_val

#' Combine band data with GCVI
USA_Ne3_LE07_dat <- bind_cols(USA_Ne3_LE07_GCVI, USA_Ne3_LE07_blue$refl, USA_Ne3_LE07_green$refl,
                              USA_Ne3_LE07_red$refl, USA_Ne3_LE07_nir$refl, USA_Ne3_LE07_swir1$refl,
                              USA_Ne3_LE07_swir2$refl)

names(USA_Ne3_LE07_dat) <- c("date", "GCVI", "blue", "green", "red", "nir", "swir1", "swir2")

#' Add sensor
USA_Ne3_LE07_dat <- USA_Ne3_LE07_dat %>% 
  mutate(sensor = "Ls7")


######################## Repeat for Landsat 8 #################################
USA_Ne3_LC08_dir_list <- list.dirs("E:\\Omaha\\PySEBAL_data\\Satellite_data\\")

#' Get metadata to convert to TOA reflectance
USA_Ne3_LC08_MTL_list <- grep(pattern="LC08", x=USA_Ne3_LC08_dir_list, value=TRUE)
USA_Ne3_LC08_MTL_list <- list.files(USA_Ne3_LC08_MTL_list, pattern="MTL", all.files=TRUE, full.names = TRUE)

#' Write a for loop to store the conversion values for each date (row) in a column for each band
USA_Ne3_LC08_conv_vals <- data.frame(
  B2_MULT_val=rep(NA,length(USA_Ne3_LC08_MTL_list)), B2_ADD_val=rep(NA,length(USA_Ne3_LC08_MTL_list)),
  B3_MULT_val=rep(NA,length(USA_Ne3_LC08_MTL_list)), B3_ADD_val=rep(NA,length(USA_Ne3_LC08_MTL_list)),
  B4_MULT_val=rep(NA,length(USA_Ne3_LC08_MTL_list)), B4_ADD_val=rep(NA,length(USA_Ne3_LC08_MTL_list)),
  B5_MULT_val=rep(NA,length(USA_Ne3_LC08_MTL_list)), B5_ADD_val=rep(NA,length(USA_Ne3_LC08_MTL_list)),
  B6_MULT_val=rep(NA,length(USA_Ne3_LC08_MTL_list)), B6_ADD_val=rep(NA,length(USA_Ne3_LC08_MTL_list)),
  B7_MULT_val=rep(NA,length(USA_Ne3_LC08_MTL_list)), B7_ADD_val=rep(NA,length(USA_Ne3_LC08_MTL_list)))

for (i in 1:length(USA_Ne3_LC08_MTL_list)) {
  f <- readLines(USA_Ne3_LC08_MTL_list[i])
  B2_MULT <- grep("REFLECTANCE_MULT_BAND_2", f, value=TRUE)
  USA_Ne3_LC08_conv_vals$B2_MULT_val[i] <- as.numeric(sub(".*=", "", B2_MULT))
  B2_ADD <- grep("REFLECTANCE_ADD_BAND_2", f, value=TRUE)
  USA_Ne3_LC08_conv_vals$B2_ADD_val[i] <- as.numeric(sub(".*=", "", B2_ADD))
  
  B3_MULT <- grep("REFLECTANCE_MULT_BAND_3", f, value=TRUE)
  USA_Ne3_LC08_conv_vals$B3_MULT_val[i] <- as.numeric(sub(".*=", "", B3_MULT))
  B3_ADD <- grep("REFLECTANCE_ADD_BAND_3", f, value=TRUE)
  USA_Ne3_LC08_conv_vals$B3_ADD_val[i] <- as.numeric(sub(".*=", "", B3_ADD))
  
  B4_MULT <- grep("REFLECTANCE_MULT_BAND_4", f, value=TRUE)
  USA_Ne3_LC08_conv_vals$B4_MULT_val[i] <- as.numeric(sub(".*=", "", B4_MULT))
  B4_ADD <- grep("REFLECTANCE_ADD_BAND_4", f, value=TRUE)
  USA_Ne3_LC08_conv_vals$B4_ADD_val[i] <- as.numeric(sub(".*=", "", B4_ADD))
  
  B5_MULT <- grep("REFLECTANCE_MULT_BAND_5", f, value=TRUE)
  USA_Ne3_LC08_conv_vals$B5_MULT_val[i] <- as.numeric(sub(".*=", "", B5_MULT))
  B5_ADD <- grep("REFLECTANCE_ADD_BAND_5", f, value=TRUE)
  USA_Ne3_LC08_conv_vals$B5_ADD_val[i] <- as.numeric(sub(".*=", "", B5_ADD))
  
  B6_MULT <- grep("REFLECTANCE_MULT_BAND_6", f, value=TRUE)
  USA_Ne3_LC08_conv_vals$B6_MULT_val[i] <- as.numeric(sub(".*=", "", B6_MULT))
  B6_ADD <- grep("REFLECTANCE_ADD_BAND_6", f, value=TRUE)
  USA_Ne3_LC08_conv_vals$B6_ADD_val[i] <- as.numeric(sub(".*=", "", B6_ADD))
  
  B7_MULT <- grep("REFLECTANCE_MULT_BAND_7", f, value=TRUE)
  USA_Ne3_LC08_conv_vals$B7_MULT_val[i] <- as.numeric(sub(".*=", "", B7_MULT))
  B7_ADD <- grep("REFLECTANCE_ADD_BAND_7", f, value=TRUE)
  USA_Ne3_LC08_conv_vals$B7_ADD_val[i] <- as.numeric(sub(".*=", "", B7_ADD))
  
}

#' Get blue band = B2
USA_Ne3_LC08_B2_file_list <- grep(pattern="LC08", x=USA_Ne3_LC08_dir_list, value=TRUE)
USA_Ne3_LC08_B2_file_list <- list.files(USA_Ne3_LC08_B2_file_list, pattern="B2", all.files=TRUE, full.names = TRUE)
USA_Ne3_LC08_B2_rast <- lapply(USA_Ne3_LC08_B2_file_list, raster) 
USA_Ne3_LC08_B2_crop <- lapply(X= USA_Ne3_LC08_B2_rast, FUN=crop, y=USA_Ne3_FFP)
USA_Ne3_LC08_B2_stack <- stack(USA_Ne3_LC08_B2_crop)
USA_Ne3_LC08_blue <- as_tibble_col(cellStats(USA_Ne3_LC08_B2_stack, stat=mean, na.rm=TRUE))

#' Get green band = B3
USA_Ne3_LC08_B3_file_list <- grep(pattern="LC08", x=USA_Ne3_LC08_dir_list, value=TRUE)
USA_Ne3_LC08_B3_file_list <- list.files(USA_Ne3_LC08_B3_file_list, pattern="B3", all.files=TRUE, full.names = TRUE)
USA_Ne3_LC08_B3_rast <- lapply(USA_Ne3_LC08_B3_file_list, raster) 
USA_Ne3_LC08_B3_crop <- lapply(X= USA_Ne3_LC08_B3_rast, FUN=crop, y=USA_Ne3_FFP)
USA_Ne3_LC08_B3_stack <- stack(USA_Ne3_LC08_B3_crop)
USA_Ne3_LC08_green <- as_tibble_col(cellStats(USA_Ne3_LC08_B3_stack, stat=mean, na.rm=TRUE))

#' Get red band = B4
USA_Ne3_LC08_B4_file_list <- grep(pattern="LC08", x=USA_Ne3_LC08_dir_list, value=TRUE)
USA_Ne3_LC08_B4_file_list <- list.files(USA_Ne3_LC08_B4_file_list, pattern="B4", all.files=TRUE, full.names = TRUE)
USA_Ne3_LC08_B4_rast <- lapply(USA_Ne3_LC08_B4_file_list, raster) 
USA_Ne3_LC08_B4_crop <- lapply(X= USA_Ne3_LC08_B4_rast, FUN=crop, y=USA_Ne3_FFP)
USA_Ne3_LC08_B4_stack <- stack(USA_Ne3_LC08_B4_crop)
USA_Ne3_LC08_red <- as_tibble_col(cellStats(USA_Ne3_LC08_B4_stack, stat=mean, na.rm=TRUE))

#' Get nir band = B5
USA_Ne3_LC08_B5_file_list <- grep(pattern="LC08", x=USA_Ne3_LC08_dir_list, value=TRUE)
USA_Ne3_LC08_B5_file_list <- list.files(USA_Ne3_LC08_B5_file_list, pattern="B5", all.files=TRUE, full.names = TRUE)
USA_Ne3_LC08_B5_rast <- lapply(USA_Ne3_LC08_B5_file_list, raster) 
USA_Ne3_LC08_B5_crop <- lapply(X= USA_Ne3_LC08_B5_rast, FUN=crop, y=USA_Ne3_FFP)
USA_Ne3_LC08_B5_stack <- stack(USA_Ne3_LC08_B5_crop)
USA_Ne3_LC08_nir <- as_tibble_col(cellStats(USA_Ne3_LC08_B5_stack, stat=mean, na.rm=TRUE))

#' Get swir1 band = B6
USA_Ne3_LC08_B6_file_list <- grep(pattern="LC08", x=USA_Ne3_LC08_dir_list, value=TRUE)
USA_Ne3_LC08_B6_file_list <- list.files(USA_Ne3_LC08_B6_file_list, pattern="B6", all.files=TRUE, full.names = TRUE)
USA_Ne3_LC08_B6_rast <- lapply(USA_Ne3_LC08_B6_file_list, raster) 
USA_Ne3_LC08_B6_crop <- lapply(X= USA_Ne3_LC08_B6_rast, FUN=crop, y=USA_Ne3_FFP)
USA_Ne3_LC08_B6_stack <- stack(USA_Ne3_LC08_B6_crop)
USA_Ne3_LC08_swir1 <- as_tibble_col(cellStats(USA_Ne3_LC08_B6_stack, stat=mean, na.rm=TRUE))

#' Get swir2 band = B7
USA_Ne3_LC08_B7_file_list <- grep(pattern="LC08", x=USA_Ne3_LC08_dir_list, value=TRUE)
USA_Ne3_LC08_B7_file_list <- list.files(USA_Ne3_LC08_B7_file_list, pattern="B7", all.files=TRUE, full.names = TRUE)
USA_Ne3_LC08_B7_rast <- lapply(USA_Ne3_LC08_B7_file_list, raster) 
USA_Ne3_LC08_B7_crop <- lapply(X= USA_Ne3_LC08_B7_rast, FUN=crop, y=USA_Ne3_FFP)
USA_Ne3_LC08_B7_stack <- stack(USA_Ne3_LC08_B7_crop)
USA_Ne3_LC08_swir2 <- as_tibble_col(cellStats(USA_Ne3_LC08_B7_stack, stat=mean, na.rm=TRUE))

#' Convert to reflectance
USA_Ne3_LC08_blue$refl <- USA_Ne3_LC08_blue$value* USA_Ne3_LC08_conv_vals$B2_MULT_val+
  USA_Ne3_LC08_conv_vals$B2_ADD_val

USA_Ne3_LC08_green$refl <- USA_Ne3_LC08_green$value* USA_Ne3_LC08_conv_vals$B3_MULT_val+
  USA_Ne3_LC08_conv_vals$B3_ADD_val

USA_Ne3_LC08_red$refl <- USA_Ne3_LC08_red$value* USA_Ne3_LC08_conv_vals$B4_MULT_val+
  USA_Ne3_LC08_conv_vals$B4_ADD_val

USA_Ne3_LC08_nir$refl <- USA_Ne3_LC08_nir$value* USA_Ne3_LC08_conv_vals$B5_MULT_val+
  USA_Ne3_LC08_conv_vals$B5_ADD_val

USA_Ne3_LC08_swir1$refl <- USA_Ne3_LC08_swir1$value* USA_Ne3_LC08_conv_vals$B6_MULT_val+
  USA_Ne3_LC08_conv_vals$B6_ADD_val

USA_Ne3_LC08_swir2$refl <- USA_Ne3_LC08_swir2$value* USA_Ne3_LC08_conv_vals$B7_MULT_val+
  USA_Ne3_LC08_conv_vals$B7_ADD_val

USA_Ne3_LC08_GCVI <- (USA_Ne3_LC08_B4_stack/USA_Ne3_LC08_B2_stack) -1
USA_Ne3_LC08_GCVI <- as_tibble_col(cellStats(USA_Ne3_LC08_GCVI, stat=mean, na.rm=TRUE))


#' Get times for Landsat 8 dates
USA_time <- str_sub(USA_Ne3_LC08_B2_file_list,56,63)
USA_time <- as_tibble_col(USA_time)
USA_time <- as.data.frame(USA_time)
USA_time$Day_times <- as.POSIXct(USA_time$value, format=("%Y%m%d"))
USA_time$date <- as.character(as.Date(USA_time$Day_times+86400))
names(USA_time) <- c("value", "time", "date")
USA_Ne3_LC08_GCVI <- bind_cols(USA_time$date, USA_Ne3_LC08_GCVI$value)
names(USA_Ne3_LC08_GCVI) <- c("date", "GCVI")

#' Combine band data with GCVI
USA_Ne3_LC08_dat <- bind_cols(USA_Ne3_LC08_GCVI, USA_Ne3_LC08_blue$refl, USA_Ne3_LC08_green$refl,
                              USA_Ne3_LC08_red$refl, USA_Ne3_LC08_nir$refl, USA_Ne3_LC08_swir1$refl,
                              USA_Ne3_LC08_swir2$refl)

names(USA_Ne3_LC08_dat) <- c("date", "GCVI", "blue", "green", "red", "nir", "swir1", "swir2")

#' Add sensor
USA_Ne3_LC08_dat <- USA_Ne3_LC08_dat %>% 
  mutate(sensor = "Ls8")

########################## Combine two satellite datasets #######################
USA_Ne3_sat_dat <- bind_rows(USA_Ne3_LE07_dat, USA_Ne3_LC08_dat)

#' Join tower dates with Landsat dates

USA_Ne3_MLvals <- USA_Ne3_MLvals  %>% 
  left_join(USA_Ne3_sat_dat, by="date")


#' Export data
write.csv(USA_Ne3_MLvals, "MLdata/USA_Ne3_soy_MLvals.csv")

