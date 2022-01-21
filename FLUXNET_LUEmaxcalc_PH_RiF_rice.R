library(raster)
library(sp)
library(tidyverse)
library(ncdf4)
library(lubridate)
library(geosphere)

#' LosBanos_irrigated
lat <- 14.1412
lon <- 121.26537
lat_m <- 1564035.073
lon_m <- 312757.979

#' Set flux footprint extent
PH_RiF_xmin <- lon_m-150
PH_RiF_xmax <- lon_m+150
PH_RiF_ymin <- lat_m-150
PH_RiF_ymax <- lat_m+150
PH_RiF_FFP <- extent(c(PH_RiF_xmin, PH_RiF_xmax,
                          PH_RiF_ymin, PH_RiF_ymax))

#' Read NDVI files
NDVI_dir_list <- list.dirs("E:\\LosBanos\\PySEBAL_data\\SEBAL_out\\")
NDVI_dir_list <- grep(pattern="vegetation", x=NDVI_dir_list, value=TRUE)
NDVI_list <- list.files(NDVI_dir_list, pattern="_NDVI", all.files=TRUE, full.names = TRUE)
NDVI_stack <- stack(NDVI_list)

#' Read Tact files
Tact_dir_list <- list.dirs("E:\\LosBanos\\PySEBAL_data\\SEBAL_out\\")
Tact_dir_list <- grep(pattern="evapotranspiration", x=Tact_dir_list, value=TRUE)
Tact_list <- list.files(Tact_dir_list, pattern="_Tact", all.files=TRUE, full.names = TRUE)
Tact_stack <- stack(Tact_list)

#' Read Tpot files
Tpot_dir_list <- list.dirs("E:\\LosBanos\\PySEBAL_data\\SEBAL_out\\")
Tpot_dir_list <- grep(pattern="evapotranspiration", x=Tpot_dir_list, value=TRUE)
Tpot_list <- list.files(Tpot_dir_list, pattern="_Tpot", all.files=TRUE, full.names = TRUE)
Tpot_stack <- stack(Tpot_list)

#' Read LAI files
LAI_dir_list <- list.dirs("E:\\LosBanos\\PySEBAL_data\\SEBAL_out\\")
LAI_dir_list <- grep(pattern="vegetation", x=LAI_dir_list, value=TRUE)
LAI_list <- list.files(LAI_dir_list, pattern="_lai", all.files=TRUE, full.names = TRUE)
LAI_stack <- stack(LAI_list)

#' Read SAVI files
SAVI_dir_list <- list.dirs("E:\\LosBanos\\PySEBAL_data\\SEBAL_out\\")
SAVI_dir_list <- grep(pattern="vegetation", x=LAI_dir_list, value=TRUE)
SAVI_list <- list.files(SAVI_dir_list, pattern="_SAVI", all.files=TRUE, full.names = TRUE)
SAVI_stack <- stack(SAVI_list)

#' Read biomass files
biomass_dir_list <- list.dirs("E:\\LosBanos\\PySEBAL_data\\SEBAL_out\\")
biomass_dir_list <- grep(pattern="biomass", x=biomass_dir_list, value=TRUE)
biomass_list <- list.files(biomass_dir_list, pattern="_Biomass_production", all.files=TRUE, full.names = TRUE)
biomass_stack <- stack(biomass_list)

#' Extract NDVI vals from stack
NDVI_soy_FFP <- crop(NDVI_stack, PH_RiF_FFP)
NDVI_vals <- cellStats(NDVI_soy_FFP, stat=mean, na.rm=TRUE)


#' Calculate fAPAR
fPAR_vals <- -0.161 + 1.257*NDVI_vals
fPAR_vals <- as_tibble_col(matrix(unlist(fPAR_vals)))

#' Calculate moisture stress biomass
Tact_soy_FFP <- crop(Tact_stack, PH_RiF_FFP)
Tact_vals <- cellStats(Tact_soy_FFP, stat=mean, na.rm=TRUE)
Tpot_soy_FFP <- crop(Tpot_stack, PH_RiF_FFP)
Tpot_vals <- cellStats(Tpot_soy_FFP, stat=mean, na.rm=TRUE)

Moisture_stress_vals <- Tact_vals/Tpot_vals
Moisture_stress_vals <- as_tibble_col(matrix(unlist(Moisture_stress_vals)))

#' Get tower data and supp. RH data
PH_RiF_daily_data <- read.csv("data/PH_RiF.csv")

Rh_dir_list <- list.dirs("E:\\LosBanos\\PySEBAL_data\\Meteo\\")
Rh_list <- list.files(Rh_dir_list, pattern="_Rh_avg", all.files=TRUE, full.names = TRUE)
Rh_stack <- lapply(Rh_list, raster)
Rh_stack <- stack(Rh_stack)
Rh_stack <- projectRaster(Rh_stack, crs=crs(NDVI_stack))
Rh_stack <- crop(Rh_stack, PH_RiF_FFP)

Rh_vals <- cellStats(Rh_stack, stat=mean, na.rm=TRUE)
Rh_dates <- as.character(as.POSIXct(str_sub(names(Rh_stack),2,10), format=("%Y_%m%d")))

Rh_df <- bind_cols(Rh_dates, Rh_vals)
names(Rh_df) <- c("date", "Rh")

#' Get tower values
PH_RiF_daily_data$Date <- as.character(as.Date(as.character(PH_RiF_daily_data$TIMESTAMP), format="%Y%m%d"))
PH_RiF_daily_df <- PH_RiF_daily_data %>% 
  dplyr::select(Date, TA_F, SW_IN_F, LW_IN_F, GPP_DT) #Add RH if available

names(PH_RiF_daily_df) <- c("date", "daily_Ta", "daily_Fsd", "daily_Fld",
                          "daily_GPP")

########################## Calculate thermal time ##############################
PH_RiF_daily_df$dateD <- as.Date(PH_RiF_daily_df$date)
PH_RiF_daily_df$DOY <-  as.numeric(strftime(PH_RiF_daily_df$dateD, format = "%j"))

PH_RiF_daily_GDD <- PH_RiF_daily_df %>% 
  mutate_at(vars(dateD), funs(year, month, day)) %>% 
  mutate(Season = if_else(year==2012 & DOY < 91,1,
                  if_else(year==2012 & DOY >= 91 & DOY < 288, 2,
                  if_else(year==2012 & DOY >= 288, 3,
                  if_else(year==2013 & DOY < 91, 3,
                  if_else(year==2013 & DOY >= 91 & DOY < 288, 4,
                  if_else(year==2013 & DOY >= 288, 5,
                  if_else(year==2014 & DOY < 91, 5,
                  if_else(year==2014 & DOY >= 91 & DOY < 288, 6,
                  if_else(year==2014 & DOY >= 288, 7, 0)))))))))) %>% 
  group_by(Season) %>% 
  mutate(GDD = cumsum(daily_Ta-10), 0) %>% 
  mutate(DAS = (cumsum(Season)/Season)) %>% 
  ungroup()

PH_RiF_daily_GDD <- PH_RiF_daily_GDD %>% 
  dplyr::select("dateD", "year", "month", "day", "DAS", "GDD")

PH_RiF_daily_df <- left_join(PH_RiF_daily_df, PH_RiF_daily_GDD, by="dateD")

#' Get times for Landsat dates
PH_RiF_time <- str_sub(names(Tpot_stack),16,25)
PH_RiF_time <- as_tibble_col(PH_RiF_time)
PH_RiF_time <- as.data.frame(PH_RiF_time)
PH_RiF_time$Day_times <- as.POSIXct(PH_RiF_time$value, format=("%Y_%m_%d"))
PH_RiF_time$date <- as.character(as.Date(PH_RiF_time$Day_times+86400))
names(PH_RiF_time) <- c("value", "time", "date")

#' Join tower dates with Landsat dates

PH_RiF_towervals <- PH_RiF_time  %>% 
  left_join(PH_RiF_daily_df, by="date") %>% 
  left_join(Rh_df, by="date")

#' Define constants
Th <- 35
Kt <- 23
Tl <- 0
rl <- 130
Jarvis_coeff <- (Th-Kt)/(Kt-Tl)

#' Calculate PAR
PH_RiF_towervals$PAR <- PH_RiF_towervals$daily_Fsd*0.48*0.0864
PH_RiF_towervals$fPAR <- fPAR_vals$value
PH_RiF_towervals$moisture_stress <- Moisture_stress_vals$value

#' Calculate vapor stress biomass
PH_RiF_towervals$esat <- 0.6108*exp(17.27*PH_RiF_towervals$daily_Ta/
                               (PH_RiF_towervals$daily_Ta+237.3))
PH_RiF_towervals$eact <- PH_RiF_towervals$Rh*PH_RiF_towervals$esat/100

PH_RiF_towervals$vapor_stress <- if (0.88-0.183*log(PH_RiF_towervals$esat-PH_RiF_towervals$eact) >1){
  1
} else {
  0.88-0.183*log(PH_RiF_towervals$esat-PH_RiF_towervals$eact)
}

#' Calculate heat stress biomass
PH_RiF_towervals$heat_stress <- ((PH_RiF_towervals$daily_Ta-Tl)*
                            (Th-PH_RiF_towervals$daily_Ta)^Jarvis_coeff)/
  ((Kt-Tl)*((Th-Kt)^Jarvis_coeff))

#' Calculate LUE
PH_RiF_towervals$APAR <- PH_RiF_towervals$fPAR*PH_RiF_towervals$PAR
PH_RiF_towervals$LUE <- PH_RiF_towervals$daily_GPP/
  (PH_RiF_towervals$APAR)
PH_RiF_towervals$LUEmax <- PH_RiF_towervals$LUE/
  (PH_RiF_towervals$heat_stress*PH_RiF_towervals$vapor_stress*PH_RiF_towervals$moisture_stress)

#' Get SEBAL GPP
biomass_FFP <- crop(biomass_stack, PH_RiF_FFP)
biomass_vals <- cellStats(biomass_FFP, stat=mean, na.rm=TRUE)
biomass_vals <- as_tibble_col(biomass_vals)

LAI_FFP <- crop(LAI_stack, PH_RiF_FFP)
LAI_vals <- cellStats(LAI_FFP, stat=mean, na.rm=TRUE)
LAI_vals <- as_tibble_col(LAI_vals)

SAVI_FPP <- crop(SAVI_stack, PH_RiF_FFP)
SAVI_vals <- cellStats(SAVI_FPP, stat=mean, na.rm=TRUE)
SAVI_vals <- as_tibble_col(SAVI_vals)

#DOY <- as_tibble_col(str_sub(names(NDVI_stack),24,26))
Site <- as_tibble_col("PH_RiF")

PH_RiF_MLvals <- bind_cols(Site, PH_RiF_towervals, biomass_vals,
                              NDVI_vals, LAI_vals, SAVI_vals)
names(PH_RiF_MLvals)[1:2] <- c("Site", "origin_date")
names(PH_RiF_MLvals)[27:30] <- c("biomass", "NDVI_vals", "LAI_vals",
                                    "SAVI_vals")
#PH_RiF_MLvals$DOY <- as.numeric(PH_RiF_MLvals$DOY)

#'Add latitude
PH_RiF_MLvals <- PH_RiF_MLvals %>% 
  mutate(Lat = lat) %>% 
  mutate(Dlength = daylength(Lat, DOY))

PH_RiF_MLvals$GPP_pred <- PH_RiF_MLvals$biomass/10

ggplot(data= PH_RiF_MLvals, aes(x=daily_GPP, y=GPP_pred))+geom_point()+
  geom_abline(intercept = 0, slope=1)
PH_RiF_daily_df$dateD <- as.Date(PH_RiF_daily_df$date)
ggplot(data=PH_RiF_daily_df, aes(x=dateD, y=daily_GPP))+geom_point()+
  scale_x_date(date_labels="%Y%m", date_breaks = "6 months")
ggplot(data=PH_RiF_MLvals, aes(x=date, y=daily_GPP))+geom_point()
PH_RiF_MLvals <- PH_RiF_MLvals %>% 
 mutate_at(vars(date), funs(year, month, day)) #%>% 
# filter(year %in% c(2001,2002,2003,2004,2005,2006,2007,2008,2009, 2010, 2011, 2012))
#PH_RiF_MLvals <- PH_RiF_MLvals %>%   
# mutate(Crop = if_else(year %in% c("2001","2003","2005","2007","2009", "2011"),
#                        "Maize", "Soy"))

#ggplot(data= PH_RiF_MLvals, aes(x=daily_GPP, y=GPP_pred, col=Crop))+geom_point()+
#  geom_abline(intercept = 0, slope=1)

#PH_RiF_MLvals <- PH_RiF_MLvals %>%   
#  filter(Crop == "Soy")
PH_RiF_MLvals <- PH_RiF_MLvals %>%   
  mutate(Crop = "Rice")

ggplot(data= PH_RiF_MLvals, aes(x=daily_GPP, y=GPP_pred, col=Crop))+geom_point()+
  geom_abline(intercept = 0, slope=1)

############################## Landsat 7 ######################################

#' Calculate and add GCVI = (NIR/Green) - 1
#' LE07 = (B4 / B2 ) -1
#' LC08 = (B5 / B3) -1
PH_RiF_LE07_dir_list <- list.dirs("E:\\LosBanos\\PySEBAL_data\\Satellite_data\\gapfilled")

#' Get metadata to convert to TOA reflectance
PH_RiF_LE07_MTL_list <- grep(pattern="LE07", x=PH_RiF_LE07_dir_list, value=TRUE)
PH_RiF_LE07_MTL_list <- list.files(PH_RiF_LE07_MTL_list, pattern="MTL", all.files=TRUE, full.names = TRUE)

#' Write a for loop to store the conversion values for each date (row) in a column for each band
PH_RiF_LE07_conv_vals <- data.frame(
  B1_MULT_val=rep(NA,length(PH_RiF_LE07_MTL_list)), B1_ADD_val=rep(NA,length(PH_RiF_LE07_MTL_list)),
  B2_MULT_val=rep(NA,length(PH_RiF_LE07_MTL_list)), B2_ADD_val=rep(NA,length(PH_RiF_LE07_MTL_list)),
  B3_MULT_val=rep(NA,length(PH_RiF_LE07_MTL_list)), B3_ADD_val=rep(NA,length(PH_RiF_LE07_MTL_list)),
  B4_MULT_val=rep(NA,length(PH_RiF_LE07_MTL_list)), B4_ADD_val=rep(NA,length(PH_RiF_LE07_MTL_list)),
  B5_MULT_val=rep(NA,length(PH_RiF_LE07_MTL_list)), B5_ADD_val=rep(NA,length(PH_RiF_LE07_MTL_list)),
  B7_MULT_val=rep(NA,length(PH_RiF_LE07_MTL_list)), B7_ADD_val=rep(NA,length(PH_RiF_LE07_MTL_list)))

for (i in 1:length(PH_RiF_LE07_MTL_list)) {
  f <- readLines(PH_RiF_LE07_MTL_list[i])
  
  B1_MULT <- grep("REFLECTANCE_MULT_BAND_1", f, value=TRUE)
  PH_RiF_LE07_conv_vals$B1_MULT_val[i] <- as.numeric(sub(".*=", "", B1_MULT))
  B1_ADD <- grep("REFLECTANCE_ADD_BAND_1", f, value=TRUE)
  PH_RiF_LE07_conv_vals$B1_ADD_val[i] <- as.numeric(sub(".*=", "", B1_ADD))
  
  B2_MULT <- grep("REFLECTANCE_MULT_BAND_2", f, value=TRUE)
  PH_RiF_LE07_conv_vals$B2_MULT_val[i] <- as.numeric(sub(".*=", "", B2_MULT))
  B2_ADD <- grep("REFLECTANCE_ADD_BAND_2", f, value=TRUE)
  PH_RiF_LE07_conv_vals$B2_ADD_val[i] <- as.numeric(sub(".*=", "", B2_ADD))
  
  B3_MULT <- grep("REFLECTANCE_MULT_BAND_3", f, value=TRUE)
  PH_RiF_LE07_conv_vals$B3_MULT_val[i] <- as.numeric(sub(".*=", "", B3_MULT))
  B3_ADD <- grep("REFLECTANCE_ADD_BAND_3", f, value=TRUE)
  PH_RiF_LE07_conv_vals$B3_ADD_val[i] <- as.numeric(sub(".*=", "", B3_ADD))
  
  B4_MULT <- grep("REFLECTANCE_MULT_BAND_4", f, value=TRUE)
  PH_RiF_LE07_conv_vals$B4_MULT_val[i] <- as.numeric(sub(".*=", "", B4_MULT))
  B4_ADD <- grep("REFLECTANCE_ADD_BAND_4", f, value=TRUE)
  PH_RiF_LE07_conv_vals$B4_ADD_val[i] <- as.numeric(sub(".*=", "", B4_ADD))
  
  B5_MULT <- grep("REFLECTANCE_MULT_BAND_5", f, value=TRUE)
  PH_RiF_LE07_conv_vals$B5_MULT_val[i] <- as.numeric(sub(".*=", "", B5_MULT))
  B5_ADD <- grep("REFLECTANCE_ADD_BAND_5", f, value=TRUE)
  PH_RiF_LE07_conv_vals$B5_ADD_val[i] <- as.numeric(sub(".*=", "", B5_ADD))
  
  B7_MULT <- grep("REFLECTANCE_MULT_BAND_7", f, value=TRUE)
  PH_RiF_LE07_conv_vals$B7_MULT_val[i] <- as.numeric(sub(".*=", "", B7_MULT))
  B7_ADD <- grep("REFLECTANCE_ADD_BAND_7", f, value=TRUE)
  PH_RiF_LE07_conv_vals$B7_ADD_val[i] <- as.numeric(sub(".*=", "", B7_ADD))
  
}

PH_RiF_LE07_B4_file_list <- list.files(PH_RiF_LE07_dir_list, pattern="_B4", all.files=TRUE, full.names = TRUE)
PH_RiF_LE07_B4_rast <- lapply(PH_RiF_LE07_B4_file_list, raster) 
PH_RiF_LE07_B4_crop <- lapply(X= PH_RiF_LE07_B4_rast, FUN=crop, y=PH_RiF_FFP)
PH_RiF_LE07_B4_stack <- stack(PH_RiF_LE07_B4_crop)

PH_RiF_LE07_B2_file_list <- list.files(PH_RiF_LE07_dir_list, pattern="_B2", all.files=TRUE, full.names = TRUE)
PH_RiF_LE07_B2_rast <- lapply(PH_RiF_LE07_B2_file_list, raster) 
PH_RiF_LE07_B2_crop <- lapply(X= PH_RiF_LE07_B2_rast, FUN=crop, y=PH_RiF_FFP)
PH_RiF_LE07_B2_stack <- stack(PH_RiF_LE07_B2_crop)

PH_RiF_LE07_GCVI <- (PH_RiF_LE07_B4_stack/PH_RiF_LE07_B2_stack) -1
PH_RiF_LE07_GCVI <- as_tibble_col(cellStats(PH_RiF_LE07_GCVI, stat=mean, na.rm=TRUE))

#' Get times for Landsat 7 dates
USA_time <- str_sub(PH_RiF_LE07_B2_file_list,68,75)
USA_time <- as_tibble_col(USA_time)
USA_time <- as.data.frame(USA_time)
USA_time$Day_times <- as.POSIXct(USA_time$value, format=("%Y%m%d"))
USA_time$date <- as.character(as.Date(USA_time$Day_times+86400))
names(USA_time) <- c("value", "time", "date")
PH_RiF_LE07_GCVI <- bind_cols(USA_time$date, PH_RiF_LE07_GCVI$value)
names(PH_RiF_LE07_GCVI) <- c("date", "GCVI")

#' Get blue band = B1
PH_RiF_LE07_B1_file_list <- list.files(PH_RiF_LE07_dir_list, pattern="_B1", all.files=TRUE, full.names = TRUE)
PH_RiF_LE07_B1_rast <- lapply(PH_RiF_LE07_B1_file_list, raster) 
PH_RiF_LE07_B1_crop <- lapply(X= PH_RiF_LE07_B1_rast, FUN=crop, y=PH_RiF_FFP)
PH_RiF_LE07_B1_stack <- stack(PH_RiF_LE07_B1_crop)
PH_RiF_LE07_blue <- as_tibble_col(cellStats(PH_RiF_LE07_B1_stack, stat=mean, na.rm=TRUE))

#' Get green band = B2
PH_RiF_LE07_B2_file_list <- list.files(PH_RiF_LE07_dir_list, pattern="_B2", all.files=TRUE, full.names = TRUE)
PH_RiF_LE07_B2_rast <- lapply(PH_RiF_LE07_B2_file_list, raster) 
PH_RiF_LE07_B2_crop <- lapply(X= PH_RiF_LE07_B2_rast, FUN=crop, y=PH_RiF_FFP)
PH_RiF_LE07_B2_stack <- stack(PH_RiF_LE07_B2_crop)
PH_RiF_LE07_green <- as_tibble_col(cellStats(PH_RiF_LE07_B2_stack, stat=mean, na.rm=TRUE))

#' Get red band = B3
PH_RiF_LE07_B3_file_list <- list.files(PH_RiF_LE07_dir_list, pattern="_B3", all.files=TRUE, full.names = TRUE)
PH_RiF_LE07_B3_rast <- lapply(PH_RiF_LE07_B3_file_list, raster) 
PH_RiF_LE07_B3_crop <- lapply(X= PH_RiF_LE07_B3_rast, FUN=crop, y=PH_RiF_FFP)
PH_RiF_LE07_B3_stack <- stack(PH_RiF_LE07_B3_crop)
PH_RiF_LE07_red <- as_tibble_col(cellStats(PH_RiF_LE07_B3_stack, stat=mean, na.rm=TRUE))

#' Get nir band = B4
PH_RiF_LE07_B4_file_list <- list.files(PH_RiF_LE07_dir_list, pattern="_B4", all.files=TRUE, full.names = TRUE)
PH_RiF_LE07_B4_rast <- lapply(PH_RiF_LE07_B4_file_list, raster) 
PH_RiF_LE07_B4_crop <- lapply(X= PH_RiF_LE07_B4_rast, FUN=crop, y=PH_RiF_FFP)
PH_RiF_LE07_B4_stack <- stack(PH_RiF_LE07_B4_crop)
PH_RiF_LE07_nir <- as_tibble_col(cellStats(PH_RiF_LE07_B4_stack, stat=mean, na.rm=TRUE))

#' Get swir1 band = B5
PH_RiF_LE07_B5_file_list <- list.files(PH_RiF_LE07_dir_list, pattern="_B5", all.files=TRUE, full.names = TRUE)
PH_RiF_LE07_B5_rast <- lapply(PH_RiF_LE07_B5_file_list, raster) 
PH_RiF_LE07_B5_crop <- lapply(X= PH_RiF_LE07_B5_rast, FUN=crop, y=PH_RiF_FFP)
PH_RiF_LE07_B5_stack <- stack(PH_RiF_LE07_B5_crop)
PH_RiF_LE07_swir1 <- as_tibble_col(cellStats(PH_RiF_LE07_B5_stack, stat=mean, na.rm=TRUE))

#' Get swir2 band = B7
PH_RiF_LE07_B7_file_list <- list.files(PH_RiF_LE07_dir_list, pattern="_B7", all.files=TRUE, full.names = TRUE)
PH_RiF_LE07_B7_rast <- lapply(PH_RiF_LE07_B7_file_list, raster) 
PH_RiF_LE07_B7_crop <- lapply(X= PH_RiF_LE07_B7_rast, FUN=crop, y=PH_RiF_FFP)
PH_RiF_LE07_B7_stack <- stack(PH_RiF_LE07_B7_crop)
PH_RiF_LE07_swir2 <- as_tibble_col(cellStats(PH_RiF_LE07_B7_stack, stat=mean, na.rm=TRUE))

#' Convert to reflectance
PH_RiF_LE07_blue$refl <- PH_RiF_LE07_blue$value* PH_RiF_LE07_conv_vals$B1_MULT_val+
  PH_RiF_LE07_conv_vals$B1_ADD_val

PH_RiF_LE07_green$refl <- PH_RiF_LE07_green$value* PH_RiF_LE07_conv_vals$B2_MULT_val+
  PH_RiF_LE07_conv_vals$B2_ADD_val

PH_RiF_LE07_red$refl <- PH_RiF_LE07_red$value* PH_RiF_LE07_conv_vals$B3_MULT_val+
  PH_RiF_LE07_conv_vals$B3_ADD_val

PH_RiF_LE07_nir$refl <- PH_RiF_LE07_nir$value* PH_RiF_LE07_conv_vals$B4_MULT_val+
  PH_RiF_LE07_conv_vals$B4_ADD_val

PH_RiF_LE07_swir1$refl <- PH_RiF_LE07_swir1$value* PH_RiF_LE07_conv_vals$B5_MULT_val+
  PH_RiF_LE07_conv_vals$B5_ADD_val

PH_RiF_LE07_swir2$refl <- PH_RiF_LE07_swir2$value* PH_RiF_LE07_conv_vals$B7_MULT_val+
  PH_RiF_LE07_conv_vals$B7_ADD_val

#' Combine band data with GCVI
PH_RiF_LE07_dat <- bind_cols(PH_RiF_LE07_GCVI, PH_RiF_LE07_blue$refl, PH_RiF_LE07_green$refl,
                             PH_RiF_LE07_red$refl, PH_RiF_LE07_nir$refl, PH_RiF_LE07_swir1$refl,
                             PH_RiF_LE07_swir2$refl)

names(PH_RiF_LE07_dat) <- c("date", "GCVI", "blue", "green", "red", "nir", "swir1", "swir2")

#' Add sensor
PH_RiF_LE07_dat <- PH_RiF_LE07_dat %>% 
  mutate(sensor = "Ls7")


######################## Repeat for Landsat 8 #################################
PH_RiF_LC08_dir_list <- list.dirs("E:\\LosBanos\\PySEBAL_data\\Satellite_data\\")

#' Get metadata to convert to TOA reflectance
PH_RiF_LC08_MTL_list <- grep(pattern="LC08", x=PH_RiF_LC08_dir_list, value=TRUE)
PH_RiF_LC08_MTL_list <- list.files(PH_RiF_LC08_MTL_list, pattern="MTL", all.files=TRUE, full.names = TRUE)

#' Write a for loop to store the conversion values for each date (row) in a column for each band
PH_RiF_LC08_conv_vals <- data.frame(
  B2_MULT_val=rep(NA,length(PH_RiF_LC08_MTL_list)), B2_ADD_val=rep(NA,length(PH_RiF_LC08_MTL_list)),
  B3_MULT_val=rep(NA,length(PH_RiF_LC08_MTL_list)), B3_ADD_val=rep(NA,length(PH_RiF_LC08_MTL_list)),
  B4_MULT_val=rep(NA,length(PH_RiF_LC08_MTL_list)), B4_ADD_val=rep(NA,length(PH_RiF_LC08_MTL_list)),
  B5_MULT_val=rep(NA,length(PH_RiF_LC08_MTL_list)), B5_ADD_val=rep(NA,length(PH_RiF_LC08_MTL_list)),
  B6_MULT_val=rep(NA,length(PH_RiF_LC08_MTL_list)), B6_ADD_val=rep(NA,length(PH_RiF_LC08_MTL_list)),
  B7_MULT_val=rep(NA,length(PH_RiF_LC08_MTL_list)), B7_ADD_val=rep(NA,length(PH_RiF_LC08_MTL_list)))

for (i in 1:length(PH_RiF_LC08_MTL_list)) {
  f <- readLines(PH_RiF_LC08_MTL_list[i])
  B2_MULT <- grep("REFLECTANCE_MULT_BAND_2", f, value=TRUE)
  PH_RiF_LC08_conv_vals$B2_MULT_val[i] <- as.numeric(sub(".*=", "", B2_MULT))
  B2_ADD <- grep("REFLECTANCE_ADD_BAND_2", f, value=TRUE)
  PH_RiF_LC08_conv_vals$B2_ADD_val[i] <- as.numeric(sub(".*=", "", B2_ADD))
  
  B3_MULT <- grep("REFLECTANCE_MULT_BAND_3", f, value=TRUE)
  PH_RiF_LC08_conv_vals$B3_MULT_val[i] <- as.numeric(sub(".*=", "", B3_MULT))
  B3_ADD <- grep("REFLECTANCE_ADD_BAND_3", f, value=TRUE)
  PH_RiF_LC08_conv_vals$B3_ADD_val[i] <- as.numeric(sub(".*=", "", B3_ADD))
  
  B4_MULT <- grep("REFLECTANCE_MULT_BAND_4", f, value=TRUE)
  PH_RiF_LC08_conv_vals$B4_MULT_val[i] <- as.numeric(sub(".*=", "", B4_MULT))
  B4_ADD <- grep("REFLECTANCE_ADD_BAND_4", f, value=TRUE)
  PH_RiF_LC08_conv_vals$B4_ADD_val[i] <- as.numeric(sub(".*=", "", B4_ADD))
  
  B5_MULT <- grep("REFLECTANCE_MULT_BAND_5", f, value=TRUE)
  PH_RiF_LC08_conv_vals$B5_MULT_val[i] <- as.numeric(sub(".*=", "", B5_MULT))
  B5_ADD <- grep("REFLECTANCE_ADD_BAND_5", f, value=TRUE)
  PH_RiF_LC08_conv_vals$B5_ADD_val[i] <- as.numeric(sub(".*=", "", B5_ADD))
  
  B6_MULT <- grep("REFLECTANCE_MULT_BAND_6", f, value=TRUE)
  PH_RiF_LC08_conv_vals$B6_MULT_val[i] <- as.numeric(sub(".*=", "", B6_MULT))
  B6_ADD <- grep("REFLECTANCE_ADD_BAND_6", f, value=TRUE)
  PH_RiF_LC08_conv_vals$B6_ADD_val[i] <- as.numeric(sub(".*=", "", B6_ADD))
  
  B7_MULT <- grep("REFLECTANCE_MULT_BAND_7", f, value=TRUE)
  PH_RiF_LC08_conv_vals$B7_MULT_val[i] <- as.numeric(sub(".*=", "", B7_MULT))
  B7_ADD <- grep("REFLECTANCE_ADD_BAND_7", f, value=TRUE)
  PH_RiF_LC08_conv_vals$B7_ADD_val[i] <- as.numeric(sub(".*=", "", B7_ADD))
  
}

#' Get blue band = B2
PH_RiF_LC08_B2_file_list <- grep(pattern="LC08", x=PH_RiF_LC08_dir_list, value=TRUE)
PH_RiF_LC08_B2_file_list <- list.files(PH_RiF_LC08_B2_file_list, pattern="B2", all.files=TRUE, full.names = TRUE)
PH_RiF_LC08_B2_rast <- lapply(PH_RiF_LC08_B2_file_list, raster) 
PH_RiF_LC08_B2_crop <- lapply(X= PH_RiF_LC08_B2_rast, FUN=crop, y=PH_RiF_FFP)
PH_RiF_LC08_B2_stack <- stack(PH_RiF_LC08_B2_crop)
PH_RiF_LC08_blue <- as_tibble_col(cellStats(PH_RiF_LC08_B2_stack, stat=mean, na.rm=TRUE))

#' Get green band = B3
PH_RiF_LC08_B3_file_list <- grep(pattern="LC08", x=PH_RiF_LC08_dir_list, value=TRUE)
PH_RiF_LC08_B3_file_list <- list.files(PH_RiF_LC08_B3_file_list, pattern="B3", all.files=TRUE, full.names = TRUE)
PH_RiF_LC08_B3_rast <- lapply(PH_RiF_LC08_B3_file_list, raster) 
PH_RiF_LC08_B3_crop <- lapply(X= PH_RiF_LC08_B3_rast, FUN=crop, y=PH_RiF_FFP)
PH_RiF_LC08_B3_stack <- stack(PH_RiF_LC08_B3_crop)
PH_RiF_LC08_green <- as_tibble_col(cellStats(PH_RiF_LC08_B3_stack, stat=mean, na.rm=TRUE))

#' Get red band = B4
PH_RiF_LC08_B4_file_list <- grep(pattern="LC08", x=PH_RiF_LC08_dir_list, value=TRUE)
PH_RiF_LC08_B4_file_list <- list.files(PH_RiF_LC08_B4_file_list, pattern="B4", all.files=TRUE, full.names = TRUE)
PH_RiF_LC08_B4_rast <- lapply(PH_RiF_LC08_B4_file_list, raster) 
PH_RiF_LC08_B4_crop <- lapply(X= PH_RiF_LC08_B4_rast, FUN=crop, y=PH_RiF_FFP)
PH_RiF_LC08_B4_stack <- stack(PH_RiF_LC08_B4_crop)
PH_RiF_LC08_red <- as_tibble_col(cellStats(PH_RiF_LC08_B4_stack, stat=mean, na.rm=TRUE))

#' Get nir band = B5
PH_RiF_LC08_B5_file_list <- grep(pattern="LC08", x=PH_RiF_LC08_dir_list, value=TRUE)
PH_RiF_LC08_B5_file_list <- list.files(PH_RiF_LC08_B5_file_list, pattern="B5", all.files=TRUE, full.names = TRUE)
PH_RiF_LC08_B5_rast <- lapply(PH_RiF_LC08_B5_file_list, raster) 
PH_RiF_LC08_B5_crop <- lapply(X= PH_RiF_LC08_B5_rast, FUN=crop, y=PH_RiF_FFP)
PH_RiF_LC08_B5_stack <- stack(PH_RiF_LC08_B5_crop)
PH_RiF_LC08_nir <- as_tibble_col(cellStats(PH_RiF_LC08_B5_stack, stat=mean, na.rm=TRUE))

#' Get swir1 band = B6
PH_RiF_LC08_B6_file_list <- grep(pattern="LC08", x=PH_RiF_LC08_dir_list, value=TRUE)
PH_RiF_LC08_B6_file_list <- list.files(PH_RiF_LC08_B6_file_list, pattern="B6", all.files=TRUE, full.names = TRUE)
PH_RiF_LC08_B6_rast <- lapply(PH_RiF_LC08_B6_file_list, raster) 
PH_RiF_LC08_B6_crop <- lapply(X= PH_RiF_LC08_B6_rast, FUN=crop, y=PH_RiF_FFP)
PH_RiF_LC08_B6_stack <- stack(PH_RiF_LC08_B6_crop)
PH_RiF_LC08_swir1 <- as_tibble_col(cellStats(PH_RiF_LC08_B6_stack, stat=mean, na.rm=TRUE))

#' Get swir2 band = B7
PH_RiF_LC08_B7_file_list <- grep(pattern="LC08", x=PH_RiF_LC08_dir_list, value=TRUE)
PH_RiF_LC08_B7_file_list <- list.files(PH_RiF_LC08_B7_file_list, pattern="B7", all.files=TRUE, full.names = TRUE)
PH_RiF_LC08_B7_rast <- lapply(PH_RiF_LC08_B7_file_list, raster) 
PH_RiF_LC08_B7_crop <- lapply(X= PH_RiF_LC08_B7_rast, FUN=crop, y=PH_RiF_FFP)
PH_RiF_LC08_B7_stack <- stack(PH_RiF_LC08_B7_crop)
PH_RiF_LC08_swir2 <- as_tibble_col(cellStats(PH_RiF_LC08_B7_stack, stat=mean, na.rm=TRUE))

#' Convert to reflectance
PH_RiF_LC08_blue$refl <- PH_RiF_LC08_blue$value* PH_RiF_LC08_conv_vals$B2_MULT_val+
  PH_RiF_LC08_conv_vals$B2_ADD_val

PH_RiF_LC08_green$refl <- PH_RiF_LC08_green$value* PH_RiF_LC08_conv_vals$B3_MULT_val+
  PH_RiF_LC08_conv_vals$B3_ADD_val

PH_RiF_LC08_red$refl <- PH_RiF_LC08_red$value* PH_RiF_LC08_conv_vals$B4_MULT_val+
  PH_RiF_LC08_conv_vals$B4_ADD_val

PH_RiF_LC08_nir$refl <- PH_RiF_LC08_nir$value* PH_RiF_LC08_conv_vals$B5_MULT_val+
  PH_RiF_LC08_conv_vals$B5_ADD_val

PH_RiF_LC08_swir1$refl <- PH_RiF_LC08_swir1$value* PH_RiF_LC08_conv_vals$B6_MULT_val+
  PH_RiF_LC08_conv_vals$B6_ADD_val

PH_RiF_LC08_swir2$refl <- PH_RiF_LC08_swir2$value* PH_RiF_LC08_conv_vals$B7_MULT_val+
  PH_RiF_LC08_conv_vals$B7_ADD_val

PH_RiF_LC08_GCVI <- (PH_RiF_LC08_B4_stack/PH_RiF_LC08_B2_stack) -1
PH_RiF_LC08_GCVI <- as_tibble_col(cellStats(PH_RiF_LC08_GCVI, stat=mean, na.rm=TRUE))


#' Get times for Landsat 8 dates
USA_time <- str_sub(PH_RiF_LC08_B2_file_list,59,66)
USA_time <- as_tibble_col(USA_time)
USA_time <- as.data.frame(USA_time)
USA_time$Day_times <- as.POSIXct(USA_time$value, format=("%Y%m%d"))
USA_time$date <- as.character(as.Date(USA_time$Day_times+86400))
names(USA_time) <- c("value", "time", "date")
PH_RiF_LC08_GCVI <- bind_cols(USA_time$date, PH_RiF_LC08_GCVI$value)
names(PH_RiF_LC08_GCVI) <- c("date", "GCVI")

#' Combine band data with GCVI
PH_RiF_LC08_dat <- bind_cols(PH_RiF_LC08_GCVI, PH_RiF_LC08_blue$refl, PH_RiF_LC08_green$refl,
                             PH_RiF_LC08_red$refl, PH_RiF_LC08_nir$refl, PH_RiF_LC08_swir1$refl,
                             PH_RiF_LC08_swir2$refl)

names(PH_RiF_LC08_dat) <- c("date", "GCVI", "blue", "green", "red", "nir", "swir1", "swir2")

#' Add sensor
PH_RiF_LC08_dat <- PH_RiF_LC08_dat %>% 
  mutate(sensor = "Ls8")

########################## Combine two satellite datasets #######################
PH_RiF_sat_dat <- bind_rows(PH_RiF_LE07_dat, PH_RiF_LC08_dat)

#' Join tower dates with Landsat dates

PH_RiF_MLvals <- PH_RiF_MLvals  %>% 
  left_join(PH_RiF_sat_dat, by="date")


#' Export data
write.csv(PH_RiF_MLvals, "MLdata/PH_RiF_MLvals.csv")
