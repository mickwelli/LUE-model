
############################### RF models for LUEmax ##########################
library(tidyverse)
library(randomForest)
library(gridExtra)
library(devtools)
library(surfin)
library(randomForestCI)
library(scales)

#' Load data
LUE_ml_df <- read.csv('MLdata/LUE_ml_df.csv')
str(LUE_ml_df)
#' The data consists of climatic info from the flux tower site e.g. Ta (temp),
#' in-situ GPP from the tower (daily_GPP) from which in-situ LUEmax is derived,
#' remotely sensed GPP from SEBAL (GPP_pred), Site, Crop, and Date info, 
#' and vegetation indices (VIs) which I will use to train a model for LUEmax.
 
#' Clip LUEmax and GPP, noting there are some Inf values for LUEmax due to 
#' issues with calculation.
summary(LUE_ml_df$LUEmax)
LUE_ml_df$LUEmax[LUE_ml_df$LUEmax > 6] <- NA
LUE_ml_df$LUEmax[LUE_ml_df$LUEmax < 0.1] <-  NA
LUE_ml_df <- LUE_ml_df[!is.na(LUE_ml_df$LUEmax),]
LUE_ml_df <- LUE_ml_df %>% 
  filter(GPP_pred >0.1) %>% 
  filter(daily_GPP>0.1)
str(LUE_ml_df)
LUE_ml_df$Site <- as.factor(LUE_ml_df$Site)
LUE_ml_df$Crop <- as.factor(LUE_ml_df$Crop)
summary(LUE_ml_df$Site)

#' Add fixed LUEmax
LUE_ml_df <- LUE_ml_df %>% 
mutate(LUEmax_fixed = if_else(Crop %in% c("Maize"),
                      4.5, 2.5))

#' Still quite a lot of low values for GPP
ggplot(LUE_ml_df, aes(x=daily_GPP))+geom_density()

#'Look at phenology
ggplot(LUE_ml_df, aes(x=GDD, y=daily_GPP))+geom_point()+facet_grid(Site~Crop, scales = "free")
ggplot(LUE_ml_df, aes(x=DAS, y=LUE))+geom_point(aes(col=Site))+facet_grid(Crop~., scales="free")+geom_smooth()+
  theme_bw()+ theme(text=element_text(size=16,  family="serif"))+xlab("Days after sowing")+ylab("LUE (gC/MJ)")
ggsave("LUEphen.jpg", units="in", width= 5, height=4, dpi=900)
ggplot(LUE_ml_df, aes(x=DAS, y=daily_GPP))+geom_point(aes(col=Site))+facet_grid(Crop~., scales="free")+geom_smooth()+
  theme_bw()+ theme(text=element_text(size=16,  family="serif"))+xlab("Days after sowing")+ylab(expression(GPP~(gC/m^2)))
ggsave("GPPphen.jpg", units="in", width= 5, height=4, dpi=900)
ggplot(LUE_ml_df, aes(x=daily_Ta, y=LUE))+geom_point()+facet_grid(Site~Crop)

#' Add crop growth stage
#LUE_ml_df <- LUE_ml_df %>% 
#  mutate(growth_stage = if_else((Crop == "Maize"&((month==4 & day >= 15)|month ==5)), "emergence",
#                        if_else((Crop == "Maize"&((month==6|month ==7|(month==8 & day<=15)))), "veg_growth",
#                        if_else((Crop == "Maize"&(((month==8 & day>15)|month==9|month ==10))), "grain_fill", 
#                        if_else((Crop == "Soy"&((month==5 & day >= 15)|(month==6 & day <= 15))), "emergence",
#                        if_else((Crop == "Soy"&((month==6 & day > 15)|month==7|month==8|(month==9 & day<=15))), "veg_growth",
#                        if_else((Crop == "Soy"&((month==9 & day > 15)|month==10)), "veg_growth","2")))))))

ggplot(LUE_ml_df, aes(x=growth_stage, y=LUE))+geom_boxplot()+facet_grid(Crop~.)
summary(factor(LUE_ml_df$growth_stage))

#' Comparing GPP from SEBAL to in-situ measures from flux towers shows that
#' SEBAL tends to overpredict with fixed LUEmax values
all_rmse <- sqrt(mean(LUE_ml_df$GPP_pred- LUE_ml_df$daily_GPP, na.rm = TRUE)^2) # Starting error
all_rmse_LUEmax <- sqrt(mean(LUE_ml_df$LUEmax- LUE_ml_df$LUEmax_fixed, na.rm = TRUE)^2)
all_me <- mean( LUE_ml_df$daily_GPP-LUE_ml_df$GPP_pred, na.rm = TRUE) # Starting error
all_me
{
allplot <- ggplot(data=LUE_ml_df, aes(x=GPP_pred, y=daily_GPP, col=Crop))+
  geom_point(size=3)+scale_y_continuous(limits=c(-5,35),oob=squish)+
    xlim(-5,35)+geom_abline(intercept = 0, slope=1)+
  annotate("text", label = "1:1 line",
           x = 25, y = 30, size = 5, colour = "black")+ 
  xlab(expression(GPP~(gC/m^2)~predicted))+ ylab(expression(GPP~(gC/m^2)~'in'-situ))+
  annotate("text",x=0,y=25, label=paste("RMSE=", round(all_rmse,2)),
  colour="black", size=5)+#+ggtitle("Fixed LUEmax")+
  annotate("text",x=0,y=20, label=paste("ME=", round(all_me,2)),
  colour="black", size=5)
}
allplot
ggsave("allplot.jpg", units="in", width= 5, height=4, dpi=900)

{
  allplot_LUEmax <- ggplot(data=LUE_ml_df, aes(x=LUEmax_fixed, y=LUEmax, col=Crop))+
    geom_point(size=1.5)+scale_y_continuous(limits=c(0,10),oob=squish)+
    xlim(0,10)+geom_abline(intercept = 0, slope=1)+
    annotate("text", label = "1:1 line",x = 8, y = 6, size = 5, colour = "black")+
    xlab("LUE max (gC/MJ) fixed")+ylab("LUEmax (gC/MJ) in-situ")+
    annotate("text",x=2,y=8, label=paste("RMSE=", round(all_rmse_LUEmax,2)),
             colour="black", size=5)+ggtitle("LUEmax")
  
}
allplot_LUEmax

#' Train RF1 with all vegetation indices (VIs), Crop type, Day of year, and Site
set.seed(1) # Keep seed consistent- results are quite variable because of small dataset?
tuneRF(LUE_ml_df[,c(34,24,23,33,25,28,27,7,11,37,36,35,38,39,40,41)],LUE_ml_df[,c(20)], stepFactor = 1.5) 
names(LUE_ml_df)
#' Include mtry=1 based on tuneRF
set.seed(1)
rf1 <- randomForest(LUE ~ GCVI + LAI_vals + NDVI_vals + SAVI_vals + Crop + 
                      Dlength + Lat + daily_Ta + Rh + red + green + blue + nir +
                      swir1 + swir2  + DAS + sensor,
                    data= LUE_ml_df, importance=TRUE, ntree=500, keep.inbag=TRUE, mtry=5)

varImpPlot(rf1, type=1)
varImpPlot(rf1, type=2)
rf1
partialPlot(rf1, LUE_ml_df, LAI_vals)

set.seed(1)
rf1_max <- randomForest(LUEmax ~ GCVI + LAI_vals + NDVI_vals + SAVI_vals + Crop + 
                      Dlength + Lat + daily_Ta + Rh + red + green + blue + nir +
                      swir1 + swir2  + DAS + sensor,
                    data= LUE_ml_df, importance=TRUE, ntree=500, keep.inbag=TRUE, mtry=5)

rf1_max

#' Plot for VarImp
rf1_VIP <- as.data.frame(varImpPlot(rf1, type=1))
rf1_VIP$varnames <- rownames(rf1_VIP)
rownames(rf1_VIP) <- NULL
names(rf1_VIP) <- c("IncMSE", "varnames")
rf1_VIP$varnames <- as.factor(rf1_VIP$varnames)
levels(rf1_VIP$varnames)
levels(rf1_VIP$varnames) <- c("blue", "Crop", "Daily av. temp",
                              "Days after sowing","Daylength","GCVI", 
                              "green", "LAI", "Latitude","NDVI", "NIR", "red",
                              "Daily av. relative humidity", "SAVI", "Sensor",
                              "SWIR 1", "SWIR 2")

rf1_VIP_gg <- ggplot(rf1_VIP, aes(x=reorder(varnames, IncMSE), y=IncMSE))+
  geom_point()+geom_segment(aes(x=varnames,xend=varnames,y=0,yend=IncMSE), size=2)+
  xlab("Variable")+ylab("% Increase in MSE")+coord_flip()+
  theme(text = element_text(size=15))
rf1_VIP_gg
ggsave("rf1_VIP_gg.jpg", units="in", width= 6, height=7, dpi=900)

#' Evaluate prediction
LUE_ml_df$pred_rf1 <- rf1$predicted # 'predicted' LUE in rf object are oob predictions
#' LUEmax to GPP
LUE_ml_df$pred_GPP_rf1 <- LUE_ml_df$pred_rf1 *LUE_ml_df$APAR #* LUE_ml_df$heat_stress *
  #LUE_ml_df$vapor_stress * LUE_ml_df$moisture_stress *  
b_rmse <- sqrt(mean(LUE_ml_df$pred_GPP_rf1- LUE_ml_df$daily_GPP, na.rm = TRUE)^2)
b_rmse
b_me <- mean(LUE_ml_df$daily_GPP- LUE_ml_df$pred_GPP_rf1, na.rm = TRUE)
b_me

#' Get uncertainty for rf1
rf1_ustat <- randomForestInfJack(rf1, LUE_ml_df, calibrate=TRUE) # Gives LUE prediction and variance
summary(rf1_ustat$y.hat) # LUEmax prediction
rf1_uplot <- ggplot(data=rf1_ustat, aes(x=y.hat, y=sqrt(var.hat)))+geom_point()+geom_smooth()
rf1_uplot # SE increases with greater values of LUEmax
#' Translate LUEmax SE to GPP SE- is this right?
LUE_ml_df$rf1_se <- sqrt(rf1_ustat$var.hat)* LUE_ml_df$APAR #* LUE_ml_df$heat_stress *
# LUE_ml_df$vapor_stress * LUE_ml_df$moisture_stress 
LUE_ml_df$rf1_var <- rf1_ustat$var.hat

################ Make faceted partial plot####################################
rf1_GCVI_pp <- data.frame(partialPlot(rf1, LUE_ml_df, GCVI)) %>% 
  mutate(Var= "GCVI")
rf1_LAI_pp <- data.frame(partialPlot(rf1, LUE_ml_df, LAI_vals)) %>% 
  mutate(Var = "LAI")
rf1_NDVI_pp <- data.frame(partialPlot(rf1, LUE_ml_df, NDVI_vals)) %>% 
  mutate(Var = "NDVI")
rf1_SAVI_pp <- data.frame(partialPlot(rf1, LUE_ml_df, SAVI_vals)) %>% 
  mutate(Var = "SAVI")
rf1_DOY_pp <- data.frame(partialPlot(rf1, LUE_ml_df, DOY)) %>% 
  mutate(Var = "DOY")
rf1_Dlength_pp <- data.frame(partialPlot(rf1, LUE_ml_df, Dlength)) %>% 
  mutate(Var = "Daylength")
rf1_Crop_pp <- data.frame(partialPlot(rf1, LUE_ml_df, Crop)) %>% 
  mutate(Var = "Crop")
rf1_DAS_pp <- data.frame(partialPlot(rf1, LUE_ml_df, DAS)) %>% 
  mutate(Var = "DAS")
rf1_temp_pp <- data.frame(partialPlot(rf1, LUE_ml_df, daily_Ta)) %>% 
  mutate(Var = "temp")
rf1_rh_pp <- data.frame(partialPlot(rf1, LUE_ml_df, Rh)) %>% 
  mutate(Var = "rh")

rf1_pp_df <- bind_rows(rf1_GCVI_pp, rf1_LAI_pp, rf1_NDVI_pp,
                        rf1_SAVI_pp) %>% 
  mutate(Model = "All predictors")
rf1_pp_ind <- ggplot(data = rf1_pp_df)+geom_line(aes(x=x, y=y), size=1.5)+
  ylab('LUE (gC/MJ)') + xlab('Index value')+
  facet_grid(~Var, scales = 'free')+ theme(text = element_text(size=13))
rf1_pp_ind
ggsave("rf1_pp_ind.jpg", units="in", width= 7, height=4, dpi=900)

rf1_pp_agr_df <- bind_rows(rf1_Dlength_pp, rf1_DAS_pp, rf1_temp_pp, rf1_rh_pp)
rf1_pp_agr_df$Var <- as.factor(rf1_pp_agr_df$Var)
levels(rf1_pp_agr_df$Var)
levels(rf1_pp_agr_df$Var) <- c("Days after sowing", "Daylength (hours)", "Relative humidity (%)",
                               "Temperature (°C)")
rf1_pp_agr <- ggplot(data = rf1_pp_agr_df)+geom_line(aes(x=x, y=y), size=1.5)+
  ylab('LUE (gC/MJ)') + xlab('Predictor value')+
  facet_grid(~Var, scales = 'free')+ theme(text = element_text(size=13))
rf1_pp_agr
ggsave("rf1_pp_arg.jpg", units="in", width= 7, height=4, dpi=900)

rf1$mse

#' Plot
{
  bplot_LUEmax <- ggplot(data=LUE_ml_df, aes(x=LUEmax, y=pred_rf1, col=Crop))+
    geom_point(size=1.5)+scale_y_continuous(limits=c(0,10),oob=squish)+
    xlim(0,10)+geom_abline(intercept = 0, slope=1)+
    annotate("text", label = "1:1 line",x = 8, y = 6, size = 5, colour = "black")+
    xlab("LUE max (gC/MJ) fixed")+ylab("LUEmax (gC/MJ) in-situ")+
    annotate("text",x=2,y=8, label=paste("RMSE=", round(sqrt(rf1$mse[500]),2)),
             colour="black", size=5)+ggtitle("LUEmax")
  
}
bplot_LUEmax
grid.arrange(allplot_LUEmax, bplot_LUEmax, ncol=2)

{
bplot <- ggplot(data=LUE_ml_df, aes(x=pred_GPP_rf1, y=daily_GPP, col=Crop, 
                                    ymin=daily_GPP-(rf1_se),
                                    ymax=daily_GPP+(rf1_se)))+
  geom_point(size=1.5)+scale_y_continuous(limits=c(-5,35),oob=squish)+
    xlim(-5,35)+geom_abline(intercept = 0, slope=1)+
  annotate("text", label = "1:1 line",x = 25, y = 30, size = 5, colour = "black")+
  xlab(expression(GPP~(gC/m^2)~predicted))+ylab(expression(GPP~(gC/m^2)~'in'-situ))+
  annotate("text",x=3,y=25, label=paste("RMSE=", round(b_rmse,3)),
           colour="black", size=5)+
  annotate("text",x=3,y=20, label=paste("ME=", round(b_me,3)),
           colour="black", size=5)+
  geom_errorbar(width=0.4)#+ggtitle("Adjusted LUEmax")
  }
bplot
ggsave("compplota.jpg", units="in", width= 5, height=4, dpi=900)

#' Now train with just VIs
tuneRF(LUE_ml_df[,c(37,31,30,32,34, 33,7,18,40,39,38,41,42,43,16,44)],LUE_ml_df[,c(27)], stepFactor = 1.5) # Not many variables
set.seed(1)
rf2 <- randomForest(LUE ~ GCVI + LAI_vals + NDVI_vals + SAVI_vals +  
                      Dlength + Lat + daily_Ta + Rh + red + green + blue + nir +
                      swir1 + swir2  + DAS + sensor,
                    data= LUE_ml_df, importance=TRUE, ntree=500, mytry=7, keep.inbag=TRUE)
varImpPlot(rf2, type=1)
varImpPlot(rf2, type=2)
rf2

#' Var Imp plot
rf2_VIP <- as.data.frame(varImpPlot(rf2, type=1))
rf2_VIP$varnames <- rownames(rf2_VIP)
rownames(rf2_VIP) <- NULL
names(rf2_VIP) <- c("IncMSE", "varnames")
rf2_VIP$varnames <- as.factor(rf2_VIP$varnames)
levels(rf2_VIP$varnames) <- c("Daylength", "Day of Year", "GCVI", "LAI",
                              "Latitude", "NDVI", "SAVI")
rf2_VIP_gg <- ggplot(rf2_VIP, aes(x=reorder(varnames, IncMSE), y=IncMSE))+
  geom_point()+geom_segment(aes(x=varnames,xend=varnames,y=0,yend=IncMSE), size=2)+
  xlab("Variable")+ylab("% Increase in MSE")+coord_flip()+
  theme(text = element_text(size=15))
rf2_VIP_gg
VIPplot <- grid.arrange(rf1_VIP_gg, rf2_VIP_gg, nrow=2)
VIPplot
ggsave("rf2_VIPplot.jpg", units="in", width= 6, height=4, dpi=900)

# Partial plots
rf2_GCVI_pp <- data.frame(partialPlot(rf2, LUE_ml_df, GCVI)) %>% 
  mutate(Var= "GCVI")
rf2_LAI_pp <- data.frame(partialPlot(rf2, LUE_ml_df, LAI_vals)) %>% 
  mutate(Var = "LAI")
rf2_NDVI_pp <- data.frame(partialPlot(rf2, LUE_ml_df, NDVI_vals)) %>% 
  mutate(Var = "NDVI")
rf2_SAVI_pp <- data.frame(partialPlot(rf2, LUE_ml_df, SAVI_vals)) %>% 
  mutate(Var = "SAVI")
rf2_DOY_pp <- data.frame(partialPlot(rf2, LUE_ml_df, DOY)) %>% 
  mutate(Var = "DOY")
rf2_Dlength_pp <- data.frame(partialPlot(rf2, LUE_ml_df, Dlength)) %>% 
  mutate(Var = "Daylength")
rf2_Crop_pp <- data.frame(partialPlot(rf2, LUE_ml_df, Crop)) %>% 
  mutate(Var = "Crop")

rf2_pp_df <- bind_rows(rf2_GCVI_pp, rf2_LAI_pp, rf2_NDVI_pp,
                       rf2_SAVI_pp)%>% 
  mutate(Model = "No 'crop' predictor")
rf_pp_df <- bind_rows(rf1_pp_df, rf2_pp_df)

rf_pp <- ggplot(data = rf_pp_df)+geom_line(aes(x=x, y=y), size=1.5)+
  ylab('LUEmax') + xlab('Index value')+
  facet_grid(Model~Var, scales = 'free_x')+ theme(text = element_text(size=13))
rf_pp

grid.arrange(rf1_pp, rf2_pp, nrow=2)


#' Look at influence of VIs
partialPlot(rf2, LUE_ml_df, GCVI)
partialPlot(rf2, LUE_ml_df, LAI_vals)
partialPlot(rf2, LUE_ml_df, NDVI_vals)
partialPlot(rf2, LUE_ml_df, SAVI_vals)

#' Evaluate prediction
LUE_ml_df$pred_rf2 <- rf2$predicted
LUE_ml_df$pred_GPP_rf2 <- LUE_ml_df$pred_rf2 * LUE_ml_df$heat_stress *
  LUE_ml_df$vapor_stress * LUE_ml_df$moisture_stress * LUE_ml_df$APAR
c_rmse <- sqrt(mean(LUE_ml_df$pred_GPP_rf2- LUE_ml_df$daily_GPP, na.rm = TRUE)^2)
c_rmse
c_me <- mean(LUE_ml_df$daily_GPP- LUE_ml_df$pred_GPP_rf2, na.rm = TRUE)
c_me


#' Get uncertainty for rf2
rf2_ustat <- randomForestInfJack(rf2, LUE_ml_df, calibrate=TRUE)
rf2_uplot <- ggplot(data=rf2_ustat, aes(x=y.hat, y=sqrt(var.hat)))+geom_point()+geom_smooth()
rf2_uplot # SE increases with greater values of LUEmax
#' Translate LUEmax SE to GPP SE- is this right?
LUE_ml_df$rf2_se <- sqrt(rf2_ustat$var.hat) * LUE_ml_df$heat_stress *
  LUE_ml_df$vapor_stress * LUE_ml_df$moisture_stress * LUE_ml_df$APAR
LUE_ml_df$rf2_var <- rf1_ustat$var.hat
{
cplot <- ggplot(data=LUE_ml_df, aes(x=pred_GPP_rf2, y=daily_GPP, col=Crop,
                                      ymin=daily_GPP-(1.96*rf2_se), 
                                      ymax=daily_GPP+(1.96*rf2_se)))+
  geom_point(size=1.5)+scale_y_continuous(limits=c(-5,35),oob=squish)+
    xlim(-5,35)+geom_abline(intercept = 0, slope=1)+
  annotate("text", label = "1:1 line",
           x = 25, y = 30, size = 5, colour = "black")+
  xlab(expression(GPP~(gC/m^2)~predicted))+ylab(expression(GPP~(gC/m^2)~'in'-situ))+
  annotate("text",x=3,y=25, label=paste("RMSE=", round(c_rmse,3)),
  colour="black", size=5)+
  annotate("text",x=3,y=20, label=paste("ME=", round(c_me,3)),
  colour="black", size=5)+
  geom_errorbar(width=0.4)#+ggtitle("Adjusted LUEmax- no crop")
}
cplot
ggsave("compplotb.jpg", units="in", width= 5, height=4, dpi=900)

grid.arrange(bplot, cplot, ncol=2)
#ggsave("compplota.jpg", units="in", width= 5, height=4, dpi=900)
grid.arrange(allplot,bplot, cplot, ncol=2)

#' It looks like rice does not improve much, so the fixed value could be
#' accurate.
#' Adjusting LUEmax using the RF model leads to underestimation of GPP for maize 
#' at higher GPP values-- this is consistent with the literature.

lines <- data.frame(lines=c(4.5,2.5,2.5), Crop=c(1, 2, 3))
geom_boxplot(data=LUE_ml_df, aes(Crop, LUEmax_fixed, linetype="Fixed"))

cropplot <- ggplot(data=LUE_ml_df, aes(Crop, LUE, col=Crop))+geom_boxplot()+
  geom_point(aes(x=1, y=4.5, fill="Fixed"), shape="diamond", size= 5, col="black")+
  geom_point(aes(x=2, y=2.5, fill="Fixed"), shape="diamond", size= 5, col="black")+
  geom_point(aes(x=3, y=2.5, fill="Fixed"), shape="diamond", size= 5, col="black")+
  ylab("LUEmax (gC/MJ)")+
  theme(legend.title= element_blank())
cropplot

LUE_cropplot_df <- LUE_ml_df %>% 
  pivot_longer(cols=c(LUEmax, LUE),names_to = "LUE_c", values_to="LUE_val") %>% 
  dplyr::select(LUE_c, LUE_val, Crop)
cropplot_b <- ggplot(data=LUE_cropplot_df, aes(Crop, LUE_val, col=Crop))+geom_boxplot()+
  ylab("LUE (gC/MJ)")+
  theme(legend.title= element_blank())+#facet_grid(~LUE_c)
  theme(text=element_text(size=16,  family="serif"))
cropplot_b
#grid.arrange(allplot,bplot, cplot, cropplot, ncol=2)
 
ggsave("cropboxplot.jpg", units="in", width= 5, height=4, dpi=900)
#### Playing with phenology etc. ####

phen_plot <- LUE_ml_df %>% 
  subset(Crop=="Maize")
ggplot(data=phen_plot, aes(x=NDVI_vals, y=LUEmax)) + geom_point()+geom_smooth()
ggplot(data=phen_plot, aes(x=DOY, y=daily_GPP)) + geom_point()+geom_smooth()

################################# LUE plot #####################################
LUE_ml_df$LUE_CASA <- LUE_ml_df$LUEmax_fixed * LUE_ml_df$heat_stress * LUE_ml_df$moisture_stress *
  LUE_ml_df$vapor_stress

LUE_CASA_me <- mean(LUE_ml_df$LUE - LUE_ml_df$LUE_CASA, na.rm=TRUE)
LUE_CASA_me
LUE_CASA_rmse <- sqrt(mean(LUE_ml_df$LUE - LUE_ml_df$LUE_CASA, na.rm=TRUE)^2)
LUE_CASA_rmse

LUE_pred_me <- mean(LUE_ml_df$LUE - LUE_ml_df$pred_rf1, na.rm=TRUE)
LUE_pred_me
LUE_pred_rmse <- sqrt(mean(LUE_ml_df$LUE - LUE_ml_df$pred_rf1, na.rm=TRUE)^2)
LUE_pred_rmse

LUE_plot_df <- LUE_ml_df %>% 
  pivot_longer(cols=c(LUE_CASA, LUE, pred_rf1),names_to = "LUE_t", values_to="LUE_val") %>% 
  dplyr::select(LUE_t, LUE_val, Crop)
LUEdotplot <- ggplot(data=LUE_plot_df, aes(y=LUE_val, x=LUE_t, col=Crop))+geom_jitter()+
  stat_summary(fun=median, show.legend = FALSE, geom="crossbar", color="black", size=0.5)+
  ylim(c(0,6))+theme_bw()+ theme(text=element_text(size=16,  family="serif"))+
  geom_point(aes(x=2, y=4.5, shape="LUEmax C4"), size= 5, col="black")+
  geom_point(aes(x=2, y=2.5, shape="LUEmax C3"), size= 5, col="black")+
  scale_x_discrete(labels=c(expression('LUE'['in-situ']), expression('LUE'['CASA']),
                            expression(hat('LUE'))))+
  ylab(expression(gC/MJ^2)) + xlab(element_blank())+
  annotate("text",x=2,y=5.5, label=paste("RMSE=", round(LUE_CASA_rmse,2)))+
  annotate("text",x=2,y=5, label=paste("ME=", round(LUE_CASA_me,2)))+
  annotate("text",x=3,y=5.5, label=paste("RMSE=", round(LUE_pred_rmse,2)))+
  annotate("text",x=3,y=5, label=paste("ME=", round(LUE_pred_me,2)))+
  theme(legend.title= element_blank())
LUEdotplot
ggsave("LUEdotplot.jpg", units="in", width= 6, height=4, dpi=900)
