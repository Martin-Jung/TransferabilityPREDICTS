library(greenbrown)
library(zoo)
library(xts)
library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
library(tidyr)
#
library(marfunky)
library(roquefort)
library(gamm4)
library(ggplot2)
library(ggthemr)
library(scales)
library(data.table)
myLog <- function(...) {
  cat(paste0("[Analysis] ", Sys.time(), " | ", ..., "\n"))
}

sites <- readRDS("sites_center.rds")

# Get a plant species study from 2010 allowing 10 years of lookback
s <- sites %>% filter(Study_common_taxon == "Plantae", startyear == 2010)
# VB1_2012__Carpenter

stop("Better do the following on cluster")
myLog(unique(s$Source_ID), " | Nr.Sites: ",nrow(s)," | Start: ",unique(s$Sample_start_earliest))
b1 <- as.data.table(readRDS("Extracts/PREDICTS_center_MCD43A4_1.rds"))
b2 <- as.data.table(readRDS("Extracts/PREDICTS_center_MCD43A4_2.rds"))
b3 <- as.data.table(readRDS("Extracts/PREDICTS_center_MCD43A4_3.rds"))
b4 <- as.data.table(readRDS("Extracts/PREDICTS_center_MCD43A4_4.rds"))
b5 <- as.data.table(readRDS("Extracts/PREDICTS_center_MCD43A4_5.rds"))
b6 <- as.data.table(readRDS("Extracts/PREDICTS_center_MCD43A4_6.rds"))
b7 <- as.data.table(readRDS("Extracts/PREDICTS_center_MCD43A4_7.rds"))

i <- interval("2000-01-01",unique(s$Sample_start_earliest))

# Subset bands to respective ids
sub_b1 <- subset(b1,SSBS %in% s$SSBS)
sub_b2 <- subset(b2,SSBS %in% s$SSBS)
sub_b3 <- subset(b3,SSBS %in% s$SSBS)
sub_b4 <- subset(b4,SSBS %in% s$SSBS)
sub_b5 <- subset(b5,SSBS %in% s$SSBS)
sub_b6 <- subset(b6,SSBS %in% s$SSBS)
sub_b7 <- subset(b7,SSBS %in% s$SSBS)

# ... that fail
sub_b1 <- sub_b1[which(ymd(sub_b1$date) %within% i),]
sub_b2 <- sub_b2[which(ymd(sub_b2$date) %within% i),]
sub_b3 <- sub_b3[which(ymd(sub_b3$date) %within% i),]
sub_b4 <- sub_b4[which(ymd(sub_b4$date) %within% i),]
sub_b5 <- sub_b5[which(ymd(sub_b5$date) %within% i),]
sub_b6 <- sub_b6[which(ymd(sub_b6$date) %within% i),]
sub_b7 <- sub_b7[which(ymd(sub_b7$date) %within% i),]

# Rbind 
sub <- rbind(sub_b1,sub_b2)
sub <- rbind(sub,sub_b3)
sub <- rbind(sub,sub_b4)
sub <- rbind(sub,sub_b5)
sub <- rbind(sub,sub_b6)
sub <- rbind(sub,sub_b7)

# Save
saveRDS(sub,paste0(unique(s$Source_ID),"_BRDF",".rds"))

#### Load and prepare for display ####
sub <- readRDS(paste0(unique(s$Source_ID),"_BRDF",".rds"))
sub$value <- sub$value * 0.0001# Correct
sub$Band <- as.factor(sub$Band)
sub$date <- as.Date(ymd(sub$date),origin="2000-01-01")
sub <- subset(sub,SSBS == "VB1_2012__Carpenter 7  1")

b1 = sub %>% filter(Band ==1) %>% dplyr::select(value,date)
b1 = apply.quarterly(xts(zoo(b1$value,b1$date)),FUN = function(x) mean(x,na.rm=T))
d <- index(b1);b1 <- as.data.frame(b1);b1$date <- as.character(d);b1$Band <- 1; 
b2 = sub %>% filter(Band ==2) %>% dplyr::select(value,date)
b2 = apply.quarterly(xts(zoo(b2$value,b2$date)),FUN = function(x) mean(x,na.rm=T))
d <- index(b2);b2 <- as.data.frame(b2);b2$date <- as.character(d);b2$Band <- 2; 
b3 = sub %>% filter(Band ==3) %>% dplyr::select(value,date)
b3 = apply.quarterly(xts(zoo(b3$value,b3$date)),FUN = function(x) mean(x,na.rm=T))
d <- index(b3);b3 <- as.data.frame(b3);b3$date <- as.character(d);b3$Band <- 3; 
b4 = sub %>% filter(Band ==4) %>% dplyr::select(value,date)
b4 = apply.quarterly(xts(zoo(b4$value,b4$date)),FUN = function(x) mean(x,na.rm=T))
d <- index(b4);b4 <- as.data.frame(b4);b4$date <- as.character(d);b4$Band <- 4; 
b5 = sub %>% filter(Band ==5) %>% dplyr::select(value,date)
b5 = apply.quarterly(xts(zoo(b5$value,b5$date)),FUN = function(x) mean(x,na.rm=T))
d <- index(b5);b5 <- as.data.frame(b5);b5$date <- as.character(d);b5$Band <- 5; 
b6 = sub %>% filter(Band ==6) %>% dplyr::select(value,date)
b6 = apply.quarterly(xts(zoo(b6$value,b6$date)),FUN = function(x) mean(x,na.rm=T))
d <- index(b6);b6 <- as.data.frame(b6);b6$date <- as.character(d);b6$Band <- 6 
b7 = sub %>% filter(Band ==7) %>% dplyr::select(value,date)
b7 = apply.quarterly(xts(zoo(b7$value,b7$date)),FUN = function(x) mean(x,na.rm=T))
d <- index(b7);b7 <- as.data.frame(b7);b7$date <- as.character(d);b7$Band <- 7; 
b <- rbind(b1,b2,b3,b4,b5,b6,b7) %>% mutate(date = as.Date(ymd(date)),Band = factor(Band)) %>% rename(value = x)
# band 1 - 620–670  = Red
# band 2 - 841–876 = NIR
# band 3 - 459–479 = Blue
# band 4 - 545–565 = Green
# band 5 - 1230–1250 = Infrared
# band 6 - 1628–1652 = Infrared
# band 7 - 2105–2155 = infrared
levels(sub$Band) <- c("B1 = Red (620–670)",
                      "B2 = NIR (841–876)",
                      "B3 = Blue (459–479)",
                      "B4 = Green (545–565)",
                      "B5 = MIDIR (1230–1250)",
                      "B6 = SWIR (1628–1652)",
                      "B7 = SWIR (2105–2155)")

# GGplot
g <- ggplot(data=sub,aes(x = date, y = value)) + theme_light() +
  geom_line(alpha = .5) +
  facet_wrap(~Band,ncol=1) +
  labs(x ="",y="Spectral reflectance",title="VB1_2012__Carpenter 7 1 \n BRDF MODIS bands") +
  scale_x_date(date_breaks = "1 year",date_labels = "%Y")
g;ggsave("VB1_2012__Carpenter_BRDF-RAW.png",plot=g)

levels(sub$Band) <- c(1,2,3,4,5,6,7)

a = sub %>% dplyr::select(date,Band,value) %>% mutate(Band = paste0("B",Band)) %>% 
  spread(Band,value) %>% 
  mutate(
    #NDVI
    # (2-1) / (2+1)
    NDVI = (B2 - B1) / (B2 + B1),
    # EVI = G * (NIR – RED)/(NIR + C1*RED - C2*BLUE + L))
    #G – Gain factor
    #L – Factor for canopy background adjustment
    #C1, C2: Coefficients for correcting aerosol influences from RED using BLUE
    #MODIS EVI algorithm: L = 1, G = 2.5, C1 = 6, C2 = 7.5
    EVI = 2.5 * ((B2 - B1) / (B2 + 6.0 * B1 - 7.5 * B3 + 1.0)),
    # EVI2
    # Using only two bands without blue
    EVI2 = 2.5 * ((B2 - B1) / (B2 + B1 + 1)),
    # SAVI
    # SAVI = (1 + L) * (NIR – RED)/(NIR + RED + L)
    SAVI = (1 + 0.5) * (B2 - B1) / (B2 + B1 + 0.5),
    # NDWI
    NDWI = (B2 - B5) / (B2 + B5),
    # GVMI 
    # Global Vegetation Moisture Index 
    #( 5 + 0.1 ) - ( 7 + 0.02 ) ( 5 + 0.1 ) + ( 7 + 0.02 )
    GVMI = (B5 + 0.1) - (B7 + 0.02) / (B5 + 0.1) + (B7 + 0.02),
    # NBRI -  Normalised Burn Ratio Index 
    # NDSWIR
    #  (nir - swir2)/(nir + swir2) 
    NBRI = (B2 - B7) / (B2 + B7)
    
  ) %>% dplyr::select(date,NDVI:NBRI) %>% 
  reshape2::melt(id.vars="date") %>% mutate(date = as.Date(date))

ff <- function(x,id){
  b = x %>% filter(variable == id) %>% dplyr::select(value,date)
  bb = apply.quarterly(xts(zoo(b$value,b$date)),FUN = function(x) mean(x,na.rm=T))
  d <- index(bb);b <- as.data.frame(bb);b$date <- as.character(d);b$variable<- id
  return(b)
}
l <- list()
for(i in unique(a$variable)) l[[i]] <- ff(a,i)
ll = rbind_all(l) %>% rename(value = x) %>% mutate(date = as.Date(date))

g2 <- ggplot(data=ll,aes(x = date, y = value)) + theme_light()+
  geom_line(alpha = .5) +
  facet_wrap(~variable,ncol=1) +
  scale_y_continuous(breaks=pretty_breaks()) + 
  labs(x="",y="Vegetation index value",title="VB1_2012__Carpenter 7 \n (quarterly aggregated mean)")+
  scale_x_date(date_breaks = "1 year",date_labels = "%Y")
g2;ggsave("VB1_2012__Carpenter_VI_quarterly.png",plot=g2)

b = ll %>% dplyr::filter(variable == "NDVI")
z <- ts(zoo(b$value, b$date),start=c(2000,1),frequency=4)
d1 <- (stl(z,s.window = "periodic"))

library(ggfortify)
g <- autoplot(d1,ts.colour = 'black') + theme_light() + labs("Seasonal decomposition (LOESS)")
ggsave("VB1_2012__Carpenter_NDVI-STL.png",plot=g)
autoplot(acf(z,plot=FALSE),conf.int.fill = '#0000FF', conf.int.value = 0.8, conf.int.type = 'ma')

b = a %>% dplyr::filter(variable == "NDVI")
plot(bfastts(b$value,b$date,type="16-day"))

z <- xts(zoo(b$value,b$date))
( z <- ts(z,start=c(2000,01),frequency = 365 / 16) ) 
d1 = TsPP(z,fpg = NULL,tsgf = TSGFlinear) # Interpolate linear
(dd = Trend(z, method = "SeasonalAdjusted",funAnnual=mean))
plot(dd)
(p <- Phenology(z,"Trs",interpolate = F))

p1 <- autoplot(as.zoo(p$sos)) + theme_light() + scale_x_yearqtr(format="%Y",breaks=pretty_breaks()) +
  labs(x ="", y = "Day of year",title="Annual start of season\n (NDVI)")
p2 <- autoplot(as.zoo(p$eos)) + theme_light() + scale_x_yearqtr(format="%Y",breaks=pretty_breaks()) +
  labs(x ="", y = "Day of year",title="Annual end of season\n (NDVI)")
p3 <- autoplot(as.zoo(p$peak)) + theme_light() + scale_x_yearqtr(format="%Y",breaks=pretty_breaks()) +
  labs(x ="", y = "NDVI",title="Annual peak value\n (NDVI)")
p4 <- autoplot(as.zoo(p$mgs)) + theme_light() + scale_x_yearqtr(format="%Y",breaks=pretty_breaks()) +
  labs(x ="", y = "NDVI",title="Annual mean growing season values\n (NDVI)")
cowplot::plot_grid(p1,p2,p3,p4)

# http://things-about-r.tumblr.com/post/106806522699/change-point-detection-in-time-series-with-r-and
#http://a-little-book-of-r-for-time-series.readthedocs.org/en/latest/src/timeseries.html
