# Package loading
library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
library(tidyr)
library(data.table)
source("000_HelperFunction.R")
myLog <- function(...) {
  cat(paste0("[Spectral] ", Sys.time(), " | ", ..., "\n"))
}

#### BRDF Spectral data ####
myLog("Starting loading extractions")
# Bind them together
b1 <- readRDS("Extracts/MCD43A4_Band1.rds")
b2 <- readRDS("Extracts/MCD43A4_Band2.rds")
b3 <- readRDS("Extracts/MCD43A4_Band3.rds")
b4 <- readRDS("Extracts/MCD43A4_Band4.rds")
b5 <- readRDS("Extracts/MCD43A4_Band5.rds")
b6 <- readRDS("Extracts/MCD43A4_Band6.rds")
b7 <- readRDS("Extracts/MCD43A4_Band7.rds")
sites <- readRDS("sites_diversity.rds")

# Get average of band values within times of sampling
results <- data.frame(SSBS = character(0),
                      Sample_midpoint = character(0), # From sites
                      # For period use 2 different types
                      # year of midyear of sampling | +/- 1 year around midyyear 
                      period = character(0), 
                      propNA = numeric(0), # Missing data
                      MODIS_version = character(0), # MODIS MCD32A4 version
                      # Average BRDF and vegetation index measurements
                      BRDF_Band1_mean = numeric(0), BRDF_Band2_mean = numeric(0), BRDF_Band3_mean = numeric(0),
                      BRDF_Band4_mean = numeric(0), BRDF_Band5_mean = numeric(0), BRDF_Band6_mean = numeric(0),BRDF_Band7_mean = numeric(0), 
                      NDVI_mean = numeric(0), NDVI_min = numeric(0),  NDVI_max = numeric(0), NDVI_cv = numeric(0),
                      EVI_mean = numeric(0), EVI_min = numeric(0),  EVI_max = numeric(0), EVI_cv = numeric(0),
                      SAVI_mean = numeric(0), SAVI_min = numeric(0),  SAVI_max = numeric(0), SAVI_cv = numeric(0),
                      NDWI_mean = numeric(0), NDWI_min = numeric(0),  NDWI_max = numeric(0), NDWI_cv = numeric(0),
                      # Spectral heterogeneity
                      PCA_BRDF_centroid = numeric(0), # Construct a PCA, calculate centroid within all (1-2) axes
                      # Spectral variability was then calculated as the mean of the Euclidean distances from the centroid of all principal components for each plot. Oldeland
                      # Also scale before hand to assess impact
                      H_sd = numeric(0), H_cdis = numeric(0)
                      )


for(siteid in unique(sites$SSBS)){
  myLog(siteid)
  sub <- subset(sites,SSBS == siteid)
  # Create interval
  i <- interval(sub$Sample_start_earliest, sub$Sample_end_latest)
  # Subset bands to respective id
  sub_b1 <- subset(b1,SSBS == siteid)
  sub_b2 <- subset(b2,SSBS == siteid)
  sub_b3 <- subset(b3,SSBS == siteid)
  sub_b4 <- subset(b4,SSBS == siteid)
  sub_b5 <- subset(b5,SSBS == siteid)
  sub_b6 <- subset(b6,SSBS == siteid)
  sub_b7 <- subset(b7,SSBS == siteid)
  
  # ... that fail
  sub_b1 <- sub_b1[which(ymd(sub_b1$date) %within% i),]
  sub_b2 <- sub_b2[which(ymd(sub_b2$date) %within% i),]
  sub_b3 <- sub_b3[which(ymd(sub_b3$date) %within% i),]
  sub_b4 <- sub_b4[which(ymd(sub_b4$date) %within% i),]
  sub_b5 <- sub_b5[which(ymd(sub_b5$date) %within% i),]
  sub_b6 <- sub_b6[which(ymd(sub_b6$date) %within% i),]
  sub_b7 <- sub_b7[which(ymd(sub_b7$date) %within% i),]
  
  # Save output metrics
  # average bands
  results$BRDF_Band1[which(results$SSBS==siteid)] <- mean(sub_b1$value,na.rm = T) * 0.0001
  results$BRDF_Band2[which(results$SSBS==siteid)] <- mean(sub_b2$value,na.rm = T) * 0.0001
  results$BRDF_Band3[which(results$SSBS==siteid)] <- mean(sub_b3$value,na.rm = T) * 0.0001
  results$BRDF_Band4[which(results$SSBS==siteid)] <- mean(sub_b4$value,na.rm = T) * 0.0001
  results$BRDF_Band5[which(results$SSBS==siteid)] <- mean(sub_b5$value,na.rm = T) * 0.0001
  results$BRDF_Band6[which(results$SSBS==siteid)] <- mean(sub_b6$value,na.rm = T) * 0.0001
  results$BRDF_Band7[which(results$SSBS==siteid)] <- mean(sub_b7$value,na.rm = T) * 0.0001
  
  # Missing values
  results$propNA[which(results$SSBS==siteid)] <- length(which(is.na(sub_b1$value))) / nrow(sub_b1)
  
  ## Metrics
  # Combine all bands
  
  sub_b1 %>% rename(Band1 = value) %>% select(SSBS,variable,Band1)
  
  bands <- sub_b1 %>% rename(Band1 = value) %>% select(SSBS,variable,Band1) %>% 
    left_join(.,(sub_b2 %>% rename(Band2 = value) %>% select(variable,Band2) ),by= c("variable")) %>% 
    left_join(.,(sub_b3 %>% rename(Band3 = value) %>% select(variable,Band3) ),by= c("variable")) %>% 
    left_join(.,(sub_b4 %>% rename(Band4 = value) %>% select(variable,Band4) ),by= c("variable")) %>% 
    left_join(.,(sub_b5 %>% rename(Band5 = value) %>% select(variable,Band5) ),by= c("variable")) %>% 
    left_join(.,(sub_b6 %>% rename(Band6 = value) %>% select(variable,Band6) ),by= c("variable")) %>% 
    left_join(.,(sub_b7 %>% rename(Band7 = value) %>% select(variable,Band7) ),by= c("variable"))
  # Correct values
  bands[,3:ncol(bands)] <- bands[,3:ncol(bands)] * 0.0001
  
  # NDVI
  # (5-4) / (5+4)
  bands$NDVI <- (bands$Band2 - bands$Band1) / (bands$Band2 + bands$Band1)
  results$NDVI[which(results$SSBS==siteid)] <- mean(bands$NDVI,na.rm = T)
  
  # EVI
  # EVI = G * (NIR – RED)/(NIR + C1*RED - C2*BLUE + L))
  #G – Gain factor
  #L – Factor for canopy background adjustment
  #C1, C2: Coefficients for correcting aerosol influences from RED using BLUE
  #MODIS EVI algorithm: L = 1, G = 2.5, C1 = 6, C2 = 7.5
  bands$EVI <- 2.5 * ((bands$Band2 - bands$Band1) / (bands$Band2 + 6.0 * bands$Band1 - 7.5 * bands$Band3 + 1.0))
  results$EVI[which(results$SSBS==siteid)] <- mean(bands$EVI,na.rm = T)
  
  # EVI2
  # Using only two bands without blue
  bands$EVI2 <- 2.5 * ((bands$Band2 - bands$Band1) / (bands$Band2 + bands$Band1 + 1))
  results$EVI2[which(results$SSBS==siteid)] <- mean(bands$EVI2,na.rm = T)
  
  # SAVI
  # SAVI = (1 + L) * (NIR – RED)/(NIR + RED + L)
  bands$SAVI <- (1 + 0.5) * (bands$Band2 - bands$Band1) / (bands$Band2 + bands$Band1 + 0.5)
  results$SAVI[which(results$SSBS==siteid)] <- mean(bands$SAVI,na.rm = T)
  
  # NDWI
  bands$NDWI <- (bands$Band2 - bands$Band5) / (bands$Band2 + bands$Band5)
  results$NDWI[which(results$SSBS==siteid)] <- mean(bands$NDWI,na.rm = T)
  
  # H_sd
  # Standard deviation of the first axis of a PCA using all bands
  mod <- select(bands,Band1:Band7) %>% subset(.,complete.cases(.)) 
  mod <- try(prcomp(mod),silent = T)
  if(class(mod)!="try-error") results$H_sd[which(results$SSBS==siteid)] <- mod$sdev[1]
  
  # H_cdis
  # Calculate average distance from the mean of all values
  bands$mean <- rowMeans(select(bands,Band1:Band7),na.rm = T)
  results$H_cdis[which(results$SSBS==siteid)] <- mean(rowMeans(abs(select(bands,Band1:Band7) - bands$mean),na.rm = T),na.rm = T)
  
  rm(sub_b1,sub_b2,sub_b3,sub_b4,sub_b5,sub_b6,sub_b7,sub,i,bands,mod) # clean up
}

saveRDS(results,"MCD43A4_BRDF_center_computed.rds")

