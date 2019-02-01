# Package loading
library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
library(tidyr)
library(data.table)
myLog <- function(...) {
  cat(paste0("[Spectral] ", Sys.time(), " | ", ..., "\n"))
}

#### BRDF Spectral data ####
myLog("Starting loading extractions")
# Bind them together
b1 <- as.data.table(readRDS("../DONE_P4_TimeSeriesDissim/Extracts/PREDICTS_center_MCD43A4_1.rds"))
b2 <- as.data.table(readRDS("../DONE_P4_TimeSeriesDissim/Extracts/PREDICTS_center_MCD43A4_2.rds"))
b3 <- as.data.table(readRDS("../DONE_P4_TimeSeriesDissim/Extracts/PREDICTS_center_MCD43A4_3.rds"))
b4 <- as.data.table(readRDS("../DONE_P4_TimeSeriesDissim/Extracts/PREDICTS_center_MCD43A4_4.rds"))
b5 <- as.data.table(readRDS("../DONE_P4_TimeSeriesDissim/Extracts/PREDICTS_center_MCD43A4_5.rds"))
b6 <- as.data.table(readRDS("../DONE_P4_TimeSeriesDissim/Extracts/PREDICTS_center_MCD43A4_6.rds"))
b7 <- as.data.table(readRDS("../DONE_P4_TimeSeriesDissim/Extracts/PREDICTS_center_MCD43A4_7.rds"))

sites <- readRDS("sites_center.rds")

# Get average of band values within times of sampling
results <- data.frame(SSBS = sites$SSBS,propNA = NA,
                      BRDF_Band1 = NA,BRDF_Band2 = NA,BRDF_Band3 = NA,BRDF_Band4 = NA,BRDF_Band5 = NA,BRDF_Band6 = NA,BRDF_Band7 = NA,
                      NDVI = NA, EVI = NA, EVI2 = NA, SAVI = NA,NDWI = NA,
                      H_sd = NA, H_cdis = NA)
for(siteid in results$SSBS){
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

