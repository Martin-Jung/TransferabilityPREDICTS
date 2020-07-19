# Package loading
library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
library(zoo)
library(imputeTS)
library(tidyr)
library(tidyverse)
library(data.table)
library(parallel)
library(doParallel)
source("000_HelperFunction.R")
myLog <- function(...) {
  cat(paste0("[Spectral] ", Sys.time(), " | ", ..., "\n"))
}
trapezoid = function(x, y) 
{ # computes the integral of y with respect to x using trapezoidal integration. 
  idx = 2:length(x)
  return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
}
# Path to the MODIS BRDF extracts
path = "/lustre/scratch/lifesci/mj291/PREDICTS_MCDv6"

cores = 6  # Number of cores to use
tp = c("midyear","yearbefore")[2] #interval
aggregate2month = TRUE # Mean aggregate to months
gapfill = TRUE # Should missing data be filled?

# ---------------------------- #
#### BRDF Spectral data ####
myLog("Starting loading extractions")
cl <- parallel::makeCluster( cores, outfile = paste0( round( as.numeric(now()) ) ,"_logfile.txt") )
doParallel::registerDoParallel(cl)

# Bind them together
b1 <- as.data.table( readRDS(paste0(path,"/MCD43A4_Band1.rds")) )
b2 <- as.data.table( readRDS(paste0(path,"/MCD43A4_Band2.rds")) )
b3 <- as.data.table( readRDS(paste0(path,"/MCD43A4_Band3.rds")) )
b4 <- as.data.table( readRDS(paste0(path,"/MCD43A4_Band4.rds")) )
b5 <- as.data.table( readRDS(paste0(path,"/MCD43A4_Band5.rds")) )
b6 <- as.data.table( readRDS(paste0(path,"/MCD43A4_Band6.rds")) )
b7 <- as.data.table( readRDS(paste0(path,"/MCD43A4_Band7.rds")) )

sites <- readRDS("resSaves/sites_diversity.rds") %>% 
  dplyr::filter( year(Sample_start_earliest) >= 2001 ) # Remove sites from before 2001 as there will be no data from MODIS for this period

# Initial subsets
myLog("subsetting to only target sites")
b1 <- subset(b1,SSBS %in% sites$SSBS) # These subsets use the data.table structure first
b2 <- subset(b2,SSBS %in% sites$SSBS)
b3 <- subset(b3,SSBS %in% sites$SSBS)
b4 <- subset(b4,SSBS %in% sites$SSBS)
b5 <- subset(b5,SSBS %in% sites$SSBS)
b6 <- subset(b6,SSBS %in% sites$SSBS)
b7 <- subset(b7,SSBS %in% sites$SSBS)

# --- #
# Eventual gap filling
# Mean aggregate all time series to monthly steps
if(aggregate2month){
  # Theoretical MODIS monthly measurements
  modis.full <- zooreg(data=NA,order.by =  as.yearmon(2000 + seq(0, (2016-2000)*12)/12),frequency = 12 ) %>% 
    as.data.frame() %>% add_rownames('date') %>% rename('value' = '.')
  # Aggregation function
  agg <- function(df,modis.full){ df %>% 
      mutate(date = as.character(as.yearmon(date)) ) %>% 
      # Those date not present, append
      bind_rows(., crossing(modis.full,SSBS= df$SSBS, Band = df$Band)) %>% 
      group_by(SSBS,Band,date) %>% dplyr::summarise(value = mean(value,na.rm = T)) %>% ungroup() %>% 
      mutate(date = as.yearmon(date)) %>% arrange(date)
  }
  # Now aggregate
  b1 <- agg(b1,modis.full)
  b2 <- agg(b2,modis.full)
  b3 <- agg(b3,modis.full)
  b4 <- agg(b4,modis.full)
  b5 <- agg(b5,modis.full)
  b6 <- agg(b6,modis.full)
  b7 <- agg(b7,modis.full)
}
# Use gap fill those
if(gapfill){
  # For each site and band, run an arima model on the time series
  doFill <- function(b){
    foreach(siteid =  unique(b$SSBS),
            .combine = bind_rows,
            .multicombine = FALSE,
            .errorhandling = 'pass',
            .packages = c('dplyr','lubridate','imputeTS','zoo'),
            .verbose = FALSE) %dopar% {
              sub <- subset(b,SSBS == siteid)
              left_join(sub %>%  mutate(date = as.character(date)),
                        na_kalman(zoo(x = sub$value,order.by = sub$date),model = 'auto.arima') %>% 
                          as.data.frame() %>% add_rownames('date') %>% rename('value.gf' = '.') %>% 
                          mutate(date = as.character(date))
                        )
            }
  }
  # Gap fill all bands
  b1 <- doFill(b1) %>% rename(value.nongf = value,value = value.gf )
  b2 <- doFill(b2) %>% rename(value.nongf = value,value = value.gf )
  b3 <- doFill(b3) %>% rename(value.nongf = value,value = value.gf )
  b4 <- doFill(b4) %>% rename(value.nongf = value,value = value.gf )
  b5 <- doFill(b5) %>% rename(value.nongf = value,value = value.gf )
  b6 <- doFill(b6) %>% rename(value.nongf = value,value = value.gf )
  b7 <- doFill(b7) %>% rename(value.nongf = value,value = value.gf )
}

# --- #
# Create interval
interv <- data.frame(SSBS = sites$SSBS)
if(tp == "midyear"){
  #interv$target_period <- interval( lubridate::floor_date(sites$Sample_midpoint,"year"), lubridate::ceiling_date(sites$Sample_midpoint,"year") )
  interv$start <- lubridate::floor_date(sites$Sample_midpoint,"year")
  interv$end <- lubridate::ceiling_date(sites$Sample_midpoint,"year")
  interv$int <- interval(interv$start,interv$end)
} else {
  #interv$target_period <- interval(sites$Sample_start_earliest - years(1), sites$Sample_start_earliest)
  interv$start <- sites$Sample_start_earliest - years(1)
  interv$end <- sites$Sample_start_earliest
  interv$int <- interval(interv$start,interv$end)
}
# Convert back to date
if(aggregate2month){b1$date <- as.Date(as.yearmon(b1$date));b2$date <- as.Date(as.yearmon(b2$date));b3$date <- as.Date(as.yearmon(b3$date));b4$date <- as.Date(as.yearmon(b4$date));b5$date <- as.Date(as.yearmon(b5$date));b6$date <- as.Date(as.yearmon(b6$date));b7$date <- as.Date(as.yearmon(b7$date))}

myLog("Interval subsetting per site")
# Subset bands to respective interval
sub_b1 <- merge(b1, interv, by = c("SSBS"))
sub_b1 <- subset(sub_b1,date %within% int)
#sub_b1 <- as.data.frame( sub_b1[date>= start & date<=end] ) # Old method
sub_b2 <- merge(b2, interv, by = c("SSBS"))
sub_b2 <- subset(sub_b2,date %within% int)
sub_b3 <- merge(b3, interv, by = c("SSBS"))
sub_b3 <- subset(sub_b3,date %within% int)
sub_b4 <- merge(b4, interv, by = c("SSBS"))
sub_b4 <- subset(sub_b4,date %within% int)
sub_b5 <- merge(b5, interv, by = c("SSBS"))
sub_b5 <- subset(sub_b5,date %within% int)
sub_b6 <- merge(b6, interv, by = c("SSBS"))
sub_b6 <- subset(sub_b6,date %within% int)
sub_b7 <- merge(b7, interv, by = c("SSBS"))
sub_b7 <- subset(sub_b7,date %within% int)
rm(b1,b2,b3,b4,b5,b6,b7) # clean up

# Get average of band values within times of sampling
results <- data.frame(SSBS = character(0),
                      # For period use 2 different types
                      # year of midyear of sampling | startyear
                      timeperiod = character(0), 
                      MODIS_version = character(0), # MODIS MCD32A4 version
                      propNA = numeric(0), # Missing data
                      # Average BRDF and vegetation index measurements
                      BRDF_Band1_mean = numeric(0), BRDF_Band2_mean = numeric(0), BRDF_Band3_mean = numeric(0),
                      BRDF_Band4_mean = numeric(0), BRDF_Band5_mean = numeric(0), BRDF_Band6_mean = numeric(0),BRDF_Band7_mean = numeric(0), 
                      NDVI_mean = numeric(0), NDVI_min = numeric(0),  NDVI_max = numeric(0), NDVI_cv = numeric(0),
                      NDVI_AUC = numeric(0),
                      EVI_mean = numeric(0), EVI_min = numeric(0),  EVI_max = numeric(0), EVI_cv = numeric(0),
                      EVI_AUC = numeric(0),
                      EVI2_AUC = numeric(0),
                      SAVI_mean = numeric(0), SAVI_min = numeric(0),  SAVI_max = numeric(0), SAVI_cv = numeric(0),
                      NDWI_mean = numeric(0), NDWI_min = numeric(0),  NDWI_max = numeric(0), NDWI_cv = numeric(0),
                      # Spectral heterogeneity
                      PCA_BRDF_variance12 = numeric(0),
                      PCA_BRDF_meancentroid = numeric(0) # Construct a PCA, calculate centroid within all (1-2) axes
                      # Spectral variability was then calculated as the mean of the Euclidean distances from the centroid of all principal components for each plot. Oldeland
                      # Also scale before hand to assess impact?
                      )

# Execute in parallel if possible
#for(siteid in unique(sites$SSBS)){
results2 <- foreach(siteid =  unique(sites$SSBS),
               .combine = rbind,
               .multicombine = FALSE,
               .errorhandling = 'pass',
               .packages = c('dplyr','lubridate','data.table'),
               .export = c('trapezoid','gapfill'),
               .verbose = FALSE) %dopar% {
  myLog(siteid)
  sub <- subset(sites,SSBS == siteid)
  
  # Subset by id
  subb_b1 <- subset(sub_b1,SSBS == siteid) # These subsets use the data.table structure first
  subb_b2 <- subset(sub_b2,SSBS == siteid)
  subb_b3 <- subset(sub_b3,SSBS == siteid)
  subb_b4 <- subset(sub_b4,SSBS == siteid)
  subb_b5 <- subset(sub_b5,SSBS == siteid)
  subb_b6 <- subset(sub_b6,SSBS == siteid)
  subb_b7 <- subset(sub_b7,SSBS == siteid)
  
  print(paste0("--> Processing following interval: ",tp))
  out <- data.frame(SSBS = as.character( siteid ), timeperiod = tp, MODIS_version = "006") # Output data.frame
  
  # ----------------------------- #
  # Collect and save output metrics

  # Missing values
  if(gapfill){
    # Adapted for non-gapfilled data
    out$propNA <- length(which(is.na(subb_b1$value.nongf))) / nrow(subb_b1)
  } else {
    out$propNA <- length(which(is.na(subb_b1$value))) / nrow(subb_b1)
  }
  
  # Average bands
  out$BRDF_Band1_mean <- mean(subb_b1$value,na.rm = T) * 0.0001
  out$BRDF_Band2_mean <- mean(subb_b2$value,na.rm = T) * 0.0001
  out$BRDF_Band3_mean <- mean(subb_b3$value,na.rm = T) * 0.0001
  out$BRDF_Band4_mean <- mean(subb_b4$value,na.rm = T) * 0.0001
  out$BRDF_Band5_mean <- mean(subb_b5$value,na.rm = T) * 0.0001
  out$BRDF_Band6_mean <- mean(subb_b6$value,na.rm = T) * 0.0001
  out$BRDF_Band7_mean <- mean(subb_b7$value,na.rm = T) * 0.0001
    
  ## Metrics
  # Combine all bands
  bands <- subb_b1 %>% rename(Band1 = value) %>% select(SSBS,date,Band1) %>% 
    left_join(.,(subb_b2 %>% rename(Band2 = value) %>% select(date,Band2) ),by= c("date")) %>% 
    left_join(.,(subb_b3 %>% rename(Band3 = value) %>% select(date,Band3) ),by= c("date")) %>% 
    left_join(.,(subb_b4 %>% rename(Band4 = value) %>% select(date,Band4) ),by= c("date")) %>% 
    left_join(.,(subb_b5 %>% rename(Band5 = value) %>% select(date,Band5) ),by= c("date")) %>% 
    left_join(.,(subb_b6 %>% rename(Band6 = value) %>% select(date,Band6) ),by= c("date")) %>% 
    left_join(.,(subb_b7 %>% rename(Band7 = value) %>% select(date,Band7) ),by= c("date"))
    
  # Correct values
  bands <- bands %>% mutate_at(vars(matches('Band')), function(x) x * 0.0001)
  
  # NDVI
  # (5-4) / (5+4)
  bands$NDVI <- (bands$Band2 - bands$Band1) / (bands$Band2 + bands$Band1)
  # -> Calc results
  out$NDVI_mean = mean(bands$NDVI, na.rm = T); out$NDVI_min = min(bands$NDVI, na.rm = T);  out$NDVI_max = max(bands$NDVI,na.rm = T); out$NDVI_cv = co.var(bands$NDVI)
  # Area under the curve
  out$NDVI_AUC = trapezoid(x = bands$date,y = bands$NDVI)
  # EVI
  # EVI = G * (NIR – RED)/(NIR + C1*RED - C2*BLUE + L))
  #G – Gain factor
  #L – Factor for canopy background adjustment
  #C1, C2: Coefficients for correcting aerosol influences from RED using BLUE
  #MODIS EVI algorithm: L = 1, G = 2.5, C1 = 6, C2 = 7.5
  bands$EVI <- 2.5 * ((bands$Band2 - bands$Band1) / (bands$Band2 + 6.0 * bands$Band1 - 7.5 * bands$Band3 + 1.0))
  # -> Calc results
  out$EVI_mean = mean(bands$EVI, na.rm = T); out$EVI_min = min(bands$EVI, na.rm = T);  out$EVI_max = max(bands$EVI,na.rm = T); out$EVI_cv = co.var(bands$EVI)
  # Area under the curve
  out$EVI_AUC = trapezoid(x = bands$date,y = bands$EVI)
  
  # EVI2
  # Using only two bands without blue
  bands$EVI2 <- 2.5 * ((bands$Band2 - bands$Band1) / (bands$Band2 + bands$Band1 + 1))
  # -> Calc results
  out$EVI2_mean = mean(bands$EVI2, na.rm = T); out$EVI2_min = min(bands$EVI2, na.rm = T);  out$EVI2_max = max(bands$EVI2,na.rm = T); out$EVI2_cv = co.var(bands$EVI2)
  # Area under the curve
  out$EVI2_AUC = trapezoid(x = bands$date,y = bands$EVI2)
  
  # SAVI
  # SAVI = (1 + L) * (NIR – RED)/(NIR + RED + L)
  bands$SAVI <- (1 + 0.5) * (bands$Band2 - bands$Band1) / (bands$Band2 + bands$Band1 + 0.5)
  # -> Calc results
  out$SAVI_mean = mean(bands$SAVI, na.rm = T); out$SAVI_min = min(bands$SAVI, na.rm = T);  out$SAVI_max = max(bands$SAVI,na.rm = T); out$SAVI_cv = co.var(bands$SAVI)
  
  # NDWI
  bands$NDWI <- (bands$Band2 - bands$Band5) / (bands$Band2 + bands$Band5)
  # -> Calc results
  out$NDWI_mean = mean(bands$NDWI, na.rm = T); out$NDWI_min = min(bands$NDWI, na.rm = T);  out$NDWI_max = max(bands$NDWI,na.rm = T); out$NDWI_cv = co.var(bands$NDWI)
  
  # H_sd
  # Standard deviation of the first axis of a PCA using all bands
  mod <- dplyr::select(bands,Band1:Band7) %>% subset(.,complete.cases(.)) 
  mod <- try(prcomp(mod,scale. = T),silent = T)
  if(class(mod)!="try-error") {
    # Calculate centroid 
    cent <- cbind(cent.PC1 = mean(mod$x[,"PC1"]), cent.PC2 =  mean(mod$x[,"PC2"]) )
    # Save the centroid value
    out$PCA_BRDF_centroid <- cent
    
    # Calculate pairwise euclidean distance matrix 
    d <- as.matrix( dist(rbind(cent, cbind(mod$x[,"PC1"],mod$x[,"PC2"]) ), method = "euclidean") )
    
    # Explained variance of first two axes
    out$PCA_BRDF_variance12 <- sum(summary(mod)$importance[2,1:2])
    out$PCA_BRDF_meancentroid <- mean( as.vector( as.matrix(d)[2:nrow(d),1] ) )
    
    rm(d,cent)
  } else{
    out$PCA_BRDF_variance12 <- NA; out$PCA_BRDF_meancentroid <- NA
  }
  #results <- rbind(results, out)
  return(out)
}
saveRDS(results2,paste0("MCD43A4_BRDF_center_computed_",tp,".rds"))
stopCluster(cl)
stop('Done')