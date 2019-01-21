# DEPRECATED
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
b1 <- as.data.table(readRDS("Extracts/PREDICTS_center_MCD43A4_1.rds"))
b2 <- as.data.table(readRDS("Extracts/PREDICTS_center_MCD43A4_2.rds"))
b3 <- as.data.table(readRDS("Extracts/PREDICTS_center_MCD43A4_3.rds"))
b4 <- as.data.table(readRDS("Extracts/PREDICTS_center_MCD43A4_4.rds"))
b5 <- as.data.table(readRDS("Extracts/PREDICTS_center_MCD43A4_5.rds"))
b6 <- as.data.table(readRDS("Extracts/PREDICTS_center_MCD43A4_6.rds"))
b7 <- as.data.table(readRDS("Extracts/PREDICTS_center_MCD43A4_7.rds"))

sites <- readRDS("sites_center.rds")

#### Burnt Area Product ####
r <-  as.data.frame(data.table::fread("Extracts/MCD45A1_BurntArea_center.csv",showProgress = T))[,-1] %>% dplyr::select(-.geo,-latitude, -longitude) %>% 
  reshape2::melt(., id.vars = "SSBS") %>%
  distinct() %>% mutate(variable = as.character(variable)) %>% 
  mutate(year = as.numeric( str_sub(variable,1,4)), # Get year out of file name
         month = as.numeric( str_sub(variable,6,7) ), # get month
         day = as.numeric( str_sub(variable,9,10) )) %>%  # get day
  # Make a date column
  mutate(date = ymd(paste(year,month,day,sep="-")))
# Kickout all those that are equal or bigger than 900 (no or not possible detection)
r <- dplyr::filter(r,value < 900)

# Create a product of 
sites <- readRDS("sites_center.rds")

# Get average of band values within times of sampling
results <- data.frame(SSBS = sites$SSBS,
                      TimeSinceLastFire = NA,Fire_nr = NA, Fire_freq_avg = NA, Fire_freq_sd = NA)

for(siteid in results$SSBS){
  myLog(siteid)
  sub <- subset(sites,SSBS == siteid)
  sub_f <- subset(r,SSBS == siteid)
  if(nrow(sub_f)==0) next() # If only invalid dates are within the period
  
  # Subset to those one months before the start of sampling
  i <- interval("2000-01-01", sub$Sample_start_earliest-months(1),tzone = "UTC")
  sub_f <- sub_f[which(ymd(sub_f$date) %within% i),]
  
  # Number of fires
  # Simply calculate the number of fires greater than zero
  results$Fire_nr[which(results$SSBS==siteid)] <- sum(sub_f$value>0)
  
  # The following metrics are only useful if fires are actually detected
  if(sum(sub_f$value>0)>0){
    
    # Time since last fire
    # Calculate the time before sampling start since a fire was detected
    d = max(which(sub_f$value>0)) # Get the last entry= Values are sorted by variable
    lf = as.Date(sub_f$value[d] - 1, origin = paste0(sub_f$year[d],"-01-01"))
    results$TimeSinceLastFire[which(results$SSBS==siteid)] <- 
      abs(as.numeric(gsub("//D",replacement = "",difftime(lf,sub$Sample_start_earliest)))) # Geta abs. difference in days
    
    if(sum(sub_f$value>0)>1){
      
      # Frequency of fires average
      # If more than two fires are detected, calculate average distance in days between them
      d = (which(sub_f$value>0)) # Get the last entry= Values are sorted by variable
      lf = as.Date(sub_f$value[d] - 1, origin = paste0(sub_f$year[d],"-01-01"))
      val <- vector()
      # Terrible inefficient code
      for(i in 1:length(lf)){
        val <- c(val,
                 abs(as.numeric(gsub("//D",replacement = "",difftime(lf[i],lf[i+1])))) 
        )
      }
      # Then calcualte the mean and sd
      results$Fire_freq_avg[which(results$SSBS==siteid)] <- mean(val,na.rm = T)
      # Frequency of fires sd
      results$Fire_freq_sd[which(results$SSBS==siteid)] <- sd(val,na.rm = T)
      # If more than two fires are detected, calculate standard deviation of distance in days between them
    }
  }
}
saveRDS(results,"MODIS_BurnEstimates_center.rds")
