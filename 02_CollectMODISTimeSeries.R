library(lubridate)
library(tidyverse)
library(jsonlite)
library(stringr)
library(data.table)
library(assertthat)
library(jsonlite)
source("000_HelperFunction.R")
options(dplyr.show_progress = T)
# Path to the extracted time seris from Google Earth Engine
#path = "C:/LocalStorage/LocalFolder/PREDICTS_MCDv6"
path = "/lustre/scratch/lifesci/mj291/PREDICTS_MCDv6"

# The combinations how they were processed in GEE
year_combinations = c("2000_2006","2006_2012","2012_2014")

# Load in al MODIS BRDF Bands
readInFormatv2 <- function(x,idv = "SSBS"){
  val =  as.data.frame(fromJSON(x,flatten = T)) %>% 
    dplyr::select(-features.properties.longitude, -features.properties.latitude,-features.geometry,-features.id,-type,-features.type) %>% 
    reshape2::melt(., id.vars = paste0("features.properties.",idv)) %>%
    distinct() %>% mutate(variable = as.character(variable)) %>% 
    dplyr::rename("SSBS" = paste0("features.properties.",idv) ) %>% 
    mutate(year = as.numeric( str_sub(variable,21,24)), # Get year out of file name
           month = as.numeric( str_sub(variable,26,27) ), # get month
           day = as.numeric( str_sub(variable,29,30) )) %>%  # get day
    # Make a date column
    mutate(date = ymd(paste(year,month,day,sep="-"))) %>% 
    # Remove variable
    dplyr::select(-variable)
  # Reorder
  val <- val[order(val$date,decreasing = F),]
  return(val)
}

for(b in seq(1,7)){
  myLog("Get all BRDF values for band ", b)
  
  out <- data.frame()
  
  r <- readInFormatv2(paste0(path,"/PREDICTS_center_MCD43A4_Band",b,"_",year_combinations[1],".geojson.json"),idv = "SSBS") %>% mutate(Band = b)
  out <- rbind(out, r)
  rm(r);gc() # Clean
  r <- readInFormatv2(paste0(path,"/PREDICTS_center_MCD43A4_Band",b,"_",year_combinations[2],".geojson.json"),idv = "SSBS") %>% mutate(Band = b)
  out <- rbind(out, r)
  rm(r);gc() # Clean
  r <- readInFormatv2(paste0(path,"/PREDICTS_center_MCD43A4_Band",b,"_",year_combinations[3],".geojson.json"),idv = "SSBS") %>% mutate(Band = b)
  out <- rbind(out, r)
  rm(r);gc() # Clean
  out <- unique( out )
  
  saveRDS(out,paste0(path,"/MCD43A4_Band",b,".rds"))
  rm(out)
  myLog("DONE")
}
stop("Done!")
