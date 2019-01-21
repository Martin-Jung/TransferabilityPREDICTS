# Package loading
library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
library(tidyr)
# Custom packages
library(roquefort)
library(marfunky)
myLog <- function(...) {
  cat(paste0("[Spectral] ", Sys.time(), " | ", ..., "\n"))
}
#### Prepare PREDICTS data to be uploaded as G-Fusion table ####
# Load
r <- readRDS("../../Data/diversity-2016-02-03-03-37-46.rds") 
r <- DropInvalidMetricsAndMethods(r)
r <- CorrectSamplingEffort(r)
sites <-SiteMetrics(diversity=r,
                        extra.cols=c("SSB","SSBS","Longitude","Latitude","Sample_start_earliest","Sample_end_latest",
                                     "Ecoregion","Biome","Country","UN_subregion","Site_name",
                                     "Sampling_method","Study_common_taxon","Max_linear_extent","Coordinates_precision_metres"
                                     ))
rm(r)
# Get suitable studies
sites <- sites %>% 
  mutate(startyear = year(ymd(Sample_start_earliest)), # Which year did sampling commense
         endyear = year(ymd(Sample_end_latest))) %>% # When did it end?
  dplyr::filter(startyear >= 2000) %>%  #Kickout everything that started before 2000
  distinct() # Get Unique fields
# How about coordinate precision ?
summary(sites$Coordinates_precision_metres)

# Get distinct 500m grid
library(rgdal)
library(sp)
library(raster)
library("snow")
library("parallel")
sp = subset(sites,select = c("SSBS","Longitude","Latitude"))
# Make spatial file
coordinates(sp) <- ~Longitude+Latitude
proj4string(sp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# read in MCD12Q1 for fishnet proxy
ref <- raster("R:/ecocon_d/shared/GIS_LandCover/MODIS_MCD12Q1/MCD12Q1_2001.tif")
# Harmonize projections
if( proj4string(sp) != proj4string(ref) ){
  warning("Projection of input files is not equal. Will try to reproject point file")
  if(is.na(proj4string(ref))) stop("CRS of RasterLayer must be set!")
  sp <- spTransform(sp,CRSobj = CRS(proj4string(ref)))
}

### Make cluster object
beginCluster( detectCores()-2 ) # leave two core for background processes
#extraxt point with df
df <- raster::extract(ref,sp,method="simple",cellnumbers=T,df=T)
endCluster()
R.utils::detachPackage("raster")
# Append to data
sites$cells <- df$cells 
sites$SameCellokay <- FALSE
# Now loop through all studies 
for( study in unique(sites$SS)) {
  print(study)
  sub <- subset(sites,SS==study)
  # Filter those studies out which fall into only one MODIS 500m cell
  check = colSums(with(sub,table(SSS,cells)))
  if(length(check)>1){
    sites$SameCellokay[which(sites$SS==study & sites$cells %in% names(check))] <- TRUE
  } else {
    print(paste("All sites of",study,"fall within one 500m cell -> Exclude"))
  }
}

# Then filter out those which only fall in one cell
sites <- sites %>% dplyr::filter(SameCellokay ==TRUE)

# Grouping correction
sites$TGrouping <- as.character(sites$Study_common_taxon)
sites$TGrouping[grep("Hymenoptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Insecta",x = sites$TGrouping,ignore.case = T)]<- "Invertebrates"
sites$TGrouping[grep("Chordata",x = sites$TGrouping,ignore.case = T)] <- "Other"
sites$TGrouping[grep("Animalia",x = sites$TGrouping,ignore.case = T)] <- "Other"
sites$TGrouping[grep("Formicidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Scarabaeidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Tracheophyta",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Strigiformes",x = sites$TGrouping,ignore.case = T)] <- "Aves"
sites$TGrouping[grep("Isoptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Coleoptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Anogeissus",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Poaceae",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Colubridae",x = sites$TGrouping,ignore.case = T)] <- "Reptilia"
sites$TGrouping[grep("Chiroptera",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Ascomycota",x = sites$TGrouping,ignore.case = T)] <- "Fungi"
sites$TGrouping[grep("Lepidoptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Bryophyta",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Sarcoptiformes",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[which(str_length(sites$TGrouping)==0)] <- "Other"
sites$TGrouping[grep("Bombus",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Apidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Arthropoda",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Drosophilidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Colletes floralis",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Phasianidae",x = sites$TGrouping,ignore.case = T)] <- "Aves"
sites$TGrouping[grep("Lophophorus impejanus",x = sites$TGrouping,ignore.case = T)] <- "Aves"
sites$TGrouping[grep("Gastropoda",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Araneae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Arachnida",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Clitellata",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Anura",x = sites$TGrouping,ignore.case = T)] <- "Amphibia"
sites$TGrouping[grep("Carabidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Hemiptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Isopoda",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Collembola",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Agaricomycetes",x = sites$TGrouping,ignore.case = T)] <- "Fungi"
sites$TGrouping[grep("Curculionidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Pongo pygmaeus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Squamata",x = sites$TGrouping,ignore.case = T)] <- "Reptilia"
sites$TGrouping[grep("Culicidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Phyllostomidae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Maerua subcordata",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Oryctolagus cuniculus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Arecaceae",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Pteropus tonganus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Nymphalidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Diptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Staphylinidae",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Opiliones",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Orthoptera",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Swietenia macrophylla",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Aenictus",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Dorylus",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Vespertilionidae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Primates",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Panthera pardus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Odocoileus virginianus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Cephalophus",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Geometridae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Rodentia",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Magnoliopsida",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Sciomyzidae",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
sites$TGrouping[grep("Liolaemus",x = sites$TGrouping,ignore.case = T)] <- "Reptilia"
sites$TGrouping[grep("Dolichopus",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Muridae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Soricidae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Lumbricidae",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Lecanoromycetes",x = sites$TGrouping,ignore.case = T)] <- "Fungi"
sites$TGrouping[grep("Clethrionomys gapperi",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Passeriformes",x = sites$TGrouping,ignore.case = T)] <- "Aves"
sites$TGrouping[grep("Dipteryx oleifera",x = sites$TGrouping,ignore.case = T)] <- "Plantae"
unique(sites$TGrouping)

# - # 
d <- sites
temp <- d %>% dplyr::select(Max_linear_extent,Sampling_method,Study_common_taxon) %>% 
  group_by(Study_common_taxon,Sampling_method) %>% 
  summarise(MLE_avg = mean(Max_linear_extent,na.rm=T))
d <- left_join(d,temp) # Join back
# Insert where empty
d$Max_linear_extent <- ifelse(is.na(d$Max_linear_extent),d$MLE_avg,d$Max_linear_extent)
d$MLE_avg <- NULL # Kickout the previously calc. average
# Then check again this time with Higher Taxa at the remaining empty values
temp <- d %>% dplyr::select(Max_linear_extent,Sampling_method,TGrouping) %>% 
  group_by(TGrouping,Sampling_method) %>% 
  summarise(MLE_avg = mean(Max_linear_extent,na.rm=T))
d <- left_join(d,temp)
d$Max_linear_extent <- ifelse(is.na(d$Max_linear_extent),d$MLE_avg,d$Max_linear_extent)
d$MLE_avg <- NULL # Kickout the previously calc. average
# Last try
temp <- d %>% dplyr::select(Max_linear_extent,Sampling_method,Study_common_taxon) %>% 
  group_by(Study_common_taxon) %>% 
  summarise(MLE_avg = mean(Max_linear_extent,na.rm=T))
d <- left_join(d,temp)
d$Max_linear_extent <- ifelse(is.na(d$Max_linear_extent),d$MLE_avg,d$Max_linear_extent)
d$MLE_avg <- NULL # Kickout the previously calc. average
rm(temp)

d <- subset(d,!is.na(d$Max_linear_extent))
# Kickout extreme outliers
# 95%  1414.214
d <- subset(d,Max_linear_extent <= quantile(d$Max_linear_extent,probs = .95,na.rm = T))

# Save finall site scores
saveRDS(d,"sites_center.rds")

