# Package loading
library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
library(tidyr)
source("00_HelperFunction.R")

#### Prepare PREDICTS data to be uploaded as G-Fusion table ####
# Load
sites <- readRDS("../../../PhD/Data/PREDICTS_v1/sites.rds")
database <- readRDS("../../../PhD/Data/PREDICTS_v1/database.rds") 
sites <- sites %>% dplyr::filter(SSS %in%  unique(sites$SSS)[which(unique(sites$SSS) %in% unique(database$SSS))] )

# Calculate sitebased biodiversity metrics
diversity <- database %>% dplyr::filter(!is.na(Sampling_effort)) %>%  # Remove sites with no recorded sampling effort
  mutate(SSS = droplevels(SSS)) %>% 
  # Use the effort corrected abundance as measurement following Newbold et al. (2014)
  dplyr::select(-Measurement) %>% dplyr::rename(Measurement = "Effort_corrected_measurement") %>% 
  SiteMetrics(diversity = ., extra.cols=c("SSBS"), sites.are.unique = T, srEstimators = T)

# Sorensen dissimilarity
dis.sor <- CompDissim2(database,"SorVeg",binary = T)
dis.bc <- database %>% mutate(Measurement = Effort_corrected_measurement) %>% CompDissim2(.,"BCVeg",binary = F)

# Calculate PIE
database$Is_abundance <- "Abundance" == database$Diversity_metric_type
site.abundance <- tapply(database$Is_abundance, database$SSS, unique)
sites$PIE <- SiteHurlbertsPie(database, site.abundance)

# ---- # 
# Taxonomic grouping and recategorization
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
sites$TGrouping[grep("Diprotodontia",x = sites$TGrouping,ignore.case = T)] <- "Mammalia"
sites$TGrouping[grep("Nematoda",x = sites$TGrouping,ignore.case = T)] <- "Invertebrates"
sites$TGrouping[grep("Glomeromycetes",x = sites$TGrouping,ignore.case = T)] <- "Fungi"
sites$TGrouping[grep("Strabomantidae",x = sites$TGrouping,ignore.case = T)] <- "Amphibia"
# Assign those in other
#unique(sites$Higher_taxa[which(sites$TGrouping == "Other")])
sites$TGrouping[ intersect( grep("Mammalia",x = sites$Higher_taxa,ignore.case = T), which(sites$TGrouping=="Other") )] <- "Mammalia"
sites$TGrouping[ intersect( grep("Coleoptera",x = sites$Higher_taxa,ignore.case = T), which(sites$TGrouping=="Other") )] <- "Invertebrates"
sites$TGrouping[ intersect( grep("Nematoda",x = sites$Higher_taxa,ignore.case = T), which(sites$TGrouping=="Other") )] <- "Invertebrates"
sites$TGrouping[ intersect( grep("Annelida",x = sites$Higher_taxa,ignore.case = T), which(sites$TGrouping=="Other") )] <- "Invertebrates"
sites$TGrouping[ intersect( grep("Aves",x = sites$Higher_taxa,ignore.case = T), which(sites$TGrouping=="Other") )] <- "Aves"
sites$TGrouping[ intersect( grep("Arachnida",x = sites$Higher_taxa,ignore.case = T), which(sites$TGrouping=="Other") )] <- "Invertebrates"
sites$TGrouping[ intersect( grep("Liliopsida,Magnoliopsida,Polypodiopsida",x = sites$Higher_taxa,ignore.case = T), which(sites$TGrouping=="Other") )] <- "Plantae"
sites$TGrouping[ intersect( grep("Liliopsida,Lycopodiopsida,Magnoliopsida,Polypodiopsida",x = sites$Higher_taxa,ignore.case = T), which(sites$TGrouping=="Other") )] <- "Plantae"
sites$TGrouping[ intersect( grep("Bryophyta,Liliopsida,Magnoliopsida,Pinopsida,Polypodiopsida",x = sites$Higher_taxa,ignore.case = T), which(sites$TGrouping=="Other") )] <- "Plantae"
sites$TGrouping[ intersect( grep("Equisetopsida,Liliopsida,Magnoliopsida,Pinopsida,Polypodiopsida",x = sites$Higher_taxa,ignore.case = T), which(sites$TGrouping=="Other") )] <- "Plantae"
sites$TGrouping[ intersect( grep(",Ascomycota,Basidiomycota,Mycetozoa",x = sites$Higher_taxa,ignore.case = T), which(sites$TGrouping=="Other") )] <- "Fungi"
#table(sites$TGrouping)

# The rest do not have a clear association to a single higher group and will be removed
sites <- sites %>% dplyr::filter(TGrouping != "Other")

#### Sampling extent approximation ####
# - # 
d <- sites
d$Max_linear_extent_metres_predicted <- is.na(d$Max_linear_extent_metres) # template column
temp <- d %>% dplyr::select(Max_linear_extent_metres,Sampling_method,Study_common_taxon) %>% 
  group_by(Study_common_taxon,Sampling_method) %>% 
  summarise(MLE_avg = mean(Max_linear_extent_metres,na.rm=T))
d <- left_join(d,temp) # Join back
# Insert where empty
d$Max_linear_extent_metres <- ifelse(is.na(d$Max_linear_extent_metres),d$MLE_avg,d$Max_linear_extent_metres)
d$MLE_avg <- NULL # Kickout the previously calc. average
# Then check again this time with Higher Taxa at the remaining empty values
temp <- d %>% dplyr::select(Max_linear_extent_metres,Sampling_method,TGrouping) %>% 
  group_by(TGrouping,Sampling_method) %>% 
  summarise(MLE_avg = mean(Max_linear_extent_metres,na.rm=T))
d <- left_join(d,temp)
d$Max_linear_extent_metres <- ifelse(is.na(d$Max_linear_extent_metres),d$MLE_avg,d$Max_linear_extent_metres)
d$MLE_avg <- NULL # Kickout the previously calc. average
# Last try
temp <- d %>% dplyr::select(Max_linear_extent_metres,Sampling_method,Study_common_taxon) %>% 
  group_by(Study_common_taxon) %>% 
  summarise(MLE_avg = mean(Max_linear_extent_metres,na.rm=T))
d <- left_join(d,temp)
d$Max_linear_extent_metres <- ifelse(is.na(d$Max_linear_extent_metres),d$MLE_avg,d$Max_linear_extent_metres)
d$MLE_avg <- NULL # Kickout the previously calc. average
rm(temp)

#d <- subset(d,!is.na(d$Max_linear_extent))
# Kickout extreme outliers
# 95%  1414.214
#d <- subset(d,Max_linear_extent <= quantile(d$Max_linear_extent,probs = .95,na.rm = T))

#### Create and join in new grouping ####
#source('01a_MethodologyIntercept.R')
stopifnot(file.exists('resSaves/SMTable1.csv'))
sm <- read_delim('resSaves/SMTable1.csv',delim = ',',na = 'NA')
# Assert that columns are present
assertthat::assert_that(
  assertthat::has_name(sm,'Sampling_grouping'),assertthat::has_name(sm,'Grouping_unit') 
)
# Formatting
sm$Sampling_grouping <- factor(sm$Sampling_grouping,levels = c('FixedPlot','Transect','Survey','Trap','Netting','Fogging') )
sm$Grouping_unit <- factor(sm$Grouping_unit, levels = c('Space','Time','Plot'))
# Join in
sites <- dplyr::left_join(sites, sm, by = c('TGrouping','Sampling_method','Sampling_effort_unit'))
# New grouping
sites$TransferGrouping <- paste0(sites$TGrouping,'_',sites$Sampling_grouping,'_',sites$Grouping_unit)
saveRDS(sites,"resSaves/sites_diversity.rds")

#### Join with biodiversity estimates ####

# Remove columns present in both dataframes
d <- d %>% dplyr::select(-Study_number,-Study_name,-Diversity_metric_type,-Site_number,-Site_name,-Block,-Use_intensity,
                         -Sample_start_earliest,-Sample_end_latest, -Longitude, -Latitude)
d <- left_join(d,diversity,by = c("Source_ID","SS","SSS","SSBS"))

# Remove those without coordinates
d <- d[-which(is.na(d$Longitude)),]

# Save final site scores
saveRDS(d,"resSaves/sites_diversity.rds")
saveRDS(dis.sor,"resSaves/sites_pairwise_sorensen.rds")
saveRDS(dis.bc,"resSaves/sites_pairwise_bc.rds")

# Save a spatial file
d <- readRDS("resSaves/sites_diversity.rds") %>% dplyr::select(SS,SSBS,Longitude,Latitude)
library(sp);library(plotKML);library(rgdal)
coordinates(d) <- ~Longitude + Latitude
proj4string(d) <- "+proj=longlat +datum=WGS84"
rgdal::writeOGR(d,"sites_diversity.kml",layer = "sites_diversity",driver = "KML")


# --- #
# Calculate distance between pairs of sites in long format
library(geosphere)
o <- data.frame()
for(study in unique(d$SS)){
  sub <- subset(d,SS == study)
  if(nrow(sub)<=2) next()
  # Distance in km
  dd <- distm(cbind(sub$Longitude,sub$Latitude), fun = distHaversine)/1000
  dd[upper.tri(dd)] <- NA;diag(dd) <- NA
  rownames(dd) <- sub$SSBS; colnames(dd) <- sub$SSBS
  dd <- reshape2::melt(dd) %>% tidyr::drop_na() %>% dplyr::rename(distance = value)
  dd <- dd %>% dplyr::mutate(SS = study)
  o <- bind_rows(o, dd)
}

# Save
saveRDS(o,paste0('resSaves/pairwiseSiteDistance.rds'))

# --- #
# Average distance to study centroid
cent <- d %>% dplyr::group_by(SS) %>% dplyr::summarise(clo = mean(Longitude),cla = mean(Latitude))

library(geosphere)
o <- data.frame()
for(study in unique(d$SS)){
  sub <- subset(d,SS == study)
  sub.c <- subset(cent, SS == study)
  if(nrow(sub)<=2) next()
  # Distance in km
  dd <- distm(cbind(sub$Longitude,sub$Latitude), sub.c[,c('clo','cla')], fun = distHaversine)/1000
  dd <- as.data.frame(dd) %>% dplyr::rename(distance = 'V1')
  dd$distance <- normalize(dd$distance)
  dd$SSBS <- sub$SSBS; dd$SS <- study
  o <- bind_rows(o, dd)
}

# Save
saveRDS(o,paste0('resSaves/centroidDistance.rds'))

