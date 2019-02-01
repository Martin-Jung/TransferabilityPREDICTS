# Package loading
library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
library(tidyr)
# Custom packages
library(roquefort)
library(yarg)
source("000_HelperFunction.R")
myLog <- function(...) {
  cat(paste0("[Spectral] ", Sys.time(), " | ", ..., "\n"))
}
#### Prepare PREDICTS data to be uploaded as G-Fusion table ####
# Load
database <- readRDS("../../Data/PREDICTS_v1/database.rds") 
sites <- readRDS("../../Data/PREDICTS_v1/sites.rds")

# Calculate sitebased biodiversity metrics
diversity <- database %>% group_by(SS,SSBS) %>% 
  summarise(
    Species.richness = n_distinct(Best_guess_binomial),
    Total.abundance = sum(Measurement, na.rm = TRUE ),
    Total.abundance.effcor = sum(Effort_corrected_measurement, na.rm = TRUE )
  )

# Sorensen dissimilarity
dis.sor <- CompDissim2(database,"SorVeg",binary = T)
dis.bc <- database %>% mutate(Measurement = Effort_corrected_measurement) %>% CompDissim2(.,"BCVeg",binary = F)

sites %>% 

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

