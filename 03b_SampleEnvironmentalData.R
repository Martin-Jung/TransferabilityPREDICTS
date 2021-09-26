library(raster)
library(exactextractr)
library(sf)
library(dplyr)
# --- #
output_path <- 'resSaves'
dir.create(output_path,showWarnings = FALSE)
# --- #

# Path to PREDICTS data
sites <- readRDS('resSaves/sites_diversity.rds') %>% 
  dplyr::select(SSBS, Longitude, Latitude)

# Accessibility to cities layer
# Rational -> Species more disturbance afinne, earlier surveyed
acc <- raster('acc_50k.tif')

# Topographic ruggedness index from Earthenv
# Rational --> More rugged habitat, more likely species did not get spotted
tri <- raster('tri_1KMmn_GMTEDmd.tif')

# Mean annual Cloud clover from earth env
# Rational --> Likely more degraded earth env sites
cc <- raster('MODCF_meanannual.tif')


ex1 <- raster::extract(acc, sites[,c('Longitude','Latitude')] )

ex2 <- raster::extract(tri, sites[,c('Longitude','Latitude')] )

ex3 <- raster::extract(cc, sites[,c('Longitude','Latitude')] )

# Write data
sites$accessibility <- ex1
sites$tri <- ex2
sites$cloud <- ex3

sites <- sites <- sites %>% dplyr::select(-Longitude, Latitude)

saveRDS(sites, paste0(output_path,'biasvariables.rds'))
