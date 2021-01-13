# Idea:
# Lump studies together based on their taxonomic group and wider intercept
library(tidyverse)

# Load dataset
sites <- readRDS('resSaves/sites_diversity.rds') %>% 
  # Filter
  dplyr::select(TGrouping, Sampling_method, Sampling_effort_unit) %>% distinct() %>% 
  dplyr::filter(TGrouping != 'Other')
# Refactorize and arrange
sites$TGrouping <- factor(sites$TGrouping, levels = c('Plantae','Fungi','Invertebrates','Amphibia','Reptilia','Aves','Mammalia'))
sites <- sites %>% arrange(TGrouping)

table(sites$Sampling_method,sites$Sampling_effort_unit)

# Save
if(!file.exists('resSaves/SMTable1.csv')){
  write.csv2(sites, 'resSaves/SMTable1.csv',row.names = FALSE)
}
# Multiplier
# Time to day
# Space to ha
# If not already an area estimate, format to area
# Circle
0.0001*(pi*x^2)
