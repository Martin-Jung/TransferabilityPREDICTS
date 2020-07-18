# Data science 
library(tidyverse)
library(reshape2)
library(stringr)
library(lubridate)
library(assertthat)
# Visualization
library(ggplot2)
library(ggthemes)
library(scales)
library(colorspace)
# Modelling
library(brms)
library(tidybayes)
library(bayesplot)
library(shinystan)
source('00_HelperFunction.R')
# Hacks to make BRMS faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# ------------------------------------ #
# Parameters 

# Path to PREDICTS data
sites <- readRDS('resSaves/sites_diversity.rds')

# Load remote-sensing data 
rs <- readRDS('resSaves/MCD43A4_BRDF_center_computed_yearbefore.rds') %>% 
  drop_na() %>% 
  # Remove sites with too much missing data
  filter(propNA < 0.5)
mean(rs$propNA);sd(rs$propNA)

# Figure path
figure_path <- 'figures'
# Intermediate results
output_path <- 'resSaves'

# ------------------------------------ #
#### Testing ####

# Test with bird point counts
bird_points <- sites %>% 
  dplyr::filter(Sampling_method == 'point counts',TGrouping == 'Aves') %>% # Filter
  mutate(logabundance = log10(Total_abundance+1)) %>% 
  left_join(., rs, by = 'SSBS')

plot(bird_points$logabundance~bird_points$NDVI_mean,col = bird_points$Biome)
subset_map(bird_points)

fit <- lme4::glmer(logabundance ~ EVI2_mean + (SS|Biome),data = bird_points, family = 'gaussian')
MuMIn::r.squaredGLMM(fit)


