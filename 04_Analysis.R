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
library(patchwork)
library(colorspace)
# Modelling
library(brms);library(rstan)
library(tidybayes)
library(bayesplot)
library(shinystan)
source('000_HelperFunction.R')
# Hacks to make BRMS faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)
# ------------------------------------ #
# Parameters 

# Path to PREDICTS data
sites <- readRDS('resSaves/sites_diversity.rds')

# Load remote-sensing data 
rs <- readRDS('resSaves/MCD43A4_BRDF_center_computed_yearbefore.rds') %>% 
  # Remove sites with too much missing data
  filter(propNA <= 0.5) %>% 
  # Remove those with no positive EVI, which likely fall in water, or where the AUC is smaller than 0
  filter(EVI2_mean >0, EVI2_AUC > 0)
# Missing data in non-gap filled estimates
mean(rs$propNA);sd(rs$propNA)
# Variance explained of PCA axes 1 & 2
mean(rs$PCA_BRDF_variance12,na.rm = T); sd(rs$PCA_BRDF_variance12,na.rm = T)
# Correlation between metrics
cor.test(rs$EVI2_mean,rs$PCA_BRDF_meancentroid,method = 'pear')

# Figure path
figure_path <- 'figures'
# Intermediate results
output_path <- 'resSaves'
dir.create(output_path,showWarnings = FALSE)
dir.create(figure_path,showWarnings = FALSE)

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


# ------------------------------------ #
#### SI Figure 1 - Variance

df <- rs %>% inner_join(sites, ., by = 'SSBS') %>% 
  dplyr::select(Biome,EVI2_mean,PCA_BRDF_meancentroid,PCA_BRDF_variance12) %>% drop_na()
biome.cols = c(
               "Boreal Forests/Taiga" = "#00ffd0",
               "Temperate Broadleaf & Mixed Forests" = "#97B669",
               "Temperate Conifer Forests" = "#75A95E",
               "Temperate Grasslands, Savannas & Shrublands" = "#FCD57A",
               "Mediterranean Forests, Woodlands & Scrub" = "#D16E3F",
               "Montane Grasslands & Shrublands" = "#c0edec",
               "Tropical & Subtropical Moist Broadleaf Forests" = "#317A22",
               "Tropical & Subtropical Grasslands, Savannas & Shrublands" = "#d1ffad",
               "Tropical & Subtropical Dry Broadleaf Forests" = "#A09700",
               "Tropical & Subtropical Coniferous Forests" = "#ccf2a0",
               "Deserts & Xeric Shrublands" =  "#DCBB50",
               "Mangroves" = "#e09dc2",
               "Flooded Grasslands & Savannas" = "#ebf0d5",
               "Tundra" = "#C1E1DD")
df$Biome <- droplevels(df$Biome)
df$Biome <- factor(df$Biome,levels = names(biome.cols))

g <- ggplot(df, aes(x = PCA_BRDF_meancentroid, y = EVI2_mean, fill = Biome)) +
  theme_few(base_size = 20) +
  geom_point(size = 2, alpha = .8, shape  = 21,colour = "grey20",stroke = 1) +
  scale_fill_manual(values = biome.cols,guide = guide_legend(ncol = 2,override.aes = list(size = 4))) + 
  theme(legend.position = 'bottom',legend.text = element_text(size = 10)) +
  scale_x_continuous(breaks = pretty_breaks(5)) +
  scale_y_continuous(breaks = pretty_breaks(5)) +
  labs(x = 'Spectral heterogeneity', y = 'Photosynthetic activity')
ggsave(filename = paste0(figure_path,'/','SIFigure1.png'),plot = g,width = 10,height = 10,dpi = 400)
