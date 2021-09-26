library(tidyverse)
library(ggplot2)
library(extrafont)
extrafont::fonts()
library(assertthat)
library(scico)
library(ggthemes)
library(scales)
library(lubridate)
library(stringr)
source('00_HelperFunction.R')
# Figure path
figure_path <- 'figures'
# Intermediate results
output_path <- 'resSaves'
dir.create(output_path,showWarnings = FALSE)
dir.create(figure_path,showWarnings = FALSE)
# --- #

# Path to PREDICTS data
sites <- readRDS('resSaves/sites_diversity.rds') %>% 
  dplyr::mutate(logabund = log10(Total_abundance+1),
                asinPIE = asin(sqrt(PIE)),
                startyear = year(Sample_start_earliest),
                sampling_duration = as.vector(
                  difftime(Sample_end_latest,Sample_start_earliest,units = 'days')
                )) %>% 
  dplyr::mutate(Predominant_land_use = fct_collapse(Predominant_land_use,
                'Primary vegetation' = 'Primary vegetation',
                'Secondary vegetation' = c('Mature secondary vegetation','Intermediate secondary vegetation',
                                           'Young secondary vegetation','Secondary vegetation (indeterminate age)'),
                'Plantation forest' = 'Plantation forest',
                'Pasture' = 'Pasture', 'Cropland' = 'Cropland',
                'Urban' = 'Urban', 'Cannot decide' = 'Cannot decide'
                                                          ))

# Add bias variable results
sites <- sites %>% dplyr::left_join(.,
                                    readRDS(paste0(output_path,'biasvariables.rds'))
                                    )

# Predictability results
results_predict1 <- readRDS('resSaves/results_predicatability_nopooling.rds') %>% 
  dplyr::mutate(term = fct_collapse(term,
    "Photosynthetic activity" = c("EVI2_mean","EVIdis"),
    "Spectral variability" = c("PCA_BRDF_meancentroid","PCAcent")
  )) %>% 
  dplyr::left_join(., sites %>% dplyr::select(SS, TGrouping, Biome,Predominant_land_use) %>% distinct(), by = 'SS' )

results_predict2 <- readRDS('resSaves/results_predicatability_pooling.rds') %>% 
  dplyr::mutate(term = fct_collapse(term,
                                    "Photosynthetic activity" = c("EVIdis"),
                                    "Spectral variability" = c("PCAcent")
  ))

results_transfer <- readRDS('resSaves/results_transferability.rds') %>% 
  dplyr::filter(dataset == 'test') %>% 
  dplyr::group_by(SS,metric,term,grouping) %>% 
    dplyr::summarise(
      mape.avg = mean(mape,na.rm = TRUE), smape.avg = mean(smape,na.rm = TRUE),
      mape.sd = sd(mape,na.rm = TRUE), smape.sd = sd(smape,na.rm = TRUE),
    ) %>% ungroup() %>% 
  dplyr::mutate(term = fct_collapse(term,
                                    "Photosynthetic activity" = c("EVIdis"),
                                    "Spectral variability" = c("PCAcent")
  )) %>% 
  dplyr::left_join(., sites %>% dplyr::select(SS, TGrouping, Biome,Predominant_land_use) %>% distinct(), by = 'SS' )

# For calculating average remote sensing measures
rs <- readRDS('resSaves/MCD43A4_BRDF_center_computed_yearbefore.rds') %>% 
  # Remove sites with too much missing data
  filter(propNA <= 0.5) %>% 
  # Remove those with no positive EVI, which likely fall in water, or where the AUC is smaller than 0
  filter(EVI2_mean > 0, EVI2_AUC > 0) %>% left_join(., sites %>% dplyr::select(SS,SSBS)) %>% 
  group_by(SS) %>% 
    dplyr::summarise(EVI2_mean = mean(EVI_mean, na.rm = TRUE),
                     SpecHetero_mean = mean(PCA_BRDF_meancentroid,na.rm = TRUE),
                     propNA_mean = mean(propNA, na.rm = TRUE)) %>% ungroup()

# ------------------------- #
#### Figure 1 ####
library(GGally)
library(MASS)
# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
# Spectral vs evi on axes
# 4 scatter plots of R2 with Abline
xx <- results_predict1 %>% dplyr::select(term, metric, SS,TGrouping, full_r2) %>% distinct() %>% tidyr::drop_na(full_r2) %>%
  mutate( metric = factor(metric,levels = c('SR','LA','PIE','SOR'),
                              labels = c('Species richness', 'Species abundance', 'Assemblage evenness', 'Assemblage composition'))
  ) %>% 
  # Also remove infinite values
  dplyr::filter(!is.infinite(full_r2))
# Calculate average per metric
xx %>% dplyr::group_by(term) %>% dplyr::summarise(m = mean(full_r2,na.rm = T), s = sd(full_r2, na.rm=T))
xx %>% dplyr::group_by(term,metric) %>% dplyr::summarise(m = mean(full_r2,na.rm = T), s = sd(full_r2, na.rm=T))
xx %>% dplyr::group_by(metric) %>% summarise(quantile(full_r2,.25))
# Split
xx <- xx %>% 
  # Spread R2 out per metric  
  tidyr::pivot_wider(id_cols = c('SS','metric','TGrouping'),names_from = 'term',values_from = 'full_r2') %>% 
  # Join in number of sites per study
  left_join(., sites %>% dplyr::group_by(SS) %>% summarise(N = n())) %>% tidyr::drop_na()
# Correlation strength
xx %>% dplyr::group_by(metric) %>% summarise(cor(`Spectral variability`, `Photosynthetic activity`))

# Calculate point density for colouring density of points within each group
xx <- xx %>% dplyr::group_by(metric) %>% 
  mutate(density = get_density(`Spectral variability`, `Photosynthetic activity`, n = 100))

# Build plot
g <- ggplot(xx, aes(x = `Spectral variability`,y = `Photosynthetic activity`, size = N, colour = density)) +
  theme_classic(base_size = 18) +
  geom_abline(slope = 1) +
  geom_point(alpha = .5) +
    scale_size_binned(range = c(.5,5), trans = 'log10',guide = guide_bins(title = 'N (log)')) +
    scale_colour_viridis_c(option = 'E',guide = guide_colourbar(title = 'Density')) +
  facet_wrap(~metric) + 
    theme(panel.spacing = unit(1, "lines"), axis.text.x.bottom = element_text(angle = 90, vjust = 0.5)) +
    scale_x_continuous(expand = c(0,0)) +
  labs( x = expression(paste(R^{2}," for Spectral variability")),
        y = expression(paste(R^{2}," for Photosynthetic activity"))
  )
# g
ggsave(plot =g ,filename = 'figures/Figure2_R2.png',width = 12,height = 6,dpi = 300)

# Also map the R2 for the SI
xx <- results_predict1 %>% # Also remove infinite values
  dplyr::filter(!is.infinite(full_r2)) %>% tidyr::drop_na(full_r2) %>% 
  dplyr::group_by(SS,metric) %>% dplyr::summarise(full_r2 = mean(full_r2)) %>% 
  mutate( metric = factor(metric,levels = c('SR','LA','PIE','SOR'),
                          labels = c('Species richness', 'Species abundance', 'Assemblage evenness', 'Assemblage composition'))
  ) %>% 
  left_join(., sites %>% dplyr::group_by(SS) %>% dplyr::summarise(Longitude = mean(Longitude,na.rm=T),
                                                                  Latitude = mean(Latitude,na.rm=T)) %>% distinct())
# Now show R2 on a map
#### SI Figure 2 - map ####
library(tmap)
data("World")
# Load NE
ne <- sf::st_read('C:/Users/Martin/Downloads/Ecoregions2017/Ecoregions2017.shp') %>% 
  dplyr::filter(REALM != 'Antarctica')
# Convert to point
xx <- st_as_sf(xx, coords = c('Longitude', 'Latitude'),crs = st_crs(World))

gg <- tm_shape(World) +
    tm_polygons(col = 'grey90', alpha = .25) +
  # tm_polygons(col = 'COLOR_BIO',border.col = NA, alpha = .25) +
tm_shape(xx) +
    tm_bubbles(size = 'full_r2',
               col = 'full_r2',
               border.col = "black", border.alpha = .5, 
               palette = 'cividis',
               # style="fixed", breaks=c(-Inf, seq(0, 6, by=2), Inf),
               title.size= expression(R^2), 
               title.col=""
    ) +
  tm_facets(by = 'metric') +
tm_format("World") +
tm_layout(panel.label.size = 2,legend.text.size = 1,legend.title.size = 1.5)
# gg
tmap::tmap_save(tm = gg,filename = 'figures/SIFigure2.png',
                # width=16, height=8, units = 'in',
                dpi = 300, asp = 0)

# -- #
# pairs plot
# r2 of model for each measure
# colour by metric
# # Combine estimates
# xx <- results_predict1 %>% dplyr::select(term, metric,SS, full_r2) %>% distinct() %>% tidyr::drop_na(full_r2) %>%
#   # Spread R2 out per metric  
#   tidyr::pivot_wider(id_cols = c('SS','term'),names_from = 'metric',values_from = 'full_r2')# %>%
#   # mutate(across(where(is.numeric), function(x) asin(sqrt(x)) ))
# 
# pm <- ggpairs(xx, 3:6,
#               mapping = ggplot2::aes(color = term),
#               diag = list(continuous = 'barDiag'),
#               upper = list(continuous = wrap(ggally_cor, displayGrid = FALSE)),
#               lower = list(continuous = "points", combo = "dot_no_facet")
#               ) + 
#   theme_classic(base_size = 18,base_family = 'Arial') +
#     scale_colour_manual(values = c('darkgreen','#D2691E'),guide = guide_legend('')) +
#     scale_fill_manual(values = c('darkgreen','#D2691E'),guide = guide_legend(''))
# pm


#### Figure 2 ####
# Show MAPE error across metrics and terms for predictability and transferability

# Combine estimates
xx <- bind_rows(
  results_predict1 %>% dplyr::select(term, metric, cv_smape) %>% dplyr::rename(smape = cv_smape) %>% 
    dplyr::mutate(method = 'Predictability') %>% distinct(),
  results_transfer %>% dplyr::select(SS,metric, term,smape.avg) %>% dplyr::rename(smape = smape.avg) %>% 
    dplyr::mutate(method = 'Transferability') %>% distinct()
#  results_predict2 %>% dplyr::select(term,metric, smape) %>% 
#    dplyr::mutate(method = 'Pooling')
)
xx$metric <- factor(xx$metric,levels = c('SR','LA','PIE','SOR'),
                    labels = c('Species richness', 'Species abundance', 'Assemblage evenness', 'Assemblage composition'))

g <- ggplot(xx, aes(x = fct_rev(metric), y = smape, fill = method)) +
  theme_classic(base_size = 18,base_family = 'Arial') +
    theme(panel.grid.major =  element_line(size = .5,colour='grey90')) +
#  theme_base(base_size = 18) +
  coord_flip() +
  geom_violin(position = position_dodge(width = .5),scale = "width",color=NA,alpha=.6,size=1) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange",
               position = position_dodge(width = .5),size=1.25, color = "black",show.legend = FALSE) +
  scale_fill_manual(values = c('darkgreen','#D2691E'),guide = guide_legend('')) +
    theme(legend.position = 'bottom') +
  facet_wrap(~term) +
  scale_y_continuous(breaks = pretty_breaks(5)) +
  labs(x = '', y = 'sMAPE')
g
ggsave(plot = g,filename = paste0(figure_path,'/Figure2.png'),width = 10,height = 8)

#### Figure 3 - Taxonomic groups and predictability/transferability ####
# Idea :
# Predictability on x-axis, Transferability on y-axis
# Points for biodiversity measure (shape) and taxonomic group (colour/phylopic)
library(cowplot)

xx <- bind_rows(
  results_predict1 %>% dplyr::select(SS,term,metric, cv_smape, TGrouping) %>%  dplyr::rename(smape = cv_smape) %>% 
    dplyr::mutate(method = 'Predictability') %>% distinct(),
  results_transfer %>% dplyr::select(SS,term,metric, smape.avg, TGrouping) %>% dplyr::rename(smape = smape.avg) %>% 
    dplyr::mutate(method = 'Transferability') %>% distinct()
) %>% 
  # Now convert to wider
  tidyr::pivot_wider(id_cols = c('SS','metric','TGrouping','term'), names_from = 'method', values_from = 'smape'  )

# Summarize
xx <- xx %>% group_by(metric, term, TGrouping) %>% 
  dplyr::summarise(
    Predictability.low = mean(Predictability,na.rm = TRUE) - sd(Predictability,na.rm = TRUE),
    Predictability.high = mean(Predictability,na.rm = TRUE) + sd(Predictability,na.rm = TRUE),
    Transferability.low = mean(Transferability,na.rm = TRUE) - sd(Transferability,na.rm = TRUE),
    Transferability.high = mean(Transferability,na.rm = TRUE) + sd(Transferability,na.rm = TRUE),
    Predictability = mean(Predictability,na.rm = TRUE),
    Transferability = mean(Transferability,na.rm = TRUE)
  ) %>% dplyr::filter(TGrouping != 'Other')

xx$TGrouping <- factor(xx$TGrouping,levels = c('Plantae','Fungi','Invertebrates','Amphibia','Reptilia','Aves','Mammalia'))
xx$metric <- factor(xx$metric,levels = c('SR','LA','PIE','SOR'),
                      labels = c('Species richness', 'Species abundance', 'Assemblage evenness', 'Assemblage composition'))
# Colours for taxonomic grouping
cols <- c("#7fbf7b","violet","#FF7F0E","#ffc425","#5D959A","tomato2","tan4","grey50")

# Plot
g <- ggplot(xx %>% dplyr::filter(term == 'Photosynthetic activity'),
            aes(x = Predictability, y = Transferability, colour = TGrouping, shape = metric)) +
  theme_classic(base_size = 20,base_family = 'Arial') +
    theme(panel.grid.major =  element_line(size = .5,colour='grey90')) +
  geom_abline(slope = 1) +
  geom_point(size = 7) +
  geom_errorbarh(aes(xmin = Predictability.low, xmax = Predictability.high),size = .5) +
  geom_errorbar(aes(ymin = Transferability.low, ymax = Transferability.high),size = .5) +
  scale_colour_manual(values = cols, guide = guide_legend(title = 'Taxonomic group')) +
  scale_shape_manual(values = c("\u25CF","\u25B2","\u25A0","\u25BC"),guide = guide_legend(title = 'Measure')) + 
    theme(legend.position = c(.85,.25),legend.background = element_blank()) +
  scale_x_continuous(breaks = pretty_breaks(7)) + scale_y_continuous(breaks = pretty_breaks(7)) +
  xlim(c(-5,100)) + ylim(c(-5,100)) +
  labs(x = 'Predictability error (%)', y = 'Transferability error (%)')
g
ggsave(plot = g,filename = paste0(figure_path,'/Figure3.png'),width = 10,height = 10)

# boxplots above and to the right?
gt1 <- ggplot(xx %>% dplyr::filter(term == 'Photosynthetic activity'),
            aes(x = Predictability, fill = TGrouping)) + 
  theme_void() + 
  geom_boxplot(colour = 'white',outlier.shape = NA) +
  xlim(c(-5,100)) +
  scale_fill_manual(values = cols, guide = guide_legend(title = 'Taxonomic group')) + guides(fill = 'none')
gt2 <- ggplot(xx %>% dplyr::filter(term == 'Photosynthetic activity'),
             aes(y = Transferability, fill = TGrouping)) + 
  theme_void() +
  geom_boxplot(colour = 'white',outlier.shape = NA) +
  ylim(c(-5,100)) +
  scale_fill_manual(values = cols, guide = guide_legend(title = 'Taxonomic group')) + guides(fill = 'none')

gg <- cowplot::plot_grid(gt1,NULL,g + theme(legend.position = 'bottom'),gt2,
                         align = 'hv',
                   nrow = 2,rel_heights = c(0.1,0,1,0.1),rel_widths = c(0.1,0,1,0.1))
cowplot::ggsave2(plot = gg,filename = paste0(figure_path,'/Figure3comb.png'),width = 10,height=10)

# ------------------------------------ #
#### Figure 4 - Other variables explaining prediction error  ####
# Idea:
# use a regression tree to identify driving study-wide factors
# Plot for both predictability and transferability

mode <- function(codes){ 
  x <- names( which.max(table(codes)) ) 
  if(is.null(x)) return(NA) else x
}

# Other variables necessary for prediction
ss <- sites %>% dplyr::select(SS,startyear, sampling_duration, Max_linear_extent_metres,
                              Rescaled_sampling_effort, Predominant_land_use,Use_intensity,
                              Sampling_grouping, Grouping_unit, Effort_multiplier,
                              TransferGrouping, N_samples,
                              accessibility,tri, cloud) %>% 
  dplyr::group_by(SS) %>% 
  dplyr::summarise(startyear = median(startyear,na.rm = TRUE),
                   Max_linear_extent_metres = median(Max_linear_extent_metres,na.rm = TRUE),
                   Rescaled_sampling_effort = mean(Rescaled_sampling_effort,na.rm = TRUE),
                   sampling_duration = median(sampling_duration,na.rm = TRUE),
                   Predominant_land_use = mode(Predominant_land_use),
                   Use_intensity = mode(Use_intensity),
                   Sampling_grouping = mode(Sampling_grouping), Grouping_unit = mode(Grouping_unit), TransferGrouping = mode(TransferGrouping),
                   Effort_multiplier = mode(Effort_multiplier),
                   accessibility = mean(accessibility,na.rm = TRUE),
                   tri = mean(tri, na.rm = TRUE),
                   cloud = mean(cloud, na.rm = TRUE),
                   NSites = n(),
                   N_samples = mean(N_samples)
  ) %>% ungroup()

df1 <- results_predict1 %>% dplyr::select(SS,term,metric, cv_smape) %>% 
  dplyr::rename(smape = cv_smape) %>% tidyr::drop_na(smape) %>% 
  dplyr::mutate(method = 'Predictability') %>% 
  dplyr::filter(term == 'Photosynthetic activity') %>% distinct() 

df2 <- results_predict2 %>% dplyr::group_by(SS,metric, term) %>% 
  dplyr::summarise(
    smape = mean(smape,na.rm = TRUE)
  ) %>% ungroup() %>% tidyr::drop_na(smape) %>%
  dplyr::filter(term == 'Photosynthetic activity') %>% 
  dplyr::mutate(method = 'Transferability')

df <- bind_rows(df1,df2)
# df1$TGrouping <- factor(df1$TGrouping,levels = c('Plantae','Fungi','Invertebrates','Amphibia','Reptilia','Aves','Mammalia'))

#Join in
df <- df %>% dplyr::left_join(.,  ss,by = 'SS') %>% distinct() %>% 
  # Also join in average remote sensing measures
  dplyr::left_join(., rs)

# Format the effort multiplier
df$Effort_multiplier[is.na(df$Effort_multiplier)] <- 0
df$newEffort[df$Effort_multiplier=='(1/24)'] <- 1/24
df$newEffort[df$Effort_multiplier=='1'] <- 1
df$newEffort[df$Effort_multiplier=='(1/1440)'] <- (1/1440)
df$newEffort[df$Effort_multiplier=='100*(pi*x^2)'] <- 100*(pi*df$Max_linear_extent_metres[which(df$Effort_multiplier=='100*(pi*x^2)')]^2)
df$newEffort[df$Effort_multiplier=='0.0001*(pi*x^2)'] <- 0.0001*(pi*df$Max_linear_extent_metres[df$Effort_multiplier=='0.0001*(pi*x^2)']^2) 
df$newEffort[df$Effort_multiplier=='7'] <- 7
df$newEffort[df$Effort_multiplier=='0.0001'] <-+ 0.0001
df$newEffort[df$Effort_multiplier=='0.0001*(x^2)'] <- 0.0001*(df$Max_linear_extent_metres[df$Effort_multiplier=='0.0001*(x^2)']^2)
df$newEffort[df$Effort_multiplier=='pi*x^2'] <- pi*df$Max_linear_extent_metres[df$Effort_multiplier=='pi*x^2']^2
df$newEffort[df$Effort_multiplier=='100'] <- 100
df$newEffort[is.na(df$newEffort)] <- 1
# New time unit and formatting
df$TimeEffort <- ifelse(df$Grouping_unit=='Time',df$newEffort,1)
df$TransferGrouping <- factor(df$TransferGrouping)
df$Grouping_unit <- factor(df$Grouping_unit)
df$Predominant_land_use <- factor(df$Predominant_land_use,levels = levels(sites$Predominant_land_use))
df$Use_intensity <- factor(df$Use_intensity,levels = levels(sites$Use_intensity))
df$metric <- factor(df$metric,levels = c('SR','LA','PIE','SOR'),
                    labels = c('Species richness', 'Species abundance', 'Assemblage evenness', 'Assemblage composition'))
# ------- #
# library(party)
# library(partykit)
# library(ggparty)
# # Build conditional inference tree
# mod_ct <- party::ctree(smape ~ TimeEffort + Max_linear_extent_metres + Biome + TGrouping + 
#                          sampling_duration + newEffort * Grouping_unit + NSites +
#                          EVI2_mean + SpecHetero_mean,
#                        data = df %>% dplyr::filter(metric ==levels(df$metric)[1]),
#                        controls = party::ctree_control()
#                        )
# # recursive partitioning of a generalized regression model
# mod_ct <- glmtree(smape ~ 0 + metric | TimeEffort + Max_linear_extent_metres + Biome + 
#                     newEffort * Grouping_unit + NSites + #Predominant_land_use +
#                     EVI2_mean + sampling_duration ,#TransferGrouping,
#                      data = df,
#                     family = gaussian())
# 
# plot(mod_ct)

# Do a different analysis and simply assess best predictors
library(lme4)
library(MuMIn)
df2 <- df
df2[,c('TimeEffort','Max_linear_extent_metres','newEffort','NSites','N_samples',
       'EVI2_mean','SpecHetero_mean',
       'sampling_duration','propNA_mean','accessibility','tri','cloud')] <- 
  apply(df2[,c('TimeEffort','Max_linear_extent_metres','newEffort','NSites','N_samples','sampling_duration',
               'EVI2_mean','SpecHetero_mean',
               'propNA_mean','accessibility','tri','cloud')],
        2, scale)

# Build formula
f <- formula(
  smape ~ Max_linear_extent_metres + #Biome +
    # newEffort
#    Max_linear_extent_metres * TimeEffort +
     TimeEffort + NSites +  N_samples + #Predominant_land_use +
                       propNA_mean + accessibility + tri + #cloud +
                        sampling_duration + # TGrouping + 
              (1| TransferGrouping) 
)

options(na.action = "na.fail")
# Method 1 - Predictability
results1 <- data.frame()
for(val in 1:4){
  x = df2 %>% dplyr::filter(metric ==levels(df$metric)[val]) %>% drop_na() %>% 
    dplyr::filter(method == 'Predictability')
  fit <- glmer(f,data = x,
               family = gaussian) # Condition for gamma family is >0
  dd <- MuMIn::dredge(fit, evaluate = TRUE,trace = TRUE,beta = 'sd')
  #or as a 95% confidence set:
  new <- model.avg(dd, subset = cumsum(weight) <= .95,fit = TRUE) # get averaged coefficients
  o <- bind_cols(
    data.frame( estimate = new$coefficients[2,] ) %>% tibble::rownames_to_column('variable'),
    as.data.frame(confint(new))
  ) %>% dplyr::filter(variable != '(Intercept)') %>% tibble::remove_rownames() %>% 
    dplyr::mutate(metric = levels(df$metric)[val],
                  method = 'Predictability',
                  full_r2_cond = performance::r2(fit)$R2_conditional,
                  full_r2_marg = performance::r2(fit)$R2_marginal )
  results1 <- bind_rows(results1, o)
}

# Method 2 - Transferability
results2 <- data.frame()
for(val in 1:4){
  x = df2 %>% dplyr::filter(metric ==levels(df$metric)[val]) %>% drop_na() %>% 
    dplyr::filter(method == 'Transferability')
  fit <- glmer(f,data = x,
               family = gaussian) # Condition for gamma family is >0
  dd <- MuMIn::dredge(fit, evaluate = TRUE,trace = TRUE,beta = 'sd')
  #or as a 95% confidence set:
  new <- model.avg(dd, subset = cumsum(weight) <= .95,fit = TRUE) # get averaged coefficients
  o <- bind_cols(
    data.frame( estimate = new$coefficients[2,] ) %>% tibble::rownames_to_column('variable'),
    as.data.frame(confint(new))
  ) %>% dplyr::filter(variable != '(Intercept)') %>% tibble::remove_rownames() %>% 
    dplyr::mutate(metric = levels(df$metric)[val],
                  method = 'Transferability',
                  full_r2_cond = performance::r2(fit)$R2_conditional,
                  full_r2_marg = performance::r2(fit)$R2_marginal )
  results2 <- bind_rows(results2, o)
}
# sjPlot::plot_model(fit,type = 'est')
# o <- sjPlot::plot_model(fit,type = 'eff')
# Save results for later #
results <- bind_rows(results1,results2)
saveRDS(results, 'resSaves/Figure4_AverageEnsembleModelStatistics.rds')

# --- #
results <- readRDS('resSaves/Figure4_AverageEnsembleModelStatistics.rds')
## Average explained by the models
# Predictability 
results %>% dplyr::filter(method == "Predictability") %>% summarize(mm = mean(full_r2_marg,na.rm=T),
                                                                    mc = mean(full_r2_cond,na.rm=T))
# Transferability
results %>% dplyr::filter(method == "Transferability") %>% summarize(mm = mean(full_r2_marg,na.rm=T),
                                                                    mc = mean(full_r2_cond,na.rm=T))

# Format metric order
results$metric <- factor(results$metric,levels = c('Species richness', 'Species abundance', 'Assemblage evenness', 'Assemblage composition'))
# Format variable names
results$variable <- fct_recode(results$variable,
                       'Sample grain' = 'Max_linear_extent_metres',
                       'Sample duration' = 'sampling_duration',
                       'Effort (sites)' = 'NSites',
                       'Effort (samples)' = 'N_samples',
                       'Effort (time)' = 'TimeEffort',
                       'Missing data (satellite)'= 'propNA_mean',
                       'Topographic ruggedness' = 'tri',
                       'Accessibility' = 'accessibility'
                                 )
results$variable <- factor(results$variable,levels = c(
  'Sample grain','Sample duration', 'Effort (sites)', 'Effort (samples)', 'Effort (time)',
  'Accessibility', 'Topographic ruggedness', 'Missing data (satellite)'
))

# Build the plot
g <- ggplot(results %>% dplyr::filter(method == "Predictability")#dplyr::filter(method == "Transferability") 
            %>% distinct(),
            aes(x = fct_rev(variable), y = estimate, ymin = `2.5 %`, ymax = `97.5 %`,
                    shape = fct_rev(metric) )) +
  theme_classic(base_size = 20,base_family = 'Arial') +
  theme(panel.grid.major =  element_line(size = .5,colour='grey90')) +
  coord_flip() +
  geom_hline(yintercept = 0,linetype = 'dotted', size = 1.25) +
  geom_pointrange(position = position_dodge(width = .5),
                  size = 1.5) +
#  facet_wrap(~method) +
  scale_shape_manual(values = c("\u25CF","\u25B2","\u25A0","\u25BC"),guide = guide_legend(title = 'Measure',nrow = 2)) + 
  scale_y_continuous(breaks = pretty_breaks(7)) +
  theme(legend.position = 'bottom',legend.background = element_blank()) +
  labs(x = '', y = 'Effect on error (Standardized coefficient)')
g
ggsave(plot = g,filename = paste0(figure_path,'/Figure4_transferability.png'),width = 10,height = 10)
# ggsave(plot = g,filename = paste0(figure_path,'/Figure4_predictability.png'),width = 10,height = 10)

# Formula extraction
library(equatiomatic)
equatiomatic::extract_eq(fit)

# ------------------------------------ #
#### Figure 5 - Transferability ####
# Idea: Map on a continuoues surface the 
# predictive horizon in 2d space, highlighting the area for 
# robust predictions can be made
# X - axis (photosynthetic avail)
# y - axis (spectral heterogeneity)
# z

# ------------------------------------ #
#### SI Table 1 ####

r <- readRDS('resSaves/results_transferability.rds') %>% 
  filter(!is.na(smape), dataset == 'test') %>% dplyr::select(SS, grouping) %>% distinct()

n_distinct(r$SS) / n_distinct(sites$SS)

# List of all groupings
sites %>% 
  dplyr::filter(SS %in% r$SS) %>% 
  dplyr::select(TGrouping, Sampling_grouping, Grouping_unit) %>% 
  group_by(TGrouping, Sampling_grouping, Grouping_unit) %>% 
  summarise(n = n()) %>% 
  write.csv(x = ., file = 'SITable.csv',row.names = FALSE)

# ------------------------------------ #
#### SI Figure 1 - Predictor comparison plot ####

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


# ------------------------------------ #
#### SI Figure 3 - Effect size - Error ####
# Plot the coefficient against the error
# Show boxplot for each taxonomic group and biodiversity measure

# Remove greatest outliers
results_predict <- results_predict1
results_predict <- results_predict %>% dplyr::filter(between(estimate,-4,4)) %>% 
  dplyr::mutate(significant = p.value <= 0.05) %>% drop_na(p.value) %>% 
  dplyr::select(term,metric,estimate,cv_smape,SS) %>% distinct()

results_predict$metric <- factor(results_predict$metric,levels = c('SR','LA','PIE','SOR'),
                    labels = c('Species richness', 'Species abundance', 'Assemblage evenness', 'Assemblage composition'))

# Summarize for error bar
ss <- results_predict %>% dplyr::group_by(metric, term) %>% 
  dplyr::summarise(
    estimate.avg = mean(estimate,na.rm=T),
    estimate.low = mean(estimate,na.rm=T) - sd(estimate,na.rm=T),estimate.high = mean(estimate,na.rm=T) + sd(estimate,na.rm=T)
  ) %>% ungroup() 

g <- ggplot(results_predict, aes(x = estimate)) +
  theme_few(base_size = 20) +
  geom_vline(xintercept = 0,size = 0.5, alpha = .5, colour = 'grey20') +
  geom_histogram(bins = 50,alpha = .8) +
  geom_pointrange(data = ss, aes(x = estimate.avg, xmin = estimate.low, xmax = estimate.high, y = -1),
                  fatten = T, colour = 'red', size = 1.5) +
  scale_x_continuous(breaks = pretty_breaks(5)) +
  facet_grid(metric~term,scales = 'free') +
  labs(x = 'Regression coefficient', y = 'Number of studies')
g
ggsave(filename = paste0(figure_path,'/','SIFigure3.png'),plot = g,width = 10,height = 12,dpi = 400)

# ------------------------------------ #
#### SI Figure 3 - Observed vs predicted ####
# Idea: 
# Make a comparison plot of observed vs predicted for each response and term
# Facet wrap with two rows
library(lme4)

# Run models from predictability section of analysis script!
assert_that(exists('fit.full.sr1'),
            exists( 'fit.full.la1'),
            exists( 'fit.full.pie1'),
            exists( 'fit.full.dis1'),
            exists('df.full1'), exists('df'))

getObsPred <- function(model){
  response <- all.vars(model@call$formula)[1]
  o <- model@frame
  n <- lme4:::predict.merMod(model,re.form = NA,newdata = o,type = 'response')
  return(
    data.frame(response,
               observed = o[,response],
               predicted = as.vector(n))
  )
}

x <- getObsPred(fit.full.la1) %>% dplyr::mutate(type='Photosynthetic activity')
plot(x$observed~x$predicted)
x <- getObsPred(fit.full.la2) %>% dplyr::mutate(type='Spectral variability')
plot(x$observed~x$predicted)


