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
    dplyr::summarise(EVI2_mean = mean(EVI_mean,na.rm = TRUE),
                     SpecHetero_mean = mean(PCA_BRDF_meancentroid,na.rm = TRUE)) %>% ungroup()


# ------------------------- #
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
library(party)
library(partykit)
library(ggparty)

mode <- function(codes){ 
  x <- names( which.max(table(codes)) ) 
  if(is.null(x)) return(NA)
  x
  }

df <- results_predict1 %>% dplyr::select(SS,term,metric, cv_smape, TGrouping, Biome) %>% 
  dplyr::rename(smape = cv_smape) %>% tidyr::drop_na(smape) %>% 
  dplyr::mutate(method = 'Predictability') %>% 
  dplyr::filter(term == 'Photosynthetic activity') %>% distinct() 
df$TGrouping <- factor(df$TGrouping,levels = c('Plantae','Fungi','Invertebrates','Amphibia','Reptilia','Aves','Mammalia'))
# Join in other variables necessary for prediction
ss <- sites %>% dplyr::select(SS,startyear, sampling_duration, Max_linear_extent_metres,
                        Rescaled_sampling_effort, Predominant_land_use,Use_intensity,
                        Sampling_grouping, Grouping_unit, Effort_multiplier,
                        TransferGrouping) %>% 
  dplyr::group_by(SS) %>% 
  dplyr::summarise(startyear = median(startyear,na.rm = TRUE),
                   Max_linear_extent_metres = median(Max_linear_extent_metres,na.rm = TRUE),
                   Rescaled_sampling_effort = mean(Rescaled_sampling_effort,na.rm = TRUE),
                   sampling_duration = median(sampling_duration,na.rm = TRUE),
                   Predominant_land_use = mode(Predominant_land_use),
                   Use_intensity = mode(Use_intensity),
                   Sampling_grouping = mode(Sampling_grouping), Grouping_unit = mode(Grouping_unit), TransferGrouping = mode(TransferGrouping),
                   Effort_multiplier = mode(Effort_multiplier),
                   NSites = n()
                   ) %>% ungroup()
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
# Build conditional inference tree
mod_ct <- party::ctree(smape ~ TimeEffort + Max_linear_extent_metres + Biome + TGrouping + 
                         sampling_duration + newEffort * Grouping_unit + NSites +
                         EVI2_mean + SpecHetero_mean,
                       data = df %>% dplyr::filter(metric ==levels(df$metric)[3]),
                       controls = party::ctree_control()
                       )
plot(mod_ct)

# recursive partitioning of a generalized regression model
mod_ct <- glmtree(smape ~ 0 + metric | TimeEffort + Max_linear_extent_metres + Biome + TGrouping + 
                    newEffort * Grouping_unit + NSites + #Predominant_land_use +
                    EVI2_mean + sampling_duration ,#TransferGrouping,
                     data = df, family = gaussian())

plot(mod_ct)
plot(mod_ct, tp_args = list(cdplot = TRUE))
plot(mod_ct, terminal_panel = NULL)

# Need to format the tree
# https://jtr13.github.io/cc19/introduction-to-package-ggparty.html
ggparty(mod_ct) +
  geom_edge() +
  geom_edge_label() +
  geom_node_label(aes(label = splitvar),
                  ids = "inner") +
  geom_node_label(aes(label = info),
                  ids = "terminal")

# Plot with ggparty
ggparty(mod_ct,
       terminal_space = 0.5,
       add_vars = list(p.value = "$node$info$p.value")) +
  geom_edge(size = 1.5) +
  geom_edge_label(colour = "grey", size = 6) +
  geom_node_plot(gglist = list(geom_point(aes(x = beauty,
                                              y = eval,
                                              col = tenure,
                                              shape = minority),
                                          alpha = 0.8),
                               theme_bw(base_size = 15)),
                 scales = "fixed",
                 id = "terminal",
                 shared_axis_labels = T,
                 shared_legend = T,
                 legend_separator = T,
                 predict = "beauty",
                 predict_gpar = list(col = "blue",
                                     size = 1.2)
  ) +
  geom_node_label(aes(col = splitvar),
                  line_list = list(aes(label = paste("Node", id)),
                                   aes(label = splitvar),
                                   aes(label = paste("p =", formatC(p.value, format = "e", digits = 2)))),
                  line_gpar = list(list(size = 12, col = "black", fontface = "bold"),
                                   list(size = 20),
                                   list(size = 12)),
                  ids = "inner") +
  geom_node_label(aes(label = paste0("Node ", id, ", N = ", nodesize)),
                  fontface = "bold",
                  ids = "terminal",
                  size = 5, 
                  nudge_y = 0.01) +
  theme(legend.position = "none")


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
#### SI Figure 2 - Effect size - Error ####
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
ggsave(filename = paste0(figure_path,'/','SIFigure2.png'),plot = g,width = 10,height = 12,dpi = 400)

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


