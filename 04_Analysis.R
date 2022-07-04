# Data science 
library(tidyverse)
library(reshape2)
library(stringr)
library(lubridate)
library(assertthat);library(assertr)
# Visualization
library(ggplot2)
library(ggthemes)
library(scales)
library(patchwork)
library(colorspace)
library(rsample) # tidy resampling
# Modelling
library(lme4)
library(mgcv)
library(broom)
source('00_HelperFunction.R')
# Hacks to make BRMS faster
options(mc.cores = parallel::detectCores()-1)
options(dplyr.summarise.inform = FALSE)
# ------------------------------------ #

# Figure path
figure_path <- 'figures'
# Intermediate results
output_path <- 'resSaves'
dir.create(output_path,showWarnings = FALSE)
dir.create(figure_path,showWarnings = FALSE)
# --- #

# Path to PREDICTS data
sites <- readRDS('resSaves/sites_diversity.rds')

# Calculate biodiversity metrics
sites <- sites %>% 
  dplyr::mutate(logabund = log10(Total_abundance+1),
                asinPIE = asin(sqrt(PIE)))

# Pairwise sorensen. Need to further subset to target studies 
sites.dis <- readRDS('resSaves/sites_pairwise_sorensen.rds')
sites.dis[names(sites.dis) %notin% unique(sites$SS)] <- NULL
# Furthermore check and remove those that are all 1 or NA.
x <- which( sapply(sites.dis, function(x){ length( unique( as.vector(x) ) ) <= 2 }) )
sites.dis[names(x)] <- NULL

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

# Load pairwise distance
pdis <- readRDS('resSaves/pairwiseSiteDistance.rds') %>% 
  dplyr::mutate(distance = log10(distance+1))

# Distance to centroid
cdis <- readRDS(paste0('resSaves/centroidDistance.rds')) %>% 
  dplyr::mutate(distance = distance + 0.5)

# ------------------------------------ #
#### Figure 1 - Predictability with and without partial pooling ####
# Idea:
# Take site based measures and assess predictability 
# within (a) each study and (b) within each study but without partial pooling

results <- data.frame()

# Now loop through each study
for(study in unique(as.character(sites$SS)) ){
  myLog(study)
  
  # Get subset
  sites.sub <- subset(sites, SS == study)
  sites.sub <- dplyr::left_join(sites.sub, cdis, by = c('SS','SSBS'))
  if(study %in% names(sites.dis)) {sites.dis.sub <- sites.dis[[study]]} else {sites.dis.sub <- NULL}
  
  # Combine with rs data
  rs.sub <- rs %>% dplyr::filter(SSBS %in% sites.sub$SSBS)
  if(nrow(rs.sub)==0) {next()}
  
  sites.sub <- sites.sub %>% left_join(., rs, by = 'SSBS')
  
  # Formulate simple GLMs for each
  fit.full.sr1_lin <- gam(Species_richness ~ EVI2_mean, family = poisson, data = sites.sub)
  fit.full.sr2_lin <- gam(Species_richness ~ PCA_BRDF_meancentroid, family = poisson, data = sites.sub)
  if(nrow(sites.sub)>10){
    fit.full.sr1_spl <- gam(Species_richness ~ s(EVI2_mean), family = poisson, data = sites.sub)
    fit.full.sr2_spl <- gam(Species_richness ~ s(PCA_BRDF_meancentroid), family = poisson, data = sites.sub)
  } else {
    fit.full.sr1_spl <- NULL; fit.full.sr2_spl <- NULL
  }
  # fit.full.sr1 <- glm(Species_richness ~ EVI2_mean, family = poisson, data = sites.sub)
  # fit.full.sr2 <- glm(Species_richness ~ PCA_BRDF_meancentroid, family = poisson, data = sites.sub)
  
  if(!all(is.na(sites.sub$logabund))){
    # fit.full.la1 <- glm(logabund ~ EVI2_mean, family = gaussian, data = sites.sub)
    # fit.full.la2 <- glm(logabund ~ PCA_BRDF_meancentroid, family = gaussian, data = sites.sub)
    fit.full.la1_lin <- gam(logabund ~ EVI2_mean, family = gaussian, data = sites.sub)
    fit.full.la2_lin <- gam(logabund ~ PCA_BRDF_meancentroid, family = gaussian, data = sites.sub)
    if(nrow(sites.sub)>10){
      fit.full.la2_spl <- gam(logabund ~ s(EVI2_mean), family = gaussian, data = sites.sub)
      fit.full.la2_spl <- gam(logabund ~ s(PCA_BRDF_meancentroid), family = gaussian, data = sites.sub)
    } else {
      fit.full.la2_spl <- fit.full.la2_spl <- NULL
    }
  } else {fit.full.la1_lin <- fit.full.la1_lin <- fit.full.la2_spl <- fit.full.la2_spl <- NULL}
  
  if(!all(is.na(sites.sub$asinPIE))){
    # Ocassionally PIE can't be calculated
    fit.full.pie1_lin <- gam(asinPIE ~ EVI2_mean, family = gaussian, data = sites.sub)
    fit.full.pie2_lin <- gam(asinPIE ~ PCA_BRDF_meancentroid, family = gaussian, data = sites.sub)
    if(nrow(sites.sub)>10){
     fit.full.pie1_spl <- gam(asinPIE ~ s(EVI2_mean), family = gaussian, data = sites.sub)
      fit.full.pie2_spl <- gam(asinPIE ~ s(PCA_BRDF_meancentroid), family = gaussian, data = sites.sub)
    } else {
      fit.full.pie1_spl <- fit.full.pie2_spl <- NULL
    }
    if(is.na(coef(fit.full.pie1_lin)[2])){ fit.full.pie1_lin <- fit.full.pie2_lin <- NULL }
  } else {
    fit.full.pie1_lin <- fit.full.pie2_lin <- fit.full.pie1_spl <- fit.full.pie2_spl <- NULL
  }
  
  if(!is.null(sites.dis.sub)){
    require(reshape2)
    # For dissimilarity, first format the data
    sites.dis.sub[upper.tri(sites.dis.sub)] <- NA # Need only one half
    df.dis <- reshape2::melt(sites.dis.sub) %>% tidyr::drop_na() %>% dplyr::rename(sr = value)
    # Transform y estimates so that those are between 0 and 1
    # https://stats.stackexchange.com/questions/48028/beta-regression-of-proportion-data-including-1-and-0
    #df.full$sr.n <- (df.full$sr * (nrow(df.full) - 1) + 0.5) / nrow(df.full)
    # Logit transform
    df.dis$sr.n <- car::logit(df.dis$sr,adjust = 0.01)
    #df.dis$sr.n <- asin(sqrt(df.dis$sr))
    
    # Calculate for each pairing the respective distance between centroids
    # Make data.frame and loop through each pairing
    # And for mean EVI 2
    rs.dis1 <- as.matrix( dist(rs.sub$EVI2_mean,method = 'manhattan',diag=F,upper=F) )
    diag(rs.dis1) <- NA; rs.dis1[upper.tri(rs.dis1)] <- NA
    rownames(rs.dis1) <- as.character(rs.sub$SSBS); colnames(rs.dis1) <- as.character(rs.sub$SSBS)
    rs.dis1 <- reshape2::melt(rs.dis1) %>% tidyr::drop_na() %>% dplyr::rename(rsdis = value)

    rs.dis2 <- as.matrix( dist(rs.sub$PCA_BRDF_centroid,method = 'manhattan',diag=F,upper=F) )
    diag(rs.dis2) <- NA; rs.dis2[upper.tri(rs.dis2)] <- NA
    rownames(rs.dis2) <- as.character(rs.sub$SSBS); colnames(rs.dis2) <- as.character(rs.sub$SSBS)
    rs.dis2 <- reshape2::melt(rs.dis2) %>% tidyr::drop_na() %>% dplyr::rename(rsdis = value)
        
    if(nrow(rs.dis1)==0 | nrow(rs.dis2) == 0){
      fit.full.dis1_1 <- fit.full.dis2_1 <- fit.full.dis2_2 <- fit.full.dis1_2 <- NULL
    } else {
      # Join
      df.full1 <- full_join(df.dis, rs.dis1, by = c('Var1','Var2')) %>% tidyr::drop_na() #EVI
      df.full2 <- full_join(df.dis, rs.dis2, by = c('Var1','Var2')) %>% tidyr::drop_na() # BRDF
      # Join in distance
      df.full1 <- dplyr::left_join(df.full1, pdis %>% dplyr::filter(SS == study) %>% 
                                     dplyr::mutate(distance = normalize(distance)), by = c('Var1','Var2'))
      df.full2 <- dplyr::left_join(df.full2, pdis %>% dplyr::filter(SS == study) %>% 
                                     dplyr::mutate(distance = normalize(distance)), by = c('Var1','Var2'))
      assert_that(nrow(df.full1)>0,
                  all( between(df.full1$sr,0,1) ))
      # Fit
      fit.full.dis1_lin <- gam(sr.n ~ distance + rsdis,data = df.full1,family = gaussian())
      fit.full.dis2_lin <- gam(sr.n ~ distance + rsdis,data = df.full2,family = gaussian())
      if(nrow(df.full1)>10){
        fit.full.dis1_spl <- gam(sr.n ~ distance + s(rsdis),data = df.full1,family = gaussian())
        fit.full.dis2_spl <- gam(sr.n ~ distance + s(rsdis),data = df.full2,family = gaussian())
      } else {
        fit.full.dis1_spl <- fit.full.dis2_spl <- NULL
      }
    }
  }
  
  # Get all estimates for any given model
  assessPredictability <- function(model){
    if(is.null(model)) {
      return(data.frame())
    }
    # Check whether one or models are NULL, if so ignore
    if(all(sapply(model, is.null))) return(data.frame())
    if(any(sapply(model, is.null))){
      model <- model[which(!sapply(model, is.null))]
    } 
    # Simple cross validation script for gam
    cvfit <- function(mod,n=10){
      o <- data.frame()
      for(i in 1:n){
        # Do k-fold Cv
        x <- mod$model[,all.vars(mod$formula)] %>% tidyr::drop_na()
        
        # Select at random, but weighted by distance
        try( train <- x[sample(1:nrow(x),size = round(nrow(x)*0.66),prob = mod$data$distance[as.numeric(rownames(x))] ),],silent = T )
        if(!exists('train')) { train <- x[sample(1:nrow(x),size = round(nrow(x)*0.66) ),]  }
        test <- x[which(rownames(x) %notin% rownames(train)),]
        
        assert_that(nrow(test)+nrow(train)==nrow(x))
        if(nrow(train)<3) {next()}
        
        # Retrain model(s)
        new <- try({ gam(mod$formula, data = train, family = mod$family)},silent = TRUE)
        if(inherits(new, "try-error")) next()
        
        # Nonlinear?
        is_smooth <- ifelse(length(grep("s\\(", names(coef(new))))>0, TRUE, FALSE)
        
        o <- bind_rows(o,
                       data.frame(i = i,type = 'mape',is_smooth = is_smooth,
                                  value = mape(observed = purrr::discard(test[, all.vars(mod$formula)[1]],is.na),
                                  predicted = purrr::discard(predict(new,newdata = test,type = 'response'),is.na) )
                                  ),
                       data.frame(i = i,type = 'smape',is_smooth = is_smooth,
                                  value = mape(observed = purrr::discard(test[, all.vars(mod$formula)[1]],is.na),
                                  predicted = purrr::discard(predict(new,newdata = test,type = 'response'),is.na),type = 'sym') )
        )
        rm(train,test)
      }
     
     if(nrow(o)>0){
       # Average
       o %>% group_by(type,is_smooth) %>% summarise(avg = mean(value), sd = sd(value)) %>% ungroup()
     } else { data.frame()}
    }
    cv <- data.frame()
    for(i in 1:length(model)){
      cv <- rbind(cv, cvfit(model[[i]]) )
    }
    # Get best model out of the two supplied ones
    check <- cv$is_smooth[which(cv$type=="smape")][which.min(cv$avg[which(cv$type=="smape")])]
    if(check) bestmodel <- model[[2]] else bestmodel <- model[[1]]
    cv <- subset(cv, is_smooth == check)
    
    # Get model coefficients without intercept
    if(check){
      cc <- broom::tidy(bestmodel)[-1,]
    } else {
      cc <- broom::tidy(bestmodel,parametric = TRUE)[-1,]
    }
    cc$is_smooth <- check
    cc$full_r2 <- performance::r2_nagelkerke(bestmodel)
    td <- data.frame(observed =  bestmodel$model[, all.vars(bestmodel$formula)[1]],
                     predicted =  predict(bestmodel,newdata = bestmodel$model,type = 'response')) %>% tidyr::drop_na()
    
    cc$full_mape <- mape(td$observed,td$predicted)
    cc$full_smape <- mape(td$observed,td$predicted,type = 'symetric')
    # --- #
    # Add cross-validated data to results output
    cc$cv_mape <- cv$avg[cv$type=='mape']
    cc$cv_smape <- cv$avg[cv$type=='smape']
    try({rm(cv,td, bestmodel)},silent = T)
    return(cc)
  }
  
  results <- dplyr::bind_rows(
    results,
    bind_rows(
      assessPredictability(list(fit.full.sr1_lin,fit.full.sr1_spl)) %>% dplyr::mutate(SS = study,metric = 'SR'),
      assessPredictability(list(fit.full.sr2_lin,fit.full.sr2_spl)) %>% dplyr::mutate(SS = study,metric = 'SR'),
      assessPredictability(list(fit.full.la1_lin,fit.full.la2_spl)) %>% dplyr::mutate(SS = study,metric = 'LA'),
      assessPredictability(list(fit.full.la2_lin,fit.full.la2_spl)) %>% dplyr::mutate(SS = study,metric = 'LA'),
      assessPredictability(list(fit.full.pie1_lin,fit.full.pie1_spl)) %>% dplyr::mutate(SS = study,metric = 'PIE'),
      assessPredictability(list(fit.full.pie2_lin,fit.full.pie2_spl)) %>% dplyr::mutate(SS = study,metric = 'PIE'),
      assessPredictability(list(fit.full.dis1_lin,fit.full.dis1_spl)) %>% dplyr::mutate(SS = study,metric = 'SOR') %>% mutate(term = 'EVIdis') ,
      assessPredictability(list(fit.full.dis2_lin,fit.full.dis2_spl)) %>% dplyr::mutate(SS = study,metric = 'SOR') %>% mutate(term = 'PCAcent')
    )
  )
  # Clean up
  rm(rs.sub,sites.dis.sub,sites.sub)
}

saveRDS(results,paste0(output_path, '/results_predicatability_nopooling.rds'))

# --- #
# Now build hierachical models but with pooling
# Build 10 sets of cross-validated studies
# Predict the overall hold out data
# Calculate sMAPE as before

results <- data.frame()

# Combine with rs data
df <- sites %>% left_join(., rs, by = 'SSBS') %>% 
  dplyr::left_join(.,cdis, by = c('SS','SSBS'))
rownames(df) <- df$SSBS

#### Predictability Model section here ####
# Formulate GLMERs for each
fit.full.sr1 <- glmer(Species_richness ~ EVI2_mean + (1|SS), family = poisson, data = df)
fit.full.sr2 <- glmer(Species_richness ~ PCA_BRDF_meancentroid + (1|SS), family = poisson, data = df)

fit.full.la1 <- glmer(logabund ~ EVI2_mean + (1|SS), family = gaussian, data = df)
fit.full.la2 <- glmer(logabund ~ PCA_BRDF_meancentroid + (1|SS), family = gaussian, data = df)

fit.full.pie1 <- glmer(asinPIE ~ EVI2_mean + (1|SS), family = gaussian, data = df)
fit.full.pie2 <- glmer(asinPIE ~ PCA_BRDF_meancentroid + (1|SS), family = gaussian, data = df)

# For Sor turnover
df.sr <- bind_rows(
    lapply(names(sites.dis), function(x) {
    d <- sites.dis[[x]]
    d[upper.tri(d)] <- NA # Need only one half
    reshape2::melt(d) %>% tidyr::drop_na() %>% 
      dplyr::rename(sr = value) %>% 
      dplyr::mutate(SS = x)
  })
)
df.sr$sr.n <- car::logit(df.sr$sr,adjust = 0.01) # Logit transform with small constant
#df.sr$sr.n <- asin(sqrt(df.sr$sr))

# Same for the remote sensing dissimilarity metrics
df.rs <- data.frame()
for(study in unique(sites$SS)){
  if(study %in% names(sites.dis)) {sites.dis.sub <- sites.dis[[study]]} else {next()}
  
  rs.sub <- rs %>% dplyr::filter(SSBS %in% union(rownames(sites.dis.sub),colnames(sites.dis.sub)) )
  if(nrow(rs.sub)==0) next()
  # Calculate for each pairing the respective distance between centroids
  # Make data.frame and loop through each pairing
  rs.dis1 <- as.matrix( dist(rs.sub$PCA_BRDF_centroid,method = 'manhattan',diag=F,upper=F) )
  diag(rs.dis1) <- NA; rs.dis1[upper.tri(rs.dis1)] <- NA
  rownames(rs.dis1) <- as.character(rs.sub$SSBS); colnames(rs.dis1) <- as.character(rs.sub$SSBS)
  rs.dis1 <- reshape2::melt(rs.dis1) %>% tidyr::drop_na() %>% dplyr::rename(rsdis = value) %>% 
    mutate(term = 'PCAcent')
  # And for mean EVI 2
  rs.dis2 <- as.matrix( dist(rs.sub$EVI2_mean,method = 'manhattan',diag=F,upper=F) )
  diag(rs.dis2) <- NA; rs.dis2[upper.tri(rs.dis2)] <- NA
  rownames(rs.dis2) <- as.character(rs.sub$SSBS); colnames(rs.dis2) <- as.character(rs.sub$SSBS)
  rs.dis2 <- reshape2::melt(rs.dis2) %>% tidyr::drop_na() %>% dplyr::rename(rsdis = value) %>% 
    mutate(term = 'EVIdis')
  
  df.rs <- bind_rows(df.rs, rs.dis1,rs.dis2)
}
# Join
df.full1 <- full_join(df.sr, df.rs %>% dplyr::filter(term == 'EVIdis'), by = c('Var1','Var2')) %>% tidyr::drop_na()
df.full2 <- full_join(df.sr, df.rs %>% dplyr::filter(term == 'PCAcent'), by = c('Var1','Var2')) %>% tidyr::drop_na()
# Join in distance
df.full1 <- dplyr::left_join(df.full1, pdis, by = c('SS','Var1','Var2'))
df.full2 <- dplyr::left_join(df.full2, pdis, by = c('SS','Var1','Var2'))

assert_that(nrow(df.full1)>0,
            all( between(df.full1$sr,0,1) ))
# Fit
fit.full.dis1 <- glmer(sr.n ~ rsdis + distance + (1|SS),data = df.full1,family = gaussian())
fit.full.dis2 <- glmer(sr.n ~ rsdis + distance + (1|SS),data = df.full2,family = gaussian())

# Function call to assess within study predictability with partial pooling
assessPredictabilityGLMER <- function(model){
  if(is.null(model)) {
    return(data.frame())
  }
  # Simple cross validation script for glm
  cvfit <- function(model,n=10){
    o <- data.frame()
    for(i in 1:n){
      # Do k-fold Cv
      x <- model@frame[,all.vars(model@call$formula)] %>% tidyr::drop_na() %>% 
        tibble::rownames_to_column()
      # Join in distance if not already in there
      if(!has_name(x,'distance')){x <- left_join(x, df %>% dplyr::select(SSBS,distance), by = c('rowname' = 'SSBS'))  %>% tidyr::drop_na()}
      # Select at random within group
      train <- x %>% group_by(SS) %>% sample_frac(.66,weight = (distance+0.5)) %>% ungroup()
      test <- x %>% dplyr::filter(rowname %notin% train$rowname)
      
      # Retrain model for training dataset
      if('family' %in% names(model@call)){
        fam <- unlist( ifelse(model@call$family == 'poisson',poisson(),gaussian()) )
      } else { fam <- gaussian()}
      new <- glmer(model@call$formula, data = train, family = fam)
      
      pb <- progress::progress_bar$new(total = length(unique(test$SS)))
      for(study in unique(test$SS)){
        sub <- test %>% dplyr::filter(SS == study)
        o <- bind_rows(
          o,
          data.frame(SS = study, i = i,type = 'mape',value = mape(observed = sub[,all.vars(model@call$formula)[1]],
                                                                  predicted = lme4:::predict.merMod(new,newdata = sub,type = 'response') )
          ),
          data.frame(SS = study, i = i,type = 'smape',value = mape(observed = sub[,all.vars(model@call$formula)[1]],
                                                                  predicted = lme4:::predict.merMod(new,newdata = sub,type = 'response'),type = 'sym')
          )
        )
        pb$tick()
      } # study-wide MAPE prediction
      pb$terminate()
    } # N for loop
    
    if(nrow(o)>0){
      # Average
      o %>% group_by(SS,type) %>% summarise(avg = mean(value), sd = sd(value)) %>% ungroup()
    } else { data.frame(type = NA, avg = NA, sd = NA)}
  }
  cv <- cvfit(model) %>% 
    dplyr::select(SS:avg) %>% tidyr::pivot_wider(names_from = 'type',values_from = 'avg')
  
  cc <- coef(model)$SS %>% tibble::rownames_to_column(var = 'SS')
  if(length(names(cc)) == 3){
    names(cc) <- c('SS','intercept','value')
  } else { names(cc) <- c('SS','intercept','value','distance'); cc <- cc %>% dplyr::select(-distance)}
  cc$term <- all.vars(model@call$formula)[2]
  cc <- dplyr::inner_join(cv,cc, by = 'SS')
  # cc <- data.frame(t(summary(model)$coefficients[2,])) %>%  # Get model coefficients without intercept
  #   dplyr::mutate(term = all.vars(model@call$formula)[2] )
  # cc$full_r2c <- performance::r2_nakagawa(model)$R2_conditional
  # cc$full_r2m <- performance::r2_nakagawa(model)$R2_marginal
  # td <- data.frame(observed =  model@frame[, all.vars(model@call$formula)[1]],
  #                  predicted =  predict(model,newdata = model@frame,type = 'response')) %>% tidyr::drop_na()
  # 
  # cc$full_mape <- mape(td$observed,td$predicted )
  # cc$full_smape <- mape(td$observed,td$predicted,type = 'symetric')
  # --- #
  return(cc)
}

results <- dplyr::bind_rows(
    assessPredictabilityGLMER(fit.full.sr1) %>% dplyr::mutate(metric = 'SR') %>% mutate(term = 'EVIdis'),
    assessPredictabilityGLMER(fit.full.sr2) %>% dplyr::mutate(metric = 'SR') %>% mutate(term = 'PCAcent'),
    assessPredictabilityGLMER(fit.full.la1) %>% dplyr::mutate(metric = 'LA') %>% mutate(term = 'EVIdis'),
    assessPredictabilityGLMER(fit.full.la2) %>% dplyr::mutate(metric = 'LA') %>% mutate(term = 'PCAcent'),
    assessPredictabilityGLMER(fit.full.pie1) %>% dplyr::mutate(metric = 'PIE') %>% mutate(term = 'EVIdis'),
    assessPredictabilityGLMER(fit.full.pie2) %>% dplyr::mutate(metric = 'PIE') %>% mutate(term = 'PCAcent'),
    assessPredictabilityGLMER(fit.full.dis1) %>% dplyr::mutate(metric = 'SOR') %>% mutate(term = 'EVIdis'),
    assessPredictabilityGLMER(fit.full.dis2) %>% dplyr::mutate(metric = 'SOR') %>% mutate(term = 'PCAcent')
)
# Save
saveRDS(results,paste0(output_path, '/results_predicatability_pooling.rds'))

# ------------------------- #
#### Transferability computation within groups ####
# Plan:
# Remove single studies within the same grouping
# Fit and predict the model on the missing estimates

# Remove any that have NA in method or unit!

sites.transf <- sites %>% tidyr::drop_na(Sampling_grouping, Grouping_unit)
n_distinct(sites.transf$TransferGrouping)

results <- data.frame()
# Now pooling
for(gg in unique(sites.transf$TransferGrouping)){
  
  # Get subset
  sites.sub <- subset(sites.transf,TransferGrouping == gg)
  sites.sub <- dplyr::left_join(sites.sub, cdis, by = c('SS','SSBS'))
  # Get studies
  sites.dis.sub <- sites.dis[ which(names(sites.dis) %in% unique(sites.sub$SS)) ]
  myLog('Processing ', n_distinct(sites.sub$SS), ' studies in this grouping')
  
  # Combine with rs data
  rs.sub <- rs %>% dplyr::filter(SSBS %in% sites.sub$SSBS) %>% dplyr::left_join(., sites.transf %>% dplyr::select(SS,SSBS), by = 'SSBS')
  if(nrow(rs.sub)==0 | n_distinct(sites.sub$SS)<=1) {next()}
  sites.sub <- sites.sub %>% left_join(., rs.sub, by = c('SS','SSBS'))
  
  assertthat::assert_that(
    all(rs.sub$SS %in% sites.sub$SS),
    all(     names(sites.dis.sub) %in% unique(sites.sub$SS) )
  )
  # Build subsets
  out <- data.frame()
  ss <- rep(1,length(unique(sites.sub$SS)));names(ss) <- unique(sites.sub$SS)
  for(k in 1:10){
    # --- #
    ss.train <- sample(names(ss),size = length(ss)*.66,prob = ss)
    ss.test <- names(ss)[names(ss) %notin% ss.train]
    if(length(ss.test)==0 || length(ss.train) == 0 ) next()
    # Change weighting for next iteration
    ss[ss.train] <- ss[ss.train] -0.1;ss[ss.test] <- ss[ss.test] + 0.1
    # Now fit for each response
    df.train <- subset(sites.sub,SS %in% ss.train)
    df.test <- subset(sites.sub, SS %in% ss.test)
    
    if(all(is.na(df.train$EVI2_mean)) || length(which(!is.na(df.train$EVI2_mean))) <= 2) next()
    
    fit.full.sr1 <- glm(Species_richness ~ EVI2_mean, family = poisson, data = df.train)
    fit.full.sr2 <- glm(Species_richness ~ PCA_BRDF_meancentroid, family = poisson, data = df.train)
    
    sub <- df.train %>% tidyr::drop_na(logabund,EVI2_mean)
    if(!all(is.na(df.train$logabund)) & nrow(sub) >2){
      fit.full.la1 <- glm(logabund ~ EVI2_mean, family = gaussian, data = sub)
      fit.full.la2 <- glm(logabund ~ PCA_BRDF_meancentroid, family = gaussian, data = sub)
      rm(sub)
    } else {fit.full.la1 <- NULL; fit.full.la2 <- NULL}
    
    sub <- df.train %>% tidyr::drop_na(asinPIE,EVI2_mean)
    if(!all(is.na(df.train$asinPIE)) & nrow(sub)>2){
      # Occasionally PIE can't be calculated
      fit.full.pie1 <- glm(asinPIE ~ EVI2_mean, family = gaussian, data = df.train)
      fit.full.pie2 <- glm(asinPIE ~ PCA_BRDF_meancentroid, family = gaussian, data = df.train)
      if(is.na(coef(fit.full.pie1)[2])){ fit.full.pie1 <- NULL; fit.full.pie2 <- NULL }
    } else {
      fit.full.pie1 <- NULL; fit.full.pie2 <- NULL
    }
    
    if(!is.null(sites.dis.sub)){
      require(reshape2)
      # For dissimilarity, first format the data
      cp <- sites.dis.sub
      # TODO: Loop through and insert NA
      for(n in names(cp)){
        cp[[n]][upper.tri(cp[[n]])] <- NA
      }
      df.dis <- reshape2::melt(cp) %>% tidyr::drop_na() %>% dplyr::rename(sr = value,SS = L1)
      # Transform y estimates so that those are between 0 and 1
      # https://stats.stackexchange.com/questions/48028/beta-regression-of-proportion-data-including-1-and-0
      #df.full$sr.n <- (df.full$sr * (nrow(df.full) - 1) + 0.5) / nrow(df.full)
      # Logit transform
      df.dis$sr.n <- car::logit(df.dis$sr,adjust = 0.01,percents = FALSE)
      #df.dis$sr.n <- asin(sqrt(df.dis$sr))
      
      # Make data.frame and loop through each pairing
      # And for mean EVI 2
      rs.dis1 <- data.frame()
      for(n in unique(rs.sub$SS)){
        x <- as.matrix( dist(rs.sub$EVI2_mean[which(rs.sub$SS==n)],method = 'manhattan',diag=F,upper=F) )
        diag(x) <- NA; x[upper.tri(x)] <- NA
        rownames(x) <- as.character(rs.sub$SSBS[rs.sub$SS == n]); colnames(x) <- as.character(as.character(rs.sub$SSBS[rs.sub$SS == n]))
        rs.dis1 <- bind_rows(rs.dis1,
                             reshape2::melt(x) %>% tidyr::drop_na() %>% dplyr::rename(rsdis = value) %>% 
                               dplyr::mutate(SS = n)
        )
      }
      # Calculate for each pairing the respective distance between centroids
      rs.dis2 <- data.frame()
      for(n in unique(rs.sub$SS)){
        x <- as.matrix( dist(rs.sub$PCA_BRDF_centroid[which(rs.sub$SS==n)],method = 'manhattan',diag=F,upper=F) )
        diag(x) <- NA; x[upper.tri(x)] <- NA
        rownames(x) <- as.character(rs.sub$SSBS[rs.sub$SS == n]); colnames(x) <- as.character(as.character(rs.sub$SSBS[rs.sub$SS == n]))
        rs.dis2 <- bind_rows(rs.dis2,
                             reshape2::melt(x) %>% tidyr::drop_na() %>% dplyr::rename(rsdis = value) %>% 
                               dplyr::mutate(SS = n)
        )
      }
      
      if(nrow(rs.dis1)==0 | nrow(rs.dis2) == 0){
        fit.full.dis1 <- NULL; fit.full.dis2 <- NULL
      } else {
        # Join
        df.full1 <- full_join(df.dis, rs.dis1, by = c('Var1','Var2','SS')) %>% tidyr::drop_na()
        df.full2 <- full_join(df.dis, rs.dis2, by = c('Var1','Var2','SS')) %>% tidyr::drop_na()
        # Join in distance
        df.full1 <- dplyr::left_join(df.full1, pdis %>% 
                                       dplyr::mutate(distance = normalize(distance)), by = c('Var1','Var2','SS'))
        df.full2 <- dplyr::left_join(df.full2, pdis %>% 
                                       dplyr::mutate(distance = normalize(distance)), by = c('Var1','Var2','SS'))
        assert_that(nrow(df.full1)>0,
                    all( between(df.full1$sr,0,1) ),
                    has_name(df.full1,'SS'),has_name(df.full2,'SS'))
        # Filter to training SS
        df.train1 <- df.full1 %>% dplyr::filter(SS %in% ss.train)
        df.train2 <- df.full2 %>% dplyr::filter(SS %in% ss.train)
        if(nrow(df.train1)<=2 | nrow(df.train2) <= 2){
          fit.full.dis1 <- NULL; fit.full.dis2 <- NULL
        } else {
          # Fit
          fit.full.dis1 <- glm(sr.n ~ distance + rsdis,data = df.train1,family = gaussian())
          fit.full.dis2 <- glm(sr.n ~ distance + rsdis,data = df.train2,family = gaussian())
        }
      }
    }
    
    # Coefficient and MAPE extacting function
    cvfit <- function(df, model, testingSS){
      if(is.null(model) | nrow(df) <= 1) return( data.frame() )
      assertthat::assert_that(has_name(df, 'SS'),
                              is.vector(testingSS))
      #  Predicts and extract stats
      df$new <- predict(model,newdata = df,type = 'response')
      # Rename observed
      names(df)[which(names(df) == all.vars(model$formula)[1])] <- 'observed'
      # Now extract and predict
      out <- df %>% tidyr::drop_na(new) %>% # Remove sites with NA prediction
        dplyr::group_by(SS) %>% 
        dplyr::summarise(
          mape = mape(observed = observed,predicted = new),
          smape = mape(observed = observed,predicted = new,type = 'smape')
        ) %>% ungroup()
      out$dataset <- ifelse(out$SS %in% testingSS,'test','train')
      return(out %>% dplyr::filter(!is.nan(mape)) )
    }
    # Now extract for each
    out <- bind_rows(
      out,
      cvfit(sites.sub, fit.full.sr1, ss.test) %>% dplyr::mutate(metric = 'SR', term = 'EVIdis',iter = k),
      cvfit(sites.sub, fit.full.sr2, ss.test) %>% dplyr::mutate(metric = 'SR', term = 'PCAcent',iter = k),
      cvfit(sites.sub, fit.full.la1, ss.test) %>% dplyr::mutate(metric = 'LA', term = 'EVIdis',iter = k),
      cvfit(sites.sub, fit.full.la2, ss.test) %>% dplyr::mutate(metric = 'LA', term = 'PCAcent',iter = k),
      cvfit(sites.sub, fit.full.pie1, ss.test) %>% dplyr::mutate(metric = 'PIE', term = 'EVIdis',iter = k),
      cvfit(sites.sub, fit.full.pie2, ss.test) %>% dplyr::mutate(metric = 'PIE', term = 'PCAcent',iter = k),
      cvfit(df.full1, fit.full.dis1, ss.test) %>% dplyr::mutate(metric = 'SOR', term = 'EVIdis',iter = k),
      cvfit(df.full2, fit.full.dis2, ss.test) %>% dplyr::mutate(metric = 'SOR', term = 'PCAcent',iter = k)
    )
  } # End of k iteration
  
  results <- bind_rows(results, out %>% dplyr::mutate(grouping = gg))
  rm(ss, out, sites.sub,sites.dis.sub, rs.sub)
}

# Save results
saveRDS(results,paste0(output_path, '/results_transferability.rds'))

#### Figure for reviewers ####
# Assess correlations between the several 

sub <- subset(rs, between(EVI2_mean, 0,1) )
plot(EVI2_min~EVI2_mean, data = sub)

cor.test(sub$EVI2_mean, sub$NDWI_mean)

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


