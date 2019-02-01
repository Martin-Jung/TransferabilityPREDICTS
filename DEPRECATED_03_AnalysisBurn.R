library(dplyr)
library(reshape2)
library(stringr)
library(lubridate)
#
library(marfunky)
library(roquefort)
library(gamm4)
library(ggplot2)
library(ggthemr)
library(scales)
myLog <- function(...) {
  cat(paste0("[Analysis] ", Sys.time(), " | ", ..., "\n"))
}

sites <- readRDS("sites_center.rds")
results <- readRDS("MODIS_BurnEstimates_center.rds")
sites <- left_join(sites,results,by="SSBS")

## Subset to within study contrast
# Get only those studies which have at least two sites with differing fire history
sites$Check <- FALSE
for(studyid in unique(sites$SS)){
  myLog(studyid)
  sub <- subset(sites,SS == studyid)
  # Check if there is at least one fire
  if(sum(sub$Fire_nr,na.rm = T)>0) sites$Check[which(sites$SS==studyid)] <- TRUE 
}
sites.sub <- subset(sites,Check == T)
sites.sub$SameCellokay <- NULL;sites.sub$Check <- NULL
sites.sub$Predominant_habitat <- factor(sites.sub$Predominant_habitat,levels = ordered(c("Primary forest","Primary non-forest","Mature secondary vegetation","Intermediate secondary vegetation",
                                                                                 "Young secondary vegetation","Secondary vegetation (indeterminate age)","Plantation forest",
                                                                                 "Pasture","Cropland","Urban","Cannot decide")))
sites.sub$Predominant_habitat<-relevel(sites.sub$Predominant_habitat,ref="Primary forest")
sites.sub$logabund <- log10(sites.sub$Total_abundance+1)
sites.sub$logsimp <- log10(sites.sub$Simpson_diversity)
#
sites.sub$TGrouping <- as.character(sites.sub$Study_common_taxon)
sites.sub$TGrouping[grep("Hymenoptera",x = sites.sub$TGrouping,ignore.case = T)] <- "Invertebrates"
sites.sub$TGrouping[grep("Insecta",x = sites.sub$TGrouping,ignore.case = T)]<- "Invertebrates"
sites.sub$TGrouping[grep("Chordata",x = sites.sub$TGrouping,ignore.case = T)] <- "Other"
sites.sub$TGrouping[grep("Animalia",x = sites.sub$TGrouping,ignore.case = T)] <- "Other"
sites.sub$TGrouping[grep("Formicidae",x = sites.sub$TGrouping,ignore.case = T)] <- "Invertebrates"
sites.sub$TGrouping[grep("Scarabaeidae",x = sites.sub$TGrouping,ignore.case = T)] <- "Invertebrates"
sites.sub$TGrouping[grep("Tracheophyta",x = sites.sub$TGrouping,ignore.case = T)] <- "Plantae"
sites.sub$TGrouping[grep("Strigiformes",x = sites.sub$TGrouping,ignore.case = T)] <- "Aves"
sites.sub$TGrouping[grep("Isoptera",x = sites.sub$TGrouping,ignore.case = T)] <- "Invertebrates"
sites.sub$TGrouping[grep("Coleoptera",x = sites.sub$TGrouping,ignore.case = T)] <- "Invertebrates"
sites.sub$TGrouping[grep("Anogeissus",x = sites.sub$TGrouping,ignore.case = T)] <- "Plantae"
sites.sub$TGrouping[grep("Poaceae",x = sites.sub$TGrouping,ignore.case = T)] <- "Plantae"
sites.sub$TGrouping[grep("Colubridae",x = sites.sub$TGrouping,ignore.case = T)] <- "Reptiles"
sites.sub$TGrouping[grep("Chiroptera",x = sites.sub$TGrouping,ignore.case = T)] <- "Mammalia"
sites.sub$TGrouping[grep("Ascomycota",x = sites.sub$TGrouping,ignore.case = T)] <- "Fungi"
sites.sub$TGrouping[grep("Lepidoptera",x = sites.sub$TGrouping,ignore.case = T)] <- "Invertebrates"
sites.sub$TGrouping[grep("Bryophyta",x = sites.sub$TGrouping,ignore.case = T)] <- "Plantae"
sites.sub$TGrouping[grep("Sarcoptiformes",x = sites.sub$TGrouping,ignore.case = T)] <- "Invertebrates"
sites.sub$TGrouping[which(str_length(sites.sub$TGrouping)==0)] <- "Other"


#----------------------------------------#
#####
# Some Figures
qplot(TimeSinceLastFire,logabund,color=TGrouping,data=sites.sub,facets=~TGrouping)

# How many with fires
a = sites.sub %>% filter(Fire_nr >0) %>% group_by(Predominant_habitat) %>% 
  summarise(N = n(),
            FN = mean(Fire_nr,na.rm=T),
            FN_se = sem(Fire_nr))
levels(a$Predominant_habitat) <- c("PV","PNV","MSV","ISV","YSV","SVI","PL","PA","CL","UR","CA")
g1 <- ggplot(a,aes(x=Predominant_habitat,y=FN,ymin=FN-FN_se,ymax=FN+FN_se)) +
  theme_minimal() +
  geom_linerange(size=1.5) +
  geom_label(aes(x=Predominant_habitat,y=FN,label=N))+
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(x="",y = "Number of Fires \n (2000 - Sampling start)")+
  theme(axis.text.x = element_blank())
g1 # Average number of fires per category

# Average timespan since last fire per category
a = sites.sub %>% filter(TimeSinceLastFire > 0) %>% group_by(Predominant_habitat) %>% 
  summarise(N = n(),
            FN = mean(TimeSinceLastFire,na.rm=T),
            FN_se = sem(TimeSinceLastFire))
levels(a$Predominant_habitat) <- c("PV","PNV","MSV","ISV","YSV","SVI","PL","PA","CL","UR","CA")
g2 <- ggplot(a,aes(x=Predominant_habitat,y=FN,ymin=FN-FN_se,ymax=FN+FN_se)) +
  theme_minimal() +
  geom_linerange(size=1.5) +
  geom_label(aes(x=Predominant_habitat,y=FN,label=N))+
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(x="",y = "Average days since last fire \n (before sampling start)")+
  theme(axis.text.x = element_text(size = 12,angle=0,vjust = 0.5))
g2 # Average number of days since last fire

library(cowplot)
g <- plot_grid(g1,g2,nrow=2,label_size = 16)
save_plot("Burnoverview.png",plot=g,base_width = 14,base_height = 8)

#### Test: Do occurence of past fires make any difference ####
sites.sub <- filter(sites.sub,Predominant_habitat != "Cannot decide")
sites.sub$Fire_nr2 <- (as.numeric(sites.sub$Fire_nr>0))
sites.sub$LU <- as.character(sites.sub$Predominant_habitat) # Reclassify
sites.sub$LU[grep("Secondary",sites.sub$LU,ignore.case = T)] <- "Secondary vegetation"
sites.sub$LU <- factor(sites.sub$LU,levels = c("Primary forest","Primary non-forest","Secondary vegetation",
                                               "Plantation forest","Pasture","Cropland","Urban"))
sites.sub$LU<-relevel(sites.sub$LU,ref="Primary forest")

res <- data.frame()
f <- paste("logabund ~","Fire_nr2 + (LU|SS) + (1|SSB)")
mod_a <- glmer(as.formula(f),data=sites.sub,family = "gaussian")
summary(mod_a)
MuMIn::r.squaredGLMM(mod_a)[1]
plot(effects::allEffects(mod_a))
res <- rbind(res,data.frame(type= "logabund",eff = fixef(mod_a)[2],
                            se = se.fixef(mod_a)[2],
                            N_obs = mod_a@devcomp$dims[1])
             )

# Fit SR
f <- paste("Species_richness ~","Fire_nr2 + (LU|SS) + (1|SSB)")
mod_s <- glmer(as.formula(f),data=sites.sub,family = "poisson")
MuMIn::r.squaredGLMM(mod_s)[1]
plot(effects::allEffects(mod_s))
res <- rbind(res,data.frame(type= "Species_richness",eff = fixef(mod_s)[2],
                            se = se.fixef(mod_s)[2],
                            N_obs = mod_s@devcomp$dims[1])
)

f <- paste("logsimp ~","Fire_nr2 + (LU|SS) + (1|SSB)")
mod_si <- glmer(as.formula(f),data=sites.sub,family = "gaussian")
MuMIn::r.squaredGLMM(mod_si)[1]
plot(effects::allEffects(mod_si))
res <- rbind(res,data.frame(type= "logsimp",eff = fixef(mod_si)[2],
                            se = se.fixef(mod_si)[2],
                            N_obs = mod_si@devcomp$dims[1])
)

# Make a quick Plot
g = ggplot(res,aes(x=type,y=eff,ymin=eff-se,ymax=eff+se,label=paste("N=",N_obs)))+
  geom_hline(yintercept=0)+
  geom_pointrange() + geom_text(vjust=1.5)+
  coord_flip() + labs(x="",y="Fixed effect",title = "Presence/Absence of fire before sampling \n (year 2000 - sampling start)")+
  scale_y_continuous(limits = c(-0.1,0.125),breaks = pretty_breaks())
ggsave("BurnPresenceBiodiversity.png",plot=g)

#### Test: Does the time since fire has an effect on Biodiversity #########
# Actual contrast would be between site with different time since disturbance
sub <- subset(sites.sub,!is.na(TimeSinceLastFire))

# Fit abundance
f <- paste("logabund ~","TimeSinceLastFire  + TGrouping +(LU|SS) + (1|SSB)")
mod_a <- glmer(as.formula(f),data=sub,family = "gaussian")
MuMIn::r.squaredGLMM(mod_a)[1]
plot(effects::allEffects(mod_a))
f <- paste("logabund ~","s(TimeSinceLastFire) + TGrouping")
mod_ga <- gamm4(as.formula(f),data=sub,random = ~ (LU|SS)+(1|SSB),family = "gaussian")
plot(mod_ga$gam, pages=1)
vis.gam(mod_ga$gam, view=c("TimeSinceLastFire"),plot.type="persp",theta=-35,color="heat")

## Get optimal randoms
# Fit SR
f <- paste("Species_richness ~","TimeSinceLastFire + (LU|SS) + (1|SSB)")
mod_s <- glmer(as.formula(f),data=sub,family = "poisson")
MuMIn::r.squaredGLMM(mod_s)[1]
plot(effects::allEffects(mod_s))
f <- paste("Species_richness ~","s(TimeSinceLastFire)")
mod_gs <- gamm4(as.formula(f),data=sub,random = ~ (LU|SS)+(1|SSB),family = "poisson")
plot(mod_gs$gam, pages=1)
vis.gam(mod_gs$gam, view=c("TimeSinceLastFire","Fire_nr"),plot.type="persp",theta=-125,color="heat")

# Simpson
rc.si <- compare_randoms(dataset = sub,responseVar = "logsimp",fitFamily = "gaussian",siteRandom = F,
                        fixedFactors = c("Predominant_habitat"))
# Fit SR
f <- paste("logsimp ~","TimeSinceLastFire + Fire_nr+",rc.si$best.random)
mod_si <- glmer(as.formula(f),data=sub,family = "gaussian")
MuMIn::r.squaredGLMM(mod_si)[1]
plot(effects::allEffects(mod_si))
f <- paste("logsimp ~","s(TimeSinceLastFire)")
mod_gi <- gamm4(as.formula(f),data=sub,random = ~ (LU|SS)+(1|SSB),family = "gaussian")
plot(mod_gi$gam, pages=1)
vis.gam(mod_gi$gam, view=c("TimeSinceLastFire","Fire_nr"),plot.type="persp",theta=-125,color="heat")
