### Title: MCMCglmm SWUS models
### Author: Alyssa Phillips
### Date: 5/12/2025

# install.packages("MCMCglmm")
library(MCMCglmm)
library(dplyr)
library(ggplot2)
library(stringr)
library(viridis)
library(cowplot)

# Resources
# https://devillemereuil.legtux.org/wp-content/uploads/2021/09/tuto_en.pdf


# (0) Load data ----
## Phenotype data
sv_data <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/phenotypes/final_leaf-level_spectra_2022-2023_03-14-24_sg_greenhouse-joined.csv")
sv_data[, c(1:5,7:9,12,16)] <- lapply(sv_data[, c(1:5,7:9,12,16)], as.factor)
str(sv_data)
dim(sv_data)

## Kinship matrix
K <- read.table("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gusrelate/WUSGgrm.filt.2025-05-05.csv", 
                row.names = 1, sep = ",", header = T)
K[1:5,1:5]
dim(K)

### Change names to genotype IDs
key <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/metadata/spectra_genotypes.csv")

genos <- colnames(K)
geno_ID_df <- merge(x = data.frame(genos), y = key, by.x = "genos", by.y = "Accession")
geno_ID_df$genos == genos # check if ordered correctly

colnames(K) <- geno_ID_df$ID_genotype_JGI
rownames(K) <- geno_ID_df$ID_genotype_JGI
K[1:5,1:5]

## Subset phenotype data to those in kinship matrix
unique(sv_data$ID_genotype_JGI) %>% length()
dim(sv_data)

sv_sub <- sv_data[sv_data$ID_genotype_JGI %in% colnames(K),]
dim(sv_sub)
unique(sv_sub$ID_genotype_JGI) %>% length()

## Amend df with new ploidy calls
gbs2ploidy_df <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/spectral_aspen.ploidycalls.2025-04-15.csv")
geno_ID_df$genos %in% gbs2ploidy_df$sample

pl_temp <- merge(x = gbs2ploidy_df, y = geno_ID_df, by.x = "sample", by.y =  "genos")

sv_sub_pl <- merge(x = pl_temp, y = sv_sub, by = "ID_genotype_JGI" )

sv_sub_pl %>% select(ID_genotype_JGI, id, ploidy, )

# (1) Visualize data ----
n_wavelengths <- colnames(sv_sub_pl) %>% str_count("w") %>% sum()
n_wavelengths

# Assess normality of the data
## One wavelength
trait_means <- aggregate( w2479.5 ~ ID_genotype_JGI + ploidy_call , sv_sub_pl, FUN = mean)
ggplot(trait_means, aes_string(x = "w2479.5")) +
  geom_density() +
  theme_bw()

## All the wavelengths
sv_long <- sv_sub_pl%>%
  tidyr::pivot_longer(
  cols = `w346.7`:`w2502.6`, 
  names_to = "wavelength",
  values_to = "value"
)
sv_long$wavelength <- as.factor(sv_long$wavelength)

ggplot(sv_long, aes(x=value, color = wavelength)) +
  geom_density(alpha=0.6) +
  theme_bw() +
  theme(
    legend.position="none"
  ) +
  ylab("density") +
  xlab("% reflectance")

# Reflectance vs state
sv_long %>%
  filter(wavelength %in% sample(levels(sv_long$wavelength), 9)) %>% # plot 9 random wavelength
  ggplot(aes(y = value, x = state)) +
    geom_boxplot() +
    theme_bw() +
  facet_wrap(~wavelength, scale="free_y") +
    theme(
      legend.position="none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    ) +
    ylab("% Reflectance") +
    xlab("State")

# Reflectance vs ploidy
sv_long %>%
  filter(wavelength %in% sample(levels(sv_long$wavelength), 9)) %>% # plot 9 random wavelength
  ggplot(aes(y = value, x = ploidy_call)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~wavelength, scale="free_y") +
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  ylab("% Reflectance") +
  xlab("Ploidy")

# Reflectance vs year
sv_long %>%
  filter(wavelength %in% sample(levels(sv_long$wavelength), 9)) %>% # plot 9 random wavelength
  ggplot(aes(y = value, x = year.x)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~wavelength, scale="free_y") +
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  ylab("% Reflectance") +
  xlab("Measurement year")

# Reflectance vs elevation
sv_long %>%
  filter(wavelength %in% sample(levels(sv_long$wavelength), 9)) %>% # plot 9 random wavelength
  ggplot(aes(y = value, x = Elevation_m)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~wavelength, scale="free_y") +
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  ylab("% Reflectance") +
  xlab("Elevation (m)")

# Reflectance geographic structure 
sv_long %>%
  filter(wavelength %in% sample(levels(sv_long$wavelength), 9)) %>% # plot 9 random wavelength
  ggplot(aes(y = Latitude, x = Longitude, col = value)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~wavelength, scale="free_y") +
  theme(
    # legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) + 
  scale_color_viridis_c() # really need to have unique legends and scales for each plot

# Reflectance vs age at spectral measurement
sv_long %>%
  filter(wavelength %in% sample(levels(sv_long$wavelength), 9)) %>% # plot 9 random wavelength
  ggplot(aes(y = value, x = age_at_spectral_measurement)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~wavelength, scale="free_y") +
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  ylab("% Reflectance") +
  xlab("Age at measurement")

# (2) Build model ----

# > Prepare dataframe ----
str(sv_sub_pl)
sv_sub_pl$ploidy_call <- as.factor(sv_sub_pl$ploidy_call)
# sv_sub_pl$sample <- as.factor(sv_sub_pl$sample)
sv_sub_pl$ID_genotype_JGI <- as.factor(sv_sub_pl$ID_genotype_JGI)

sv_sub_pl$age_scaled <- scale(sv_sub_pl$age_at_spectral_measurement)
# hist(sv_sub_pl$age_at_spectral_measurement)
# hist(sv_sub_pl$age_scaled)

# > Specify priors ----
# with sufficient sample size, the prior shouldn't matter much
# MCMCglmm uses an inverse-gamma distribution for residual and random priors 
# fixed effects have a normal prior

prior <- list(R = list(V = 1, nu = 0.002), # priors on Residual
              G = list(G1 = list(V = 1, nu = 0.002), # priors on Random effects
                       G2 = list(V = 1, nu = 0.002),
                       G3 = list(V = 1, nu = 0.002),
                       G4 = list(V = 1, nu = 0.002)
                       ))

# > Run model ----
## id is nested within ID_genotype_JGI -- id must have unique labels
mod <- MCMCglmm(w2388.3 ~ ploidy_call + age_scaled,  
                   random = ~ ID_genotype_JGI + ID_genotype_JGI:id + state + year.x,
                   ginverse = list(ID_genotype_JGI = K), # link kinship matrix
                   data = sv_sub_pl,
                   family = "gaussian",
                   prior = prior,
                   nitt = 100000, # number of iterations
                   burnin = 10000, # iterations to drop
                   thin = 10, # iterations to keep (every 10)
                   pl = T, pr = T) # pr = T means the posterior distribution of the random
#plot(ft_mod$VCV)


# Evaluate chain convergence and autocorrelation

# Extract efect sizes ----
###extract the breeding values for flowering time and plot them against true value
solutions<-data.frame(ft_mod$Sol)
subset_cols<- select(solutions, starts_with("X.intercept."),starts_with("id."), starts_with("gxe."), starts_with("trmt."))
# interc<-posterior.mode(mod1$Sol[,1])
subset_cols_modes<-posterior.mode(subset_cols)
bv<-subset_cols_modes
class(subset_cols_modes)
names(bv)
head(subset_cols_modes)

d6tmp<-
  expand.grid( unique(dataset1_accession_selection$id), unique(dataset3_all_replicates_and_treatments$trmt) )#update these to reference the .rda file
d6tmp <- d6tmp %>% rename(id=Var1, trmt = Var2) %>%
  mutate(gxe=paste(id,trmt,sep=".")) %>%
  mutate(id_name=paste0("id.",id))%>%
  mutate(gxe_name=paste0("gxe.",gxe))%>%
  mutate(trmt_name=paste0("trmt.",trmt))
head(d6tmp)

#unique(d6tmp$gxe_name)
d6tmp$id_value = bv[d6tmp$id_name]
d6tmp$gxe_value = bv[d6tmp$gxe_name]#"gxe.9775.Daily_min"
#  mutate(across(where(is.character), ~ str_trim(.)))
d6tmp$trmt_value = bv[d6tmp$trmt_name]
d6tmp$Intercept = bv["X.Intercept."]
head(d6tmp)
d6tmp<-d6tmp %>% mutate(breeding_value_flowering_time = (id_value + gxe_value + trmt_value + Intercept))

d6<-merge(d6tmp,dataset5, by = c("id","trmt"))
d6$ft<-as.numeric(d6$ft)
d6$trmt<-as.character(d6$trmt)

bv_ft<-d6
# d6 <- d6 %>%
#   mutate(trmt = factor(trmt, levels = treatment_order))
ggplot(data = d6, aes(x = ft, y = breeding_value_flowering_time, color = trmt))+
  geom_point()

colnames(d6)
d6$breeding_value_flowering_time[which(is.nan(d6$breeding_value_flowering_time))] = NA
d6$breeding_value_flowering_time[which(d6$breeding_value_flowering_time==Inf)] = NA
d6$ft[which(is.nan(d6$ft))] = NA
d6$ft[which(d6$ft==Inf)] = NA

lm(d6$ft~ d6$breeding_value_flowering_time)#still not working, look into this later

#extract heritability of flowering time overall for the field
#ft_mod <- MCMCglmm(ft ~ 1,  random = ~ id + trmt + gxe + tray + zone, data = dataset5[!is.na(dataset5$ft),] , pl = T, pr = T)
ft_herit<-ft_mod$VCV[,"id"] / (ft_mod$VCV[,"id"] + ft_mod$VCV[,"trmt"] + ft_mod$VCV[,"gxe"] + ft_mod$VCV[,"tray"] + ft_mod$VCV[,"zone"]+ ft_mod$VCV[,"units"])
ft_herit
effectiveSize(ft_herit)
mean(ft_herit)
HPDinterval(ft_herit)

ft_herit<-(ft_mod$VCV[,"id"] + ft_mod$VCV[,"gxe"]) / (ft_mod$VCV[,"id"] + ft_mod$VCV[,"trmt"] + ft_mod$VCV[,"gxe"] + ft_mod$VCV[,"tray"] + ft_mod$VCV[,"zone"]+ ft_mod$VCV[,"units"])



#######silique number
#####
hist(dataset5$sil)
temp<-dataset5 %>% filter(sil >0) %>%
  mutate(log_sil_number = log(sil))
hist(temp$log_sil_number)
sum(is.na(temp$log_sil_number))
log_sil_mod <- MCMCglmm(log_sil_number ~ 1,  random = ~ id + trmt + gxe + tray + zone, data = temp[!is.na(temp$log_sil_number),] , pl = T, pr = T)

###extract the breeding values for flowering time and plot them against true value
solutions<-data.frame(log_sil_mod$Sol)
subset_cols<- select(solutions, starts_with("X.intercept."),starts_with("id."), starts_with("gxe."), starts_with("trmt."))
# interc<-posterior.mode(mod1$Sol[,1])
subset_cols_modes<-posterior.mode(subset_cols)
bv<-subset_cols_modes
class(subset_cols_modes)
names(bv)
head(subset_cols_modes)

d6tmp<-
  expand.grid( unique(dataset1_accession_selection$id), unique(dataset3_all_replicates_and_treatments$trmt) )
d6tmp <- d6tmp %>% rename(id=Var1, trmt = Var2) %>%
  mutate(gxe=paste(id,trmt,sep=".")) %>%
  mutate(id_name=paste0("id.",id))%>%
  mutate(gxe_name=paste0("gxe.",gxe))%>%
  mutate(trmt_name=paste0("trmt.",trmt))
head(d6tmp)


#unique(d6tmp$gxe_name)
d6tmp$id_value = bv[d6tmp$id_name]
d6tmp$gxe_value = bv[d6tmp$gxe_name]#"gxe.9775.Daily_min"
#  mutate(across(where(is.character), ~ str_trim(.)))
d6tmp$trmt_value = bv[d6tmp$trmt_name]
d6tmp$Intercept = bv["X.Intercept."]
head(d6tmp)
d6tmp<-d6tmp %>% mutate(bv_log_sil_number = (id_value + gxe_value + trmt_value + Intercept))

d6<-merge(d6tmp,temp, by = c("id","trmt"))
d6$log_sil_number<-as.numeric(d6$log_sil_number)
d6$trmt<-as.character(d6$trmt)
# d6 <- d6 %>%
#   mutate(trmt = factor(trmt, levels = treatment_order))
bv_sil<-d6
ggplot(data = d6, aes(x = log_sil_number, y = bv_log_sil_number, color = trmt))+
  geom_point()

ggplot(data = d6, aes(x = exp(log_sil_number), y = exp(bv_log_sil_number), color = trmt))+
  geom_point()


#extract heritability of flowering time overall for the field
#ft_mod <- MCMCglmm(ft ~ 1,  random = ~ id + trmt + gxe + tray + zone, data = dataset5[!is.na(dataset5$ft),] , pl = T, pr = T)
log_sil_herit<-log_sil_mod$VCV[,"id"] / (log_sil_mod$VCV[,"id"] + log_sil_mod$VCV[,"trmt"] + log_sil_mod$VCV[,"gxe"] + log_sil_mod$VCV[,"tray"] + log_sil_mod$VCV[,"zone"]+ log_sil_mod$VCV[,"units"])
log_sil_herit
effectiveSize(log_sil_herit)
mean(log_sil_herit)
HPDinterval(log_sil_herit)
#H2 of sil number is very very low? 0.04-0.09
log_sil_herit<-(log_sil_mod$VCV[,"id"] + log_sil_mod$VCV[,"gxe"]) / (log_sil_mod$VCV[,"id"] + log_sil_mod$VCV[,"trmt"] + log_sil_mod$VCV[,"gxe"] + log_sil_mod$VCV[,"tray"] + log_sil_mod$VCV[,"zone"]+ log_sil_mod$VCV[,"units"])
##########

######survival
range(dataset5$surv)
hist(dataset5$surv, breaks = 4)
#dataset5$surv_tf<-ifelse(dataset5$surv == 1, T, F)
summary(dataset5$surv)
dataset5$surv_cat<-ifelse(dataset5$surv == 1, "yes", "no")
str(dataset5$surv_cat)
dataset5$surv<-as.factor(dataset5$surv)
surv_mod <- MCMCglmm(surv_cat ~ 1,  random = ~ id + trmt + gxe + tray + zone, data = dataset5[!is.na(dataset5$surv_cat),] , pl = T, pr = T, family = "categorical")#pr = T means the posterior distribution of the random

###extract the breeding values for flowering time and plot them against true value
solutions<-data.frame(surv_mod$Sol)
subset_cols<- select(solutions, starts_with("X.intercept."),starts_with("id."), starts_with("gxe."), starts_with("trmt."))
# interc<-posterior.mode(mod1$Sol[,1])
subset_cols_modes<-posterior.mode(subset_cols)
bv<-subset_cols_modes
summary(subset_cols_modes)
names(bv)
head(subset_cols_modes)

d6tmp<-
  expand.grid( unique(dataset1_accession_selection$id), unique(dataset3_all_replicates_and_treatments$trmt) )
d6tmp <- d6tmp %>% rename(id=Var1, trmt = Var2) %>%
  mutate(gxe=paste(id,trmt,sep=".")) %>%
  mutate(id_name=paste0("id.",id))%>%
  mutate(gxe_name=paste0("gxe.",gxe))%>%
  mutate(trmt_name=paste0("trmt.",trmt))
head(d6tmp)

#unique(d6tmp$gxe_name)
d6tmp$id_value = bv[d6tmp$id_name]
d6tmp$gxe_value = bv[d6tmp$gxe_name]#"gxe.9775.Daily_min"
#  mutate(across(where(is.character), ~ str_trim(.)))
d6tmp$trmt_value = bv[d6tmp$trmt_name]
d6tmp$Intercept = bv["X.Intercept."]
head(d6tmp)
d6tmp<-d6tmp %>% mutate(bv_surv = (id_value + gxe_value + trmt_value + Intercept))

d7<-merge(d6tmp,dataset5, by = c("id","trmt"))
head(d7)
bv_surv<-d7
#d7$surv<-as.numeric(d6$surv)
#d6$trmt<-as.character(d6$trmt)
# d6 <- d6 %>%
#   mutate(trmt = factor(trmt, levels = treatment_order))
dev.off()
ggplot(data = d7, aes(x = surv, y = bv_surv, color = trmt))+
  geom_point()#where are the zero values for the original surv df

range(d6$surv)


#### merge 3 breeding values
all_breeding_values<-merge(select(bv_ft,id,trmt,breeding_value_flowering_time,tray,zone), select(bv_sil,bv_log_sil_number,id,tray), by = c("id","tray"))

all_breeding_values<-merge(all_breeding_values, select(bv_surv,bv_surv,id,tray), by = c("id","tray"))#definitely go through this and make sure the merges are okay bc we go from 18,000 ro 25,000 when we add in survival... i think ones that were na in the other fields wont have values and that surivial is the only one without NA for the entire column and that is why it takes it back to the full value? IDK now thos it looks like all rows have a value for silique number and only a handfull are missing for flowering time


all_breeding_values$bv_sil_number_transformed<-exp(all_breeding_values$bv_log_sil_number)

write.csv(all_breeding_values,paste0(mypath,"data/Dataset_5_all_phenotype_breeding_values.csv"),row.names = F)




#####Now how does heritbaility of ft change over env treatments?
treatment<-as.character(unique(dataset5$trmt))
model_outputs<-c()
for (t in treatment){
  tmp<- dat %>% dplyr::filter(trmt == t)
  ft_model_loop <- MCMCglmm(ft ~ 1, random = ~ id + trmt + gxe + tray + zone, data = dataset5)
  herit<-ft_model_loop[["VCV"]][ , "id"] / rowSums(ft_model_loop$VCV)
  model_outputs_ft <- rbind(model_outputs_ft, data.frame(trmt = t, herit = mean(herit), effectiveSize = effectiveSize(herit), lower = HPDinterval(herit)[,1], upper = HPDinterval(herit)[,2]))
  
}

desired_order<-c("Daily max","Daily med","Daily min","Weekly max","Weekly med","Weekly min","Biweekly max","Biweekly med","Biweekly min","Monthly max","Monthly med","Monthly min","Natural CA","-50% CA")


model_outputs$trmt<-as.factor(model_outputs$trmt, levels = desired_order)
ggplot(data = model_outputs, aes(x =trmts, y = herit))+
  geom_point()

###sil
treatment<-as.character(unique(dataset5$trmt))
model_outputs_log_sil<-c()
for (t in treatment){
  tmp<- dat %>% dplyr::filter(trmt == t)
  log_sil_model_loop <- MCMCglmm(log_sil_number ~ 1, random = ~ id + trmt + gxe + tray + zone, data = temp)
  herit<-log_sil_model_loop[["VCV"]][ , "id"] / rowSums(log_sil_model_loop$VCV)
  model_outputs_ft <- rbind(model_outputs_ft, data.frame(trmt = t, herit = mean(herit), effectiveSize = effectiveSize(herit), lower = HPDinterval(herit)[,1], upper = HPDinterval(herit)[,2]))
  print(paste0("done with ",t))
}