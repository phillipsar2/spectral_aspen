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
library(Matrix)
library(parallel)

# Resources
# https://devillemereuil.legtux.org/wp-content/uploads/2021/09/tuto_en.pdf
# https://wildanimalmodels.org/docs/univariate/mcmcglmm/fixed_random_effects/#computing-heritability-without-accounting-for-fixed-effects


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

colnames(K) <- geno_ID_df$ID_genotype_JGI # add names that match genotype IDs
rownames(K) <- geno_ID_df$ID_genotype_JGI
K[1:5,1:5]

### Convert kinship matrix to inverted spare matrix

# ev <- eigen(K) # check if eigen values are positive
# ev$values # last one is negative

## Set to zero then try again
K[K < 0] <- 0

## Invert (by solving) and convert to sparse matrix
K <- as(K, "sparseMatrix") # must be a sparse matrix!!

# cond_number <- kappa(K)
# print(cond_number)
# 
# eps <- 100 # resolve singularity
# cond_number <- kappa(K + diag(eps, nrow(K)))
# cond_number < 10e12
# cond_number
# kinship_sparse <- Matrix(K + diag(eps, nrow(K)), sparse = TRUE)
# K_inv <- solve(kinship_sparse)# create inverted matrix
# kinship_sparse <- Matrix(K, sparse = TRUE)

k <- solve(K) # invert matrix

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

## Make sure categorical variables are factors
sv_sub_pl$ploidy_call <- as.factor(sv_sub_pl$ploidy_call)
sv_sub_pl$ID_genotype_JGI <- as.factor(sv_sub_pl$ID_genotype_JGI)
sv_sub_pl$age_scaled <- scale(sv_sub_pl$age_at_spectral_measurement)
str(sv_sub_pl)
dim(sv_sub_pl)

# (1) Visualize data ----
n_wavelengths <- colnames(sv_sub_pl) %>% str_count("w") %>% sum()
n_wavelengths

# Sample map
library(raster)
library(sf)
library(ggplot2)
library(ggspatial)
library(USAboundaries)
library(purrr)

states <- us_states()

events_sf <- unique(sv_sub_pl[,12:13]) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) 
dim(events_sf)

events_sf$ploidy <- unique(sv_sub_pl[,c(3,12:13)])$ploidy_call

ggplot() + 
  geom_sf(data = states, size = 4, color = "black", fill = NA) +
  # ggtitle("Lake Erie Outline") + 
  # coord_sf() +
  geom_sf(data = events_sf, size = 2, aes(color = ploidy), alpha = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        # axis.text = element_blank(),
        axis.line = element_blank(),
        # axis.ticks = element_blank(),
        legend.title = element_blank(),
        # legend.text = element_blank()
  ) +
  xlim(c(-123, -103)) + ylim(c(32,47)) +
  scale_color_manual(values = c("blue", "red"))

# Assess normality of the data
## One wavelength
trait_means <- aggregate( w701 ~ ID_genotype_JGI + ploidy_call , sv_sub_pl, FUN = mean)
ggplot(trait_means, aes(x = w701)) +
  geom_histogram() +
  theme_bw()

ggplot(sv_sub_pl, aes(y = w701, x = year.x)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  ylab("% Reflectance") +
  xlab("Measurement year")

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
# table(sv_sub_pl$ID_genotype_JGI, sv_sub_pl$year.x)

wv = "w701"
sv_sub_comp <- sv_sub_pl %>% 
  filter(year.x == '2023') %>% # year are non-overlapping
  select(all_of(wv), ploidy_call, age_scaled, ID_genotype_JGI, id, state, year.x, leaf_num) %>%
  filter(complete.cases(.)) %>%
  rename(year = year.x) %>%
  # group_by(ID_genotype_JGI, id) %>% 
  # slice_sample(n = 3) %>% # How to subsample the genotypes????
  as.data.frame()

table(sv_sub_pl$ID_genotype_JGI, sv_sub_pl$year.x)
# length(unique(sv_sub_comp$id))

# > Run model ----
## id is nested within ID_genotype_JGI -- id must have unique labels

## >> Full model ----

## Specify priors
# with sufficient sample size, the prior shouldn't matter much
# MCMCglmm uses an inverse-gamma distribution for residual and random priors 
# fixed effects have a normal prior
prior <- list(R = list(V = 1, nu = 0.002), # priors on Residual
              G = list(G1 = list(V = 1, nu = 0.002), # priors on Random effects
                       G2 = list(V = 1, nu = 0.002),
                       G3 = list(V = 1, nu = 0.002),
                       G4 = list(V = 1, nu = 0.002)
              ))


mod <- MCMCglmm( as.formula( paste0(wv, "~ ploidy_call + age_scaled" ) ),  
                   random = ~ ID_genotype_JGI + id + year + state , # dropped year
                   ginverse = list(ID_genotype_JGI = k), # link kinship matrix
                   data = sv_sub_comp,
                   family = "gaussian",
                   prior = prior,
                   nitt = 200000, # number of iterations
                   burnin = 10000, # iterations to drop
                   thin = 100, # iterations to keep, determined by degree of autocorrelation
                   # pl = T, # estimate BLUPs
                pr = T) # pr = T means the posterior distribution of the random effects were estimated (the BLUPs)


saveRDS(mod, "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/h2_mods/mod.rds")

## >> Model without state ----
prior.state <- list(R = list(V = 1, nu = 0.002), # priors on Residual
              G = list(G1 = list(V = 1, nu = 0.002), # priors on Random effects
                       G2 = list(V = 1, nu = 0.002),
                       G3 = list(V = 1, nu = 0.002)
              ))
mod.state <- MCMCglmm( as.formula( paste0(wv, "~ ploidy_call + age_scaled" ) ),  
                 random = ~ ID_genotype_JGI + id + year , # dropped year
                 ginverse = list(ID_genotype_JGI = k), # link kinship matrix
                 data = sv_sub_comp,
                 family = "gaussian",
                 prior = prior.state,
                 nitt = 200000, # number of iterations
                 burnin = 10000, # iterations to drop
                 thin = 100, # iterations to keep, determined by degree of autocorrelation
                 # pl = T, # estimate BLUPs
                 pr = T) # pr = T means the posterior distribution of the random effects were estimated (the BLUPs)

saveRDS(mod.state, "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/h2_mods/mod.state.rds")

# mod <- lapply(mod, function(m) m$Sol)
# mod <- do.call(mcmc.list, mod)

# assess burn-in length
# mod.burn <- MCMCglmm( as.formula( paste0(wv, "~ ploidy_call + age_scaled" ) ),
#                  random = ~ ID_genotype_JGI + id + state ,
#                  ginverse = list(ID_genotype_JGI = K), # link kinship matrix
#                  data = sv_sub_comp,
#                  family = "gaussian",
#                  prior = prior,
#                  nitt = 70000, # number of iterations
#                  burnin = 1, # iterations to drop
#                  thin = 100, # iterations to keep
#                  pr = T)
# plot(mod.burn[["VCV"]])

# saveRDS(mod, "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/h2_mods/test.mod.rds")
mod <- readRDS("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/h2_mods/mod.rds")
mod.state <- readRDS("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/h2_mods/mod.state.rds")

# (3) Evaluate chain convergence and autocorrelation ----

## > Chain mixing ----
plot(mod.state$VCV) # random effects; units are the residuals
plot(mod.state$Sol) # fixed effects, also has the BLUPs

# library(coda)
# 
# par(mfrow=c(4,2), mar=c(2,2,1,2))
# gelman.plot(mod, auto.layout=F)
# 
# ### Estimate scale reduction factors
# gelman.diag(mod) # good is close to 1
# 
# ### Look at chain mixing
# par(mfrow=c(8,2), mar=c(2, 1, 1, 1))
# plot(mod, ask=T, auto.layout=F)

## > Assess the effective sample size - the corrected samples size after ----
## removing autocorrelated samples
effectiveSize(mod$Sol[,1:3])
effectiveSize(mod[["VCV"]]) # ideally want about 1,000

autocorr.diag(mod$VCV) # assess whether the thinning is appropriate (want small numbers)
autocorr.diag(mod$Sol[,1:3])

## > Chain convergence ----
heidel.diag(mod$VCV)

## > Check estimates look like observed data ----
# solutions <- data.frame(mod$Sol) # BLUPs
# geno_blups <- select(solutions, starts_with("ID_genotype_JGI.")) %>%
#   colMeans() %>%
#   as.data.frame()
# geno_blups$ID <- str_split(rownames(geno_blups), pattern = "\\.", simplify = T)[,2]
# geno_blups <- arrange(geno_blups, ID)
# 
# geno_means <- aggregate( as.formula( paste0( wv, " ~ ID_genotype_JGI" ) ) , sv_sub_pl, FUN = mean)

# intercept <- mean(solutions$X.Intercept.)
# plot( x = geno_values$w701, y = intercept + geno_values[,3],
#       xlab = "observed", ylab = "estimated")
# abline(a = 0, b = 1)


# (4) Model summary ----
summary(mod) #  DIC: 7601.521 
summary(mod.state) #  DIC: 7600.978 

# (5) Compute heritability estimate ----
## Additive genetic variance
mean(mod.state$VCV[,"ID_genotype_JGI"])

colMeans(mod$VCV) # Variances

## h2
### Compute the variance in model predictions (without accounting for random effects) as the matrix product of predictors by parameter estimates
# mod$X # fixed effect design matrix
# mod$Sol # Posterior Distribution of MME solutions
predictions <- mod.state$X %*% t(mod.state$Sol[,1:3])
dim(mod.state$X)
dim(mod.state$Sol)

### Compute the variance of each posterior sample (columns)
fixef_variance <- apply(predictions, MARGIN = 2, var)
# median(fixef_variance)
# hist(fixef_variance)

herit <- mod.state$VCV[,"ID_genotype_JGI"] / (mod.state$VCV[,"ID_genotype_JGI"] + mod.state$VCV[,"id"] + mod.state$VCV[,"units"] + fixef_variance)

effectiveSize(herit) # check to make sure posterior sampling is over 1000
hist(herit) # posterior distribution

mean(herit)
median(herit)
HPDinterval(herit, prob = 0.95)


### Parallelize it ----
year = "2022"
n_wavelengths <- colnames(sv_sub_pl) %>% str_count("w") %>% sum() - 1 
# wv = "w351.1"

set.seed(1)
mod <- mclapply(801:n_wavelengths, function(i) {
  wv <- colnames(sv_sub_pl)[-c(1:19)][i]
  print(wv)
  
  dir <- paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/h2_mods/", year, "/", wv)
  dir.create(dir, showWarnings = T)
  
  sv_sub_comp <- sv_sub_pl %>% 
    filter(year.x == year) %>% # year are non-overlapping
    select(all_of(wv), ploidy_call, age_scaled, ID_genotype_JGI, id, year.x, leaf_num) %>%
    filter(complete.cases(.)) %>%
    rename(year = year.x) %>%
    # group_by(ID_genotype_JGI, id) %>% # Subsample replicates to 3 per genotype
    # slice_sample(n = 3) %>% 
    as.data.frame()
  
  prior <- list(R = list(V = 1, nu = 0.002), # priors on Residual
                G = list(G1 = list(V = 1, nu = 0.002), # priors on Random effects
                         G2 = list(V = 1, nu = 0.002)
                ))
  
  mod <- MCMCglmm( as.formula( paste0(wv, "~ ploidy_call + age_scaled" ) ),  
                   random = ~ ID_genotype_JGI + id, # run within a single year; dropped state
                   ginverse = list(ID_genotype_JGI = k), # link kinship matrix
                   data = sv_sub_comp,
                   family = "gaussian",
                   prior = prior,
                   nitt = 200000, # number of iterations
                   burnin = 10000, # iterations to drop
                   thin = 100) # iterations to keep, determined by degree of autocorrelation

  saveRDS(mod, paste0(dir, "/mod.", wv ,".rds"))
  
  # Plot posterior probability density distributions and trace plot
  pdf(paste0(dir, "/trace_posterior.", wv,".pdf"))
  plot(mod$VCV) # random effects; units are the residuals
  plot(mod$Sol) # fixed effects, also has the BLUPs
  dev.off()
  
  # Heritability
  predictions <- mod$X %*% t(mod$Sol) 
  fixef_variance <- apply(predictions, MARGIN = 2, var) # variance of fixed effects
  
  h2 <- mod$VCV[,"ID_genotype_JGI"] / (mod$VCV[,"ID_genotype_JGI"] + mod$VCV[,"id"] + mod$VCV[,"units"] + fixef_variance)
  
  pdf(paste0(dir, "/h2_posterior.", wv,".pdf"))
  hist(h2)
  dev.off()
  
  h2_df <- c( median(h2), mean(h2), HPDinterval(h2, prob = 0.95) ) %>% as.data.frame() %>% t()
  write.table(h2_df, paste0(dir, "/h2_estimates.", wv, ".csv"), row.names = F, col.names = F, sep = ",")
  
}, 
mc.cores = 10)

# Plot the heritability results ----

### 2023
h2_files <- Sys.glob("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/h2_mods/2023/*/*.csv")
h2_list <- lapply(h2_files, function(x) read.table(x, sep = ","))
h2_df <- do.call(rbind, h2_list)
colnames(h2_df) <- c("median", "mean", "lower_ci", "upper_ci")

h2_df$wv <- str_split(h2_files, pattern = "h2_estimates.w", simplify = T)[,2] %>%
  gsub(pattern = ".csv", replacement = "") %>%
  as.numeric()

head(h2_df) 
str(h2_df)

ggplot(h2_df, aes(x = wv, y = median)) +
  geom_line() +
  theme_bw() +
  geom_line(aes(y = lower_ci), col = "gray") +
  geom_line(aes(y = upper_ci), col = "gray") +
  ylab("Median h2") + xlab("Wavelength") + ggtitle("2023") +
  scale_x_continuous(breaks=seq(from = 300, to = 2500, by = 200))

### 2022
h2_files_22 <- Sys.glob("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/h2_mods/2022/*/*.csv")
h2_list_22 <- lapply(h2_files_22, function(x) read.table(x, sep = ","))
h2_df_22 <- do.call(rbind, h2_list_22)
colnames(h2_df_22) <- c("median", "mean", "lower_ci", "upper_ci")

h2_df_22$wv <- str_split(h2_files_22, pattern = "h2_estimates.w", simplify = T)[,2] %>%
  gsub(pattern = ".csv", replacement = "") %>%
  as.numeric()

head(h2_df_22) 
str(h2_df_22)

ggplot(h2_df_22, aes(x = wv, y = median)) +
  geom_line() +
  theme_bw() +
  geom_line(aes(y = lower_ci), col = "gray") +
  geom_line(aes(y = upper_ci), col = "gray") +
  ylab("Median h2") + xlab("Wavelength") + ggtitle("2022") +
  scale_x_continuous(breaks=seq(from = 300, to = 2500, by = 200))


# (6) Effect sizes ----
## Ploidy effect sizes
n_wavelengths <- colnames(sv_sub_pl) %>% str_count("w") %>% sum() - 1 
mod_files <- Sys.glob("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/h2_mods/2023/*/mod.*.rds")

ploidy_eff <- matrix(NA, nrow = n_wavelengths, ncol = 4)

for (j in 1:n_wavelengths){
  mod_wv <- readRDS(mod_files[j])
  fixed_eff <- median(mod_wv$Sol[,2])
  fixed_CI <- HPDinterval(mod_wv$Sol[,2], prob = 0.95) %>% as.data.frame()
  
  wv_name <- str_split(mod_files[j], pattern = "mod.w", simplify = T)[,2] %>%
    gsub(pattern = ".rds", replacement = "")
  
  ploidy_eff[j,] <- c( wv_name, fixed_eff, fixed_CI ) %>% unlist() 
}

colnames(ploidy_eff) <- c("wv", "median", "lower", "upper")
ploidy_eff_num <- apply(ploidy_eff, 2, as.numeric)

ggplot(ploidy_eff_num, aes(x = wv, y= median)) +
  geom_line() +
  theme_bw() +
  geom_line(aes(y = lower), col = "gray") +
  geom_line(aes(y = upper), col = "gray") +
  ylab("Median effect of triploidy") + xlab("Wavelength (nm)") + ggtitle("2023") +
  scale_x_continuous(breaks=seq(from = 300, to = 2500, by = 200))

# Individual models
random_eff <- apply(mod.state$VCV, 2, median)
random_CI <- HPDinterval(mod.state$VCV, prob = 0.95) %>% as.data.frame()

fixed_eff <- apply(mod.state$Sol[,1:3], 2, median)
fixed_CI <- HPDinterval(mod.state$Sol[,1:3], prob = 0.95) %>% as.data.frame()

eff_df <- rbind(random_CI, fixed_CI)
eff_df$eff <- c(random_eff, fixed_eff)
rownames(eff_df)[5] <- "Intercept"
eff_df$variable <- rownames(eff_df)

filter(eff_df, variable != "Intercept") %>%
  ggplot() +
  geom_segment( aes(x=variable, xend=variable, y=lower, yend=upper), color="grey", size = 1) +
  geom_point( aes(x=variable, y=eff), color="black", size=2 ) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Median of the posterior distribution") +
  xlab("Variable") +
  ylim(-10,40)