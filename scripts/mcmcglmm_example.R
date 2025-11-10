### Title: MCMCglmm - example script
### Author: Alyssa Phillips
### Date: 6/28/2025

library(MCMCglmm)
library(dplyr)
library(ggplot2)
library(stringr)
library(viridis)
library(cowplot)
library(Matrix)
library(parallel)

# Resources
## On general MCMCglmm usage
# https://math.umd.edu/~slud/s770/LecSlides/Lecs19to26/CourseNotesMCMCglmm.pdf

## Moving from lme4 to MCMCglmm
# https://github.com/tmalsburg/MCMCglmm-intro

## On estimating heritability and animal models
# https://wildanimalmodels.org/docs/univariate/mcmcglmm/fixed_random_effects/#computing-heritability-without-accounting-for-fixed-effects
# https://devillemereuil.legtux.org/wp-content/uploads/2021/09/tuto_en.pdf

## On properly extracting predicted values
# https://tomhouslay.com/wp-content/uploads/2017/02/indivvar_plasticity_tutorial_mcmcglmm1.pdf

## On DIC scores and model selection
# https://deepthoughtsandsilliness.blogspot.com/2007/12/focus-on-dic.html
# https://st540.wordpress.ncsu.edu/files/2020/03/compare.pdf

# (0) Load & prep data ----
## > Phenotype data ----
df <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/phenotypes/final_leaf-level_spectra_2022-2023_03-14-24_sg_greenhouse-joined.csv")
str(df)
dim(df)

## Set categorical variables to factors
df[, c(1:5,7:9,12,16)] <- lapply(df[, c(1:5,7:9,12,16)], as.factor) # specify by column numbers


## > Kinship matrix ----
K <- read.table("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gusrelate/WUSGgrm.filt.2025-05-05.csv", 
                row.names = 1, sep = ",", header = T)
K[1:5,1:5] 
dim(K)

colnames(K) %in% df$

### Change names to match dataframe (if needed)
# key <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/metadata/spectra_genotypes.csv")
# 
# genos <- colnames(K)
# geno_ID_df <- merge(x = data.frame(genos), y = key, by.x = "genos", by.y = "Accession")
# geno_ID_df$genos == genos # check if ordered correctly
# 
# colnames(K) <- geno_ID_df$ID_genotype_JGI # add names that match genotype IDs
# rownames(K) <- geno_ID_df$ID_genotype_JGI
# K[1:5,1:5]

### Convert kinship matrix to inverted sparse matrix
### Required by MCMCglmm

#### Set negative values to zero
K[K < 0] <- 0

#### Invert (by solving) and convert to sparse matrix
K <- as(K, "sparseMatrix")
k <- solve(K) # invert matrix

## Subset phenotype data genotypes in kinship matrix
unique(df$ID_genotype_JGI) %>% length()
dim(df)

df_sub <- df[df$ID_genotype_JGI %in% colnames(k),] # ID_genotype_JGI == whatever column that matches your kinship matrix IDs

length(unique(df_sub$ID_genotype_JGI)) == dim(k)[1] # check that subset worked; must be TRUE

## Add Alyssa's ploidy calls to df, if needed
gbs2ploidy_df <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/spectral_aspen.ploidycalls.2025-04-15.csv")

df_sub_pl <- merge(x = gbs2ploidy_df, y = df_sub, by.x = "sample", by.y =  "ID_genotype_JGI")

## Scale age
df_sub_pl$age_scaled <- scale(df_sub_pl$age_at_spectral_measurement)

## Make sure categorical variables are factors
str(df_sub_pl)
df_sub_pl$ploidy_call <- as.factor(df_sub_pl$ploidy_call)

## Verify all variables have the correct number of levels
str(df_sub_pl)
# table(sv_sub_pl$ID_genotype_JGI, sv_sub_pl$year.x)

## If plant is a replicate within genotype, plants need unique IDs

# (1) Visualize data ----
n_wavelengths <- colnames(df_sub_pl) %>% str_count("w") %>% sum()
n_wavelengths

# Assess normality of the data
## One wavelength
trait_means <- aggregate(w701 ~ ID_genotype_JGI + ploidy_call , df_sub_pl, FUN = mean) # aggregate means for one trait

## genotype means
ggplot(trait_means, aes(x = w701)) +
  geom_histogram() +
  theme_bw()

## all data across years
ggplot(df_sub_pl, aes(y = w701, x = year.x)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  ylab("% Reflectance") +
  xlab("Measurement year")

# trait vs all other variables 
df_sub_pl %>%
  ggplot(aes(y = w701, x = state)) +
    geom_boxplot() +
    theme_bw() +
    theme(
      legend.position="none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )

# Reflectance across space
df_sub_pl %>%
  ggplot(aes(y = Latitude, x = Longitude, col = w701)) +
  geom_point() +
  theme_bw() +
  theme(
    # legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) + 
  scale_color_viridis_c() # really need to have unique legends and scales for each plot


# (2) Build model ----

# > Prepare dataframe ----
trait = "w701" # set trait name

## Subset dataframe to desired variables and only keep complete cases
df_sub_comp <- df_sub_pl %>% 
  select(all_of(trait), ploidy_call, age_scaled, ID_genotype_JGI, id, state, year, leaf_num) %>%
  filter(complete.cases(.)) %>%
  as.data.frame()

dim(df_sub_comp)

# > Run model ----
## id is nested within ID_genotype_JGI -- id must have unique labels

## >> Full model ----

## Specify priors
# with sufficient sample size, the prior shouldn't matter much
# MCMCglmm uses an inverse-gamma distribution for residual and random priors 
# fixed effects have a normal prior (default)
prior <- list(R = list(V = 1, nu = 0.002), # priors on Residual
              G = list(G1 = list(V = 1, nu = 0.002), # priors on Random effects
                       G2 = list(V = 1, nu = 0.002),
                       G3 = list(V = 1, nu = 0.002),
                       G4 = list(V = 1, nu = 0.002)
              ))


mod <- MCMCglmm( as.formula( paste0(trait, "~ ploidy_call + age_scaled" ) ),  # specify fixed effects
                   random = ~ ID_genotype_JGI + id + year + state , # random effects
                   ginverse = list(ID_genotype_JGI = k), # link kinship matrix to genotype variable
                   data = df_sub_comp,
                   family = "gaussian", # specify normal distribution
                   prior = prior,
                   nitt = 200000, # number of iterations
                   burnin = 10000, # iterations to drop
                   thin = 100) # iterations to keep, determined by degree of autocorrelation
                   # pl = T, # save posterior distribution of latent effects 
                # pr = T) # pr = T means the posterior distribution of the random effects were estimated (the BLUPs)

## save model
saveRDS(mod, "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/h2_mods/mod.rds")

## load in model file
# mod <- readRDS("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/h2_mods/mod.rds")


# (3) Evaluate chain convergence and autocorrelation ----

## > Assess if burn-in length is sufficient ----
## just run burn-in and see if it reaches stable
mod.burn <- MCMCglmm( as.formula( paste0(trait, "~ ploidy_call + age_scaled" ) ),
                      random = ~ ID_genotype_JGI + id + state ,
                      ginverse = list(ID_genotype_JGI = k), # link kinship matrix
                      data = df_sub_comp,
                      family = "gaussian",
                      prior = prior,
                      nitt = 10000, # number of iterations
                      burnin = 1, # iterations to drop
                      thin = 100, # iterations to keep
                      pr = T)
plot(mod.burn[["VCV"]])

## > Chain mixing ----
## Evaluate if trace plots look like caterpillars
## Also see if 'units' (the residuals) posterior is normally distributed 
plot(mod$VCV) # random effects
plot(mod$Sol) # fixed effects, also has the BLUPs

## > Autocorrelation ----
# assess whether the thinning is appropriate (want small numbers)
autocorr.diag(mod$VCV) 
autocorr.diag(mod$Sol)

## > Effective sample size ----
## Check the corrected samples size after thinning (removing autocorrelated samples) is high enough
## Want above 1,000 for all variables
effectiveSize(mod$Sol)
effectiveSize(mod[["VCV"]]) 

## > Within chain convergence ----
heidel.diag(mod$VCV)

## > Convergence between chains ----
## IN DEVELOPMENT
# https://theoreticalecology.wordpress.com/2011/12/09/mcmc-chain-analysis-and-convergence-diagnostics-with-coda-in-r/

# (4) Extract effect sizes ----
summary(mod)

## > Random effects (variances) ----
alpha = 0.05
rand_eff <- apply(mod$VCV, 1, median) %>% as.data.frame()
rand_CI <- HPDinterval(mod$VCV, prob = 1 - alpha) %>% as.data.frame()

## > Fixed effects ----
## Where CI doesn't cross zero, effect is significant
fix_eff <- apply(mod$Sol, 1, median) %>% as.data.frame()
fix_CI <- HPDinterval(mod$Sol, prob = 1 - alpha) %>% as.data.frame()