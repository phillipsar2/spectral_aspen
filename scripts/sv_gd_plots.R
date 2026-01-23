### Title: SV ~ GD plots
### Author: Alyssa Phillips
### Date: 11/15/25

library(dplyr)
library(ggplot2)
library(viridis)
library(broom)

# Load data ----
sampling <-  c("southnorth", "random", "outwards", "northsouth", "inwards")
window <- c("25000", "50000")

sv <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/spectral_variation/ms_arr_spectra.csv")
outdir = '/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/sv_by_gd'

s = 'outwards'
w = 25000

for (s in sampling){
  for(w in window){
    
  gd <- read.csv(paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/",s,".",w,".div.estimates.csv"))

  sv_sub <- filter(sv, grid == s)
  # sv_sub <- sv[sv$grid == "outwards",]
  # 370*3*7 # 3 months, 7 years

  # sv_sub$id <- as.factor(sv_sub$id) # zero based
  sv_sub$index <- as.numeric(sv_sub$id) + 1

  df <- merge(sv_sub, gd, by.y = 'X', by.x = 'index')
  df$month <- factor(df$month, levels = c("jul", "aug", "sep"))
  df$year <- as.factor(df$year)

  ## A and Asq headers are swapped in all the genetic div files
  pdf(file = paste0(outdir, "/M_by_Akm2.",s,".",w,".pdf"), width = 5, height = 4)
  ggplot(df, aes(x = area_km2, y = M)) +
    geom_point() +
    theme_bw() +
    geom_smooth(method = glm, formula = y ~ log(x)) +
    ylab("Area (Km2)") +
    xlab("M")
  dev.off()
  
  dplyr::filter(df, year == '2018') %>%
  ggplot( aes(x = area_km2, y = specvar, color = month)) +
  geom_point() +
    theme_bw() +
    geom_smooth(method = glm, formula = y ~ log(x)) +
    ylab("Area (Km2)") +
    xlab("Spectral variation") +
    scale_color_discrete(palette = hex, labels = c("July", "August", "September")) 
  
  ggplot(df, aes(x = area_km2, y = M)) +
    geom_point() +
    theme_bw() +
    geom_smooth(method = glm, formula = y ~ log(x)) +
    ylab("Area (Km2)") +
    xlab("M")
  
  hex <- c("#117733", "#88CCEE", "#DDCC77")
  pdf(file = paste0(outdir, "/logM_by_logSV.",s,".",w,".pdf"), width = 8, height = 8)
  ggplot(df, aes(x = log(M), y = log(specvar), col = month)) +
    geom_point(alpha = 0.5) +
    theme_bw() +  geom_smooth(method = 'lm') +
    ylab("log(SV)") +
    xlab("log(M)") +
    facet_wrap(~year, scales = "free_y", ncol = 4) +
    scale_color_discrete(palette = hex)+
    scale_color_discrete(palette = hex, labels = c("July", "August", "September")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  dev.off()
  
  filter(df, year == '2018') %>%
  ggplot(aes(x = log(M), y = log(specvar), col = month)) +
    geom_point(alpha = 0.5) +
    theme_bw() +  geom_smooth(method = 'lm') +
    ylab("log(SV)") +
    xlab("log(M)") +
    scale_color_discrete(palette = hex)+
    scale_color_discrete(palette = hex, labels = c("July", "August", "September")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  pdf(file = paste0(outdir, "/M_by_AV.",s,".",w,".pdf"), width = 8, height = 8)
  ggplot(df, aes(x = M, y = specvar, col = month)) +
    geom_point(alpha = 0.5) +
    theme_bw() +  
    geom_smooth(method = 'lm') +
    ylab("SV") +
    xlab("M (10^3)") +
    facet_wrap(~year, scales = "free_y") +
    scale_color_discrete(palette = hex, labels = c("July", "August", "September")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_x_continuous(labels=function(x){x/1000})
  dev.off()

  # cor(df$M, df$specvar)
  
  cor_df <- df %>%
    group_by(year, month) %>%
    summarize(cor = cor(specvar, M))
  
  pdf(file = paste0(outdir, "/cor_by_month.",s,".",w,".pdf"), width = 5, height = 4)
  ggplot(cor_df, aes(x = month, y = cor, group = year, color = year)) +
    geom_point() +
    geom_line() +
    scale_color_viridis() +
    ylab("Pearson correlation coefficient") + xlab("Month") + 
    labs(color = "Year") +
    theme_bw() +
    scale_x_discrete(labels = c("July", "August", "September") )
  dev.off()
  
  model_df <- df %>%
    group_by(year, month) %>%
    do(model = lm(specvar ~ M, data = .)) 
  b_df <- model_df %>% 
    summarise(intercept = coef(model)[1], slope = coef(model)[2]) %>%
    select(slope)
  
  output_df <- cbind(cor_df, b_df)
  write.csv(output_df, paste0(outdir,"/cor_beta_df.",s,".",w,".csv" ))
  
  pdf(file = paste0(outdir, "/beta_by_month.",s,".",w,".pdf"), width = 5, height = 4)
  ggplot(output_df, aes(x = month, y = slope, group = year, color = year)) +
    geom_point() +
    geom_line() +
    scale_color_viridis() +
    ylab("SV/M") + xlab("Month") + 
    labs(color = "Year") +
    theme_bw() +
    scale_x_discrete(labels = c("July", "August", "September") )
  dev.off()
  }
}

# Mantel test ----
library(ade4)

## Spectral distance ----
# sd <- 
outdir = '/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/sv_by_gd'

s = 'outwards'
w = 25000

## Genetic distance ---
gd <- read.csv(paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/",s,".",w,".div.estimates.csv"))

mantel.rtest(station.dists, ozone.dists, nrepet = 9999)