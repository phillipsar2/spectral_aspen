### Title: Pixy output
### Author: Alyssa Phillips
### Date: 9/8/25

library(ggplot2)
library(tidyverse)
library(raster)
library(sars)

# Load pixy data ----
dir <- "/global/scratch/users/arphillips/spectral_aspen/data/pixy/09292025/"
window = 25000 # 25000 & 50000
sampling = "random" # southnorth, random, outwards, northsouth, inwards

# Calculate pi ----
## Load pi data ----
pi_files <- list.files(path = dir, pattern = paste0(sampling, ".*w", window, "_pi.txt"))
pi_list <- lapply(paste0(dir, pi_files), function(x){ read.table(x, sep = "\t", header = T) })
pi_df <- bind_rows(pi_list)
pi_df$pop_num <- gsub(x = pi_df$pop, pattern = "pop", replacement = "")
dim(pi_df)
str(pi_df)

pi_df <- arrange(pi_df, as.numeric(pop_num))

## Sort by chromosome number
chroms <- unique(pi_df$chromosome)
chrOrder <- sort(chroms)
pi_df$chrOrder <- factor(pi_df$chromosome, levels = chrOrder)

## Calculate average pi ----
hist(pi_df$no_sites/window * 100,
     main = "percent of sites in each window")

# (window 1 count_diffs + window 2 count_diffs) / (window 1 comparisons + window 2 comparisons)
# sum(inp$count_diffs) / sum(inp$count_comparisons)
avg_theta <- data.frame("pop" = unique(pi_df$pop_num), # create empty df
                     "avg_pi" = NA,
                     "avg_watterson" = NA,
                     "M" = NA)

for (p in unique(pi_df$pop_num)){
  sub <- pi_df[pi_df$pop_num == p,]
  sub <- filter(sub, no_sites != 0) # drop windows with zero sequenced sites
  avg_theta$avg_pi[avg_theta$pop == p] <- sum(sub$count_diffs) / sum(sub$count_comparisons)
}
head(avg_theta)

hist(avg_theta$avg_pi)
mean(avg_theta$avg_pi)

# Calculate Watterson's theta ----
## Load theta w data ----
w_files <- list.files(path = dir, pattern = paste0(sampling, ".*w", window, "_watterson_theta.txt"))
w_list <- lapply(paste0(dir, w_files), function(x){ read.table(x, sep = "\t", header = T) })
w_df <- bind_rows(w_list)
w_df$pop_num <- gsub(x = w_df$pop, pattern = "pop", replacement = "")
dim(w_df)
str(w_df)

str(w_df)

w_df <- arrange(w_df, as.numeric(pop_num))

## Calc avg theta 2 ----
# averaging the theta column 
for (w in unique(w_df$pop_num)){
  sub <- w_df[w_df$pop_num == w,]
  sub <- filter(sub, no_sites != 0) # drop windows with zero sequenced sites
  avg_theta$avg_watterson[avg_theta$pop == w] <- mean(sub$avg_watterson_theta)
}
head(avg_theta$avg_watterson)

hist(avg_theta$avg_watterson)

mean(avg_theta$avg_watterson)

# Count mutations (S) ----
## Calc total number of segregating sites ----
## no_var_sites = number of variable sites
for (m in unique(w_df$pop_num)){
  sub <- w_df[w_df$pop_num == m,]
  sub <- filter(sub, no_sites != 0) # drop windows with zero sequenced sites
  avg_theta$M[avg_theta$pop == m] <- sum(sub$no_var_sites)
}
head(avg_theta$M)

hist(avg_theta$M)
mean(avg_theta$M)

# Load Area data ----
mardir <- "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/"

## Load mar object
load(paste0(mardir, "mardflist_RMBL_aspen.rda"))

## Extract area
area <- mardflist[[sampling]][,c("A","Asq")]

## Extract sample size
n <- mardflist[[sampling]][,"N"]

pdf(paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/figures/", sampling, ".nvsA.pdf"))
plot(area$A, n) # n increases with A
dev.off()

# Plot ----
area_df <- data.frame(area, avg_theta, n)
area_df$n <- as.numeric(area_df$n)
area_df$pop <- as.factor(area_df$pop)
str(area_df)

# area_df <- filter(area_df, n > 5) # filter to more than 5 genos

ggplot(area_df, aes(x = Asq, y = M)) +
  # geom_point(aes(size = n), alpha = 0.5) +
  geom_point() +
  theme_bw() +
  geom_smooth() +
  ylab("M") +
  xlab("Area") +
  ggtitle(sampling)

ggplot(area_df, aes(x = log(Asq), y = log(M))) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  geom_smooth() +
  ylab("log M") +
  xlab("log Area") +
  ggtitle(sampling)

ggplot(area_df, aes(x = Asq, y = avg_watterson)) +
  # geom_point(aes( size = n), alpha = 0.5) +
  geom_point() +
  theme_bw() +
  geom_smooth() +
  ylab("Watterson's theta") +
  xlab("Area") +
  ggtitle(sampling)

ggplot(area_df, aes(x = Asq, y = avg_pi)) +
  # geom_point(aes( size = n), alpha = 0.5) +
  geom_point() +
  theme_bw() +
  geom_smooth() +
  ylab("Pi") +
  xlab("Area") +
  ggtitle(sampling)

# Write data frame ----
colnames(area_df) <- c( "A", "Asq", "pop","thetapi", "thetaw","M","n")
write.csv(area_df, 
          paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/", sampling, ".", window, ".div.estimates.csv"))

# Estimate MAR ----
# Fit power law equation
# area <- 1:10
# mutations <- c(50, 75, 100, 125, 150, 175, 200, 225, 250, 275)
# data <- data.frame(A = area, M = mutations)
# mar <- MARcalc(data, Mtype = "M", Atype = "A")

mar <- mar::MARcalc(area_df, "M", "A")
war <- mar::MARcalc(area_df, "thetaw", "A")
par <- mar::MARcalc(area_df, "thetapi", "A")

write.csv(summary(mar)$Parameters, paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/", sampling, ".", window, ".model.summary.csv"))

pdf(paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/figures/", sampling , ".", window, ".mar.pdf"),
    width = 8, height =3)
par(mfrow = c(1,3))
plot(mar, xlab = "Area", ylab = "M")
plot(war, xlab = "Area", ylab = "Watterson's theta")
plot(par, xlab = "Area", ylab = "Nucleotide diversity")
dev.off()
