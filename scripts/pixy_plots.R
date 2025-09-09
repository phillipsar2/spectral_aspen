### Title: Pixy output
### Author: Alyssa Phillips
### Date: 9/8/25

library(ggplot2)
library(tidyverse)

dir <- "/global/scratch/users/arphillips/spectral_aspen/data/pixy/"

# Provide path to input. Can be pi or Dxy.
# NOTE: this is the only line you should have to edit to run this code:
inp <- read.table(paste0(dir, "inwards.Chr01.w10000_pi.txt"), sep = "\t", header = T)
window = 10000

# Find the chromosome names and order them: first numerical order, then any non-numerical chromosomes
#   e.g., chr1, chr2, chr22, chrX
chroms <- unique(inp$chromosome)
chrOrder <- sort(chroms)
inp$chrOrder <- factor(inp$chromosome,levels=chrOrder)

# Plot pi for each population found in the input file
# Saves a copy of each plot in the working directory


# Load converting wide to long ----
pixy_to_long <- function(pixy_files){ # input file path
  pixy_df <- list()
  for(i in 1:length(pixy_files)){
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop)
      pixy_df[[i]] <- df
  }
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
}

pixy_files <- paste0(dir, "inwards.Chr01.w10000_pi.txt")
pixy_df <- pixy_to_long(pixy_files)
head(pixy_df)

# Plot along chromosome ----
pixy_df %>%
  filter(chromosome == "Chr01") %>%
  filter(statistic %in% c("avg_pi")) %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/1000000) %>%
  ggplot(aes(x = chr_position, y = value, color = statistic))+
  geom_line(size = 0.25)+
  # facet_grid(statistic ~ .,
  #            scales = "free_y", switch = "x", space = "free_x",
  #            labeller = labeller(statistic = pixy_labeller,
  #                                value = label_value))+
  xlab("Position on Chromosome (Mb)")+
  ylab("Statistic Value")+
  theme_bw()+
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_color_brewer(palette = "Set1")

# Calculate average thetas ----
# (window 1 count_diffs + window 2 count_diffs) / (window 1 comparisons + window 2 comparisons)
# sum(inp$count_diffs) / sum(inp$count_comparisons)

p = unique(inp$pop)[1]
for (p in unique(inp$pop)){
  sub <- inp[inp$pop == p,]
  # sub <- filter(sub, no_sites / window >= 0.1) # remove windows with < 10% of sites
  sub <- filter(sub, no_sites != 0)
  avg_pi <- sum(sub$count_diffs) / sum(sub$count_comparisons)
}

# avg_pi <- tapply(inp, inp$pop, function(x){sum(inp$count_diffs) / sum(inp$count_comparisons)})
# length(avg_pi)

