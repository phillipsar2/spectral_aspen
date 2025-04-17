### Interpreting gbs2ploidy output
### Alyssa Phillips
### 4/14/24

library(dplyr)
library(ggplot2)

dataset = "RMBL"

dir <- as.character("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/")

## Import propOut values ----
files <- list.files(path = dir, pattern = "*propOut.csv")
files <- files[grep("flow", files, invert = T)] # exclude flow cytometry test files
# files <- files[grep("SRR2", files, invert = T)] # exclude spectral genotypes

propOut <- lapply(paste0(dir, files), function(x) read.table(x, header = T, sep = ","))
names(propOut) <- files

# Separate samples that failed filtering criteria
failed_samples <- propOut[grep("fail", propOut)] %>% names()
passed_samples <- propOut[grep("fail", propOut, invert = T)] %>% names()

# Bind passed samples into dataframe
propOut_df <- bind_rows(propOut[passed_samples])
head(propOut_df)

## Plot results ----
# rows <- (dim(propOut_df)[1] / 15) / 5

facet_multiple <- function(plot = NULL, facets = NULL, ncol = 2, nrow = 2, scales = 'fixed') {
  
  if(is.null(plot)) {   # Check plot argument
    stop('Argument \"plot\" required')
  }
  
  if(is.null(facets)) {   # Check facets argument
    message('Argument \"facets\" not provided. Ploting single panel')
    return(plot)
  }
  
  if(!all(facets %in% colnames(plot$data))) {   # Ensure facets exists
    stop(paste('The facets:', facets, 'could not be found in the data'))
  }
  
  if(is.null(ncol) | is.null(nrow)) {   # Check ncol and nrow arguments
    stop('Arguments \"ncol\" and \"nrow\" required')
  }
  
  # Get info on layout
  n_panel_tot <- nrow(unique(plot$data[, facets, drop = FALSE]))
  n_layout    <- ncol*nrow
  n_pages     <- ceiling(n_panel_tot/n_layout)
  plot        <- plot + facet_wrap(facets = facets, ncol = ncol, scales = scales)
  
  # When no multiple page needed
  if(n_pages == 1) {
    return(plot)
  }
  
  # Extract ggplot2 data and title
  data   <- plot$data
  title  <- plot$labels$title
  
  # Work with the scales
  if(!scales %in% c('free', 'free_x') &&                             # if scale fixed on x
     is.numeric(eval(plot$mapping$x, data)) &&                       # and x is numeric
     length(grep('xmax', plot$scales$scales, fixed = TRUE)) == 0) {  # and x-scale hasn't been defined in ggplot2
    plot$coordinates$limits$x <- range(eval(plot$mapping$x, data))
  }
  
  if(!scales %in% c('free', 'free_y') &&                             # if scale fixed on y
     is.numeric(eval(plot$mapping$y, data)) &&                       # and y is numeric
     length(grep('ymax', plot$scales$scales, fixed = TRUE)) == 0) {  # and y-scale hasn't been defined in ggplot2
    plot$coordinates$limits$y <- range(eval(plot$mapping$y, data))
  }
  
  # Prepare the grouping
  data$groups <- findInterval(unclass(interaction(data[,facets])),
                              seq(from = 1, by = n_layout, length.out = n_pages)[-1])+1
  
  # Plot each page
  for (i in seq_along(1:n_pages)) {
    plot <- plot %+% data[data$groups == i,] + 
      ggtitle(label = bquote(atop(bold(.(title)), atop(italic(Page~.(i)~of~.(n_pages))))))
    
    # For last page call facet_layout
    if(i == n_pages) {
      plot <- facet_layout(plot = plot, facets = facets, ncol = ncol, nrow = nrow, scales = scales)
    }
    
    # Print plots
    if(!is.null(plot)) {
      print(plot)
    }
  } # End for loop
  
} 
facet_layout <- function(plot = NULL, facets = NULL, nrow = 2, ncol = 2, scales = 'fixed') {
  
  if(is.null(plot)) {   # Check plot argument
    stop('Argument \"plot\" required')
  }
  
  if(is.null(facets)) {   # Check facets argument
    message('Argument \"facets\" not provided. Ploting single panel')
    return(plot)
  }
  
  if(!all(facets %in% colnames(plot$data))) { # Ensure facets exists
    stop(paste('The facets:', facets, 'could not be found in the data'))
  }
  
  if(is.null(ncol) | is.null(nrow)) {   # Check ncol and nrow arguments
    stop('Arguments \"ncol\" and \"nrow\" required')
  }
  
  # Get info on layout
  n_panel_tot <- nrow(unique(plot$data[, facets, drop = FALSE]))
  n_layout    <- ncol*nrow
  
  if(n_panel_tot > n_layout) {   # Check layout
    stop('nrow * ncol >= n is not TRUE, use \"facet_multiple()\" instead')
  }
  
  n_missing  <- n_layout - n_panel_tot
  nrow_last  <- max(which(n_panel_tot >= seq(1, n_layout, by = ncol)))
  panel_last <- n_panel_tot - ncol*(nrow_last - 1)
  
  if(nrow_last == nrow & nrow != 1) {
    plot <- plot + facet_wrap(facets = facets, ncol = ncol, scales = scales)
    return(plot)
  }
  
  # Clean up factors
  plot$data[, 'ghosts'] <- factor(interaction(plot$data[, facets]))
  
  # Add ghost panels
  plot$data[, 'ghosts'] <- factor(x = plot$data[, 'ghosts'], 
                                  levels = c(levels(plot$data[, 'ghosts']), 
                                             paste0('ghost=', 1:n_missing)))
  
  plot <- plot + facet_wrap(facets = 'ghosts', ncol = ncol, drop = FALSE, scales = scales)
  
  drop_grobs <- function(sep = '', ...) {
  }
  
  # Convert ggplot to grob
  g <- ggplotGrob(plot)
  
  # Remove empty panels
  g$grobs[names(g$grobs) %in% drop_grobs()] <- NULL
  
  # Remove unwanted axes
  g$layout <- g$layout[!g$layout$name %in% drop_grobs(sep = '-'),]
  
  # Move bottom axis closer to panels
  g$layout[g$layout$name %in% paste0('axis_b-', (1:panel_last)+((nrow-1)*ncol)), 'b'] <- seq(1, 40, by = 5)[nrow_last+1]
  
  # Print the plot 
  grid::grid.newpage()
  grid::grid.draw(g)
  
}

# pdf(paste0(dir, "final_", dataset,"_predictions.allelicratio.allsnps.mind6.minsites10k.",Sys.Date(),"plots.pdf"))
pdf(paste0(dir, dataset,"_predictions.allelicratio.allsnps.mind6.minsites10k.",Sys.Date(),"plots.pdf"))
p <- propOut_df %>%
  dplyr::filter(quantile == '75%') %>%
  ggplot(aes(x = as.factor(allelic_ratio), y = proportion)) +
  geom_point() +
  # facet_wrap_paginate( ~ sample, nrow = 5, ncol = 4, page = 5) +
  theme_bw()+
  xlab("Allelic ratio") +
  ylab("75th quantile of the posterior distribuion of allelic proportions")
facet_multiple(plot = p, facets = 'sample', ncol = 4, nrow = 5)
dev.off()


## allelic proportions for the first nine individuals
# geno_names <- colnames(ad)
# pdf("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/nquack/gbs2ploidy_alleleratios.pdf")
# par(mfrow=c(3,3))
# for(i in 1:9){ # change number from 1 to 9
#   plot(propOut[[i]][,5], ylim=c(0,1), axes=FALSE,
#        xlab="allelic ratios",
#        # ylab="75th quantile of the posterior distribuion of allelic proportions",
#        ylab="proportions",
#        main = geno_names[i])
#   axis(1, at = 1:3,c("1:2","1:1","2:1"))
#   axis(2)
#   box()
#   segments(1:5, propOut[[i]][,1], 1:5,propOut[[i]][,5])
#   # title(main=paste("true ploidy =",dat[[3]][i]))
# }
# dev.off()

## What is the winning ploidy ----
n_genos <- dim(propOut_df)[1]/15
winner_prop <- propOut_df %>%
  group_by(sample, quantile) %>%
  filter(proportion == max(proportion)) %>%
  arrange(sample, quantile)

winner_num <- winner_prop %>%
  group_by(sample) %>%
  count(allelic_ratio, name = 'Count')
dim(winner_num)

winner_num$ploidy_call <- ifelse(winner_num$allelic_ratio == '0.5', yes = 'diploid', no = 'triploid')
head(winner_num)

## Measure accuracy compared to flow samples
# accuracy <- c()
# for (i in 1:dim(winner_num)[1]){
#   true_ploidy <- flow_meta$flow_ploidy[which(flow_meta$accession == winner_num$sample[i])]
#   accuracy[i] <- true_ploidy == winner_num$ploidy_call[i]
# }
# sum(accuracy)/length(accuracy) * 100

# propcalls <- merge(winner_num, flow_meta, by.x = 'sample', by.y = 'accession')
# dim(propcalls)so
# write.csv(propcalls, "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/flowcyt_predict.propOut.allsnps.mindp6.minsites10k.accuracy.csv",
#           row.names = F)

## Export ploidy calls ----
out_df <- dplyr::select(winner_num, sample, ploidy_call)
out_df <- df[df$sites >= 2000,] %>% 
  dplyr::select( sample, ploidy_call)
write.csv(out_df, paste0(dir, dataset, ".ploidycalls.", Sys.Date(), ".csv"), row.names = F, quote = F)

## Export names of failed samples ----
fail_names <- gsub(x = failed_samples, pattern = ".propOut.csv", replacement = "")
# fail_names <-  df$sample[df$sites < 2000]
write.csv(fail_names, paste0(dir, dataset,".failedsamples.", Sys.Date(), ".csv"), 
          row.names = F, quote = F)

###
### Evaluate calls vs quality
###

# winner_num <- filter(winner_num) %>%
#   filter(Count >= 3)
# qual <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/RMBL_qualitystats.csv")
# 
# df <- merge(x = winner_num, y = qual, by.x = "sample", by.y = "sample")
# head(df)
# 
# plot(df$sites, df$avg_alt_depth,
#      col = as.factor(df$ploidy_call),
#      xlab = "Number of sites",
#      ylab = "Mean alternate read depth")
# 
# plot(df$sites, df$dp,
#      col = as.factor(df$ploidy_call),
#      xlab = "Number of sites",
#      ylab = "Mean read depth")
# 
# plot(df$avg_alt_depth, df$avg_ref_depth,
#      col = as.factor(df$ploidy_call),
#      xlab = "Mean alternate read depth",
#      ylab = "Mean reference read depth")

# plot(df$avg_alt_depth, df$avg_ref_depth,
#      # pch = as.factor(df$ploidy_call),
#      col = ifelse(df$sites < 2000,'red','green'),
#      xlab = "Mean alternate read depth",
#      ylab = "Mean reference read depth")
# 
# df$sample[df$sites < 2500]
