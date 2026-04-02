#create supplementary figure showing total human autopod ATAC-seq WGCNA results 
#Alexander Okamoto
#April 2, 2026

#load packages
library(tidyverse)
library(cowplot)
library(ggpubr) #for plotting options
library(grid) # for plotting options
library(magick) # for adding images to ggplots
library(AnnotationDbi)
library(org.Hs.eg.db)
library(rGREAT)

#modules

files <- paste("~/Desktop/Autopod_Skeleton_Figures/human_total_peak_cluster_", 1:20, ".png", sep = "")
# 2. Read them all into one magick object
images <- image_read(files)

# Converts images 1:38 to ggplots and tiles them into 4 rows
human_total_atac_plot <- plot_grid(plotlist = lapply(images[1:20], image_ggplot), 
                                  ncol = 10, 
                                  labels = 1:20)

ggsave(plot = human_total_atac_plot + panel_border(color = "black", size = 1) + 
         bgcolor("white"), 
       filename = "~/Dropbox/Autopod Paper/Autopod_Paper_Figures/Fig_S6_Human_wgcna_atac_total.png", 
       device = "png", dpi = 300, height = 90, width = 180, units = "mm")

