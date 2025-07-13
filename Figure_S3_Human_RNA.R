#create supplementary figure showing humana autopod RNA-seq results 
#Alexander Okamoto
#February 14, 2025

#load packages
library(tidyverse)
library(cowplot)

#load cluster plots

plt1 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_early_RNA_cluster_1.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt2 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_early_RNA_cluster_2.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt3 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_late_RNA_cluster_1.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt4 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_late_RNA_cluster_2.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt5 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_late_RNA_cluster_3.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt6 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_late_RNA_cluster_4.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 


#make final plot
#make final plot
skel_h <- 0.33
skel_w <- 0.5

human_RNA_plot <- ggdraw() +
  draw_plot(plt1, x = 0.0, y = 0.66, width = skel_w, height = skel_h) +
  draw_plot(plt2, x = 0.5, y = 0.66, width = skel_w, height = skel_h) +
  draw_plot(plt3, x = 0.0, y = 0.33, width = skel_w, height = skel_h) +
  draw_plot(plt4, x = 0.5, y = 0.33, width = skel_w, height = skel_h) +
  draw_plot(plt5, x = 0.0, y = 0, width = skel_w, height = skel_h) +
  draw_plot(plt6, x = 0.5, y = 0, width = skel_w, height = skel_h) +
  draw_plot_label(label = c("a1", "a2", "b1", "b2", "b3", "b4"), 
                  size = 12,
                  x = c(0, 0.5, 0, 0.5, 0, 0.5), 
                  y = c(0.99, 0.99, 0.66, 0.66, 0.33, 0.33)) + 
  bgcolor("white")

human_RNA_plot 

ggsave(plot = human_RNA_plot + panel_border(color = "black", size = 1), 
       filename = "~/Dropbox/Autopod Paper/Autopod_Paper_Figures/S3_Human_RNA.png", 
       device = "png", dpi = 300, height = 180, width = 120, units = "mm")

