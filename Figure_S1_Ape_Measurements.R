#create supplementary figure showing human and mouse ATAC-seq overlaps
#Alexander Okamoto
#February 12, 2025

#load packages

library(tidyverse) #for data manipulation
library(ggpubr) #for plotting options
library(grid) # for plotting options
library(magick) # for adding images to ggplots
library(cowplot) #for manipulating graphs

#load custom functions for analysis
source("~/Desktop/Capellini_Lab/Weekly_Coding/RNA_Analysis_Functions.R") 

#load in data from Campbell Rolian
Hominoid_Limb_Stats <- read_csv("~/Desktop/Capellini_Lab/Rolian_Data/Hominoid_Limb_Stats_R_Clean_Final.csv", col_names = T)
#clearn data and convert to long format
Hominoid_Limb_Stats_clean <- Hominoid_Limb_Stats %>% pivot_longer(!Statistic, names_to = "Anatomical_Region", values_to = "Measurement") %>% pivot_wider(names_from = Statistic, values_from = Measurement)

#add information about which axis (length or width) the measurement pertains to for filtering
Hominoid_Limb_Stats_clean$Axis <- c(rep("Length", 6), rep(c("Length", "Width"), 28))
#add information on the specific tissue
Hominoid_Limb_Stats_clean$Tissue <- NA
Hominoid_Limb_Stats_clean$Tissue[1:6] <- Hominoid_Limb_Stats_clean$Anatomical_Region[1:6]
Hominoid_Limb_Stats_clean$Tissue[7:62] <- Hominoid_Limb_Stats_clean$Anatomical_Region[7:62] %>% gsub(pattern = "L", replacement = "") %>% gsub(pattern = "W",replacement = "") %>% gsub(pattern = "MC",replacement = "HM") %>% gsub(pattern = "MT",replacement = "FM") %>% gsub(pattern = "FFPP",replacement = "FPP") %>% gsub(pattern = "FPP",replacement = "FP") %>% gsub(pattern = "PP",replacement = "HP") %>% gsub(pattern = "FMP",replacement = "FP") %>% gsub(pattern = "MP",replacement = "HP") 
Hominoid_Limb_Stats_clean$Phal_Element <- NA
Hominoid_Limb_Stats_clean$Phal_Element[str_which(pattern = "MP", string =  Hominoid_Limb_Stats_clean$Anatomical_Region)] <- "Medial"
Hominoid_Limb_Stats_clean$Phal_Element[str_which(pattern = "PP", 
                                                 string = Hominoid_Limb_Stats_clean$Anatomical_Region)] <- "Proximal"

#combine proximal and medial phalangeal measurements
Hominoid_Limb_Stats_Reduced <- Hominoid_Limb_Stats_clean %>%  group_by(Tissue, Axis) %>% 
  summarize(
    `Homo_sample_size_(n)` = median(`Homo_sample_size_(n)`), 
    `Homo_mean_(mm)` = sum(`Homo_mean_(mm)`),
    `Homo_SD_(mm)` = sqrt(sum(`Homo_Variance_(mm)`)), 
    `Pan_sample_size_(n)` = median(`Pan_sample_size_(n)`), 
    `Pan_mean_(mm)` = sum(`Pan_mean_(mm)`),
    `Pan_SD_(mm)` = sqrt(sum(`Pan_Variance_(mm)`)),
    `Gorilla_sample_size_(n)` = median(`Gorilla_sample_size_(n)`), 
    `Gorilla_mean_(mm)` = sum(`Gorilla_mean_(mm)`),
    `Gorilla_SD_(mm)` = sqrt(sum(`Gorilla_Variance_(mm)`))
  )

#calculate percent change in chimpanzee
Hominoid_Limb_Stats_Reduced$Homo_Pan_Perc_Change <- (Hominoid_Limb_Stats_Reduced$`Homo_mean_(mm)`/Hominoid_Limb_Stats_Reduced$`Pan_mean_(mm)`)*100

#calculate percent change in gorilla
Hominoid_Limb_Stats_Reduced$Homo_Gorilla_Perc_Change <- (Hominoid_Limb_Stats_Reduced$`Homo_mean_(mm)`/Hominoid_Limb_Stats_Reduced$`Gorilla_mean_(mm)`)*100

#now calculate absolute value of change
Hominoid_Limb_Stats_Reduced$Homo_Pan_Dist_Change <- abs(Hominoid_Limb_Stats_Reduced$Homo_Pan_Perc_Change - 100)

#add in gorilla
#calculate percent change in chimpanzee
Hominoid_Limb_Stats_Reduced$Pan_Gorilla_Perc_Change <- (Hominoid_Limb_Stats_Reduced$`Pan_mean_(mm)`/Hominoid_Limb_Stats_Reduced$`Gorilla_mean_(mm)`)*100


#get chimp plots 
Hominoid_Limb_Stats_Reduced %>% 
  dplyr::filter(Axis == "Length", Tissue %in% c("FM1", "FM2", "FM3", "FM4", "FM5", "FP1", "FP2", "FP3", "FP4", "FP5", "HM1", "HM2", "HM3", "HM4", "HM5", "HP1", "HP2", "HP3", "HP4", "HP5")) %>% 
  dplyr::select(Tissue, Homo_Pan_Perc_Change) %>% 
  dplyr::rename(tissue = Tissue, value = Homo_Pan_Perc_Change) %>% 
  drop_na() %>% 
  plot_autopod_expression(file_name = "Homo_Pan_Perc_Change_Length2", 
                          opt_text = "Homo_Pan%-L", center = T, center_color_value = 100, min_color_value = 49, max_color_value = 121)

Hominoid_Limb_Stats_Reduced %>% 
  dplyr::filter(Axis == "Width", Tissue %in% c("FM1", "FM2", "FM3", "FM4", "FM5", "FP1", "FP2", "FP3", "FP4", "FP5", "HM1", "HM2", "HM3", "HM4", "HM5", "HP1", "HP2", "HP3", "HP4", "HP5")) %>% 
  dplyr::select(Tissue, Homo_Pan_Perc_Change) %>% dplyr::rename(tissue = Tissue, value = Homo_Pan_Perc_Change) %>% 
  drop_na() %>%
  plot_autopod_expression(file_name = "Homo_Pan_Perc_Change_Width2", opt_text = "Homo_Pan%-W", center = T, center_color_value = 100, min_color_value = 65, max_color_value = 156)

#get gorilla plots
Hominoid_Limb_Stats_Reduced %>% 
  dplyr::filter(Axis == "Length", Tissue %in% c("FM1", "FM2", "FM3", "FM4", "FM5", "FP1", "FP2", "FP3", "FP4", "FP5", "HM1", "HM2", "HM3", "HM4", "HM5", "HP1", "HP2", "HP3", "HP4", "HP5")) %>% 
  dplyr::select(Tissue, Homo_Gorilla_Perc_Change) %>% 
  dplyr::rename(tissue = Tissue, value = Homo_Gorilla_Perc_Change) %>% 
  drop_na() %>% 
  plot_autopod_expression(file_name = "Homo_Gorilla_Perc_Change_Length2", opt_text = "Homo_Gorilla%-L", center = T, center_color_value = 100, min_color_value = 49, max_color_value = 121)

Hominoid_Limb_Stats_Reduced %>% 
  dplyr::filter(Axis == "Width", Tissue %in% c("FM1", "FM2", "FM3", "FM4", "FM5", "FP1", "FP2", "FP3", "FP4", "FP5", "HM1", "HM2", "HM3", "HM4", "HM5", "HP1", "HP2", "HP3", "HP4", "HP5")) %>% 
  dplyr::select(Tissue, Homo_Gorilla_Perc_Change) %>% 
  dplyr::rename(tissue = Tissue, value = Homo_Gorilla_Perc_Change) %>% 
  drop_na() %>% 
  plot_autopod_expression(file_name = "Homo_Gorilla_Perc_Change_Width2", opt_text = "Homo_Gorilla%-W", center = T, center_color_value = 100, min_color_value = 65, max_color_value = 156)

Hominoid_Limb_Stats_Reduced %>% 
  dplyr::filter(Axis == "Length", Tissue %in% c("FM1", "FM2", "FM3", "FM4", "FM5", "FP1", "FP2", "FP3", "FP4", "FP5", "HM1", "HM2", "HM3", "HM4", "HM5", "HP1", "HP2", "HP3", "HP4", "HP5")) %>% 
  dplyr::select(Tissue, Pan_Gorilla_Perc_Change) %>% 
  dplyr::rename(tissue = Tissue, value = Pan_Gorilla_Perc_Change) %>% 
  drop_na() %>% 
  plot_autopod_expression(file_name = "Pan_Gorilla_Perc_Change_Length2", opt_text = "Pan_Gorilla%-L", center = T, center_color_value = 100, min_color_value = 49, max_color_value = 121)

Hominoid_Limb_Stats_Reduced %>% 
  dplyr::filter(Axis == "Width", Tissue %in% c("FM1", "FM2", "FM3", "FM4", "FM5", "FP1", "FP2", "FP3", "FP4", "FP5", "HM1", "HM2", "HM3", "HM4", "HM5", "HP1", "HP2", "HP3", "HP4", "HP5")) %>% 
  dplyr::select(Tissue, Pan_Gorilla_Perc_Change) %>% 
  dplyr::rename(tissue = Tissue, value = Pan_Gorilla_Perc_Change) %>% 
  drop_na() %>%  
  plot_autopod_expression(file_name = "Pan_Gorilla_Perc_Change_Width2", opt_text = "Pan_Gorilla%-W", center = T, center_color_value = 100, min_color_value = 65, max_color_value = 156)

#load plots
pan_L <- image_read("~/Desktop/Autopod_Skeleton_Figures/Homo_Pan_Perc_Change_Length2.png") %>%
  image_ggplot()
pan_W <- image_read("~/Desktop/Autopod_Skeleton_Figures/Homo_Pan_Perc_Change_Width2.png") %>%
  image_ggplot()
gorilla_L <- image_read("~/Desktop/Autopod_Skeleton_Figures/Homo_Gorilla_Perc_Change_Length2.png") %>%
  image_ggplot()
gorilla_W <- image_read("~/Desktop/Autopod_Skeleton_Figures/Homo_Gorilla_Perc_Change_Width2.png") %>%
  image_ggplot()
ape_L <- image_read("~/Desktop/Autopod_Skeleton_Figures/Pan_Gorilla_Perc_Change_Length2.png") %>%
  image_ggplot()
ape_W <- image_read("~/Desktop/Autopod_Skeleton_Figures/Pan_Gorilla_Perc_Change_Width2.png") %>%
  image_ggplot()


#make final plot

ape_measurements_plot <- ggdraw() +
  draw_plot(pan_L, x = 0.0, y = 0.67, width = 0.5, height = 0.33) +
  draw_plot(pan_W, x = 0.5, y = 0.67, width = 0.5, height = 0.33) +
  draw_plot(gorilla_L, x = 0, y = .34, width = 0.5, height = .33) +
  draw_plot(gorilla_W, x = 0.5, y = 0.34, width = 0.5, height = 0.33) +
  draw_plot(ape_L, x = 0, y = .01, width = 0.5, height = .33) +
  draw_plot(ape_W, x = 0.5, y = 0.01, width = 0.5, height = 0.33) +
  draw_plot_label(label = c("a", "b", "c", "d", "e", "f"), 
                  size = 12,
                  x = c(0, 0.5, 0, 0.5, 0, 0.5), 
                  y = c(1, 1, 0.67, 0.67, 0.33, 0.33)) + 
  bgcolor("white")

ape_measurements_plot

ggsave(plot = ape_measurements_plot + panel_border(color = "black", size = 1), 
       filename = "~/Dropbox/Autopod Paper/Autopod_Paper_Figures/S1_Ape_Measurements.png", 
       device = "png", dpi = 300, height = 6, width = 4, units = "in")





