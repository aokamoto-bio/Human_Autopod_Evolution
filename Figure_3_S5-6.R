#create figure showing overlaps with autopod ATAC data
#Alexander Okamoto
#June 6, 2025
library(magick)

#note that images read in are created in the autopod_overlap_with_other_skeleton

#create function to all for modification of all file reading simultaneously
read_image_gg <- function(file){
  gg_image <- image_read(file) %>% image_crop(geometry = "290x270+40+30") %>%
    image_ggplot() 
  return(gg_image)
}

make_overlap_figure <- function(skeletal_heatmap_list, filename){
  
  #load plots 
  plot1 <- read_image_gg(skeletal_heatmap_list[1])
  plot2 <- read_image_gg(skeletal_heatmap_list[2]) 
  plot3 <- read_image_gg(skeletal_heatmap_list[3]) 
  plot4 <- read_image_gg(skeletal_heatmap_list[4]) 
  plot5 <- read_image_gg(skeletal_heatmap_list[5]) 
  plot6 <- read_image_gg(skeletal_heatmap_list[6]) 
  plot7 <- read_image_gg(skeletal_heatmap_list[7]) 
  plot8 <- read_image_gg(skeletal_heatmap_list[8]) 
  plot9 <- read_image_gg(skeletal_heatmap_list[9]) 
  plot10 <- read_image_gg(skeletal_heatmap_list[10]) 
  plot11 <- read_image_gg(skeletal_heatmap_list[11]) 
  plot12 <- read_image_gg(skeletal_heatmap_list[12]) 
  plot13 <- read_image_gg(skeletal_heatmap_list[13]) 
  plot14 <- read_image_gg(skeletal_heatmap_list[14]) 
  plot15 <- read_image_gg(skeletal_heatmap_list[15])
  plot16 <- read_image_gg(skeletal_heatmap_list[16]) 
  plot17 <- read_image_gg(skeletal_heatmap_list[17]) 
  plot18 <- read_image_gg(skeletal_heatmap_list[18])


#make final plot
#make final plot
skel_h <- 0.16
skel_w <- 0.32
x_pos <- rep(c(0, skel_w, skel_w*2), 6) + 0.04
y_pos <- c(rep((0.95 - skel_h), 3), 
           rep((0.95 - skel_h*2), 3), 
           rep((0.95 - skel_h*3), 3), 
           rep((0.95 - skel_h*4), 3), 
           rep((0.95 - skel_h*5), 3), 
           rep((0.95 - skel_h*6), 3)) 

overlap_plot <- ggdraw() +
  draw_plot(plot1, x = x_pos[1], y = y_pos[1], width = skel_w, height = skel_h) + 
  draw_plot(plot2, x = x_pos[2], y = y_pos[2], width = skel_w, height = skel_h) + 
  draw_plot(plot3, x = x_pos[3], y = y_pos[3], width = skel_w, height = skel_h) + 
  draw_plot(plot4, x = x_pos[4], y = y_pos[4], width = skel_w, height = skel_h) + 
  draw_plot(plot5, x = x_pos[5], y = y_pos[5], width = skel_w, height = skel_h) + 
  draw_plot(plot6, x = x_pos[6], y = y_pos[6], width = skel_w, height = skel_h) + 
  draw_plot(plot7, x = x_pos[7], y = y_pos[7], width = skel_w, height = skel_h) + 
  draw_plot(plot8, x = x_pos[8], y = y_pos[8], width = skel_w, height = skel_h) + 
  draw_plot(plot9, x = x_pos[9], y = y_pos[9], width = skel_w, height = skel_h) + 
  draw_plot(plot10, x = x_pos[10], y = y_pos[10], width = skel_w, height = skel_h) + 
  draw_plot(plot11, x = x_pos[11], y = y_pos[11], width = skel_w, height = skel_h) + 
  draw_plot(plot12, x = x_pos[12], y = y_pos[12], width = skel_w, height = skel_h) + 
  draw_plot(plot13, x = x_pos[13], y = y_pos[13], width = skel_w, height = skel_h) + 
  draw_plot(plot14, x = x_pos[14], y = y_pos[14], width = skel_w, height = skel_h) + 
  draw_plot(plot15, x = x_pos[15], y = y_pos[15], width = skel_w, height = skel_h) + 
  draw_plot(plot16, x = x_pos[16], y = y_pos[16], width = skel_w, height = skel_h) + 
  draw_plot(plot17, x = x_pos[17], y = y_pos[17], width = skel_w, height = skel_h) + 
  draw_plot(plot18, x = x_pos[18], y = y_pos[18], width = skel_w, height = skel_h) + 
  draw_plot_label(label = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r"), 
      size = 12,
      x = x_pos + 0.01, 
      y = (y_pos + skel_h)) + 
  draw_plot_label(label = c("Total", "HARs", "HAQERs", "hCONDELs", "Inversions", "SDRs"), 
                  size = 12,
                  x = rep(0, 6), 
                  y = 0.88 - seq(from = 0, to = skel_h*5, by = skel_h), 
                  angle = 90, 
                  hjust = 0.5) + 
  draw_plot_label(label = c("Brain-filtered", "Autopod-specific", "Tissue-specific"), 
                  size = 12,
                  x = c(0.20, 0.52, 0.84), 
                  y = rep(0.975, 3), 
                  hjust = 0.5) + 
  bgcolor("white") + geom_rect( 
                            aes(xmin = 0, xmax = 1, 
                                ymin = 0.8, ymax = 0.95), 
                            color = NA, fill = "gray", size = 1, alpha = 0.25) + 
  geom_segment(aes(x = 0.18, y = 0.985, xend = 0.85, yend = 0.985),
                arrow = arrow(length = unit(2, "mm"), type = "closed"),
                size = 2,
                color = "gray"
  ) + 
  draw_plot_label(label = "specificity", 
                  size = 10,
                  x = 0.10, 
                  y = 0.998, 
                  hjust = 0.5,
                  color = "darkgray")


overlap_plot  

#save the plot
ggsave(plot = overlap_plot + panel_border(color = "black", size = 1), 
       filename = filename, 
       device = "png", dpi = 300, height = 220, width = 140, units = "mm")

}

#generate the plots
#main figure
combined_list <- paste("~/Desktop/Autopod_Skeleton_Figures/", rep(c("all", "hars", "haqers", "hcondels", "inversions", "sdrs"), each = 3), "_", c("all", "autopod", "specific"), "_combined.png", sep = "")
make_overlap_figure(skeletal_heatmap_list = combined_list, 
                    filename = "~/Dropbox/Autopod Paper/Autopod_Paper_Figures/figure_3_overlaps.png")

#early supplementary figure
early_list <- paste("~/Desktop/Autopod_Skeleton_Figures/", rep(c("all", "hars", "haqers", "hcondels", "inversions", "sdrs"), each = 3), "_", c("all", "autopod", "specific"), "_early.png", sep = "")
make_overlap_figure(skeletal_heatmap_list = early_list, 
                    filename = "~/Dropbox/Autopod Paper/Autopod_Paper_Figures/figure_S4_early_overlaps.png")

#late supplementary figure
late_list <- paste("~/Desktop/Autopod_Skeleton_Figures/", rep(c("all", "hars", "haqers", "hcondels", "inversions", "sdrs"), each = 3), "_", c("all", "autopod", "specific"), "_late.png", sep = "")
make_overlap_figure(skeletal_heatmap_list = late_list, 
                    filename = "~/Dropbox/Autopod Paper/Autopod_Paper_Figures/figure_S5_late_overlaps.png")



