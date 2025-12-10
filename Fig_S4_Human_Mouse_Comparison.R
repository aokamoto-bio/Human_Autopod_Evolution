#create supplementary figure showing human and mouse ATAC-seq overlaps
#Alexander Okamoto
#December 10, 2025

#load packages
library(tidyverse) #for data manipulation
library(ggpubr) #for plotting options
library(VennDiagram) #for plotting Venn Diagrams
library(grid) # for plotting options
library(scales) # for plotting options
library(magick) # for plotting options

#overall overlaps

#count human peaks 
human_total_n <- system("wc -l ~/Desktop/Capellini_Lab/Human_ATAC/Human_Autopod_IDR_peaks/Brain_Filtered/human_all_merge_IDR_0.05_brain_filtered.bed", intern=T) %>% readr::parse_number()

#if needed to recreate mouse data
#mouse_IDR_beds <- dir("~/Desktop/Mouse_NarrowPeaks_Fixed/IDR_called_peaks/Brain_Filtered", pattern = "*.bed", full.names = TRUE)
#mouse_IDR_sorted_beds <- gsub(pattern = ".bed", replacement = "_sorted_mm10.bed", x = mouse_IDR_beds)
#for(i in 1:12){
#  system(paste("sort -k 1,1 -k2,2n", mouse_IDR_beds[i], ">", mouse_IDR_sorted_beds[i]))
#}
#mouse_multiinter <- correct_multiinter(mouse_IDR_sorted_beds, save_bed = T, output_bed_name = "~/Desktop/Capellini_Lab/Mouse_ATAC/mouse_multiinter_mm10.bed")

#liftover_to_hg38("~/Desktop/Capellini_Lab/Mouse_ATAC/mouse_multiinter_simple_mm10.bed", chain = "mm10tohg38")

#count mouse peaks lifted to human
mouse_total_n <- system("wc -l ~/Desktop/Capellini_Lab/Mouse_ATAC/mouse_multiinter_simple_hg38.bed", intern=T) %>% readr::parse_number()
#count overlap
human_mouse_total_overlap <- system("bedtools intersect -a ~/Desktop/Capellini_Lab/Human_ATAC/Human_Autopod_IDR_peaks/Brain_Filtered/human_all_merge_IDR_0.05_brain_filtered.bed -b ~/Desktop/Capellini_Lab/Mouse_ATAC/mouse_multiinter_simple_hg38.bed | wc -l", intern = T) %>% parse_number()

human_mouse_total_overlap/human_total_n

##Create Venn Diagram
# move to new plotting page
grid.newpage()

# create pairwise Venn diagram
all_venn <- draw.pairwise.venn(area1 = human_total_n, 
                               area2 = mouse_total_n,
                               cross.area=human_mouse_total_overlap,
                               fill=c("red","grey"), 
                               category = c("All Human", "Mouse"), 
                               margin  = 0.05, 
                               cat.cex = 0.8,
                               cat.pos =c(0, 0))

#early overlaps
#list samples for comparison
comp_samples <- c("FM1", "FM3", "FM5", "FP1", "FP3", "FP5", "HM1", "HM3", "HM5", "HP1", "HP3", "HP5")

human_sample_list <- paste("~/Desktop/Capellini_Lab/Human_ATAC/Human_Autopod_IDR_peaks/Brain_Filtered/",comp_samples, "_early_IDR_0.05_brain_filtered.narrowPeak", sep = "", collapse = " ")
human_sample_list_early <- paste("~/Desktop/Capellini_Lab/Human_ATAC/Human_Autopod_IDR_peaks/Brain_Filtered/",comp_samples, "_early_IDR_0.05_brain_filtered.narrowPeak", sep = "", collapse = " ")

#specify files
mouse_file <- "~/Desktop/Capellini_Lab/Mouse_ATAC/mouse_multiinter_simple_hg38.bed"

#count peaks
#if reading a file use commented line
#human_n <- system(paste("wc -l ", human_file, sep = ""), intern=T) %>% readr::parse_number()

#combine desired bed files, sort, and merge to combine overlapping peaks. 
human_n_early <- system(paste("cat ", human_sample_list, " | sort -k 1,1 -k2,2n | bedtools merge -i stdin | wc -l ", sep = ""), intern=T) %>% readr::parse_number()
mouse_n <- system(paste("wc -l ", mouse_file, sep = ""), intern=T) %>% readr::parse_number()

#calculate overlap
#combine desired bed files, sort, and merge to combine overlapping peaks, finally overlap with moyse
#human_mouse_overlap <- system(paste("bedtools intersect -a ", human_file, " -b ", mouse_file,  " | wc -l ", sep = ""), intern=T) %>% readr::parse_number()
human_mouse_overlap_early <- system(paste("cat ", human_sample_list_early, " | sort -k 1,1 -k2,2n | bedtools merge -i stdin | bedtools intersect -a stdin -b ", mouse_file,  " | wc -l ", sep = ""), intern=T) %>% readr::parse_number()

##Create Venn Diagram
# move to new plotting page
grid.newpage()

# create pairwise Venn diagram
early_venn <- draw.pairwise.venn(area1 = human_n_early, 
                                 area2 = mouse_n,
                                 cross.area = human_mouse_overlap_early,
                                 fill = c("red","grey"), 
                                 category = c("Early Human", "Mouse"), 
                                 margin  = 0.05, 
                                 cat.cex = 0.8,
                                 cat.pos =c(0, 0))

#list samples for comparison
human_sample_list_late <- paste("~/Desktop/Capellini_Lab/Human_ATAC/Human_Autopod_IDR_peaks/Brain_Filtered/",comp_samples, "_late_IDR_0.05_brain_filtered.narrowPeak", sep = "", collapse = " ")

#combine desired bed files, sort, and merge to combine overlapping peaks. 
human_n_late <- system(paste("cat ", human_sample_list, " | sort -k 1,1 -k2,2n | bedtools merge -i stdin | wc -l ", sep = ""), intern=T) %>% readr::parse_number()

#calculate overlap
#combine desired bed files, sort, and merge to combine overlapping peaks, finally overlap with moyse
#human_mouse_overlap <- system(paste("bedtools intersect -a ", human_file, " -b ", mouse_file,  " | wc -l ", sep = ""), intern=T) %>% readr::parse_number()
human_mouse_overlap_late <- system(paste("cat ", human_sample_list_late, " | sort -k 1,1 -k2,2n | bedtools merge -i stdin | bedtools intersect -a stdin -b ", mouse_file,  " | wc -l ", sep = ""), intern=T) %>% readr::parse_number()

##Create Venn Diagram
# move to new plotting page
grid.newpage()
bgcolor("white")
# create pairwise Venn diagram
late_venn <- draw.pairwise.venn(area1 = human_n_late, 
                                area2=mouse_n,
                                cross.area = human_mouse_overlap_late,
                                fill=c("red","grey"), 
                                category = c("Late Human", "Mouse"), 
                                margin  = 0.1, 
                                cat.cex = 0.8, 
                                cat.pos =c(0, 0)) 

#ggsave(filename = "~/Desktop/Capellini_Lab/Human_ATAC/ATAC_human_mouse_overlap_venn_plot.png", plot = venn_plot, device = "png", height = 4, width = 3.5, units = "in")

#now for DARs

#count human peaks 
human_DAR_n <- system("wc -l /Users/alexanderokamoto/Desktop/Capellini_Lab/Human_ATAC/human_DARs_hg38.bed", intern=T) %>% readr::parse_number()

#now for DARs
liftover_to_hg38("/Users/alexanderokamoto/Desktop/Capellini_Lab/Mouse_ATAC/mouse_DARs_mm10.bed", chain = "mm10tohg38")

#count mouse peaks lifted to human
mouse_DAR_n <- system("wc -l ~/Desktop/Capellini_Lab/Mouse_ATAC/mouse_DARs_hg38.bed", intern=T) %>% readr::parse_number()
#count overlap
human_mouse_DAR_overlap <- system("bedtools intersect -a /Users/alexanderokamoto/Desktop/Capellini_Lab/Human_ATAC/human_DARs_hg38.bed -b ~/Desktop/Capellini_Lab/Mouse_ATAC/mouse_DARs_hg38.bed | wc -l", intern = T) %>% parse_number()

human_mouse_DAR_overlap/human_DAR_n

##Create Venn Diagram
# move to new plotting page
grid.newpage()

# create pairwise Venn diagram
DAR_venn <- draw.pairwise.venn(area1 = human_DAR_n, 
                               area2 = mouse_DAR_n,
                               cross.area=human_mouse_DAR_overlap,
                               fill=c("red","grey"), 
                               category = c("Human DARs", "Mouse DARs"), 
                               margin  = 0.05, 
                               cat.cex = 0.8,
                               cat.pos =c(0, 0))




#human-mouse DEG orthologs
human_DEGs <- unique(read.delim(file= "/Users/alexanderokamoto/Desktop/Capellini_Lab/Human Autopod RNA/All_Human_Bio_Pairwise_Autopod_DEGs2.tsv", sep = "\t", header = T))
mouse_DEGs <- unique(read.delim(file= "/Users/alexanderokamoto/Desktop/Capellini_Lab/Mouse Autopod RNA/All_Mouse_Bio_Pairwise_Autopod_DEGs.tsv", sep = "\t", header = T))

ortholog_table <- read.table(file = "~/Desktop/Capellini_Lab/Weekly_Coding/human_mouse_orthologs.txt", header = T, sep = "\t")

dim(ortholog_table)
human_orthologs <- length(which(ortholog_table$human_gene_id %in% human_DEGs$ensembl_gene_id))
mouse_orthologs <- length(which(ortholog_table$mouse_gene_id %in% mouse_DEGs$ensembl_gene_id))

DEG_overlap <- length(intersect(human_DEGs$ensembl_gene_id, mouse_to_human_gene(mouse_DEGs$ensembl_gene_id)))


##Create Venn Diagram
# move to new plotting page
grid.newpage()

# create pairwise Venn diagram
DEG_venn <- draw.pairwise.venn(area1 = human_orthologs, 
                                 area2 = mouse_orthologs,
                                 cross.area = DEG_overlap,
                                 fill = c("red","grey"), 
                                 category = c("Human DEGs", "Mouse DEGs"), 
                                 margin  = 0.05, 
                                 cat.cex = 0.8,
                                 cat.pos =c(0, 0))

#create function to all for modification of all file reading simultaneously
read_image_gg <- function(file){
  gg_image <- image_read(file) %>% image_crop(geometry = "290x270+40+30") %>%
    image_ggplot() 
  return(gg_image)
}

#create function to all for modification of all file reading simultaneously
unique_overlaps <- read_image_gg(file = "~/Desktop/Autopod_Skeleton_Figures/human_mouse_unique_overlaps.png")

#generate final plot
human_mouse_plot <- ggdraw() +
  draw_plot(all_venn, x = 0, y = 0.67, width = 0.5, height = 0.33) +
  draw_plot(early_venn, x = 0.5, y = 0.67, width = 0.5, height = 0.33) +
  draw_plot(late_venn, x = 0, y = 0.34, width = 0.5, height = 0.33) +
  draw_plot(DAR_venn, x = 0.5, y = 0.34, width = 0.5, height = 0.33) +
  draw_plot(DEG_venn, x = 0, y = 0.01, width = 1, height = 0.33) +
  
  draw_plot_label(label = c("a", "b", "c", "d", "e"), 
                  size = 15,
                  x = c(0, 0.5, 0, 0.5, 0), 
                  y = c(1, 1, 0.67, 0.67, 0.34)) + 
  bgcolor("white")

ggsave(plot = human_mouse_plot + panel_border(color = "black", size = 1), 
       filename = "~/Dropbox/Autopod Paper/Autopod_Paper_Figures/S4_Human_Mouse_comparison.png", 
       device = "png", dpi = 300, height = 220, width = 180, units = "mm")



