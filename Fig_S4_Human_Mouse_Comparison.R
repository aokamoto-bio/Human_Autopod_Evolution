#create supplementary figure showing human and mouse ATAC-seq overlaps
#Alexander Okamoto
#June 6, 2025

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

## PART 2
## Compare human pleiotropy with mouse overlap

multiinter_both_df <- read.table( "~/Desktop/Capellini_Lab/Human_ATAC/corrected_multiinter_autopod_both_tp.bed")
colnames(multiinter_both_df) <- c("chr", "start", "end", "num", "list", paste(c("FM1", "FM2", "FM3", "FM4", "FM5", "FP1", "FP2", "FP3", "FP4", "FP5", "HM1", "HM2", "HM3", "HM4", "HM5", "HP1", "HP2", "HP3", "HP4", "HP5"), "_Early", sep = ""), paste(c("FM1", "FM2", "FM3", "FM4", "FM5", "FP1", "FP2", "FP3", "FP4", "FP5", "HM1", "HM2", "HM3", "HM4", "HM5", "HP1", "HP2", "HP3", "HP4", "HP5"), "_Late", sep = ""))

#make function generate limbtype specific and shared sets
add_multiinter_count_col <- function(dataset, grep_term, names){
  dataset[names] <- rowSums(dataset[,grep(pattern = grep_term, x = colnames(dataset))])
  return(dataset)
}

#add columns
multiinter_both_df <- add_multiinter_count_col(dataset = multiinter_both_df,grep_term = "Early", names = "Early")
multiinter_both_df <- add_multiinter_count_col(dataset = multiinter_both_df, grep_term = "Late", names = "Late")
multiinter_both_df$tp_max <- multiinter_both_df$Early
multiinter_both_df$tp_max[which(multiinter_both_df$Early < multiinter_both_df$Late)] <- multiinter_both_df$Late[which(multiinter_both_df$Early < multiinter_both_df$Late)] 


#read.table("~/Desktop/Capellini_Lab/Mouse_ATAC/IDR_called_peaks/Brain_Filtered/corrected_multiinter_autopod_mouse.bed") %>% 
#  dplyr::select(V1:V4) %>% 
#  write_bed("~/Desktop/Capellini_Lab/Mouse_ATAC/IDR_called_peaks/Brain_Filtered/corrected_multiinter_autopod_mouse_simple_mm10.bed")
#make corrected mouse multinter file
#liftover_to_hg38(bedfile = "~/Desktop/Capellini_Lab/Mouse_ATAC/IDR_called_peaks/Brain_Filtered/corrected_multiinter_autopod_mouse_simple_mm10.bed", chain = "mm10ToHg38")


#find mouse overlaps
system("bedtools intersect -a ~/Desktop/Capellini_Lab/Human_ATAC/corrected_multiinter_autopod_both_tp.bed -b ~/Desktop/Capellini_Lab/Mouse_ATAC/mouse_multiinter_simple_hg38.bed -wa > ~/Desktop/Capellini_Lab/Human_ATAC/human_mouse_corrected_multiinter_autopod_overlap.bed")

human_mouse_multiinter <- read.table( "~/Desktop/Capellini_Lab/Human_ATAC/human_mouse_corrected_multiinter_autopod_overlap.bed")
colnames(human_mouse_multiinter) <- c("chr", "start", "end", "num", "list", paste(c("FM1", "FM2", "FM3", "FM4", "FM5", "FP1", "FP2", "FP3", "FP4", "FP5", "HM1", "HM2", "HM3", "HM4", "HM5", "HP1", "HP2", "HP3", "HP4", "HP5"), "_Early", sep = ""), paste(c("FM1", "FM2", "FM3", "FM4", "FM5", "FP1", "FP2", "FP3", "FP4", "FP5", "HM1", "HM2", "HM3", "HM4", "HM5", "HP1", "HP2", "HP3", "HP4", "HP5"), "_Late", sep = ""))

#add columns
human_mouse_multiinter<- add_multiinter_count_col(dataset = human_mouse_multiinter,grep_term = "Early", names = "Early")
human_mouse_multiinter <- add_multiinter_count_col(dataset = human_mouse_multiinter, grep_term = "Late", names = "Late")
human_mouse_multiinter$tp_max <- human_mouse_multiinter$Early
human_mouse_multiinter$tp_max[which(human_mouse_multiinter$Early < human_mouse_multiinter$Late)] <- human_mouse_multiinter$Late[which(human_mouse_multiinter$Early < human_mouse_multiinter$Late)] 

#reformat data for informative plotting
human_mouse_multiinter_sum <- human_mouse_multiinter %>% group_by(tp_max) %>% dplyr::select(tp_max) %>% summarize(Freq = n())

multiinter_both_sum <- multiinter_both_df %>% group_by(tp_max) %>% dplyr::select(tp_max) %>% summarize(total = n())

human_mouse_multiinter_sum <- merge(human_mouse_multiinter_sum, multiinter_both_sum)

human_mouse_multiinter_sum$perc <- human_mouse_multiinter_sum$Freq / human_mouse_multiinter_sum$total * 100

#visualize
human_mouse_sharing_plot <- human_mouse_multiinter_sum %>% 
  ggplot(aes(x = tp_max, y = Freq, fill = perc)) +
  geom_bar(stat="identity") + 
  scale_x_continuous(breaks = breaks_pretty()) + 
  labs(x = "n tissues with accessibility", y = "n regulatory elements", fill ="% of all \nhuman REs") + 
  theme(legend.position = "top")


#unique overlaps 
mouse_unique <- dir("~/Desktop/Capellini_Lab/Mouse_ATAC/IDR_called_peaks/Brain_Filtered/Mouse_Autopod_Liftover_hg38", 
                    pattern = "*_unique_peaks_hg38.bed", full.names = TRUE)
human_unique <- dir("~/Desktop/Capellini_Lab/Human_ATAC/Human_Autopod_IDR_peaks/Brain_Filtered/Unique_Peaks", 
                    pattern = "*_both_tp_unique_peaks.bed", full.names = TRUE)


unique_df <- data.frame(tissue = comp_samples, mouse_bed = mouse_unique, human_bed = human_unique[c(1,3,5:6, 8, 10:11, 13, 15:16, 18, 20)])
unique_df$overlaps <- NA

for (i in 1:12){
  #count overlap
  unique_df$overlaps[i] <- system(
    paste("bedtools intersect -a ", unique_df$human_bed[i],  " -b ", unique_df$mouse_bed[i], "| wc -l", sep = ""), intern = T) %>% parse_number()
  
}
#calculate median number of overlaps

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
  draw_plot(early_venn, x = 0, y = 0.34, width = 0.5, height = 0.33) +
  draw_plot(late_venn, x = 0, y = 0.01, width = 0.5, height = 0.33) +
  draw_plot(human_mouse_sharing_plot, x = 0.5, y = 0.67, width = 0.5, height = 0.33) +
  draw_plot(unique_overlaps, x = 0.5, y = 0.34, width = 0.5, height = 0.33) +
  draw_plot(DEG_venn, x = 0.5, y = 0.01, width = 0.5, height = 0.33) +
  
  draw_plot_label(label = c("a", "b", "c", "d", "e", "f"), 
                  size = 15,
                  x = c(0, 0, 0, 0.5, 0.5, 0.5), 
                  y = c(1, 0.67, 0.34, 1, 0.67, 0.34)) + 
  bgcolor("white")

ggsave(plot = human_mouse_plot + panel_border(color = "black", size = 1), 
       filename = "~/Dropbox/Autopod Paper/Autopod_Paper_Figures/S4_Human_Mouse_comparison.png", 
       device = "png", dpi = 300, height = 220, width = 180, units = "mm")

mouse_multiinter%>% <- read.table("~/Desktop/Capellini_Lab/Mouse_ATAC/mouse_multiinter_mm10.bed", header =T)
colnames(mouse_multiinter) <- c("chr", "start", "end", "num", "list", "FM1",  "FM3",  "FM5", "FP1",  "FP3",  "FP5", "HM1",  "HM3",  "HM5", "HP1",  "HP3", "HP5")
mouse_multiinter$mm_id <- paste(mouse_multiinter$chr, ":", mouse_multiinter$start+1, "-", mouse_multiinter$end, sep = "")

#load dataset for translating components
human_mouse_multiinter <- read.table( "~/Desktop/Capellini_Lab/Mouse_ATAC/mouse_multiinter_simple_hg38.bed") %>% unique()
colnames(human_mouse_multiinter) <- c("chr", "start", "end", "mm_id", "something")
human_mouse_multiinter$id <- paste(human_mouse_multiinter$chr, ":", human_mouse_multiinter$start, "-", human_mouse_multiinter$end, sep = "")

mouse_multiinter2 <- merge(mouse_multiinter, human_mouse_multiinter %>% dplyr::select(id, mm_id)) %>% 
  dplyr::filter(num == 1)

mouse_multiinter2$list <- gsub(pattern = ",13", replacement = "", mouse_multiinter2$list)
mouse_multiinter2$mm_tissue <- NA
tissues <- colnames(mouse_multiinter2)[7:18]
for(i in 12:1){
  mouse_multiinter2$mm_tissue[which(mouse_multiinter2$list == i)] <- tissues[i]

}

multiinter_both_df3 <- multiinter_both_df2 %>% 
  dplyr::filter(digit_count == 1)
multiinter_both_df3$hs_tissue <- NA
multiinter_both_df3$id <- paste(multiinter_both_df3$chr, ":", multiinter_both_df3$start, "-", multiinter_both_df3$end, sep = "")


multiinter_both_df3$hs_tissue[which(multiinter_both_df3$HP1 == 1)] <- "HP1"
multiinter_both_df3$hs_tissue[which(multiinter_both_df3$HP3 == 1)] <- "HP3"
multiinter_both_df3$hs_tissue[which(multiinter_both_df3$HP5 == 1)] <- "HP5"
multiinter_both_df3$hs_tissue[which(multiinter_both_df3$FP1 == 1)] <- "FP1"
multiinter_both_df3$hs_tissue[which(multiinter_both_df3$FP3 == 1)] <- "FP3"
multiinter_both_df3$hs_tissue[which(multiinter_both_df3$FP5 == 1)] <- "FP5"
multiinter_both_df3$hs_tissue[which(multiinter_both_df3$HM1 == 1)] <- "HM1"
multiinter_both_df3$hs_tissue[which(multiinter_both_df3$HM3 == 1)] <- "HM3"
multiinter_both_df3$hs_tissue[which(multiinter_both_df3$HM5 == 1)] <- "HM5"
multiinter_both_df3$hs_tissue[which(multiinter_both_df3$FM1 == 1)] <- "FM1"
multiinter_both_df3$hs_tissue[which(multiinter_both_df3$FM3 == 1)] <- "FM3"
multiinter_both_df3$hs_tissue[which(multiinter_both_df3$FM5 == 1)] <- "FM5"

hs_mm_tissue_df <- merge(multiinter_both_df3 %>% dplyr::select(id, hs_tissue),  mouse_multiinter2 %>% dplyr::select(id, mm_tissue))


multiinter_both_df2 <- multiinter_both_df %>% mutate(
  FM1 = pmax(FM1_Early, FM1_Late),
  FM2 = pmax(FM2_Early, FM2_Late),
  FM3 = pmax(FM3_Early, FM3_Late),
  FM4 = pmax(FM4_Early, FM4_Late),
  FM5 = pmax(FM5_Early, FM5_Late),
  FP1 = pmax(FP1_Early, FP1_Late),
  FP2 = pmax(FP2_Early, FP2_Late),
  FP3 = pmax(FP3_Early, FP3_Late),
  FP4 = pmax(FP4_Early, FP4_Late),
  FP5 = pmax(FP5_Early, FP5_Late),
  HM1 = pmax(HM1_Early, HM1_Late),
  HM2 = pmax(HM2_Early, HM2_Late),
  HM3 = pmax(HM3_Early, HM3_Late),
  HM4 = pmax(HM4_Early, HM4_Late),
  HM5 = pmax(HM5_Early, HM5_Late),
  HP1 = pmax(HP1_Early, HP1_Late),
  HP2 = pmax(HP2_Early, HP2_Late),
  HP3 = pmax(HP3_Early, HP3_Late),
  HP4 = pmax(HP4_Early, HP4_Late),
  HP5 = pmax(HP5_Early, HP5_Late),
) %>% 
  mutate(digit_count = rowSums(across(FM1:HP5)))

human_mouse_multiinter2 <- human_mouse_multiinter %>% mutate(
  

#too fancy, now something simple
  unique_sharing_df <- data.frame(tissue = tissues, value = NA)
  
for (i in tissues){
  # mouse_bed <- paste("~/Desktop/Mouse_NarrowPeaks_Fixed/IDR_called_peaks/Brain_Filtered/Mouse_Unique/mouse_", i, "_unique_mm10.bed", sep = "")
   mouse_bed_hg38 <- paste("~/Desktop/Mouse_NarrowPeaks_Fixed/IDR_called_peaks/Brain_Filtered/Mouse_Unique/mouse_", i, "_unique_hg38.bed", sep = "")
   human_bed <- paste("~/Desktop/Mouse_NarrowPeaks_Fixed/IDR_called_peaks/Brain_Filtered/Mouse_Unique/human_", i, "_unique_hg38.bed", sep = "")
  # #make mouse file
  # mouse_multiinter2 %>%
  #   dplyr::filter(num == 1) %>% 
  #   dplyr::filter(mm_tissue == i) %>% 
  #   dplyr::select(chr, start, end) %>% 
  #   write_bed(filename = mouse_bed, ncol = 3)
  # liftover_to_hg38(mouse_bed, chain = "mm10tohg38")
  # multiinter_both_df3 %>%
  #   dplyr::filter(digit_count == 1) %>% 
  #   dplyr::filter(hs_tissue == i) %>% 
  #   dplyr::select(chr, start, end) %>% 
  #   write_bed(filename = human_bed, ncol = 3)
  
  print(i)
  unique_sharing_df$value[which(unique_sharing_df$tissue == i)] <- system(paste("bedtools intersect -a", human_bed, "-b", mouse_bed_hg38, "| wc -l", sep = " "), intern=T) %>% readr::parse_number()
  
  
}

plot_autopod_expression(expr_df = unique_sharing_df, file_name = "human_mouse_unique_overlaps", center = T)



