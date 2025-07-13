#create supplementary figure showing mouse autopod ATAC-seq results 
#Alexander Okamoto
#February 12, 2025

#load packages
library(tidyverse) #for data manipulation
library(readxl) #for loading excel files
library(ChIPQC) #for reading genomic ranges
library(scales) #for plotting options
library(ggpubr) #for plotting options
library("cowplot") #for plotting options
library(tidyverse) #for data manipulation
#load custom functions for analysis
source("~/Desktop/Capellini_Lab/Weekly_Coding/RNA_Analysis_Functions.R") 


#load data
mouse_IDR_peaks <- dir("~/Desktop/Mouse_NarrowPeaks_Fixed/IDR_called_peaks/Brain_Filtered", pattern = "*_IDR_0.05_brain_filtered.narrowPeak", 
                       full.names = TRUE)

#create dataframe with metadata
#read in sample info
mouse_IDR_BF_sample_info <- read_excel("~/Desktop/Capellini_Lab/Mouse_ATAC/IDR_called_peaks/Brain_Filtered/Mouse_IDR_BF_Metadata.xlsx")

#make dataframe of ATAC peak sharing to use for PCA
mouse_IDR_Peaks <- lapply(mouse_IDR_peaks, ChIPQC:::GetGRanges, simple = TRUE)
names(mouse_IDR_Peaks) <- mouse_IDR_BF_sample_info$sample
mouseIDR_GRangesList<-GRangesList(mouse_IDR_Peaks)   
mouse_IDR_reduced <- GenomicRanges::reduce(unlist(mouseIDR_GRangesList))
mouse_IDR_consensusIDs <- paste0("consensus_", seq(1, length(mouse_IDR_reduced)))
mcols(mouse_IDR_reduced) <- do.call(cbind, lapply(mouseIDR_GRangesList, function(x) (mouse_IDR_reduced %over% x) + 0))
mouse_IDR_reducedConsensus <- mouse_IDR_reduced
mcols(mouse_IDR_reducedConsensus) <- cbind(as.data.frame(mcols(mouse_IDR_reducedConsensus)), mouse_IDR_consensusIDs)

#convert to dataframe and run PCA analysis
mouse_IDR_ATAC_PCA_matrix <- as.data.frame(elementMetadata(mouse_IDR_reducedConsensus)) %>% dplyr::select(-mouse_IDR_consensusIDs)  %>% as.matrix() %>% t()
mouse_IDR_ATAC_PCA_results <- prcomp(mouse_IDR_ATAC_PCA_matrix)

# Variance explained by PCs
mouse_IDR_pc_eigenvalues <- mouse_IDR_ATAC_PCA_results$sdev^2

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
mouse_IDR_pc_eigenvalues <- tibble(PC = factor(1:length(mouse_IDR_pc_eigenvalues)), 
                                   variance = mouse_IDR_pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

#add sample metadata for plotting
mouse_IDR_BF_sample_info$tissue_region <- c(rep(x = "Metatarsal", 3), rep(x = "Foot Phalanges", 3), rep(x = "Metacarpal", 3), rep(x = "Hand Phalanges", 3))
mouse_IDR_ATAC_PCA_results_tidy <- mouse_IDR_ATAC_PCA_results$x %>% data.frame %>% mutate(sample = rownames(.)) %>% merge(y = mouse_IDR_BF_sample_info)

#PC1: phalanges versus metapodials
#PC2: early versus late
mouse_PCA <- mouse_IDR_ATAC_PCA_results_tidy %>% 
  ggplot(aes(x = PC1, y = PC2, shape = digit, color= tissue_region, label = tissue)) + 
  geom_point() + 
  labs(x = paste0("PC1 (", round(mouse_IDR_pc_eigenvalues$pct[1], 1), "%)"), y = paste0("PC2 (", round(mouse_IDR_pc_eigenvalues$pct[2], 1), "%)  ")) + 
  labs(shape = "Digit", color = "Tissue") + 
  scale_shape_manual(values = c(16, 17, 3)) +
  theme(legend.position='none')

#plot counts of element sharing

#make corrected multiinter file
mouse_IDR_beds <- gsub(pattern = ".narrowPeak", replacement = ".sorted.bed", mouse_IDR_peaks)
#mouse_multiinter <- correct_multiinter(file_list = mouse_IDR_beds, save_bed = T, output_bed_name = "~/Desktop/Capellini_Lab/Mouse_ATAC/IDR_called_peaks/Brain_Filtered/corrected_multiinter_autopod_mouse.bed")
mouse_multiinter_df <- read.table("~/Desktop/Capellini_Lab/Mouse_ATAC/IDR_called_peaks/Brain_Filtered/corrected_multiinter_autopod_mouse.bed")
colnames(mouse_multiinter_df) <- c("chr", "start", "end", "num", "list", "FM1", "FM3",  "FM5", "FP1",  "FP3", "FP5", "HM1", "HM3", "HM5", "HP1", "HP3", "HP5")

#plot results

pleiotropy_plot <- mouse_multiinter_df %>% 
  ggplot(aes(x = num)) +
  geom_bar() + 
  scale_x_continuous(breaks = breaks_pretty()) + 
  labs(x = "Number of Tissues with Accessibility", y = "Number of \nRegulatory \nElements")



## RNA-seq 

#Metadata is group you want to look at vs individual sample names
sampleMetaData <- data.frame(sample = c("HP1A", "HP3A", "HP5A", "HM1A", "HM3A", "HM5A", "FP1A", "FP3A", "FP5A", "FM1A", "FM3A", "FM5A", "HP1D", "HP3D", "HP5D", "HM1D", "HM3D", "HM5D", "FP1D", "FP3D", "FP5D", "FM1D", "FM3D", "FM5D", "HP1F", "HP3F", "HP5F", "HM1F", "HM3F", "HM5F", "FP1F", "FP3F", "FP5F", "FM1F", "FM3F", "FM5F", "HP1G", "HP3G", "HP5G", "HM1G", "HM3G", "HM5G", "FP1G", "FP3G", "FP5G", "FM1G", "FM3G", "FM5G", "HP1H", "HP3H", "HP5H", "HM1H", "HM3H", "HM5H", "FP1H", "FP3H", "FP5H", "FM1H", "FM3H", "FM5H", "HP1I", "HP3I", "HP5I", "HM1I", "HM3I", "HM5I", "FP1I", "FP3I", "FP5I", "FM1I", "FM3I", "FM5I", "HP3D_alt", "FP3D_alt"), 
                             limbtype=c(rep(c(rep(c("Forelimb"), 6), rep(c("Hind_Limb"), 6)), 6), "Forelimb", "Hind_Limb"), 
                             
                             digit=c(rep(c("I", "III", "V"), 24), "III", "III"), 
                             tissue = c(rep(c(rep("Phalange", 3), rep("Metapodial", 3), rep("Phalange", 3), rep("Metapodial", 3)), 6), "Phalange", "Phalange"), 
                             tissue2 = c(rep(c(rep("Hand Phalanges", 3), rep("Metacarpal", 3), rep("Foot Phalanges", 3), rep("Metatarsal", 3)), 6), "Hand Phalanges", "Foot Phalanges"), 
                             
                             replicate = c(rep('A', 12), rep('D', 12), rep('F', 12), rep('G', 12), rep('H', 12), rep('I', 12), 'D', 'D'), 
                             
                             sampleType = c(rep(c("HP1", "HP3", "HP5", "HM1", "HM3", "HM5", "FP1", "FP3", "FP5", "FM1", "FM3", "FM5"), 6), "HP3", "FP3"), 
                             
                             sex = c(rep("Female", 12), rep("Male", 12), rep("Female", 12), rep("Male", 24), rep("Female", 12), "Male", "Male"), 
                             
                             prep  =  c(rep("SMART-Seq", 72), rep("KAPA", 2)))

rownames(sampleMetaData) <- sampleMetaData$sample
sampleMetaData$limbtype <- as.factor(sampleMetaData$limbtype)
sampleMetaData$replicate <- as.factor(sampleMetaData$replicate)
sampleMetaData$sex <- as.factor(sampleMetaData$sex)
#modify digit for 
sampleMetaData$digit[73:74] <- "III KAPA"

#correct path for alternate gene results
sampleMetaData$sample_file <- paste("Mouse RSEM Output Reduced/", sampleMetaData$sample, sep = "")
sampleMetaData$sample_file[which(sampleMetaData$prep == "KAPA")] <- gsub(" Reduced", replacement = "", sampleMetaData$sample_file[which(sampleMetaData$prep == "KAPA")])

#load the counts data
mouse_ct <- read_in_counts(samplenames = sampleMetaData$sample_file, filepath = "~/Desktop/Capellini_Lab/Mouse Autopod RNA/")
colnames(mouse_ct)[2:ncol(mouse_ct)] <- sampleMetaData$sample

#glm design
# Again, we can include more things beside tissue, like we could look at the effect of tissue *and* age, or we could control for sex here. 
design <- formula(~replicate + limbtype)
#reduced <- formula(~limbtype)

#build DESeqDataSet
rsem.in <- DESeqDataSetFromMatrix(countData = mouse_ct, colData = sampleMetaData, design = design, tidy = T)

#fix gene names in rsem.in
#don't have to do this--when deseq reads in gene names from rsem file it duplicates name of gene. Gets rid of second name. 
#rownames(rsem.in) <-  sapply(strsplit(c(rownames(rsem.in)), split = '_', fixed = TRUE), function(x) (x[1]))

#calculate DGE likelihood ratio test
#LRT needs reduced. Just W t

#, test = "LRT", reduced = reducedest--pairwise does not need it. Or post-hoc tests can compare between pairs. 
rsem.de <- DESeq(rsem.in)

#perform variance standardizing transformation
rsem.vst <- vst(rsem.de)
mat <- assay(rsem.vst)
mm <- model.matrix(~replicate, colData(rsem.vst))
#mat <- limma::removeBatchEffect(mat, batch=colData(rsem.vst)$replicate)
mat <- limma::removeBatchEffect(mat, batch=colData(rsem.vst)$prep)

assay(rsem.vst) <- mat

# perform a PCA on the data in assay(x) for the selected genes
mouse_sample_pca <- prcomp(t(assay(rsem.vst)))

# Variance explained by PCs
mouse_pc_eigenvalues <- mouse_sample_pca$sdev^2

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
mouse_pc_eigenvalues <- tibble(PC = factor(1:length(mouse_pc_eigenvalues)), 
                               variance = mouse_pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

#visualize
mouse_pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

# The PC scores are stored in the "x" value of the prcomp object
mouse_pc_scores <- mouse_sample_pca$x %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

mouse_pc_scores_full <- merge(x = mouse_pc_scores, by.x = 'sample', by.y = "row.names", y = sampleMetaData)


# create the plot showing batch effect correction

mouse_rna_pca_batch <- mouse_pc_scores_full %>% 
  ggplot(aes(x = PC1, y = PC2, shape = digit, color = tissue2)) +
  geom_point() + 
  labs(x = paste0("PC1 (", round(mouse_pc_eigenvalues$pct[1], 1), "%)"), 
       y = paste0("PC2 (", round(mouse_pc_eigenvalues$pct[2], 1), "%)  "), 
       color = "Tissue", shape = "Digit") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  scale_shape_manual(values = c(16, 17, 2, 3)) + 
  guides(color = guide_legend(ncol = 1, nrow = 4), shape = guide_legend(ncol = 1, nrow = 4))

#get just the legend for plotting separately
pca_legend <- get_legend(mouse_rna_pca_batch)

#remove legend from figure 
mouse_rna_pca_batch <- mouse_rna_pca_batch + theme(legend.position='none')

mouse_rna_pca_batch

#repeat without KAPA samples
#just running everything to avoid any tiny errors

#remove kappa samples 
sampleMetaData2 <- sampleMetaData %>% dplyr::filter(prep != "KAPA")


#load the counts data
mouse_ct2 <- read_in_counts(samplenames = sampleMetaData2$sample_file, filepath = "~/Desktop/Capellini_Lab/Mouse Autopod RNA/")
colnames(mouse_ct2)[2:ncol(mouse_ct2)] <- sampleMetaData2$sample

#glm design
# Again, we can include more things beside tissue, like we could look at the effect of tissue *and* age, or we could control for sex here. 
design <- formula(~replicate + limbtype)
#reduced <- formula(~limbtype)

#build DESeqDataSet
rsem.in <- DESeqDataSetFromMatrix(countData = mouse_ct2, colData = sampleMetaData2, design = design, tidy = T)

#fix gene names in rsem.in
#don't have to do this--when deseq reads in gene names from rsem file it duplicates name of gene. Gets rid of second name. 
#rownames(rsem.in) <-  sapply(strsplit(c(rownames(rsem.in)), split = '_', fixed = TRUE), function(x) (x[1]))

#calculate DGE likelihood ratio test
#LRT needs reduced. Just W t

#, test = "LRT", reduced = reducedest--pairwise does not need it. Or post-hoc tests can compare between pairs. 
rsem.de <- DESeq(rsem.in)

#perform variance standardizing transformation
rsem.vst <- vst(rsem.de)
mat <- assay(rsem.vst)
mm <- model.matrix(~replicate, colData(rsem.vst))
mat <- limma::removeBatchEffect(mat, batch=colData(rsem.vst)$replicate)


assay(rsem.vst) <- mat

# perform a PCA on the data in assay(x) for the selected genes
mouse_sample_pca <- prcomp(t(assay(rsem.vst)))

# Variance explained by PCs
mouse_pc_eigenvalues2 <- mouse_sample_pca$sdev^2

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
mouse_pc_eigenvalues2 <- tibble(PC = factor(1:length(mouse_pc_eigenvalues2)), 
                               variance = mouse_pc_eigenvalues2) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

#visualize
mouse_pc_eigenvalues2 %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

# The PC scores are stored in the "x" value of the prcomp object
mouse_pc_scores <- mouse_sample_pca$x %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

mouse_pc_scores_full2 <- merge(x = mouse_pc_scores, by.x = 'sample', by.y = "row.names", y = sampleMetaData)


# create the plot showing batch effect correction

mouse_rna_pca <- mouse_pc_scores_full2 %>% 
  ggplot(aes(x = PC1, y = PC2, shape = digit, color = tissue2)) +
  geom_point() + 
  labs(x = paste0("PC1 (", round(mouse_pc_eigenvalues2$pct[1], 1), "%)"), y = paste0("PC2 (", round(mouse_pc_eigenvalues2$pct[2], 1), "%)  "), color = "Limb Type", shape = "Tissue") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_shape_manual(values = c(16, 17, 3)) + 
  theme(legend.position="none")

mouse_rna_pca

#load cluster plots

plt1 <- image_read("~/Desktop/Autopod_Skeleton_Figures/mouse_RNA_cluster_1.png") %>% image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt2 <- image_read("~/Desktop/Autopod_Skeleton_Figures/mouse_RNA_cluster_2.png") %>% image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt3 <- image_read("~/Desktop/Autopod_Skeleton_Figures/mouse_RNA_cluster_3.png") %>% image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt4 <- image_read("~/Desktop/Autopod_Skeleton_Figures/mouse_RNA_cluster_4.png") %>% image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt5 <- image_read("~/Desktop/Autopod_Skeleton_Figures/mouse_RNA_cluster_5.png") %>% image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt6 <- image_read("~/Desktop/Autopod_Skeleton_Figures/mouse_RNA_cluster_6.png") %>% image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt7 <- image_read("~/Desktop/Autopod_Skeleton_Figures/mouse_RNA_cluster_7.png") %>% image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt8 <- image_read("~/Desktop/Autopod_Skeleton_Figures/mouse_RNA_cluster_8.png") %>% image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt9 <- image_read("~/Desktop/Autopod_Skeleton_Figures/mouse_RNA_cluster_9.png") %>% image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt10 <- image_read("~/Desktop/Autopod_Skeleton_Figures/mouse_RNA_cluster_10.png") %>% 
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt11 <- image_read("~/Desktop/Autopod_Skeleton_Figures/mouse_RNA_cluster_11.png") %>% image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt12 <- image_read("~/Desktop/Autopod_Skeleton_Figures/mouse_RNA_cluster_12.png") %>% image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
plt12


#make final plot
skel_h <- 0.2
skel_w <- 0.2

mouse_plot <- ggdraw() +
  draw_plot(mouse_PCA, x = 0, y = 0.8, width = 0.5, height = .2) +
  draw_plot(pleiotropy_plot, x = 0.5, y = 0.8, width = 0.5, height = .2) +
  draw_plot(pca_legend, x = 0, y = 0.4, width = 0.6, height = .2) +
  draw_plot(mouse_rna_pca_batch, x = 0, y = 0.6, width = 0.5, height = .2) +
  draw_plot(mouse_rna_pca, x = 0.5, y = 0.6, width = 0.5, height = .2) +
  draw_plot(plt1, x = 0.6, y = 0.4, width = skel_w, height = skel_h) +
  draw_plot(plt2, x = 0.8, y = 0.4, width = skel_w, height = skel_h) +
  draw_plot(plt3, x = 0, y = 0.2, width = skel_w, height = skel_h) +
  draw_plot(plt4, x = 0.2, y = 0.2, width = skel_w, height = skel_h) +
  draw_plot(plt5, x = 0.4, y = 0.2, width = skel_w, height = skel_h) +
  draw_plot(plt6, x = 0.6, y = 0.2, width = skel_w, height = skel_h) +
  draw_plot(plt7, x = 0.8, y = 0.2, width = skel_w, height = skel_h) +
  draw_plot(plt8, x = 0, y = 0, width = skel_w, height = skel_h) +
  draw_plot(plt9, x = 0.2, y = 0, width = skel_w, height = skel_h) +
  draw_plot(plt10, x = 0.4, y = 0, width = skel_w, height = skel_h) +
  draw_plot(plt11, x = 0.6, y = 0, width = skel_w, height = skel_h) +
  draw_plot(plt12, x = 0.8, y = 0, width = skel_w, height = skel_h) +
  
  draw_plot_label(label = c("a", "b", "c", "d", "e1", "e2", "e3", "e4", "e5", "e6","e7", "e8", "e9", "e10","e11", "e12"), 
                  size = 15,
                  x = c(0, 0.5, 0, 0.5, 0.6, 0.8, 0, skel_w, skel_w*2, skel_w*3, skel_w*4, 0, skel_w, skel_w*2, skel_w*3, skel_w*4), 
                  y = c(1, 1, 0.8, 0.8, 0.6, 0.6, 0.4, 0.4, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.2, 0.2)) + 
  bgcolor("white")

ggsave(plot = mouse_plot + panel_border(color = "black", size = 1), 
       filename = "~/Dropbox/Autopod Paper/Autopod_Paper_Figures/S2_Mouse_ATAC_RNA.png", 
       device = "png", dpi = 300, height = 220, width = 180, units = "mm")
