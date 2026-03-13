#create supplementary figure showing mouse autopod RNA and ATAC-seq results 
#Alexander Okamoto
#December 10, 2025

#load packages
library(tidyverse) #for data manipulation
library(readxl) #for loading excel files
library(ChIPQC) #for reading genomic ranges
library(scales) #for plotting options
library(ggpubr) #for plotting options
library("cowplot") #for plotting options
library(tidyverse) #for data manipulation
library(DESeq2)
#load custom functions for analysis
source("~/Desktop/Capellini_Lab/Weekly_Coding/RNA_Analysis_Functions.R") 

#load read count information
mouse_ATAC_peak_counts <- read.delim(file = "~/Desktop/Capellini_Lab/Mouse_ATAC/mouse_ATAC_peaks_countMatrix.txt", header = T)
#clean up column names
colnames(mouse_ATAC_peak_counts) <- gsub(
  pattern = "X.n.netscratch.capellini_lab.Lab.Mouse_noMT.", 
  replacement = "", x = colnames(mouse_ATAC_peak_counts) )
colnames(mouse_ATAC_peak_counts) <- gsub(
  pattern = ".noMT.bam", 
  replacement = "", x = colnames(mouse_ATAC_peak_counts) )
colnames(mouse_ATAC_peak_counts)[1] <- "peak"


#Metadata is group you want to look at vs individual sample names
mouse_ATAC_sampleMetaData <- data.frame(sample = colnames(mouse_ATAC_peak_counts)[7:42], 
                                        limbtype=c(rep(c(rep(c("Hind_limb"), 6), rep(c("Forelimb"), 6)), 3)), 
                                        
                                        digit=rep(c("I", "III", "V"), 12), 
                                        tissue = rep(c(rep("Metapodial", 3), rep("Phalange", 3), rep("Metapodial", 3), rep("Phalange", 3)), 3), 
                                        tissue2 = rep(c(rep("Metatarsal", 3), rep("HL_phalange", 3), rep("Metacarpal", 3), rep("FL_phalange", 3)), 3), 
                                        replicate = c(rep('rep1', 12), rep('rep2', 12), rep('rep3', 12)),
                                        sampleType = rep(c("FM1", "FM3", "FM5", "FP1", "FP3", "FP5","HM1", "HM3", "HM5", "HP1", "HP3", "HP5"), 3)
)

rownames(mouse_ATAC_sampleMetaData) <- mouse_ATAC_sampleMetaData$sample
mouse_ATAC_sampleMetaData$limbtype <- as.factor(mouse_ATAC_sampleMetaData$limbtype)
mouse_ATAC_sampleMetaData$replicate <- as.factor(mouse_ATAC_sampleMetaData$replicate)

mouse_ATAC_peak_ct <- mouse_ATAC_peak_counts[,c(1, 7:42)]

#glm design
# Again, we can include more things beside tissue, like we could look at the effect of tissue *and* age, or we could control for sex here. 
design <- formula(~replicate + limbtype)
#reduced <- formula(~limbtype)

#build DESeqDataSet
rsem.in <- DESeqDataSetFromMatrix(countData = mouse_ATAC_peak_ct, colData = mouse_ATAC_sampleMetaData, design = design, tidy = T)

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

mouse_pc_scores_full <- merge(x = mouse_pc_scores, by.x = 'sample', by.y = "row.names", y = mouse_ATAC_sampleMetaData)

#add plotting metadata
mouse_pc_scores_full$tissue2 <- gsub(pattern = "FL_phalange", replacement = "Hand Phalanges", x = mouse_pc_scores_full$tissue2)
mouse_pc_scores_full$tissue2 <- gsub(pattern = "HL_phalange", replacement = "Foot Phalanges", x = mouse_pc_scores_full$tissue2)
mouse_pc_scores_full$tissue2 <- factor(mouse_pc_scores_full$tissue2, levels = c("Metatarsal", "Metacarpal", "Hand Phalanges", "Foot Phalanges"))

#PC1: phalanges versus metapodials
#PC2: early versus late
mouse_PCA <- mouse_pc_scores_full %>% 
  ggplot(aes(x = PC1, y = PC2, shape = digit, color= tissue2, label = tissue)) + 
  geom_point() + 
  labs(x = paste0("PC1 (", round(mouse_pc_eigenvalues$pct[1], 1), "%)"), y = paste0("PC2 (", round(mouse_pc_eigenvalues$pct[2], 1), "%)  ")) + 
  labs(shape = "Digit", color = "Tissue") + 
  scale_shape_manual(values = c(16, 17, 3)) +
  theme(legend.position='none') + 
  scale_color_manual(values = c("#E69F00", "#009E73", "#D55E00", "#0072B2"))

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
sampleMetaData$sample_file <- paste("Mouse RSEM Output Reduced Corrected/", sampleMetaData$sample, sep = "")

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
  guides(color = guide_legend(ncol = 2, nrow = 2), 
         shape = guide_legend(ncol = 2, nrow = 2)) + 
  scale_color_manual(values = c("#E69F00", "#009E73", "#D55E00", "#0072B2"))

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


#investigate most associated genes 
refgenes <- read.table("~/Desktop/Capellini_Lab/Capellini_Lab_Peak_Sets/Other_Peak_sets/mm10_refGene.txt")

loadings <- as.data.frame(mouse_sample_pca$rotation)
top_genes_pc1 <- rownames(loadings[order(abs(loadings$PC1), decreasing = TRUE), ])[1:50]
top_genes_pc2 <- rownames(loadings[order(abs(loadings$PC2), decreasing = TRUE), ])[1:50]

top_genes_pc1_annot <- mouse_annot %>% dplyr::filter(ensembl_gene_id %in% top_genes_pc1)
top_genes_pc2_annot <- mouse_annot %>% dplyr::filter(ensembl_gene_id %in% top_genes_pc2)

top_genes_pc1_GO <- simple_go_analysis(top_genes_pc1_annot$mgi_symbol, universe = refgenes$V13, database = org.Mm.eg.db )
top_genes_pc2_GO <- simple_go_analysis(top_genes_pc2_annot$mgi_symbol, universe = refgenes$V13, database = org.Mm.eg.db )
top_genes_pc1_GO$GeneRatio <- paste("'", top_genes_pc1_GO$GeneRatio, "'", sep ="")
top_genes_pc2_GO$GeneRatio <- paste("'", top_genes_pc2_GO$GeneRatio, "'", sep ="")

write_tsv(x = top_genes_pc1_annot, file = "~/Desktop/Capellini_Lab/Mouse Autopod RNA/mouse_top_genes_pc1.tsv")
write_tsv(x = top_genes_pc2_annot, file = "~/Desktop/Capellini_Lab/Mouse Autopod RNA/mouse_top_genes_pc2.tsv")
write_tsv(x = top_genes_pc1_GO, file = "~/Desktop/Capellini_Lab/Mouse Autopod RNA/mouse_top_genes_pc1_GO.tsv")
write_tsv(x = top_genes_pc2_GO, file = "~/Desktop/Capellini_Lab/Mouse Autopod RNA/mouse_top_genes_pc2_GO.tsv")

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
  theme(legend.position="none") + 
  scale_color_manual(values = c("#E69F00", "#009E73", "#D55E00", "#0072B2"))

#make final plot
mouse_plot <- ggdraw() +
  draw_plot(mouse_PCA, x = 0, y = 0.5, width = 0.5, height = .5) +
  draw_plot(mouse_rna_pca_batch, x = 0.5, y = 0.5, width = 0.5, height = .5) +
  draw_plot(mouse_rna_pca, x = 0, y = 0, width = 0.5, height = .5) +
  draw_plot(pca_legend, x = 0.5, y = 0, width = 0.5, height = .5) +
  draw_plot_label(label = c("a", "b", "c"), 
                  size = 15,
                  x = c(0, 0.5, 0), 
                  y = c(1, 1, 0.5)) + 
  bgcolor("white")

mouse_plot <- ggarrange(plotlist = c(mouse_PCA, mouse_rna_pca_batch, mouse_rna_pca),
          ncol = 1, legend.grob = pca_legend, 
          legend = "bottom", 
          labels = c("a", "b", "c"))


ggsave(plot = mouse_plot + theme(plot.background = element_rect(fill = "white")) + panel_border(color = "black", size = 1), 
       filename = "~/Dropbox/Autopod Paper/Autopod_Paper_Figures/S2_Mouse_ATAC_RNA_PCAs.png", 
       device = "png", dpi = 300, height = 240, width = 180, units = "mm")
