#create PCA/Volcano plot figure showing human RNA and ATAC-seq information
#Alexander Okamoto
#December 10, 2025

#load packages
library(readxl)
library(GenomicRanges)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(DESeq2)
library(ggrepel)
library(ggbreak) #for breaking axes on plots

#load read count information
human_ATAC_peak_counts <- read.delim(file = "~/Desktop/Capellini_Lab/Human_ATAC/human_ATAC_peaks_countMatrix.txt", header = T)
#clean up column names
colnames(human_ATAC_peak_counts) <- gsub(
  pattern = "X.n.netscratch.capellini_lab.Lab.aokamoto.Human_noMT.", 
  replacement = "", x = colnames(human_ATAC_peak_counts) )
colnames(human_ATAC_peak_counts) <- gsub(
  pattern = ".noMT.bam", 
  replacement = "", x = colnames(human_ATAC_peak_counts) )
colnames(human_ATAC_peak_counts)[1] <- "peak"
sample_info <- read.table("~/Desktop/Capellini_Lab/Human_ATAC/ATAC_files2.txt", header=TRUE)

sample_info$limbtype <- as.factor(sample_info$limbtype)
sample_info$replicate <- as.factor(sample_info$H_ID)

human_ATAC_peak_ct <- human_ATAC_peak_counts[,c(1, 7:146)]

#remove sample with too few reads
human_ATAC_peak_ct$`28597_FP5` <- NULL
sample_info <- sample_info %>% dplyr::filter(sample != "28597_FP5")

## Run DESeq

#glm design
# Again, we can include more things beside tissue, like we could look at the effect of tissue *and* age, or we could control for sex here. 
design <- formula(~replicate + limbtype)
#reduced <- formula(~limbtype)

#build DESeqDataSet
rsem.in <- DESeqDataSetFromMatrix(countData = human_ATAC_peak_ct, colData = sample_info, design = design, tidy = T)

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
human_sample_pca <- prcomp(t(assay(rsem.vst)))

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
# Variance explained by PCs
human_atac_pc_eigenvalues <- human_sample_pca$sdev^2

human_atac_pc_eigenvalues <- tibble(PC = factor(1:length(human_atac_pc_eigenvalues)), 
                               variance = human_atac_pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# The PC scores are stored in the "x" value of the prcomp object
human_pc_scores <- human_sample_pca$x %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

human_pc_scores_full <- merge(x = human_pc_scores, by = 'sample', y = sample_info)


#PC1: phalanges versus metapodials
#PC2: early versus late

#add plotting metadata
human_pc_scores_full$tissue_type <- gsub(pattern = "FL_phalanges", replacement = "Hand Phalanges", x = human_pc_scores_full$tissue_type)
human_pc_scores_full$tissue_type <- gsub(pattern = "HL_phalanges", replacement = "Foot Phalanges", x = human_pc_scores_full$tissue_type)
human_pc_scores_full$tissue_type <- factor(human_pc_scores_full$tissue_type, levels = c("Metatarsal", "Metacarpal", "Hand Phalanges", "Foot Phalanges"))


#make the plot
human_atac_pca <- human_pc_scores_full %>% 
  ggplot(aes(x = PC1, y = PC2, shape = digit, color = tissue_type, fill = tissue_type, alpha = as.character(timepoint))) + 
  geom_point() + 
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = c("#E69F00", "#009E73", "#D55E00", "#0072B2")) + 
  scale_fill_manual(values = c("#E69F00", "#009E73", "#D55E00", "#0072B2")) +
  #scale_alpha_manual(values=c("Early", "Late")) + 
  labs(x = paste0("PC1 (", round(human_atac_pc_eigenvalues$pct[1], 1), "%)"), 
       y = paste0("PC2 (", round(human_atac_pc_eigenvalues$pct[2], 1), "%)"), 
       #color = "Tissue Type", 
       shape = "Digit", 
       alpha = "Stage", 
       fill = "Tissue",
       title = NULL ) +
  theme(legend.position = "bottom", 
        legend.box = "horizontal", 
        legend.margin = margin(), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 8), 
        ) + 
  geom_point(alpha = 1, fill = NA) 

timepoint_DARs <- runDESEQ2_df_hs_atac(design_in = formula(~timepoint), count_data = human_ATAC_peak_ct, col_data =  sample_info)
timepoint_DARs$expressed <- "Not Significant"
timepoint_DARs$expressed[timepoint_DARs$log2FoldChange < -0.58 & timepoint_DARs$padj < 0.01] <- "Upregulated early"
timepoint_DARs$expressed[timepoint_DARs$log2FoldChange > 0.58 & timepoint_DARs$padj < 0.01] <- "Upregulated late"

atac_stage_volcano <- ggplot(data = timepoint_DARs, 
                        aes( x = log2FoldChange, y =-log10(padj), colour = expressed)) + 
  geom_point(size = 0.1) +
  scale_color_manual(values = c("grey", "#00AFBB", "brown3"))+ 
  geom_hline(yintercept = -log10(0.01), 
             col = "gray", 
             linetype = 'dashed') +
  geom_vline(xintercept = c(-0.58, 0.58), col = "gray", linetype = 'dashed') +
  labs(x = expression("log"[2]*"FC"), 
       y = expression("-log"[10]*" adj. p-value"), 
       color = NULL) + 
  theme_classic() +
  theme(legend.position="bottom", 
        legend.box = "vertical", 
        legend.margin=margin(), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10), 
        legend.key.size = unit(0.1, "cm"))



#now the RNA plot

#Determine the file names. I change this so much that I like to hardcode it, but you can also just input a text file or table.
sampleNames <- read.table(file = "~/Desktop/Capellini_Lab/Human Autopod RNA/Human RSEM Output CORRECTED/human_autopod_RNA_samples.txt", header = F)
#sort sample names alphabetically 
sampleNames <- sort(sampleNames$V1)

#make DESeq object. Again, you could hardcode or have in a file. I like to hardcode here just so that I can change it more easily, but if 
# I wanted to add additional variables (like sex for the human data for example), it would be easier to have a metadata file.
#"tissue" is part of your formula design, so we can change it or add to it as necessary. For example, if I also wanted to look at
# the effect of timepoint, I could code that here as well (again though I would want to have it in a file)

#Metadata is group you want to look at vs individual sample names
sampleMetaData <- data.frame(sample = sampleNames, limbtype=c(rep(c("Hind_limb"), 10), rep(c("Forelimb"), 10)), digit=c("I", "II", "III", "IV", "V"), tissue = c(rep("Metapodial", 5), rep("Phalange", 5), rep("Metapodial", 5), rep("Phalange", 5)), tissue2 = c(rep("Metatarsal", 5), rep("Hand Phalanges", 5), rep("Metacarpal", 5), rep("Foot Phalange", 5)), replicate = c(rep('H28662', 20), rep('H28670', 20), rep('H28693A', 20), rep('H28693B', 20), rep('H28743', 20), rep('H28769', 20), rep('H28841', 20), rep('H28852', 20), rep('H28853', 20), rep('H28867', 20), rep('H28906', 20), rep('H29057', 20)), sampleType = c("FM1", "FM2", "FM3", "FM4", "FM5", "FP1", "FP2", "FP3", "FP4", "FP5", "HM1", "HM2", "HM3", "HM4", "HM5", "HP1", "HP2", "HP3", "HP4", "HP5"), sex = c(rep("Male", 100), rep("Female", 20), rep("Male", 40), rep("Female", 20), rep("Male", 60)), timepoint = c(rep('Late', 100), rep('Early', 40), rep('Late', 20), rep('Early', 80)), stage = c(rep('72', 20), rep('67', 20), rep('74', 20), rep('72', 20), rep('67', 20), rep('57', 20), rep('53', 20), rep('72', 20), rep('59', 20), rep('54', 20), rep('57', 20), rep('54', 20)))

#remove problematic samples
sampleMetaData <- sampleMetaData %>% dplyr::filter(!sample %in% c("28867_HM1", "29057_FM3"))

rownames(sampleMetaData) <- sampleMetaData$sample

#change variables to factors
sampleMetaData$limbtype <- as.factor(sampleMetaData$limbtype)
sampleMetaData$replicate <- as.factor(sampleMetaData$replicate)
sampleMetaData$sex <- as.factor(sampleMetaData$sex)
sampleMetaData$timepoint <- as.factor(sampleMetaData$timepoint)
sampleMetaData$stage <- as.factor(sampleMetaData$stage)
sampleMetaData$sample_tp <- paste(sampleMetaData$sampleType, "_", sampleMetaData$timepoint, sep = "")

#load the counts data
ct <- read_in_counts(samplenames = sampleMetaData$sample, filepath = "~/Desktop/Capellini_Lab/Human Autopod RNA/Human RSEM Output CORRECTED//")

#remove version number from gene_id so it will be compatible with biomart
ct$gene_id = sapply(strsplit(ct$gene_id, ".", fixed=T), function(x) x[1])

#will sum the counts for matches to the same gene
ct <- aggregate(. ~ gene_id, data= ct, FUN=sum)

#change row names to unique values
rownames(ct) <- ct$gene_id
names(ct)[1] <- "ensembl_gene_id"

#remove samples with low read counts
ct$`28867_HM1` <- NULL
ct$`29057_FM3` <- NULL

# we can include more things beside tissue, like we could look at the effect of tissue *and* age, or we could control for sex here. 
# note the DESeq2 cannot handle perfectly correlated variables like individual and sex since one contains a full set of the other
# due to this, I only control for individual (replicate)
design <- formula(~replicate + tissue)

#build DESeqDataSet
human_rsem.in <- DESeqDataSetFromMatrix(countData = ct, colData = sampleMetaData, design = design, tidy = T)
#, test = "LRT", reduced = reducedest--pairwise does not need it. Or post-hoc tests can compare between pairs. 
human_rsem.de <- DESeq(human_rsem.in)
human_normalized_counts <- counts(human_rsem.de, normalized=TRUE)

#perform variance standardizing transformation
human_rsem.vst <- vst(human_rsem.de)
mat <- assay(human_rsem.vst)
mm <- model.matrix(~replicate, colData(human_rsem.vst))
mat <- limma::removeBatchEffect(mat, batch=colData(human_rsem.vst)$replicate)
assay(human_rsem.vst) <- mat

# perform a PCA on the data in assay(x) for the selected genes
human_sample_pca <- prcomp(t(assay(human_rsem.vst)))

#investigate most associated genes 
loadings <- as.data.frame(human_sample_pca$rotation)
top_genes_pc1 <- rownames(loadings[order(abs(loadings$PC1), decreasing = TRUE), ])[1:50]
top_genes_pc2 <- rownames(loadings[order(abs(loadings$PC2), decreasing = TRUE), ])[1:50]

top_genes_pc1_annot <- annot %>% dplyr::filter(ensembl_gene_id %in% top_genes_pc1)
top_genes_pc2_annot <- annot %>% dplyr::filter(ensembl_gene_id %in% top_genes_pc2)

top_genes_pc1_GO <- simple_go_analysis(top_genes_pc1_annot$hgnc_symbol)
top_genes_pc2_GO <- simple_go_analysis(top_genes_pc2_annot$hgnc_symbol)
top_genes_pc1_GO$GeneRatio <- paste("'", top_genes_pc1_GO$GeneRatio, "'", sep ="")
top_genes_pc2_GO$GeneRatio <- paste("'", top_genes_pc2_GO$GeneRatio, "'", sep ="")

write_tsv(x = top_genes_pc1_annot, file = "~/Desktop/Capellini_Lab/Human Autopod RNA/human_top_genes_pc1.tsv")
write_tsv(x = top_genes_pc2_annot, file = "~/Desktop/Capellini_Lab/Human Autopod RNA/human_top_genes_pc2.tsv")
write_tsv(x = top_genes_pc1_GO, file = "~/Desktop/Capellini_Lab/Human Autopod RNA/human_top_genes_pc1_GO.tsv")
write_tsv(x = top_genes_pc2_GO, file = "~/Desktop/Capellini_Lab/Human Autopod RNA/human_top_genes_pc2_GO.tsv")
                    
# Variance explained by PCs
human_pc_eigenvalues <- human_sample_pca$sdev^2

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
human_pc_eigenvalues_RNA <- tibble(PC = factor(1:length(human_pc_eigenvalues)), 
                               variance = human_pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# The PC scores are stored in the "x" value of the prcomp object
human_pc_scores <- human_sample_pca$x %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

human_pc_scores_full <- merge(x = human_pc_scores, by = 'sample', y = sampleMetaData)

human_pc_scores_full$`Tissue Type` <- human_pc_scores_full$tissue2
#human_pc_scores_full$`Tissue Type` <- gsub(pattern = "Metatarsal", replacement = "MT", x = human_pc_scores_full$`Tissue Type`)
#human_pc_scores_full$`Tissue Type` <- gsub(pattern = "Metacarpal", replacement = "MC", x = human_pc_scores_full$`Tissue Type`)
human_pc_scores_full$`Tissue Type` <- gsub(pattern = "Foot Phalange", replacement = "Foot Phalanges", human_pc_scores_full$`Tissue Type`)
#human_pc_scores_full$`Tissue Type` <- gsub(pattern = "Hand Phalanges", replacement = "Hand Phal.", x = human_pc_scores_full$`Tissue Type`)
#human_pc_scores_full$`Tissue Type` <- gsub(pattern = "Foot Phalanges", replacement = "Foot Phal.", x = human_pc_scores_full$`Tissue Type`)
human_pc_scores_full$`Tissue Type` <- factor(human_pc_scores_full$`Tissue Type`, levels = c("Metatarsal", "Metacarpal", "Hand Phalanges", "Foot Phalanges"))

#make the plot
human_rna_pca <- human_pc_scores_full %>% 
  ggplot(aes(x = PC1, y = PC2, color = `Tissue Type`, shape = digit)) + 
  geom_point(aes(fill = `Tissue Type`, alpha = as.character(timepoint))) + 
  scale_color_manual(values = c("#E69F00", "#009E73", "#D55E00", "#0072B2")) + 
  scale_fill_manual(values = c("#E69F00", "#009E73", "#D55E00", "#0072B2")) +
  scale_shape_manual(values=c(21:25)) +
  #scale_alpha_manual(values=c("Early", "Late")) + 
  labs(x = paste0("PC1 (", round(human_pc_eigenvalues_RNA$pct[1], 1), "%)"), 
       y = paste0("PC2 (", round(human_pc_eigenvalues_RNA$pct[2], 1), "%)"), 
       #color = "Tissue Type", 
       shape = "Digit", 
       alpha = "Stage", 
       color = "Tissue",
       fill = "Tissue",
       title = NULL ) +
  theme(legend.position = "bottom", 
        legend.box = "horizontal", 
        legend.margin = margin(), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.1, "cm")) +
    geom_point(alpha = 1, fill = NA)


human_rna_pca2 <- human_pc_scores_full %>% 
  ggplot(aes(x = PC3, y = PC4, color = `Tissue Type`, shape = digit)) + 
  geom_point(aes(fill = `Tissue Type`, alpha = as.character(timepoint))) + 
  scale_shape_manual(values=c(21:25)) +
  scale_color_manual(values = c("#E69F00", "#009E73", "#D55E00", "#0072B2")) + 
  scale_fill_manual(values = c("#E69F00", "#009E73", "#D55E00", "#0072B2")) +
  #scale_alpha_manual(values=c("Early", "Late")) + 
  labs(x = paste0("PC3 (", round(human_pc_eigenvalues_RNA$pct[3], 1), "%)"), 
       y = paste0("PC4 (", round(human_pc_eigenvalues_RNA$pct[4], 1), "%)"), 
       #color = "Tissue Type", 
       shape = "Digit", 
       alpha = "Stage", 
       color = "Tissue",
       fill = "Tissue",
       title = NULL ) +
  theme(legend.position = "bottom", 
        legend.box = "horizontal", 
        legend.margin = margin(), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.1, "cm")) +
  geom_point(alpha = 1, fill = NA)

human_rna_pca

#now timepoint volcano plot
#initialize function
#design_in is the DESeq model input

runDESEQ2_df_hs <- function(design_in, count_data, col_data, padj.cutoff = 0.05){
  
  #glm design
  # Again, we can include more things beside tissue, like we could look at the effect of tissue *and* age, or we could control for sex here. 
  
  #build DESeqDataSet
  rsem.in <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = design_in, tidy = T)
  #pre-filter rows with 0 or 1 read
  rsem.in <- rsem.in[ rowSums(counts(rsem.in)) >= 1, ]
  rsem_cols <-as.data.frame(row.names(rsem.in))
  colnames(rsem_cols) <- 'ensembl_gene_id'
  require(biomaRt)
  hs_ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl', host = "https://asia.ensembl.org")
  
  #make annotation data frame
  annot <- getBM(
    attributes = c(
      'hgnc_symbol',
      'ensembl_gene_id',
      'gene_biotype'),
    filters = 'ensembl_gene_id',
    values = ct$ensembl_gene_id,
    mart = hs_ensembl)
  
  rsem_annot <- merge(x= rsem_cols, y =  annot,
                      by= 'ensembl_gene_id', no.dups = TRUE, all.x=TRUE) %>% setorder(cols = "ensembl_gene_id")  
  #remove duplicate rows
  rsem_annot <- rsem_annot[!duplicated(rsem_annot$ensembl_gene_id),]
  mcols(rsem.in) <- cbind(mcols(rsem.in), rsem_annot$hgnc_symbol)
  
  
  #fix gene names in rsem.in
  #don't have to do this--when deseq reads in gene names from rsem file it duplicates name of gene. Gets rid of second name. 
  #rownames(rsem.in) <-  sapply(strsplit(c(rownames(rsem.in)), split = '_', fixed = TRUE), function(x) (x[1]))
  
  #calculate DGE likelihood ratio test
  #LRT needs reduced. Just W t
  
  #, test = "LRT", reduced = reducedest--pairwise does not need it. Or post-hoc tests can compare between pairs. 
  rsem.de <- DESeq(rsem.in)
  
  #rsem.in variable is building the dataset 
  #now model is in rsem.de
  
  #summarized experiment input
  #can set alpha level at whatever
  rsem.de.res <- results(rsem.de, alpha = 0.05)
  #
  rsem.de.res$mgi_symbol <- mcols(rsem.in)['rsem_annot$hgnc_symbol']
  
  #summary shows summary of results
  summary(rsem.de.res)
  
  #This is a way to export a table of the differentially expressed genes, ordered by p-value
  rsem.de.res.df <- rsem.de.res[order(rsem.de.res$padj),] %>% as.data.frame() %>% dplyr::filter(padj < padj.cutoff)
  
  return(rsem.de.res.df)
} #end function


timepoint_DEGs <- runDESEQ2_df_hs(design_in = formula(~timepoint), count_data = ct, col_data =  sampleMetaData, padj.cutoff = 1)
timepoint_DEGs$expressed <- "Not Significant"
timepoint_DEGs$expressed[timepoint_DEGs$log2FoldChange < -0.58 & timepoint_DEGs$padj < 0.01] <- "Upregulated early"
timepoint_DEGs$expressed[timepoint_DEGs$log2FoldChange > 0.58 & timepoint_DEGs$padj < 0.01] <- "Upregulated late"
top_DEG_n <- 10
top_upgenes <- timepoint_DEGs %>% dplyr::filter(expressed == "Upregulated late") %>% arrange(padj) %>% head(n = top_DEG_n) %>% pull(rsem_annot.hgnc_symbol)
top_downgenes <- timepoint_DEGs %>% dplyr::filter(expressed =="Upregulated early") %>% arrange(padj) %>% head(n = top_DEG_n) %>% pull(rsem_annot.hgnc_symbol)
timepoint_DEGs$delabel <- ifelse(timepoint_DEGs$rsem_annot.hgnc_symbol %in% c(top_downgenes, top_upgenes), timepoint_DEGs$rsem_annot.hgnc_symbol, NA)


stage_volcano <- ggplot(data = timepoint_DEGs, 
       aes( x = log2FoldChange, y =-log10(padj), colour = expressed, label = delabel)) + 
  geom_point(size = 0.1) +
  scale_color_manual(values = c("grey", "#00AFBB", "brown3"))+ 
  geom_hline(yintercept = -log10(0.01), 
             col = "gray", 
             linetype = 'dashed') +
  geom_vline(xintercept = c(-0.58, 0.58), col = "gray", linetype = 'dashed') +
  labs(x = expression("log"[2]*"FC"), 
       y = expression("-log"[10]*" adj. p-value"), 
       color = NULL) + 
  theme_classic() + 
geom_text_repel(max.overlaps = Inf, 
                show.legend = FALSE, 
                size = 2, 
                force = 1.5, 
                force_pull = 0.5, 
                min.segment.length = 0, color = "black") +
  theme(legend.position="bottom", 
        legend.box = "vertical", 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10), 
        legend.key.size = unit(0.1, "cm"))


#make final plot

#first arange PCR plots with shared legend
pca_part <- ggarrange(print(human_rna_pca), 
                      print(human_atac_pca),
                      nrow = 1, 
                      ncol = 2, 
                      common.legend = T, 
                      legend = "bottom" 
                      )

#first arange PCR plots with shared legend
volcano_part <- ggarrange(plotlist=list(stage_volcano, atac_stage_volcano),
                      nrow = 1, 
                      ncol = 2, 
                      common.legend = T, 
                      legend = "bottom" 
)

#combine plots to create figure
figure2 <- ggdraw() +
  draw_plot(pca_part, x = 0, y = 0.5, width = 1, height = .5) +
  draw_plot(volcano_part, x = 0, y = 0, width = 1, height = .45) +
  draw_plot_label(label = c("a", "b", "c", "d"), 
                  size = 12,
                  x = c(0, 0.5, 0, 0.5), 
                  y = c(1, 1, 0.5, 0.5)) + 
  bgcolor("white")


#save resulting figure 
ggsave(plot = figure2 + panel_border(color = "black", size = 1), 
       filename = "~/Dropbox/Autopod Paper/Autopod_Paper_Figures/Figure_2_PCA_volcano.png", 
       device = "png", 
       dpi = 300, 
       height = 150, width = 180, units = "mm")
