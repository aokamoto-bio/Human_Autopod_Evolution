#create supplementary figure showing human autopod ATAC-seq WGCNA results 
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

#early modules
DAR_plt1 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_early_DAR_cluster_1.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
DAR_plt2 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_early_DAR_cluster_2.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
DAR_plt3 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_early_DAR_cluster_3.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
DAR_plt4 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_early_DAR_cluster_4.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
DAR_plt5 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_early_DAR_cluster_5.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
DAR_plt6 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_early_DAR_cluster_6.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
DAR_plt7 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_early_DAR_cluster_7.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 

#load cluster plots
#early modules
early_plt1 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_early_RNA_cluster_1.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
early_plt2 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_early_RNA_cluster_2.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
early_plt3 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_early_RNA_cluster_3.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 
early_plt4 <- image_read("~/Desktop/Autopod_Skeleton_Figures/human_early_RNA_cluster_4.png") %>%
  image_crop(geometry = "290x270+40+30") %>%
  image_ggplot() 

#now plot the correlation
#create RNA module beds
RNA_modules <- read.table(file = "~/Desktop/Capellini_Lab/Human Autopod RNA/human_early_gene_modules.tsv", header = T, sep = "\t")


#create DAR module beds
DAR_modules <- read.table(file = "~/Desktop/Capellini_Lab/Human_ATAC/human_early_DAR_modules.tsv", header = T, sep = "\t") %>% 
  separate_wider_delim(peak, ":", names = c("chr", "range")) %>% 
  separate_wider_delim(range, "-", names = c("start", "end")) %>% 
  dplyr::select(chr, start, end, cluster)

#columns(org.Hs.eg.db) # returns list of available keytypes
get_cluster_entrez_ids <- function(cluster_df, clust) {
  results <- mapIds(org.Hs.eg.db,
                    keys = cluster_df %>% 
                      dplyr::filter(cluster == clust) %>% 
                      pull(ensembl_gene_id), #Column containing Ensembl gene ids
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
  return(results)
}

rna_modules_list <- list(
  "b1" = get_cluster_entrez_ids(cluster_df = RNA_modules, clust = 1), 
  "b2" = get_cluster_entrez_ids(cluster_df = RNA_modules, clust = 2), 
  "b3" = get_cluster_entrez_ids(cluster_df = RNA_modules, clust = 3), 
  "b4" = get_cluster_entrez_ids(cluster_df = RNA_modules, clust = 4)
)

#gr = randomRegions(genome = "hg38")
mod_enrich_df <- data.frame()

for (DAR_mod in 1:max(DAR_modules$cluster)){
  dar_bed <- rtracklayer::import(paste("~/Desktop/Capellini_Lab/Human_ATAC/human_late_DAR_cluster_", DAR_mod, "_hg38.bed", sep = ""), format = "BED")
  res = great(dar_bed, gene_sets = rna_modules_list, "hg38")
  res2 <- getEnrichmentTables(res)
  #only save results if some enrichments found. 
  if(nrow(res2 > 0)){
    res2$DAR_mod <- DAR_mod
    mod_enrich_df <-rbind(mod_enrich_df, res2)
  }
}

mod_enrich_df2 <- mod_enrich_df%>% 
  dplyr::filter(p_adjust < 0.05)

#for visualization purposes, the minimum p-value was set as 1e-10
mod_enrich_df2$p_adjust[which(mod_enrich_df2$p_adjust == 0)] <- 1e-10
#factor plotting variables
mod_enrich_df2$id <- factor(mod_enrich_df2$id, levels = paste("b", 1:4, sep = ""))
mod_enrich_df2$DAR_mod <- paste("a", mod_enrich_df2$DAR_mod, sep ="")
mod_enrich_df2$DAR_mod <- factor(mod_enrich_df2$DAR_mod, levels = paste("a", 1:7, sep = ""))

early_wgcna_corr_great <- ggplot(data = mod_enrich_df2, 
                                aes(x = DAR_mod, y = id, 
                                    color = -log10(p_adjust), 
                                    size = fold_enrichment)) + 
  geom_point() +
  theme_bw() + 
  labs(x = "DAR Modules", 
       y = "DEG Modules", 
       size = "Fold enrichment", color = 
         "-log10(Padj)") + 
  scale_x_discrete(drop=FALSE) + 
  scale_y_discrete(drop=FALSE) +
  scale_color_gradient(low = "gray", high = "darkgreen") +
  scale_size_continuous(breaks = c(1.2, 1.6, 2))

#now for TF enrichment
RNA_modules <- read.table(file = "~/Desktop/Capellini_Lab/Human Autopod RNA/human_early_gene_modules.tsv", header = T, sep = "\t")

early_TF_df2 <- read.table(sep = "\t", "~/Desktop/Capellini_Lab/Human_ATAC/early_module_TF_enrichments.tsv", header = T) %>% 
  dplyr::filter(SampleID != "all_peaks")
early_TF_df3 <- merge(early_TF_df2, y = RNA_modules, by.x = "motif.name", by.y ="hgnc_symbol" ) %>% 
  dplyr::select(SampleID, cluster) %>% 
  group_by(SampleID, cluster) %>% 
  summarize(n = n())

early_TF_per_cluster <- merge(early_TF_df2, y = RNA_modules, by.x = "motif.name", by.y ="hgnc_symbol" ) %>% 
  dplyr::select(SampleID, cluster) %>% 
  group_by(cluster) %>% 
  summarize(total_n = n())


early_TF_df4 <- merge(early_TF_df3, early_TF_per_cluster) %>% mutate(perc = n/total_n)

early_TF_df4$cluster <- paste("b", early_TF_df4$cluster, sep = "")
early_TF_df4$SampleID <- gsub(pattern = "early_", replacement = "a", early_TF_df4$SampleID)
early_tf_corr <- ggplot(data = early_TF_df4 %>% dplyr::filter(cluster != "b0"), 
                        aes(y = cluster, x = SampleID, 
                            fill = perc*100)) + 
  geom_tile() +
  theme_classic() + 
  labs(y = "DEG Modules", 
       x = "DAR Modules", 
       fill = "% TFs") + 
  scale_fill_gradient(low = "white", high = "darkred")


#make final plot
#make final plot
skel_h <- 0.25
skel_w <- 1/7

human_early_plot <- ggdraw() +
  draw_plot(DAR_plt1, x = 0.0, y = skel_h*3, width = skel_w, height = skel_h) +
  draw_plot(DAR_plt2, x = skel_w, y = skel_h*3, width = skel_w, height = skel_h) +
  draw_plot(DAR_plt3, x = skel_w*2, y = skel_h*3, width = skel_w, height = skel_h) +
  draw_plot(DAR_plt4, x = skel_w*3, y = skel_h*3, width = skel_w, height = skel_h) +
  draw_plot(DAR_plt5, x = skel_w*4, y = skel_h*3, width = skel_w, height = skel_h) +
  draw_plot(DAR_plt6, x = skel_w*5, y = skel_h*3, width = skel_w, height = skel_h) +
  draw_plot(DAR_plt7, x = skel_w*6, y = skel_h*3, width = skel_w, height = skel_h) +
  draw_plot(early_plt1, x = 0, y = skel_h*2, width = skel_w, height = skel_h) +
  draw_plot(early_plt2, x = skel_w, y = skel_h*2, width = skel_w, height = skel_h) +
  draw_plot(early_plt3, x = skel_w*2, y = skel_h*2, width = skel_w, height = skel_h) +
  draw_plot(early_plt4, x = skel_w*3, y = skel_h*2, width = skel_w, height = skel_h) +
  draw_plot(early_wgcna_corr_great, x = 0, y = 0, width = 0.6, height = skel_h*2) +
  draw_plot(early_tf_corr, x = 0.6, y = 0, width = 0.4, height = skel_h*2) +
  draw_plot_label(label = c(paste("a", 1:7, sep = ""), paste("b", 1:4, sep = ""), "c", "d"), 
                  size = 12,
                  x = c((0:6*skel_w), (0:3*skel_w), 0, 0.6), 
                  y = c(rep(1, times = 7), rep(1-skel_h, times = 4), 1-skel_h*2, 1-skel_h*2)
  ) + #now add gene names
  draw_plot_label(
    label = c("synapse org.", "ECM", "fat cell diff.", "skin dev."),
    #fontface = "italic",
    size = 8,
    x = c((0:3*(skel_w)) + 0.05),
    y = c(rep(1-skel_h, times =4)),
    hjust = 0
  ) + #add GO enrichments
  draw_plot_label(
    label = c("cartilage dev.", "limb dev.", "collagen", "nitrogen resp.", "abnorm. limb", "steroid resp.", "skeletal dev."),
    size = 8,
    x = c((0:6*(skel_w)) + 0.05),
    y = c(rep(1, times =7)),
    hjust = 0, 
  ) + 
  bgcolor("white")

human_early_plot 

ggsave(plot = human_early_plot + panel_border(color = "black", size = 1), 
       filename = "~/Dropbox/Autopod Paper/Autopod_Paper_Figures/Fig3_Human_wgcna_early.png", 
       device = "png", dpi = 300, height = 144, width = 180, units = "mm")

