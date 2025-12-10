#create plot showing enrichment of evolutionary signal in DAR modules
#Alexander Okamoto
#December 10, 2025

library("LOLA")
library(rtracklayer)
library(tidyverse)

#create lola database at "~/Desktop/Capellini_Lab/Autopod_LOLA"
system("cp ~/Desktop/Capellini_Lab/Structural_Variant_Beds/human_specific_inversions_hg38.bed ~/Desktop/Capellini_Lab/Autopod_LOLA/regions")
system("cp ~/Desktop/Capellini_Lab/Structural_Variant_Beds/human_SDRs_hg38.bed  ~/Desktop/Capellini_Lab/Autopod_LOLA/regions")
system("cp ~/Desktop/Capellini_Lab/HARs/HARs_AO_2024_merged_hg38.bed  ~/Desktop/Capellini_Lab/Autopod_LOLA/regions")
system("cp ~/Desktop/Capellini_Lab/Capellini_Lab_Peak_Sets/Other_Peak_sets/hCONDELs_Xue_McLean_hg38.bed  ~/Desktop/Capellini_Lab/Autopod_LOLA/regions")
system("cp ~/Desktop/Capellini_Lab/Structural_Variant_Beds/HAQERS_T2T_hg38.bed  ~/Desktop/Capellini_Lab/Autopod_LOLA/regions")
system("cp ~/Desktop/Capellini_Lab/Human_ATAC/corrected_multiinter_autopod_human_vars.bed  ~/Desktop/Capellini_Lab/Autopod_LOLA/regions")
system("cp ~/Desktop/Capellini_Lab/Human_ATAC/corrected_multiinter_autopod_chimp_vars.bed  ~/Desktop/Capellini_Lab/Autopod_LOLA/regions")
system("cp ~/Desktop/Capellini_Lab/Human_ATAC/corrected_multiinter_autopod_gorilla_vars.bed  ~/Desktop/Capellini_Lab/Autopod_LOLA/regions")
system("cp ~/Desktop/290k_MPRA_Cluster_Results/beds/290k_all_emvars_hg38.bed  ~/Desktop/Capellini_Lab/Autopod_LOLA/regions")
system("cp ~/Desktop/Capellini_Lab/Structural_Variant_Beds/human_specific_inversionsflank_500000_hg38.bed ~/Desktop/Capellini_Lab/Autopod_LOLA/regions")
system("cp ~/Desktop/Capellini_Lab/Structural_Variant_Beds/human_SDRsflank_500000_hg38.bed ~/Desktop/Capellini_Lab/Autopod_LOLA/regions")
system("cp ~/Desktop/Capellini_Lab/Structural_Variant_Beds/human_SDRsflank_500000_hg38.bed ~/Desktop/Capellini_Lab/Autopod_LOLA/regions")
system("bedtools intersect -a ~/Desktop/Capellini_Lab/Human_ATAC/Human_Autopod_IDR_peaks/Brain_Filtered/human_all_merge_IDR_0.05_brain_filtered.bed -b ~/Desktop/Capellini_Lab/Mouse_ATAC/mouse_multiinter_simple_hg38.bed > ~/Desktop/Capellini_Lab/Autopod_LOLA/regions/human_mouse_overlap.bed")

#load datasets for LOLA
human_evolDB = loadRegionDB(dbLocation = "~/Desktop/Capellini_Lab/", collections = "Autopod_LOLA")
all_peaks <- import("/Users/alexanderokamoto/Desktop/Capellini_Lab/Capellini_Lab_Peak_Sets/All_Capellini_Sets/human_all_merge_IDR_0.05_brain_filtered.bed")

DAR_modules_list <- list(
  "all_DARs" = import("/Users/alexanderokamoto/Desktop/Capellini_Lab/Human_ATAC/human_DARs_hg38.bed"),
  "early0" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_early_DAR_cluster_0_hg38.bed"), 
  "early1" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_early_DAR_cluster_1_hg38.bed"), 
  "early2" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_early_DAR_cluster_2_hg38.bed"), 
  "early3" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_early_DAR_cluster_3_hg38.bed"), 
  "early4" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_early_DAR_cluster_4_hg38.bed"), 
  "early5" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_early_DAR_cluster_5_hg38.bed"), 
  "early6" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_early_DAR_cluster_6_hg38.bed"), 
  "early7" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_early_DAR_cluster_7_hg38.bed"), 
  "late0" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_late_DAR_cluster_0_hg38.bed"), 
  "late1" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_late_DAR_cluster_1_hg38.bed"), 
  "late2" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_late_DAR_cluster_2_hg38.bed"), 
  "late3" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_late_DAR_cluster_3_hg38.bed"), 
  "late4" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_late_DAR_cluster_4_hg38.bed"), 
  "late5" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_late_DAR_cluster_5_hg38.bed"), 
  "late6" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_late_DAR_cluster_6_hg38.bed"), 
  "late7" = import("~/Desktop/Capellini_Lab/Human_ATAC/human_late_DAR_cluster_7_hg38.bed")
)

#run LOLA and filter results
locResults = runLOLA(userSets = DAR_modules_list, userUniverse = all_peaks, regionDB = human_evolDB, cores=1)
head(locResults)

locResults2 <- locResults %>% 
  dplyr::filter(qValue< 0.05)


#for visualization purposes, the minimum p-value was set as 1e-10
locResults2$qValue[which(locResults2$qValue == 0)] <- 1e-10
#factor plotting variables
#mod_enrich_df2$id <- factor(mod_enrich_df2$id, levels = paste("b", 1:13, sep = ""))
#mod_enrich_df2$DAR_mod <- paste("a", mod_enrich_df2$DAR_mod, sep ="")
#mod_enrich_df2$DAR_mod <- factor(mod_enrich_df2$DAR_mod, levels = paste("a", 1:10, sep = ""))
locResults2$description <-factor(locResults2$description, 
                                 levels = c("human_SNVs", "chimp_SNVs", "gorilla_SNVs", "mouse_conserved", "HARs","HAQERs", "hCONDELs",  "SDR_flank", "SDRs", "inversion_flank", "inversions", "MPRA_Diff_Act"))

human_evol_mod_enrich <- ggplot(data = locResults2, 
                                 aes(x = userSet, y = description, 
                                     color = -log10(qValue), 
                                     size = oddsRatio)) + 
  geom_point() +
  theme_bw() + 
  labs(x = "DAR Module", 
       y = "Enrichment Set", 
       size = "Odds Ratio", color = 
         expression("-log"[10]*" adj. p-value")) + 
  scale_y_discrete(drop=FALSE, 
                   labels = c("chimp_SNVs" = "Chimp SNVs",
                              "gorilla_SNVs" = "Gorilla SNVs",
                              "human_SNVs" = "Human SNVs", 
                              "MPRA_Diff_Act" = "Differentially active (MPRA)",
                              "inversion_flank" = "Inversion flanks",
                              "inversions" = "Inversions",
                              "SDR_flank" = "SDR flanks",
                              "mouse_conserved" = "Human-mouse conserved")) + 
  scale_x_discrete(drop=FALSE, 
                   labels = c("all_DARs" = "All DARs", 
                              "early2" = "Early 2", 
                              "early4" = "Early 4", 
                              "early7" = "Early 7",
                              "late2" = "Late 2", 
                              "late3" = "Late 3", 
                              "late4" = "Late 4")) +
  scale_size_continuous(breaks = c(1, 2, 3)) +
  scale_color_gradient(low = "gray", high = "darkgreen")

ggsave(plot = human_evol_mod_enrich + panel_border(color = "black", size = 1), 
       filename = "~/Dropbox/Autopod Paper/Autopod_Paper_Figures/Fig5_evol_enrichments.png", 
       device = "png", dpi = 300, height = 90, width = 180, units = "mm")
