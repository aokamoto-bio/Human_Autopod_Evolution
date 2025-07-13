#code for figure S10 and table S3 GREAT/TF enrichments

#load packages
library(readxl)
library(GenomicRanges)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(regioneR) #for reading bed files as GRanges
library(ggnewscale)

#load libraries for TF analysis
library("monaLisa")
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(mia)

#load libraries for gene ontology enrichments
library("rGREAT")


#load custom functions for analysis
source("~/Desktop/Capellini_Lab/Weekly_Coding/RNA_Analysis_Functions.R") 


multiinter_both_df <- read.table( "~/Desktop/Capellini_Lab/Human_ATAC/corrected_multiinter_autopod.bed")
colnames(multiinter_both_df) <- c("chr", "start", "end", "num", "list", paste(c("FM1", "FM2", "FM3", "FM4", "FM5", "FP1", "FP2", "FP3", "FP4", "FP5", "HM1", "HM2", "HM3", "HM4", "HM5", "HP1", "HP2", "HP3", "HP4", "HP5"), "_Early", sep = ""), paste(c("FM1", "FM2", "FM3", "FM4", "FM5", "FP1", "FP2", "FP3", "FP4", "FP5", "HM1", "HM2", "HM3", "HM4", "HM5", "HP1", "HP2", "HP3", "HP4", "HP5"), "_Late", sep = ""))
multiinter_both_df$ID <- paste(multiinter_both_df$chr, multiinter_both_df$start, multiinter_both_df$end, sep = "_")

multiinter_df_non_autopod <- read.delim(file = "~/Desktop/Capellini_Lab/Human_ATAC/corrected_multiinter_autopod_shared_non_autopod.bed", header = F, sep = "\t")
colnames(multiinter_df_non_autopod) <- c("chr", "start", "end", "num", "list", paste(c("FM1", "FM2", "FM3", "FM4", "FM5", "FP1", "FP2", "FP3", "FP4", "FP5", "HM1", "HM2", "HM3", "HM4", "HM5", "HP1", "HP2", "HP3", "HP4", "HP5"), "E", sep = ""), paste(c("FM1", "FM2", "FM3", "FM4", "FM5", "FP1", "FP2", "FP3", "FP4", "FP5", "HM1", "HM2", "HM3", "HM4", "HM5", "HP1", "HP2", "HP3", "HP4", "HP5"), "L", sep = ""))
multiinter_df_non_autopod$ID <- paste(multiinter_df_non_autopod$chr, multiinter_df_non_autopod$start, multiinter_df_non_autopod$end, sep = "_")

multiinter_autopod_only <- multiinter_both_df %>% dplyr::filter(!ID %in% multiinter_df_non_autopod$ID) 


#"biological enrichments
#create function to do 3 things:
#1 perform GREAT analysis
#2 plot elements near genes
#3 TF binding enrichments

#first load useful datasets 

#get pwm matrix for TFs
pwms <- getMatrixSet(JASPAR2020,
                     opts = list(matrixtype = "PWM",
                                 species    = "9606"))


#TF analysis function
pairwise_TF_analysis <- function(bed1, 
                                 bed2, 
                                 ids = c("bed1", "bed2")
){
  #import bed files as rtracklayer
  range1 <- rtracklayer::import(con = bed1, format = "bed")
  range2 <- rtracklayer::import(con = bed2, format = "bed")
  #extract human sequence
  seqs1 <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, range1)
  seqs2 <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, range2)
  
  #combine sequences into single object
  comp_seqs <-  c(seqs1, seqs2)
  bins <- rep(ids, c(length(seqs1), length(seqs2)))
  bins <- factor(bins)
  #perform motif enrichment
  se <- calcBinnedMotifEnrR(seqs = comp_seqs, bins = bins,
                            pwmL = pwms, BPPARAM = BiocParallel::MulticoreParam(4))
  se_long <- merge(meltAssay(se, assay.type = "negLog10Padj", add_row_data = T, ) %>% 
                     dplyr::select(-motif.pwm, -motif.pfm) %>% 
                     as.data.frame(), meltAssay(se, assay.type = "log2enr", add_row_data = T) %>% 
                     dplyr::select(-motif.pwm, -motif.pfm) %>% 
                     as.data.frame()) %>% 
    dplyr::filter(negLog10Padj > 4)
  
  return(se_long)
}


#massive function to process multiinter comparisons

#takes a multiinter dataset
#col_end for 'add_multiinter_count_col' sub function
#file name to save plot
#great_tsv for name of tsv file containing GREAT ontology enrichments
#optional variants to specify which parts of the function should be run
categories_bio_data <- function(multiinter, col_end, filename, great_tsv, run_GREAT = T, run_TF_motif = T){
  
  #add columns to multiinter needed for determining subsets
  multiinter <- add_multiinter_count_col(dataset = multiinter, grep_term = "Early", names = "Early", column_end = col_end)
  multiinter <- add_multiinter_count_col(dataset = multiinter, grep_term = "Late", names = "Late", column_end = col_end)
  
  multiinter  <- add_multiinter_count_col(dataset = multiinter, grep_term = "P", names = "Phalangeal", column_end = col_end)
  multiinter  <- add_multiinter_count_col(dataset = multiinter, grep_term = "M", names = "Metapodial", column_end = col_end)
  
  multiinter <- add_multiinter_count_col(dataset = multiinter, grep_term = "H", names = "Hand", column_end = col_end)
  multiinter <- add_multiinter_count_col(dataset = multiinter, grep_term = "F", names = "Foot", column_end = col_end)
  
  multiinter$stage_status <-"BOTH_STAGES"
  multiinter$stage_status[which(multiinter$Early > 0 & multiinter$Late == 0)] <- "EARLY"
  multiinter$stage_status[which(multiinter$Early ==0 & multiinter$Late > 0)] <- "LATE"
  
  multiinter$region_status <-"BOTH_REGIONS"
  multiinter$region_status[which(multiinter$Metapodial > 0 & multiinter$Phalangeal == 0)] <- "METAPODIAL"
  multiinter$region_status[which(multiinter$Metapodial == 0 & multiinter$Phalangeal > 0)] <- "PHALANGEAL"
  
  multiinter$limb_status <-"SHARED"
  multiinter$limb_status[which(multiinter$Foot > 0 & multiinter$Hand == 0)] <- "FOOT"
  multiinter$limb_status[which(multiinter$Foot == 0 & multiinter$Hand > 0)] <- "HAND"
  
  #create dfs to store results 
  great_df <- data.frame()
  tf_df <- data.frame()
  
  #iterate over timepoint options
  for (s in c("BOTH_STAGES", "EARLY", "LATE")){
    stage_bed <- paste("~/Desktop/TEMP_ATAC_DIR/", s, "_temp.bed", sep ="")
    multiinter %>% 
      dplyr::filter(stage_status == s) %>% 
      dplyr::select(chr, start, end) %>% 
      write_bed(filename = stage_bed, ncol = 3)
    if(run_GREAT){
      stage_great <- bed_GREAT_enrich(bedfile = stage_bed)
      stage_great$Set <- s
      #store results 
      great_df <- rbind(great_df, stage_great)
    } #end GREAT analysis
  } #end loop over stages
  
  #iterate over region options
  for (r in c("BOTH_REGIONS", "METAPODIAL", "PHALANGEAL")){
    region_bed <- paste("~/Desktop/TEMP_ATAC_DIR/", r, "_temp.bed", sep ="")
    multiinter %>% 
      dplyr::filter(region_status == r) %>% 
      dplyr::select(chr, start, end) %>% 
      write_bed(filename = region_bed, ncol = 3)
    if(run_GREAT){
      region_great <- bed_GREAT_enrich(bedfile = region_bed)
      region_great$Set <- r
      #store results 
      great_df <- rbind(great_df, region_great)
    } #end GREAT analysis
  } #end loop over regions
  
  #iterate over region options
  for (l in c("SHARED", "HAND", "FOOT")){
    limb_bed <- paste("~/Desktop/TEMP_ATAC_DIR/", l, "_temp.bed", sep ="")
    multiinter %>% 
      dplyr::filter(limb_status == l) %>% 
      dplyr::select(chr, start, end) %>% 
      write_bed(filename = limb_bed, ncol = 3)
    if(run_GREAT){
      limb_great <- bed_GREAT_enrich(bedfile = limb_bed)
      limb_great$Set <- l
      #store results 
      great_df <- rbind(great_df, limb_great)
    } #end GREAT analysis
  } #end loop over regions
  
  #now for pairwise portion
  stages <- c("BOTH_STAGES", "EARLY", "LATE")
  regions <- c("BOTH_REGIONS", "METAPODIAL", "PHALANGEAL")
  limbs <- c("SHARED", "HAND", "FOOT")
  #regions <-
  #limbs <- 
  pairwise_df <- data.frame(comp_id = 1:3, 
                            stages1 = c(stages[1], stages[1], stages[2]),
                            stages2 = c(stages[2], stages[3], stages[3]),
                            regions1 = c(regions[1], regions[1], regions[2]),
                            regions2 = c(regions[2], regions[3], regions[3]),
                            limbs1 = c(limbs[1], limbs[1], limbs[2]),
                            limbs2 = c(limbs[2], limbs[3], limbs[3]))
  
  
  for (p in pairwise_df$comp_id){
    if(run_TF_motif){
      #stage
      stage_tf <- pairwise_TF_analysis(bed1 = paste("~/Desktop/TEMP_ATAC_DIR/", pairwise_df$stages1[p], "_temp.bed", sep =""), 
                                       bed2 = paste("~/Desktop/TEMP_ATAC_DIR/", pairwise_df$stages2[p], "_temp.bed", sep = ""), 
                                       ids = c(pairwise_df$stages1[p], pairwise_df$stages2[p]))
      stage_tf$comparison <- paste(pairwise_df$stages1[p], "vs", pairwise_df$stages2[p], sep = "_")
      stage_tf$type <- "STAGE"
      tf_df <- rbind(tf_df, stage_tf)
      
      #region
      region_tf <- pairwise_TF_analysis(bed1 = paste("~/Desktop/TEMP_ATAC_DIR/", pairwise_df$regions1[p], "_temp.bed", sep =""), 
                                        bed2 = paste("~/Desktop/TEMP_ATAC_DIR/", pairwise_df$regions2[p], "_temp.bed", sep = ""), 
                                        ids = c(pairwise_df$regions1[p], pairwise_df$regions2[p]))
      region_tf$comparison <- paste(pairwise_df$regions1[p], "vs", pairwise_df$regions2[p], sep = "_")
      region_tf$type <- "REGION"
      tf_df <- rbind(tf_df, region_tf)
      
      #limb
      limb_tf <- pairwise_TF_analysis(bed1 = paste("~/Desktop/TEMP_ATAC_DIR/", pairwise_df$limbs1[p], "_temp.bed", sep =""), 
                                      bed2 = paste("~/Desktop/TEMP_ATAC_DIR/", pairwise_df$limbs2[p], "_temp.bed", sep = ""), 
                                      ids = c(pairwise_df$limbs1[p], pairwise_df$limbs2[p]))
      limb_tf$comparison <- paste(pairwise_df$limbs1[p], "vs", pairwise_df$limbs2[p], sep = "_")
      limb_tf$type <- "limb"
      tf_df <- rbind(tf_df, limb_tf)
      
    } #end if run_TF_motif
  } #end loop over pairwise stages 
  
  #make the plots
  
  #timepoint
  early_late_genes_plot <- compare_features_near_gene(feature_bed1 = "~/Desktop/TEMP_ATAC_DIR/EARLY_temp.bed", feature_bed2 = "~/Desktop/TEMP_ATAC_DIR/LATE_temp.bed",
                                                      threshold_type = "min", 
                                                      autopod_threshold = 10, 
                                                      comparison = "EARLYvLATE", 
                                                      return_opt = "plot",
                                                      xlab = "Early REs",
                                                      ylab = "Late REs")
  
  both_stages_early_genes_plot <- compare_features_near_gene(feature_bed1 = "~/Desktop/TEMP_ATAC_DIR/BOTH_STAGES_temp.bed", feature_bed2 = "~/Desktop/TEMP_ATAC_DIR/EARLY_temp.bed",
                                                             threshold_type = "min", 
                                                             autopod_threshold = 10, 
                                                             comparison = "BOTH_STAGESvEARLY", 
                                                             return_opt = "plot",
                                                             xlab = "Both stages REs",
                                                             ylab = "Early REs")
  
  both_stages_late_genes_plot <- compare_features_near_gene(feature_bed1 = "~/Desktop/TEMP_ATAC_DIR/BOTH_STAGES_temp.bed", feature_bed2 = "~/Desktop/TEMP_ATAC_DIR/LATE_temp.bed",
                                                            threshold_type = "min", 
                                                            autopod_threshold = 10, 
                                                            comparison = "BOTH_STAGESvLATE", 
                                                            return_opt = "plot",
                                                            xlab = "Both stages REs",
                                                            ylab = "Late REs")
  
  #limb
  hand_foot_genes_plot <- compare_features_near_gene(feature_bed1 = "~/Desktop/TEMP_ATAC_DIR/HAND_temp.bed", feature_bed2 = "~/Desktop/TEMP_ATAC_DIR/FOOT_temp.bed",
                                                     threshold_type = "min", 
                                                     autopod_threshold = 10, 
                                                     comparison = "HANDvFOOT", 
                                                     return_opt = "plot",
                                                     xlab = "Hand REs",
                                                     ylab = "Foot REs")
  
  shared_foot_genes_plot <- compare_features_near_gene(feature_bed1 = "~/Desktop/TEMP_ATAC_DIR/SHARED_temp.bed", feature_bed2 = "~/Desktop/TEMP_ATAC_DIR/FOOT_temp.bed",
                                                       threshold_type = "min", 
                                                       autopod_threshold = 10, 
                                                       comparison = "SHAREDvFOOT", 
                                                       return_opt = "plot",
                                                       xlab = "Shared REs",
                                                       ylab = "Foot REs")
  
  shared_hand_genes_plot <- compare_features_near_gene(feature_bed1 = "~/Desktop/TEMP_ATAC_DIR/SHARED_temp.bed", feature_bed2 = "~/Desktop/TEMP_ATAC_DIR/HAND_temp.bed",
                                                       threshold_type = "min", 
                                                       autopod_threshold = 10, 
                                                       comparison = "SHAREDvHAND",
                                                       return_opt = "plot",
                                                       xlab = "Shared REs",
                                                       ylab = "Hand REs")
  #region
  metapodial_phalangeal_genes_plot <- compare_features_near_gene(feature_bed1 = "~/Desktop/TEMP_ATAC_DIR/METAPODIAL_temp.bed", feature_bed2 = "~/Desktop/TEMP_ATAC_DIR/PHALANGEAL_temp.bed",
                                                                 threshold_type = "min", 
                                                                 autopod_threshold = 10, 
                                                                 comparison = "METAPODIALvPHALANGEAL", 
                                                                 return_opt = "plot",
                                                                 xlab = "Metapodial REs",
                                                                 ylab = "Phalangeal REs")
  
  both_regions_metapodial_genes_plot <- compare_features_near_gene(feature_bed1 = "~/Desktop/TEMP_ATAC_DIR/BOTH_REGIONS_temp.bed", feature_bed2 = "~/Desktop/TEMP_ATAC_DIR/METAPODIAL_temp.bed",
                                                                   threshold_type = "min", 
                                                                   autopod_threshold = 10, 
                                                                   comparison = "BOTH_REGIONSvMETAPODIAL", 
                                                                   return_opt = "plot",
                                                                   xlab = "Both regions REs",
                                                                   ylab = "Metapodial REs")
  
  both_regions_phalangeal_genes_plot <- compare_features_near_gene(feature_bed1 = "~/Desktop/TEMP_ATAC_DIR/BOTH_REGIONS_temp.bed", feature_bed2 = "~/Desktop/TEMP_ATAC_DIR/PHALANGEAL_temp.bed",
                                                                   threshold_type = "min", 
                                                                   autopod_threshold = 10, 
                                                                   comparison = "BOTH_STAGESvPHALANGEAL",
                                                                   return_opt = "plot",
                                                                   xlab = "Both regions REs",
                                                                   ylab = "Phalangeal REs")
  
  
  
  #make final plot
  plot_h <- 0.33
  plot_w <- 0.33
  
  gene_total_plot <- ggdraw() +
    #hand versus foot
    draw_plot(hand_foot_genes_plot, 
              x = 0, y = (0.01 + plot_h*2), 
              width = plot_w, height = plot_h) +
    draw_plot(shared_hand_genes_plot, 
              x = plot_w, y = (0.01 + plot_h*2), 
              width = plot_w, height = plot_h) +
    draw_plot(shared_foot_genes_plot, 
              x = plot_w*2, y = (0.01 + plot_h*2), 
              width = plot_w, height = plot_h) +
    #phalanges versus metapodials
    draw_plot(metapodial_phalangeal_genes_plot, 
              x = 0, y = (0.01 + plot_h), 
              width = plot_w, height = plot_h) +
    draw_plot(both_regions_metapodial_genes_plot, 
              x = plot_w, y = (0.01 + plot_h), 
              width = plot_w, height = plot_h) +
    draw_plot(both_regions_phalangeal_genes_plot, 
              x = plot_w*2, y = (0.01 + plot_h), 
              width = plot_w, height = plot_h) +
    #early versus late
    draw_plot(early_late_genes_plot, x = 0, y = 0.01, width = plot_w, height = plot_h) +
    draw_plot(both_stages_early_genes_plot, x = plot_w, y = 0.01, width = plot_w, height = plot_h) +
    draw_plot(both_stages_late_genes_plot, x = plot_w*2, y = 0.01, width = plot_w, height = plot_h) +
    
    draw_plot_label(label = c("a", "b", "c", "d", "e", "f", "g", "h","i"), 
                    size = 15,
                    x = rep(c(0, 0.33, 0.66), 3),                    y = c(rep(1, 3), rep(0.67, 3), rep(0.34, 3))) + 
    bgcolor("white")
  
  ggsave(plot = gene_total_plot + panel_border(color = "black", size = 1), 
         filename = "~/Dropbox/Autopod Paper/Autopod_Paper_Figures/fig_S10_gene_total_plot.png", 
         device = "png", dpi = 300, height = 180, width = 180, units = "mm")
  
  
  
  #final steps
  
  #save great results to file is desired 
  if(run_GREAT){
    write.table(x = great_df, file = great_tsv, col.names = T, row.names = F, sep = "\t")
  }
  if(run_TF_motif){
    tf_df %>% arrange(type, comparison) %>% write.table(col.names = T, row.names = F, sep = "\t", "~/Desktop/Capellini_Lab/Human_ATAC/categories_TF_enrichments.tsv")
    return(tf_df %>% arrange(type, comparison))
  }
}

multiinter_autopod_only_great <- categories_bio_data(multiinter = multiinter_autopod_only,
                                                     col_end = 45, 
                                                     filename = "~/Dropbox/Autopod Paper/Autopod_Paper_Figures/fig_S10_gene_total_plot.png", 
                                                     great_tsv = "~/Desktop/Capellini_Lab/Human_ATAC/autopod_only_categories_GREAT.tsv", run_GREAT = T, 
                                                     run_TF_motif = T)