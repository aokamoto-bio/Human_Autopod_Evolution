#calculate scores for genomic positions fixed along the human or chimp lineages

#load packages
require(tidyverse)
require(BSgenome.Hsapiens.UCSC.hg38)
require(BSgenome.Ptroglodytes.UCSC.panTro6)
require(BSgenome.Ggorilla.UCSC.gorGor6)
require(Biostrings)
require(GenomicRanges)

#to avoid positions with poor sequence alignment between the species, we start with the good set:
fixed_bp_table <- read.delim(file = "~/Desktop/Capellini_Lab/Human_ATAC/fixed_bps_table.txt", header = T, sep = "\t")
fixed_bp_table %>% separate(col = ID, sep = "_", into = c('chr', 'start', 'end'), remove = F) %>% dplyr::select(chr, start, end, ID) %>% write_bed(filename = "~/Desktop/TEMP_ATAC_DIR/fixed_bp_table_hg38.bed")

human_fixed_ids <- fixed_bp_table %>% dplyr::filter(human_fixed > 0) %>% pull(ID)
chimp_fixed_ids <- fixed_bp_table %>% dplyr::filter(chimp_fixed > 0) %>% pull(ID)


bedfile <- "~/Desktop/TEMP_ATAC_DIR/fixed_bp_table_hg38.bed"

#load genomes & chains
chimp_genome <- BSgenome.Ptroglodytes.UCSC.panTro6
chimp_chain <- "~/Desktop/Polymorphisms_Project/LiftOver/Chains/hg38ToPanTro6.over.chain"
gorilla_genome <- BSgenome.Ggorilla.UCSC.gorGor6
gorilla_chain <- "~/Desktop/Polymorphisms_Project/LiftOver/Chains/hg38.gorGor6.all.chain"

#remove human & chimp SNP positions/repeat elements
system(paste("bedtools subtract -a ", bedfile, " -b ~/Desktop/Capellini_Lab/Capellini_Lab_Peak_Sets/Other_Peak_sets/ENCODE_BLACKLIST_hg38.bed | ",
             "bedtools subtract -a stdin -b ~/Desktop/Capellini_Lab/Weekly_Coding/variants_dbSNP151_autopod_hg38.bed | bedtools subtract -a stdin -b ~/Desktop/Capellini_Lab/Weekly_Coding/greatape_autopod_SNPs_hg38.bed > ",
             "~/Desktop/TEMP_ATAC_DIR/temp_human_ape_bps.bed", 
             sep = ""
))

#load bedfile and start processing
bed_table <- read.delim("~/Desktop/TEMP_ATAC_DIR/temp_human_ape_bps.bed", sep = "\t", header = F)
#add required column names for easier calling later on
colnames(bed_table)[1:4] <- c("chr", "start", "end", "ID")
bed_table <- bed_table %>% 
  dplyr::filter(chr %in% c(paste("chr", 1:22, sep = ""), "chrX", "chrY"))

#get human DNA sequences
bed_table$human_seqs <- bed_table %>% 
  dplyr::select(chr, start, end) %>% 
  dplyr::filter(chr %in% c(paste("chr", 1:22, sep = ""), "chrX", "chrY")) %>% 
  makeGRangesFromDataFrame() %>% 
  getSeq(x=BSgenome.Hsapiens.UCSC.hg38) %>% 
  as.character()
#get lengths of sequence to expect
bed_table$human_length <- bed_table$end - bed_table$start + 1
#create stable ID value
bed_table$ID2 <- paste(bed_table$chr, bed_table$start, bed_table$end, sep = "")

#get ape DNA sequences

#get granges of liftover object
chimp_granges <- bed_table %>% 
  dplyr::filter(chr %in% c(paste("chr", 1:22, sep = ""), "chrX", "chrY")) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% 
  liftOver(chain = import.chain(chimp_chain)) %>% 
  unlist()

#make a dataframe to store extracted bits of sequence
chimp_info <-  chimp_granges %>% as.data.frame()
chimp_info$seq <- chimp_granges %>% getSeq(x = chimp_genome)  %>% 
  as.character()
chimp_concat_info <- chimp_info %>% 
  dplyr::select(seq, ID2, width) %>% 
  group_by(ID2) %>% 
  mutate(chimp_seq = paste0(seq, collapse = ""), chimp_length = sum(width)) %>%
  dplyr::select(ID2, chimp_seq, chimp_length) %>% 
  unique()

#get granges of liftover object
gorilla_granges <- bed_table %>% 
  dplyr::filter(chr %in% c(paste("chr", 1:22, sep = ""), "chrX", "chrY")) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% 
  liftOver(chain = import.chain(gorilla_chain)) %>% 
  unlist()

#make a dataframe to store extracted bits of sequence
gorilla_info <-  gorilla_granges %>% as.data.frame()
gorilla_info$seq <- gorilla_granges %>% getSeq(x = gorilla_genome)  %>% 
  as.character()
gorilla_concat_info <- gorilla_info %>% 
  dplyr::select(seq, ID2, width) %>% 
  group_by(ID2) %>% 
  mutate(gorilla_seq = paste0(seq, collapse = ""), gorilla_length = sum(width)) %>% 
  dplyr::select(ID2, gorilla_seq, gorilla_length) %>% 
  unique()

#combined ape and human bedtables
combined_table <- merge(bed_table, chimp_concat_info, by = 'ID2') %>% merge(gorilla_concat_info, by = 'ID2')

#subset table to regions of equal length between species and identify human and chimp fixed bps
reduced_table <- combined_table %>% 
  as_tibble() %>%
  dplyr::filter(human_length == chimp_length & human_length == gorilla_length) %>% 
  group_by(ID2) %>% 
  rowwise() %>% 
  mutate(bps_difs = compare_DNA_polarized(human_seqs, chimp_seq, gorilla_seq), chimp_bps_difs = compare_DNA_polarized(chimp_seq, human_seqs, gorilla_seq)) 

#remove segments with poor alignments
reduced_table$per_bp_difs <- reduced_table$bps_difs/reduced_table$human_length
reduced_table$chimp_per_bp_difs <- reduced_table$chimp_bps_difs/reduced_table$human_length
reduced_table2 <- reduced_table %>% dplyr::filter(per_bp_difs < 0.75) %>% dplyr::filter(chimp_per_bp_difs < 0.75) 

#get positions
#compare polarized DNA sequences and return the positions with fixed differences
compare_DNA_polarized_positions <- function(x, y, z){
  count <- 0
  x_bps <- unlist(str_split(x, pattern = ""))
  y_bps <- unlist(str_split(y, pattern = ""))
  z_bps <- unlist(str_split(z, pattern = ""))
  if(length(x_bps) == length(y_bps) && length(x_bps) == length(z_bps)){
    pos <- which(sapply(seq(length(x_bps)),
                                 function(i){
                                   x_bps[i] != y_bps[i] & y_bps[i] == z_bps[i]
                                 })) %>% paste(collapse = ",")
    
    if(nchar(pos) == 0){
      pos <- "none"
    }
    return(pos)
  }
}

#clean up large object that are no longer needed
rm(chimp_granges)
rm(gorilla_granges)
rm(chimp_concat_info)
rm(chimp_info)
rm(gorilla_concat_info)
rm(gorilla_info)
rm(bed_table)


#filter by same sequence length, remove any regions which match even less than expected by chance
#create human and chimp subsets

reduced_table_human <- reduced_table2 %>% 
  ungroup() %>% 
  dplyr::filter(ID %in% human_fixed_ids) %>% 
  rowwise() %>% 
  dplyr::mutate(bps_difs = compare_DNA_polarized_positions(human_seqs, chimp_seq, gorilla_seq)) %>% 
  dplyr::filter(bps_difs != "none")

reduced_table_chimp <- reduced_table2 %>% 
  ungroup() %>% 
  dplyr::filter(ID %in% chimp_fixed_ids) %>% 
  rowwise() %>% 
  dplyr::mutate(bps_difs = compare_DNA_polarized_positions(chimp_seq, human_seqs, gorilla_seq)) %>% 
  dplyr::filter(bps_difs != "none")

#expand to get one row per fixed bp difference
chimp_fixed_df <- reduced_table_chimp %>% 
  ungroup() %>% 
  dplyr::select(chr, start, bps_difs, ID) %>% 
  separate_longer_delim(bps_difs, delim = ",") %>% 
  mutate(end = start + as.numeric(bps_difs) - 1) %>% 
  mutate(start = end - 1) %>% 
  dplyr::select(chr, start, end, ID) 

human_fixed_df <- reduced_table_human %>% 
  ungroup() %>% 
  dplyr::select(chr, start, bps_difs, ID) %>% 
  separate_longer_delim(bps_difs, delim = ",") %>% 
  mutate(end = start + as.numeric(bps_difs) - 1) %>% 
  mutate(start = end - 1) %>% 
  dplyr::select(chr, start, end, ID) 
  
#create bed files of all fixed positions for chimp
write_bed(chimp_fixed_df, filename = "~/Desktop/Capellini_Lab/Human_ATAC/chimp_fixed_bps_hg38.bed")

#intersect positions with bed file containing per base pair phyloP scores
#this was done using a small private cluster
#sort -k1,1 -k2,2n chimp_fixed_bps_hg38.bed | bedtools intersect -a stdin -b /cold/aokamoto/241-mammalian-2020v2.bed -sorted -wa -wb  > chimp_fixed_bps_phyloP_hg38.bed 

#create bed files of all fixed positions for human 
write_bed(human_fixed_df, filename = "~/Desktop/Capellini_Lab/Human_ATAC/human_fixed_bps_hg38.bed")

#intersect positions with bed file containing per base pair phyloP scores
#sort -k1,1 -k2,2n human_fixed_bps_hg38.bed | bedtools intersect -a stdin -b /cold/aokamoto/241-mammalian-2020v2.bed -sorted -wa -wb  > human_fixed_bps_phyloP_hg38.bed 

#read in results 
human_fixed_bp_phylop <- read.delim(file = "~/Desktop/Capellini_Lab/Human_ATAC/human_fixed_bps_phyloP_hg38.bed", header = F, sep = "\t")
colnames(human_fixed_bp_phylop) <- c("chr", "start", "end", "ID", "pos_chr", "pos_start", "pos_end", "pos_id", "phylop")

chimp_fixed_bp_phylop <- read.delim(file = "~/Desktop/Capellini_Lab/Human_ATAC/chimp_fixed_bps_phyloP_hg38.bed", header = F, sep = "\t")
colnames(chimp_fixed_bp_phylop) <- c("chr", "start", "end", "ID", "pos_chr", "pos_start", "pos_end", "pos_id", "phylop")

#calculate stats for phyloP comparisons
mean(human_fixed_bp_phylop$phylop)
mean(chimp_fixed_bp_phylop$phylop)

#wilcoxon rank sum test
wilcox.test(human_fixed_bp_phylop$phylop, chimp_fixed_bp_phylop$phylop, paired = F)

#now for autopod specific
human_fixed_bp_phylop_autopod <- human_fixed_bp_phylop %>% dplyr::filter(ID %in% multiinter_autopod_only$ID)
chimp_fixed_bp_phylop_autopod <- chimp_fixed_bp_phylop %>% dplyr::filter(ID %in% multiinter_autopod_only$ID)
wilcox.test(human_fixed_bp_phylop_autopod$phylop, chimp_fixed_bp_phylop_autopod$phylop, paired = F)