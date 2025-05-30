#!/usr/bin/Rscript

commands <- commandArgs()

print(commands)
print(length(commands))

idr_calling_updated <- function(list_set, pooled_peaks, outdir, CUTOFF = NULL, makeplots=FALSE) {
  
  #####
  ##Updated to reflect new github version of IDR, as well as latest ENCODE CHIP-seq guidelines (October 2017).
  ##
  
  ##Performs IDR on specified data.
  ##Takes as input a list of .narrowPeak files (called from MACS2) for individual replicates (e.g. Prox1, Prox2, called with MACS2 separately).
  ##
  ##Note: this will require a 'pooled_peaks' file. This is generated
  ##by pooling together all of the duplicate filtered .bed data used to call individual peak sets,
  ##then calling MACS2 on this pooled .bed file.
  ##pooled_peaks should be absolute filepath.
  
  ##Note that this will require absolute pathnames, given that IDR needs to be called from its 'idrCode' folder as it assumes all the code is there.
  ##These are easily made with something like:
  ##
  ##ls $PWD/*.narrowPeak > peak_list.txt
  ##
  ##Sorting shouldn't matter.
  
  ##Starting off with batch-consistency-analysis.r 
  
  input_list <- file(list_set)
  rep_data <- readLines(input_list)
  close(input_list)
  
  #file_list_char <- paste(rep_data, collapse=" ")
  
  if(is.null(outdir)) {
    outdir <- paste(prefix, "_outdir", sep="")
  }
  
  try(dir.create(outdir)) ##If it hasn't already been made.
  
  ##Now, IDR does things in a pair-wise fashion, so we'll need to do labelling. I think
  ##we'll be using just numeric identifiers, as using filenames would be far too long.
  
  keep_track_table <- data.frame(filename=rep_data, rep_number =1:length(rep_data))
  write.csv(keep_track_table, paste(outdir, "/", "idr_repcompare_labelling_conventions", ".csv", sep=""))
  
  out_counter <- 1
  
  out_file_prefixes <- list()
  
  while(out_counter <= length(rep_data)) {
    
    int_counter <- out_counter
    
    while(int_counter <= length(rep_data)) {
      
      if(int_counter == out_counter) {
        print("ignored self-self")
      }else{
        
        ##idr --samples ${REP1_PEAK_FILE} ${REP2_PEAK_FILE} --peak-list ${POOLED_PEAK_FILE} --input-file-type narrowPeak --output-file ${IDR_OUTPUT} --rank signal.value --soft-idr-threshold ${IDR_THRESH} --plot --use-best-multisummit-IDR
        
        #sys_command <- paste0("idr --samples ", rep_data[out_counter], " ", rep_data[int_counter], " --peak-list ", pooled_peaks, " --input-file-type narrowPeak --output-file ", paste(outdir, "/", "rep_", out_counter, "_IDR_", "rep_", int_counter, sep=""), " --rank p.value --soft-idr-threshold 0.05 --plot --use-best-multisummit-IDR")
        sys_command <- paste0("idr --samples ", rep_data[out_counter], " ", rep_data[int_counter], " --peak-list ", pooled_peaks, " --input-file-type narrowPeak --output-file ", paste(outdir, "/", "rep_", out_counter, "_IDR_", "rep_", int_counter, sep=""), " --rank p.value --soft-idr-threshold 0.05 --plot")
        ##"--use-best-multisummit-IDR" is warned about, but it's mentioned in the ENCODE Chip-seq pipeline, so I'd trust it?
        
        #idr --samples PF1_dup.bed_biorep_shift_single_peaks.narrowPeak PF3_dup.bed_biorep_shift_single_peaks.narrowPeak --peak-list proximal_femur_pooled_dupremoved.bed_biorep_shift_single_peaks.narrowPeak --input-file-type narrowPeak --output-file test   --use-old-output-format --rank p.value --soft-idr-threshold 0.05 --plot
        
        
        #sys_command <- paste("Rscript batch-consistency-analysis.r ", rep_data[out_counter], " ", rep_data[int_counter],  
        #" -1 ", paste(outdir, "/", "rep_", out_counter, "_IDR_", "rep_", int_counter, sep=""), " 0 F p.value", sep="")
        
        ##^Note that this is hard-coded for MACS2, and will definitely be different for other peak callers.
        
        system(sys_command)
        
        out_file_prefixes <- c(out_file_prefixes, paste(outdir, "/", "rep_", out_counter, "_IDR_", "rep_", int_counter, sep=""))
        
        
      }
      
      int_counter <- int_counter + 1
      print("finished another comparison")
      print(length(out_file_prefixes))
    }
    
    out_counter <- out_counter + 1
  }
  
  if(makeplots == TRUE) {
    ##Just making some quality control plots with an IDR script.
    
    sys_command_plot <- paste("Rscript batch-consistency-plot.r ", length(out_file_prefixes), " ", paste(outdir, "/", "all_idr_rep_compares", "_graphs", sep=""), " ", 
                              paste(out_file_prefixes, collapse=" "), sep="")
    
    system(sys_command_plot)
    
  }
  
  ##Here would be pseudo-replicate IDR data analysis.
  
  ##Next: actual slicing.
  
  ##Need to decide what cut-off IDR value to use initially:
  
  peak_sizes <- sapply(rep_data, function(x) unlist(strsplit(system(paste("wc -l ", x, sep=""), intern=TRUE), " "))[1])
  peak_sizes <- sapply(peak_sizes, as.integer)
  ##Let's just go for the max.
  
  if(is.null(CUTOFF)) {
    
    if(max(peak_sizes) < 1e5) { ##Recommended by IDR.
      cutoff_thresh <- 0.05
      
    }else{ ##Be more stringent.
      #cutoff_thresh <- 0.01
      cutoff_thresh <- 0.05 ##FORCED.
    }
  }else{
    print("CUSTOM CUTOFF")
    print(CUTOFF)
    cutoff_thresh <- CUTOFF
    
  }
  
  print("selected IDR cutoff threshold")
  print(cutoff_thresh)
  
  ##Now to count the number of peaks in our IDR output files which have IDR values below this threshold.
  
  #( awk '$11 <= 0.05 {print $0}' distal_idr_replicates-overlapped-peaks.txt | wc -l )
  
  ##Now we just go by prefix.
  
  overlap_peak_files <- sapply(unlist(out_file_prefixes), function(x) paste(x, "-overlapped-peaks.txt", sep=""))
  ##for readability.
  
  ##Updated to $12 for new IDR output (global IDR, $11 being local IDR values).
  ##
  ##Also, note that now IDR reports p-values as log 10's
  
  trans_cutoff <- -log(cutoff_thresh)/log(10)
  
  #num_passed_peaks <- sapply(overlap_peak_files, function(x) system(paste("( awk '$11 <= ", cutoff_thresh, " {print $0}' ", x, " | wc -l )", sep=""), intern=TRUE))
  num_passed_peaks <- sapply(unlist(out_file_prefixes), function(x) system(paste("( awk '$12 >= ", trans_cutoff, " {print $0}' ", x, " | wc -l )", sep=""), intern=TRUE))
  num_passed_peaks <- sapply(num_passed_peaks, as.integer)
  
  ##And then we simply take the max.
  
  peaks_to_slice <- max(num_passed_peaks)
  
  ##For the record.
  
  idr_results <- data.frame(filename=overlap_peak_files, peaks_passed=num_passed_peaks)
  write.csv(idr_results, paste(outdir, "/", "idr_repcompare_cutoff_numbers_", cutoff_thresh, ".csv", sep=""))
  
  ##I won't be worrying about 'conservative' vs. 'optimal', as I don't know how pseudoReplicate data could be generated
  ##for single-end read data.
  
  #Finally, slicing the pooled peaks based on the maximum number of replicate-consolidated peaks (below the IDR threshold).
  
  ##Now, sorting the merged file by p-value (8th column of macs output) and cutting off the first 'N' number of peaks, where
  ##'N' is the max num of peaks passing the IDR-threshold.
  
  thresh_call <- paste("sort -k8nr,8nr ", pooled_peaks, " | head -n ", peaks_to_slice, " > ", paste(outdir, "/", gsub(".narrowPeak", paste("_replicate_consolidated_IDR_threshold", cutoff_thresh, ".narrowPeak", sep=""), unlist(strsplit(pooled_peaks, "/"))[length(unlist(strsplit(pooled_peaks, "/")))]), sep=""), sep="")
  print(cutoff_thresh)
  print(thresh_call)
  system(thresh_call)
  
}

idr_calling_updated(commands[7], commands[8], commands[9], CUTOFF=0.05)