# List SNP comparison files
snp_files <- list.files("/path/to/files/", pattern = "\\.snp$", full.names = TRUE)
extracted_numbers <- regmatches(snp_files, gregexpr("\\d+", snp_files))

# Prepare output data frame
snp_diffs <- data.frame(Sample_1 = integer(), Sample_2 = integer(), N_differences = integer(), stringsAsFactors=FALSE)

# Loop over each .snp file
for (i in 1:length(snp_files)) {
  # Define output file path
  final_compare_file <- paste0("/path/to/files/", 
                               extracted_numbers[[i]][1], "vs", extracted_numbers[[i]][2], "_final.compare")
  
  # If the final.compare file already exists, read it instead of recomputing
  if (file.exists(final_compare_file)) {
    differences_filtered <- read.delim(final_compare_file, header = TRUE)
  } else {
    # Load input files
    compare <- read.delim(snp_files[i], header = TRUE)
    samp1 <- read.table(paste0("/path/to/files/", extracted_numbers[[i]][1], ".cns"), header = TRUE)
    samp2 <- read.table(paste0("/path/to/files/", extracted_numbers[[i]][2], ".cns"), header = TRUE)

    # Normalize frequency columns
    compare[, 8] <- as.numeric(sub("%", "", compare[, 8], fixed = TRUE)) / 100
    compare[, 14] <- as.numeric(sub("%", "", compare[, 14], fixed = TRUE)) / 100

    # Recombination filtering
    thresh <- 4
    window <- 100
    differences_filtered <- compare

    for (r in 1:nrow(compare)) {
      pos <- compare[r, 1]
      if (samp1[samp1$Position == pos, "Var"] != ".") {
        background <- samp1[samp1$Position > (pos - window) & samp1$Position < (pos + window),]
        n_mutations <- nrow(background[background$Var != ".",])
        if (n_mutations >= thresh || abs(compare[r, 8] - compare[r, 14]) < 0.4) {
          differences_filtered <- differences_filtered[differences_filtered$Position != pos,]
        }
      }
      if (samp2[samp2$Position == pos, "Var"] != ".") {
        background <- samp2[samp2$Position > (pos - window) & samp2$Position < (pos + window),]
        n_mutations <- nrow(background[background$Var != ".",])
        if (n_mutations >= thresh || abs(compare[r, 8] - compare[r, 14]) < 0.4) {
          differences_filtered <- differences_filtered[differences_filtered$Position != pos,]
        }
      }
    }

    # Save filtered output
    write.table(differences_filtered, final_compare_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
  }

  # Record results
  snp_diffs[i, 1] <- extracted_numbers[[i]][1]
  snp_diffs[i, 2] <- extracted_numbers[[i]][2]
  snp_diffs[i, 3] <- nrow(differences_filtered)
}

# Save summary
write.table(snp_diffs,
            "/path/to/files/Pairwise_SNP_Differences.txt",
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
