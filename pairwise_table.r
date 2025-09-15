# === LIST final.compare FILES ===
compare_files <- list.files("/path/to/compare/files/",
                            pattern = "_final\\.compare$", full.names = TRUE)
extracted_numbers <- regmatches(compare_files, gregexpr("\\d+", compare_files))

# === PREPARE OUTPUT ===
snp_diffs <- data.frame(Sample_1 = character(),
                        Sample_2 = character(),
                        N_differences = integer(),
                        stringsAsFactors = FALSE)

# === LOOP OVER FINAL COMPARE FILES ===
for (i in seq_along(compare_files)) {
  final_compare_file <- compare_files[i]
  sample1 <- extracted_numbers[[i]][1]
  sample2 <- extracted_numbers[[i]][2]

  cat("Processing:", basename(final_compare_file), "=> Sample1:", sample1, "Sample2:", sample2, "\n")

  if (is.na(sample1) || is.na(sample2)) {
    cat("Could not extract numeric sample IDs. Skipping.\n")
    next
  }

  file_size <- file.info(final_compare_file)$size
  if (is.na(file_size)) {
    cat("File does not exist or is unreadable:", final_compare_file, "\n")
    next
  }

  if (file_size == 0) {
    cat("  📭 Empty file. Recording 0 differences.\n")
    snp_diffs[nrow(snp_diffs) + 1, ] <- c(sample1, sample2, 0)
    next
  }

  # Attempt to read the file
  compare <- tryCatch({
    read.delim(final_compare_file, header = TRUE)
  }, error = function(e) {
    cat("Failed to read file:", final_compare_file, " - ", e$message, "\n")
    return(data.frame())  # Return empty data frame
  })

  # If read returned a data frame but it is empty, record as 0
  if (nrow(compare) == 0) {
    cat("  📭 File is empty or unreadable content. Recording 0 differences.\n")
    snp_diffs[nrow(snp_diffs) + 1, ] <- c(sample1, sample2, 0)
  } else {
    n_diff <- nrow(compare)
    cat("Read success. Differences found:", n_diff, "\n")
    snp_diffs[nrow(snp_diffs) + 1, ] <- c(sample1, sample2, n_diff)
  }
}

# === WRITE SUMMARY ===
cat("Total pairs written:", nrow(snp_diffs), "\n")

write.table(snp_diffs,
            "/path/to/output/Pairwise_SNP_Differences.txt",
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
