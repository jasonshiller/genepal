# Count different gene features and genes from gff file
summarise_gff <- function(gff_file) {
  gff <- read.csv(gff_file, sep = "\t", comment.char = "#", header = FALSE)
  summary <- table(gff[[3]])
  splicing_false_count <- sum(grepl("canonical_splicing=false", gff[[9]]))
  summary_df <- as.data.frame(as.list(summary))
  summary_df$Non_canon_splice_sites <- splicing_false_count
  summary_df$ID <- sub("\\.gff3$", "", basename(gff_file))
  return(summary_df)
}

count_eggnog_rows <- function(file_path) {
  read.csv(file_path, sep = "\t", comment.char = "#", header = FALSE) %>%
    nrow()
}

extract_eggnog_stats <- function() {
  annotation_files <- list.files(
    path = "results/annotations",
    pattern = "\\.annotations$",
    recursive = TRUE,
    full.names = TRUE
  )

  ids <- gsub(
    pattern = "\\.emapper\\.annotations$",
    replacement = "",
    x = basename(annotation_files)
  )

  row_counts <- unname(sapply(annotation_files, count_eggnog_rows))

  data.frame(
    ID = ids,
    "Eggnog annotations" = row_counts,
    check.names = FALSE
  )
}

# Parse BUSCO file (short summary)
parse_busco_file <- function(file_path) {
  busco_output <- readLines(file_path)
  patterns <- list(
    Complete_BUSCOs = "\\s*Complete BUSCOs \\(C\\)\\s*",
    Single_Copy_BUSCOs = "\\s*Complete and single-copy BUSCOs \\(S\\)\\s*",
    Duplicated_BUSCOs = "\\s*Complete and duplicated BUSCOs \\(D\\)\\s*",
    Fragmented_BUSCOs = "\\s*Fragmented BUSCOs \\(F\\)\\s*",
    Missing_BUSCOs = "\\s*Missing BUSCOs \\(M\\)\\s*",
    Total_BUSCOs_Searched = "\\s*Total BUSCO groups searched\\s*"
  )
  results <- sapply(names(patterns), function(name) {
    line <- busco_output[grepl(patterns[[name]], busco_output)]
    as.numeric(str_match(line, "(\\d+)")[, 2])
  }, USE.NAMES = TRUE)

  # Calculate percentages (Its easier than parsing the annoying format)
  results <- as.list(results)
  results$Complete_BUSCOs_Percentage <-
    round((results$Complete_BUSCOs / results$Total_BUSCOs_Searched) * 100,
      digits = 1
    )
  results$Single_Copy_BUSCOs_Percentage <-
    round((results$Single_Copy_BUSCOs / results$Total_BUSCOs_Searched) * 100,
      digits = 1
    )
  results$Duplicated_BUSCOs_Percentage <-
    round((results$Duplicated_BUSCOs / results$Total_BUSCOs_Searched) * 100,
      digits = 1
    )
  results$Fragmented_BUSCOs_Percentage <-
    round((results$Fragmented_BUSCOs / results$Total_BUSCOs_Searched) * 100,
      digits = 1
    )
  results$Missing_BUSCOs_Percentage <-
    round((results$Missing_BUSCOs / results$Total_BUSCOs_Searched) * 100,
      digits = 1
    )
  return(results)
}

parse_busco_folder <- function(folder_path, col_prefix = "Genome") {
  list_of_files <- list.files(folder_path, pattern = "short_summary.specific.*.txt$", full.names = TRUE)

  if (length(list_of_files) < 1) {
    return(NULL)
  }

  list_of_files %>%
    lapply(parse_busco_file) %>%
    bind_rows() %>%
    mutate(
      ID = sapply(strsplit(basename(list_of_files), "\\."), `[`, 3),
      Lineage = sapply(strsplit(basename(list_of_files), "\\."), `[`, 4)
    ) %>%
    select(
      ID,
      Lineage,
      Complete_BUSCOs_Percentage,
      Single_Copy_BUSCOs_Percentage,
      Duplicated_BUSCOs_Percentage,
      Fragmented_BUSCOs_Percentage,
      Missing_BUSCOs_Percentage
    ) %>%
    (
      function(df) {
        colnames(df) <- c(
          "ID",
          "Lineage",
          paste(col_prefix, "level BUSCO Complete (%)"),
          paste(col_prefix, "level BUSCO Single Copy (%)"),
          paste(col_prefix, "level BUSCO Duplicated (%)"),
          paste(col_prefix, "level BUSCO Fragmented (%)"),
          paste(col_prefix, "level BUSCO Missing (%)")
        )
        df
      })
}

# Convert the orthofinder hog data into counts for heatmap and other operations
# this data comes from "Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
# file from orthofinder. Function takes path to this as input
transform_hogs <- function(hogs_path) {
  count_genes <- function(cell) {
    if (cell == "") {
      return(0)
    } else {
      return(length(strsplit(cell, ",")[[1]]))
    }
  }
  tmp_data <- list()
  hogs <- utils::read.csv(file = hogs_path, sep = "\t")
  # pull out just the required cols (HOG id column and the sample columns)
  hogs <- hogs[, c(1, 4:ncol(hogs))]
  # Iterate over the columns (excluding the HOG column)
  for (col in names(hogs)[-1]) {
    # Apply the count_genes function to each cell in the column
    tmp_data[[col]] <- sapply(hogs[[col]], count_genes)
  }
  new_df <- as.data.frame(tmp_data)
  rownames(new_df) <- hogs$HOG
  colnames(new_df) <- gsub(
    x = colnames(new_df),
    pattern = ".pep", replacement = ""
  )
  return(new_df)
}

# process the Phylogenetic_Hierarchical_Orthogroups (N0.tsv)
# preps data for heatmap
process_hog_counts <- function(file_path) {
  # Check if the file exists, return NULL silently if it doesn't
  if (!file.exists(file_path)) {
    return(NULL)
  }
  # Proceed with the rest of the code if the file exists
  hogs_counts <- transform_hogs(file_path)
  # Get a vector of the rows that contain 1 orth per species
  # (in other words uninteresting for heatmap or display purposes)
  hca_1 <- apply(hogs_counts, 1, function(x) all(x == 1))
  hca_1_or_more <- apply(hogs_counts, 1, function(x) all(x >= 1))
  # Create heatmap data and calculate core/accessory orthogroups
  heatmap_hogs_counts_df <- t(hogs_counts[!hca_1, ])
  core <- sum(hca_1_or_more)
  accessory <- nrow(hogs_counts) - core
  # Generate heatmap title
  heat_title <- paste(
    "Core Orthogroups = ",
    core, ": Accessory orthogroups = ",
    accessory
  )
  return(list(heatmap = heatmap_hogs_counts_df, title = heat_title))
}

# Summarise protein clustering outputs
process_protein_clustering <- function(file_path) {
  row_id <- NULL # Explicit declaration to avoid lintr warning
  if (!file.exists(file_path)) {
    return(NULL) # No output if the file doesn't exist
  }
  df <- read.table(
    file = file_path,
    sep = "\t",
    header = TRUE,
    nrows = 10,
    row.names = 1
  )
  df <- df %>%
    rownames_to_column(var = "row_id") %>% # Changed to snake_case
    pivot_longer(cols = -row_id, names_to = "ID", values_to = "Value") %>%
    pivot_wider(names_from = row_id, values_from = "Value")
  colnames(df) <- gsub(
    x = colnames(df),
    pattern = "genes",
    replacement = "proteins"
  )
  return(df)
}

parse_gff_compare_file <- function(file_path) {
  sensitivity_precision <-
    read.csv(file_path, comment.char = "#", header = FALSE, nrows = 6, sep = "\t") %>%
    mutate(
      Description = str_trim(str_extract(V1, "^[^:]+")),
      Sensitivity = as.numeric(str_extract(V1, "(?<=:)[ ]*([0-9.]+)")),
      Precision = as.numeric(str_extract(V1, "(?<=\\|)[ ]*([0-9.]+)"))
    ) %>%
    mutate(ID = strsplit(basename(file_path), "\\.")[[1]][1]) %>%
    select(ID, Description, Sensitivity, Precision)

  matching <-
    read.csv(file_path, header = FALSE, sep = "\t", skip = 16, nrows = 3) %>%
    mutate(
      Description = str_trim(str_extract(V1, "^[^:]+")),
      Counts = as.numeric(str_extract(V1, "(?<=:)[ ]*([0-9.]+)"))
    ) %>%
    mutate(ID = strsplit(basename(file_path), "\\.")[[1]][1]) %>%
    select(ID, Description, Counts)

  missing_novel <- read.csv(file_path, header = FALSE, sep = "\t", skip = 20, nrows = 6) %>%
    mutate(
      Type = str_trim(str_extract(V1, "^[^:]+")),
      Count = as.integer(str_extract(V1, "(?<=:)[ ]*[0-9]+")),
      Total = as.integer(str_extract(V1, "(?<=/)[0-9]+")),
      Percentage = round(x = (Count / Total) * 100, digits = 1)
    ) %>%
    mutate(ID = strsplit(basename(file_path), "\\.")[[1]][1]) %>%
    select(ID, Type, Count, Total, Percentage)

  return(list(sensitivity_precision, matching, missing_novel))
}

parse_gff_compare_folder <- function(folder_path) {
  list_of_files <- list.files(folder_path, pattern = "\\.stats$", full.names = TRUE)

  if (length(list_of_files) < 1) {
    return(NULL)
  }

  parsed_data <- list_of_files %>%
    lapply(parse_gff_compare_file)

  sensitivity_precision <- parsed_data %>%
    lapply(function(x) x[[1]]) %>%
    bind_rows()

  matching <- parsed_data %>%
    lapply(function(x) x[[2]]) %>%
    bind_rows()

  missing_novel <- parsed_data %>%
    lapply(function(x) x[[3]]) %>%
    bind_rows()

  return(list(sensitivity_precision = sensitivity_precision, matching = matching, missing_novel = missing_novel))
}
