library(tidyverse)
library(arrow)

source('/home/alegbe/repos/evidence_datasource_parsers/modules/baseline_expression/AdaTiSS_fn.R')

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) < 1) {
  cat("Usage: Rscript process_adatiss.R <data_path> [datasource_filter] [datatype_filter] [unit_filter]\n")
  cat("Example: Rscript process_adatiss.R /path/to/parquet/directory DICE 'bulk rna-seq' TPM\n")
  cat("Example: Rscript process_adatiss.R /path/to/parquet/directory GTEX 'bulk rna-seq' TPM\n")
  cat("Example: Rscript process_adatiss.R /path/to/parquet/directory PRIDE 'mass-spectrometry proteomics' 'PPB (iBAQ)'\n")
  stop("Please provide the data path as the first argument")
}

# Set default values
data_path <- args[1]
datasource_filter <- if (length(args) >= 2) args[2] else NULL
datatype_filter <- if (length(args) >= 3) args[3] else NULL
unit_filter <- if (length(args) >= 4) args[4] else NULL

# Get the last part of the path for output directory naming
data_name <- basename(data_path)

cat("Data path:", data_path, "\n")
cat("Data source filter:", datasource_filter, "\n")
cat("Data type filter:", datatype_filter, "\n")
cat("Unit filter:", unit_filter, "\n")

# Load parquet data from all files
cat("Loading parquet data...\n")
parquet_files <- list.files(data_path, pattern = "\\.parquet$", full.names = TRUE)

if (length(parquet_files) == 0) {
  stop("No parquet files found in the specified directory")
}

# Read all parquet files and combine
data_list <- lapply(parquet_files, function(f) {
  cat("Reading", basename(f), "\n")
  read_parquet(f)
})

# Combine all data
data <- bind_rows(data_list)
cat("Total records:", nrow(data), "\n")

# Apply filters if specified
if (!is.null(datasource_filter)) {
  data <- data %>% filter(datasourceId == datasource_filter)
  cat("After datasource filter:", nrow(data), "records\n")
}

if (!is.null(datatype_filter)) {
  data <- data %>% filter(datatypeId == datatype_filter)
  cat("After datatype filter:", nrow(data), "records\n")
}

if (!is.null(unit_filter)) {
  data <- data %>% filter(unit == unit_filter)
  cat("After unit filter:", nrow(data), "records\n")
}

# Check available columns and data types
cat("Available columns:", colnames(data), "\n")
cat("Data types:", unique(data$datatypeId), "\n")
cat("Units:", unique(data$unit), "\n")
cat("Data sources:", unique(data$datasourceId), "\n")

# Transform from long format to wide format matrix (genes x samples)
# Create sample identifier based on available columns
if ("celltypeBiosampleFromSource" %in% colnames(data)) {
  # For cell-type data (DICE, Tabula Sapiens)
  data <- data %>%
    mutate(sample_id = paste(donorId, celltypeBiosampleFromSource, sep = "_"))
  biosample_col <- "celltypeBiosampleFromSource"
} else if ("tissueBiosampleFromSource" %in% colnames(data)) {
  # For tissue data (GTEX, PRIDE)
  data <- data %>%
    mutate(sample_id = paste(donorId, tissueBiosampleFromSource, sep = "_"))
  biosample_col <- "tissueBiosampleFromSource"
} else {
  # Fallback to just donorId
  data <- data %>%
    mutate(sample_id = donorId)
  biosample_col <- "donorId"
}

# Create expression matrix
cat("Creating expression matrix...\n")
expr_matrix <- data %>%
  select(targetId, sample_id, expression) %>%
  pivot_wider(names_from = sample_id, values_from = expression, values_fill = 0) %>%
  column_to_rownames("targetId") %>%
  as.matrix()

cat("Expression matrix dimensions:", dim(expr_matrix), "\n")

# Create phenotype data for biosample abundance calculation
p.dat <- data %>%
  select(sample_id, all_of(biosample_col)) %>%
  distinct() %>%
  as.data.frame()

cat("Phenotype data dimensions:", dim(p.dat), "\n")
cat("Tissue/cell types:", unique(p.dat[[biosample_col]]), "\n")

# Run AdaTiSS preprocessing
cat("Running AdaTiSS preprocessing...\n")

# Determine data type for preprocessing
if (!is.null(datatype_filter)) {
  if (str_detect(datatype_filter, "mass-spectrometry")) {
    dat_type <- "intensity"
    proc_zero <- "ceiled to 1"
    filter_col_prp <- 1
    exp_thres <- 1
  } else if (str_detect(datatype_filter, "rna-seq")) {
    dat_type <- "TPM or RPKM"
    proc_zero <- "ceiled to 1"
    filter_col_prp <- 1
    exp_thres <- 1
  } else {
    dat_type <- "TPM or RPKM"  # default
    proc_zero <- "ceiled to 1"
    filter_col_prp <- 1
    exp_thres <- 1
  }
} else {
  # Default to TPM/RPKM processing
  dat_type <- "TPM or RPKM"
  proc_zero <- "ceiled to 1"
  filter_col_prp <- 1
  exp_thres <- 1
}

cat("Using data type:", dat_type, "\n")
X <- preproc.filter.fn(expr_matrix, dat.type = dat_type, proc.zero = proc_zero, 
                       filter.col.prp = filter_col_prp, exp.thres = exp_thres)

cat("Filtered expression matrix dimensions:", dim(X), "\n")

# Calculate biosample abundance
cat("Calculating biosample abundance...\n")
biosample.abd <- tiss.abd.fn(X, p.dat)

# Run AdaTiSS analysis
cat("Running AdaTiSS analysis...\n")
cat("Processing", nrow(X), "genes across", ncol(X), "samples\n")
cat("This may take a while for large datasets...\n")

# For testing, you can subset to fewer genes by uncommenting the next line:
# X <- X[1:min(1000, nrow(X)), ]  # Process only first 1000 genes for testing

result <- AdaTiSS(X, tiss.abd = biosample.abd, adjust = TRUE, adjust.opt = 0)

# Save results
# Create output directory name with data source info
if (!is.null(datasource_filter)) {
  output_name <- paste0(data_name, "_", datasource_filter)
} else {
  output_name <- data_name
}

output_dir <- paste0("/home/alegbe/results/adatiss/", output_name)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("Saving results...\n")
write.table(result$ada.s, file = file.path(output_dir, "adatiss_sample_scores.txt"), 
            sep = "\t", quote = FALSE, row.names = TRUE)
write.table(result$ada.z, file = file.path(output_dir, "adatiss_biosample_scores.txt"), 
            sep = "\t", quote = FALSE, row.names = TRUE)
write.table(result$pop.fit.mx, file = file.path(output_dir, "adatiss_population_fit.txt"), 
            sep = "\t", quote = FALSE, row.names = TRUE)

cat("AdaTiSS analysis completed successfully!\n")
cat("Results saved to:", output_dir, "\n")
cat("Sample scores dimensions:", dim(result$ada.s), "\n")
cat("Biosample scores dimensions:", dim(result$ada.z), "\n")
