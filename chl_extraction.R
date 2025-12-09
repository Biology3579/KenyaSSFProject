# combine_l3b_chlor.R
# Extracts tar archive, reads all L3b files with oceancolouR, outputs combined CSV/parquet

library(here)
library(oceancolouR)

library(data.table)
library(lubridate)
library(arrow)  # optional, for parquet

# ---- CONFIG ----
tar_file <- here::here("requested_files_1.tar")  # your tar archive
extract_dir <- here::here("raw_data")            # where to extract
out_csv <- here::here("all_chlor_a_combined.csv")
out_parquet <- here::here("all_chlor_a_combined.parquet")

# Optional: bounding box (Comoros/Kenya/Tanzania region)
bbox <- list(xmin = 38, xmax = 46, ymin = -14, ymax = 6)  # NULL = keep all

# ---- STEP 1: Extract tar if not already extracted ----
if (!dir.exists(extract_dir) || length(list.files(extract_dir, pattern = "\\.nc$")) == 0) {
  message("Extracting tar archive...")
  dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
  untar(tar_file, exdir = extract_dir)
  message("  Extracted to: ", extract_dir)
}

# ---- STEP 2: Find all .nc/.h5 files ----
files <- list.files(extract_dir, pattern = "\\.(nc|h5|hdf5)$", 
                    recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
if (length(files) == 0) stop("No .nc/.h5 files found in: ", extract_dir)
message(sprintf("Found %d files to process", length(files)))

# ---- STEP 3: Helper function to parse date from filename ----
parse_date <- function(fname) {
  # Extracts YYYYMMDD from typical NASA filename patterns
  m <- regexpr("(\\d{8})", fname)
  if (m[1] != -1) {
    d <- substring(fname, m[1], m[1] + 7)
    return(as.Date(d, format = "%Y%m%d"))
  }
  return(NA)
}

# ---- STEP 4: Read all files and combine ----
all_data <- rbindlist(lapply(seq_along(files), function(i) {
  f <- files[i]
  message(sprintf("[%d/%d] %s", i, length(files), basename(f)))
  
  # Read with oceancolouR (handles bin decoding automatically)
  dt <- tryCatch({
    l3 <- read_h5_L3b(f)
    as.data.table(l3)
  }, error = function(e) {
    message("  FAILED: ", e$message)
    return(NULL)
  })
  
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  
  # Standardize chlorophyll column name
  chl_col <- grep("chlor", names(dt), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(chl_col) && chl_col != "chlor_a") setnames(dt, chl_col, "chlor_a")
  
  # Add metadata
  dt[, `:=`(
    date = parse_date(basename(f)),
    source_file = basename(f)
  )]
  
  # Filter by bounding box if specified
  if (!is.null(bbox)) {
    dt <- dt[lon >= bbox$xmin & lon <= bbox$xmax & 
               lat >= bbox$ymin & lat <= bbox$ymax]
  }
  
  # Keep only valid chlor_a values
  dt <- dt[!is.na(chlor_a) & is.finite(chlor_a)]
  
  return(dt)
}), use.names = TRUE, fill = TRUE)

message(sprintf("\nCombined %s rows from %d files", 
                format(nrow(all_data), big.mark = ","), length(files)))

# ---- STEP 5: Write outputs ----
fwrite(all_data, out_csv)
message("✓ Wrote CSV: ", out_csv)

write_parquet(all_data, out_parquet)
message("✓ Wrote Parquet: ", out_parquet)

# ---- OPTIONAL: Quick summary ----
message("\nSummary:")
message("  Date range: ", min(all_data$date, na.rm = TRUE), " to ", 
        max(all_data$date, na.rm = TRUE))
message("  Chlor-a range: ", round(min(all_data$chlor_a, na.rm = TRUE), 3), 
        " to ", round(max(all_data$chlor_a, na.rm = TRUE), 3), " mg/m³")
message("  Mean: ", round(mean(all_data$chlor_a, na.rm = TRUE), 3), " mg/m³")