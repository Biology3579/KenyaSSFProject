## ---------------------------
##
## Script name: extracting.r
##
## Purpose of script: 
##       A file of functions to extract additional data
##
## Author: Candela Ferrer Diez
##
## Date Created: 2025-11-15
##
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


# A function to explore the dataset ----
explore_data_structure <- function(data, name) {
  cat("\n========================================\n")
  cat(toupper(name), "DATA\n")
  cat("========================================\n\n")
  
  cat("Dimensions:", paste(dim(data), collapse = " x "), "\n\n")
  
  cat("Structure:\n")
  print(glimpse(data))
  
}

# A function to clean column names ----

clean_column_names <- function(data) {
  data %>%
    # standardize and clean spaces 
    rename_with(~stringr::str_trim(.)) %>%             # Remove leading/trailing spaces
    rename_with(~stringr::str_squish(.)) %>%           # Replace multiple spaces with a single space
    rename_with(~stringr::str_replace_all(., " ", "_")) %>% # Convert remaining single spaces to underscores
    
    # 2. Use janitor for final snake case conversion and handling of other special chars
    janitor::clean_names()  # Convert column names to snakecase
}

# A function to drop specified columns ----
drop_columns <- function(data, drop_cols) {
  
  # Check if drop_cols is NULL or empty; if so, return data unchanged
  if (is.null(drop_cols) || length(drop_cols) == 0) {
    return(data)
  }
  
  # Use select() with -any_of() to remove the specified columns
  data <- data %>%
    dplyr::select(-dplyr::any_of(drop_cols))
  
  return(data)
}

# A function to standardize text ----
standardize_text <- function(data, text_cols) {
  data %>%
    mutate(across(all_of(text_cols), str_trim)) %>%
    mutate(across(all_of(text_cols), str_squish)) %>%
    mutate(across(all_of(text_cols), str_to_lower)) %>%
    mutate(across(all_of(text_cols), ~str_replace_all(., " ", "_")))
}

# A function to convert column types ----
convert_column_types <- function(data, 
                                 date_cols = NULL,
                                 integer_cols = NULL,
                                 numeric_cols = NULL,
                                 factor_cols = NULL) {
  
  if (!is.null(date_cols)) {
    data <- data %>%
      mutate(across(all_of(date_cols), dmy))
  }
  
  if (!is.null(integer_cols)) {
    data <- data %>%
      mutate(across(all_of(integer_cols), ~suppressWarnings(as.integer(.))))
  }
  
  if (!is.null(numeric_cols)) {
    data <- data %>%
      mutate(across(all_of(numeric_cols), ~suppressWarnings(as.numeric(.))))
  }
  
  if (!is.null(factor_cols)) {
    data <- data %>%
      mutate(across(all_of(factor_cols), as.factor))
  }
  
  return(data)
}

# A function to fix coordinates (if needed) ----
fix_coordinates <- function(data, 
                            lat_col = "latitude", 
                            lon_col = "longitude") {
  data %>%
    mutate(
      "{lat_col}" := if_else(abs(.data[[lat_col]]) > 90, NA_real_, .data[[lat_col]]),
      "{lon_col}" := if_else(abs(.data[[lon_col]]) > 180, NA_real_, .data[[lon_col]])
    )
}

