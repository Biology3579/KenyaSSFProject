## ---------------------------
##
## Script name: exploring_cleaning_processing.r
##
## Purpose of script: 
##      # A file of functions for cleaning the Palmer Penguins dataset
##
## Author: Candela Ferrer Diez
##
## Date Created: 2025-10-26
##
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


# A function to explore the dataset
explore_data <- function(data, name) {
  cat("\n========================================\n")
  cat(toupper(name), "DATA\n")
  cat("========================================\n\n")
  
  cat("Dimensions:", paste(dim(data), collapse = " x "), "\n\n")
  
  cat("Structure:\n")
  print(glimpse(data))
  
  cat("\n--- Key Location Names ---\n")
  if("Site" %in% names(data)) cat("Unique Sites:", length(unique(data$Site)), "\n", paste(unique(data$Site), collapse = ", "), "\n")
  if("Station" %in% names(data)) cat("Unique Stations:", length(unique(data$Station)), "\n", paste(unique(data$Station), collapse = ", "), "\n")
  cat("\n")
}


# A function to standardize site name formatting
## i.e. trim and lowercase snake_case
standardize_site_format <- function(data) {
  data %>% 
    mutate(
      # Remove leading/trailing whitespace
      site = str_trim(site),
      # Remove extra spaces between words
      site = str_squish(site),
      # Convert to lowercase
      site = str_to_lower(site),
      # Replace spaces with underscores for snake_case
      site = str_replace_all(site, " ", "_")
    )
}

# A function to convert column types
clean_data <- function(data, 
                       date_cols = NULL,
                       numeric_cols = NULL,
                       factor_cols = NULL,
                       remove_na_cols = TRUE) {
  
  cleaned <- data
  
  # Convert date columns
  if (!is.null(date_cols)) {
    cleaned <- cleaned %>% 
      mutate(across(all_of(date_cols), ~dmy(.)))
  }
  
  # Convert numeric columns
  if (!is.null(numeric_cols)) {
    cleaned <- cleaned %>% 
      mutate(across(all_of(numeric_cols), as.numeric))
  }
  
  # Convert factor columns
  if (!is.null(factor_cols)) {
    cleaned <- cleaned %>% 
      mutate(across(all_of(factor_cols), as.factor))
  }
  
  # Remove entirely NA columns (optional)
  if (remove_na_cols) {
    cleaned <- cleaned %>% 
      select(where(~!all(is.na(.))))
  }
  
  return(cleaned)
}