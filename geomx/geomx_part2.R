###############################################################################
# GeoMx RNA-seq Analysis Pipeline (Part II)
#
# Description:
# This script performs downstream analysis of GeoMx RNA-seq data, including:
#   1. Library loading and processed data import
#   2. Segment filtering
#   3. Differential expression analysis (Q3 and background normalization)
#   4. Visualization and saving results
#
# Notes:
# - File paths below are examples and should be updated for your local or HPC setup.
# - This script assumes that preprocessed data from Part I is available.
#
# Reference:
# Adapted from the Bioconductor GeoMx workflow:
# https://www.bioconductor.org/packages/release/workflows/vignettes/
# GeomxTools_RNA-NGS_Analysis.html
#
# Author: Xue Li
# Date: 2026-03-19
###############################################################################

# 1. Load Libraries and Data -------------------------------------------------
# Example paths (update as needed)
library_path <- "../software/Rlib"
working_dir <- "../geomx/output/"
RData_file <- "GeoMx_data_processed.RData"
output_dir <- paste0("../geomx/output_", Sys.Date(), "/diff/")


# Set up library paths
.libPaths(library_path)  # Adjust path as needed

# Check R version (optional)
cat("R version:", R.Version()$version.string, "\n")
#  R version: R version 4.4.1 (2024-06-14) 

# Load required libraries
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(knitr)
library(scales)

# Set working directory and define data paths
setwd(working_dir)  # Adjust path as needed
# Get today's date

dir.create(output_dir, showWarnings = FALSE)
load(paste0(working_dir, "/GeoMx_data_processed.RData"))

###############################################################################
# 2. Filter Segment -----------------------------------------------------------
table(pData(target_geoMx_data)$Segment__Name__Label_)
#    CD163      CD3 Vimentin 
#      16        5       58 

# Filter out CD3
# Filter out rows where Segment__Name__Label_ equals "CD3"
target_geoMx_data <- target_geoMx_data[, pData(target_geoMx_data)$Segment__Name__Label_ != "CD3"]

table(pData(target_geoMx_data)$Segment__Name__Label_)
# CD163 Vimentin 
#    16       58 

###############################################################################
# 3. Differential Expression Analysis -----------------------------------------

# Replace values in the Factor_Site column
pData(target_geoMx_data)$Factor_Site <- gsub("OS Lung Metastasis #1", "Lung", pData(target_geoMx_data)$Factor_Site)
pData(target_geoMx_data)$Factor_Site <- gsub("OS Renal Metastasis", "Renal", pData(target_geoMx_data)$Factor_Site)
pData(target_geoMx_data)$Factor_Site <- gsub("Primary OS", "Primary", pData(target_geoMx_data)$Factor_Site)


#############################################
# Diff Test starts here


get_results <- function(mixedOutmc) {
  # Combine results from the "lsmeans" element
  results <- do.call(rbind, mixedOutmc["lsmeans", ])
  rowname <- rownames(results)
  results <- as.data.frame(results)
  results$Contrast <- rowname
  
  
  # Add the gene information by repeating column names for each row
  results$Gene <- unlist(lapply(colnames(mixedOutmc), function(col) {
    rep(col, nrow(mixedOutmc["lsmeans", ][[1]]))
  }))
  
  # Add adjusted p-values (FDR)
  results$FDR <- p.adjust(results$`Pr(>|t|)`, method = "fdr")
  
  # Subset and rename columns for clarity
  results <- results[, c("Gene", "Contrast", "Estimate", "Pr(>|t|)", "FDR")]
  results$P <- results$`Pr(>|t|)`
  
  return(results)
}

# Applies a log transformation to an assay data element (q_norm) in the object
assayDataElement(object = target_geoMx_data, elt = "log_q") <-
    assayDataApply(target_geoMx_data, 2, FUN = log, base = 2, elt = "q_norm")
assayDataElement(object = target_geoMx_data, elt = "log_neg") <-
	assayDataApply(target_geoMx_data, 2, FUN = log, base = 2, elt = "neg_norm")
		

# Convert test variables to factors 
pData(target_geoMx_data)[["dogid"]] <- factor(pData(target_geoMx_data)[["Factor_DogID"]])
pData(target_geoMx_data)[["slides"]] <-	factor(pData(target_geoMx_data)[["slide name"]])
		

# Define a reusable function to perform mixed model tests
process_mixed_model <- function(geoMx_data, elt, contrast_column, segment_label = NULL, factorsite = NULL, nCores = 4) {
  # Prepare the contrast factor
  pData(geoMx_data)$contrast <- factor(pData(geoMx_data)[[contrast_column]])
  
  # Filter by segment and factor site if provided
  if (!is.null(segment_label) & !is.null(factorsite)) {
    ind <- pData(geoMx_data)$Segment__Name__Label_ == segment_label & 
           pData(geoMx_data)$Factor_Site == factorsite
    geoMx_data <- geoMx_data[, ind]
  } else if (!is.null(segment_label)) {
    ind <- pData(geoMx_data)$Segment__Name__Label_ == segment_label
    geoMx_data <- geoMx_data[, ind]
  } else if (!is.null(factorsite)) {
    ind <- pData(geoMx_data)$Factor_Site == factorsite
    geoMx_data <- geoMx_data[, ind]
  }
  
  if((any(pData(geoMx_data)$contrast %in% "High")&any(pData(geoMx_data)$contrast %in% "Low"))|(!grepl("Infiltration", contrast_column)&!grepl("Thres", contrast_column))){
	  # Perform mixed model differential expression
	  mixedOutmc <- mixedModelDE(
	    geoMx_data,
	    elt = elt,  # Test for the specified normalization
	    modelFormula = ~ contrast + (1 + contrast | slides),
	    groupVar = "contrast",
	    nCores = nCores,
	    multiCore = FALSE
	  )
	  cat("Test done!\n")
  
	  # Extract and format results
	  results <- get_results(mixedOutmc)
	#  results <- results[, -2]  # Remove "Subset" column

	  # Split and rank results
	  split_column <- "Contrast"
	  rank_column <- "P"
	  split_results <- split(results, results[[split_column]])
	  ranked_results <- lapply(split_results, function(df) {
	    df[order(df[[rank_column]]), ]
	  })
  
	  return(ranked_results)
  
  }else{
	  cat("Test can be applied only to factors with 2 or more levels, STOP and no results!\n")
	  return("No results!\n")
  }
}

########
# Define contrast columns
contrast_columns <- c(
  "Factor_Site",
  "Factor_DogID",
  "MonoEnrichment_HigherThreshold",
  "MonoEnrichment_LowerThreshold",
  "Tcell_Infiltration_HighT",
  "Tcell_Infiltration_LowT"

)

# Define the parameter combinations
combinations <- list(
  list(segment_label = NULL, factorsite = NULL),              # All samples
  list(segment_label = "Vimentin", factorsite = NULL),        # Vimentin
  list(segment_label = "CD163", factorsite = NULL),           # CD163
  list(segment_label = NULL, factorsite = "Lung"),            # Lung
  list(segment_label = NULL, factorsite = "Primary"),         # Primary
  list(segment_label = "Vimentin", factorsite = "Lung"),      # Vimentin & Lung
  list(segment_label = "Vimentin", factorsite = "Primary"),    # Vimentin & Primary
  list(segment_label = "CD163", factorsite = "Lung"),      # CD163 & Lung
  list(segment_label = "CD163", factorsite = "Primary")    # CD163 & Primary
)

# Run tests for both normalizations: Q3 (`log_q`) and negative normalization (`log_neg`)
# Define normalization methods
normalizations <- c("log_neg", "log_q")


# Initialize an empty list to store results
results <- list()

# Nested loops: Iterate through contrast columns, normalizations, and parameter combinations
for (contrast_column in contrast_columns) {
  for (norm in normalizations) {
    for (combo in combinations) {
		
	  if(contrast_column=="Factor_Site"&&!is.null(combo$factorsite)){ next }
	  
	  cat(paste0("Performing test for: ", contrast_column,":", norm, "\n"))
	  cat("        segment&factorsite:\n")
	  cat(paste0("        ", combo$segment_label, "&", combo$factorsite, "\n"))
      # Extract parameters
      segment_label <- combo$segment_label
      factorsite <- combo$factorsite
      
      # Generate a unique key for storing the result
      key <- paste0(
        contrast_column, "_", norm,
        if (!is.null(segment_label)) paste0("_", segment_label) else "",
        if (!is.null(factorsite)) paste0("_", factorsite) else ""
      )
      
	  key <- gsub("log_q", "Q3", key)
	  key <- gsub("log_neg", "Neg", key)
	  key <- gsub("Vimentin", "Vim", key)
	  key <- gsub("MonoEnrichment_HigherThreshold","Mono_HighThres",key)
	  key <- gsub("MonoEnrichment_LowerThreshold","Mono_LowThres",key)
	  key <- gsub("Tcell_Infiltration_HighT","Tcell_HighThres",key)
	  key <- gsub("Tcell_Infiltration_LowT","Tcell_LowThres",key)
	  key <- gsub("Factor_Site","FacSite",key)
	  key <- gsub("Factor_DogID","DogID",key)

	  
      # Run the process_mixed_model function
      results[[key]] <- process_mixed_model(
        geoMx_data = target_geoMx_data,
        elt = norm,
        contrast_column = contrast_column,
        segment_label = segment_label,
        factorsite = factorsite
      )
	 cat("Results saved to: ", key, "\n")
	 cat("=============================================\n")
    }
  }
}

results <- results[!sapply(results, is.character)]



# To write two Excel sheets, one for Q3 normalization and one for Neg normalization, we can use the openxlsx package.
# Initialize two lists to store results for Q3 and Neg
results_Q3 <- list()
results_Neg <- list()

# Separate results into Q3 and Neg based on their keys
for (key in names(results)) {
  # Extract result
  result <- results[[key]]
  
  # Check for Q3 or Neg and store appropriately
  if (grepl("_Q3", key)) {
    results_Q3[[key]] <- result
  } else if (grepl("_Neg", key)) {
    results_Neg[[key]] <- result
  }
}

#########################
# Function to reverse the order of substrings in a column
reverse_column <- function(data, column_name) {
  # Check if the column exists in the data frame
  if (!column_name %in% colnames(data)) {
    stop("Column not found in the data frame!")
  }
  
  # Reverse the substrings for the entire column
  data[[column_name]] <- sapply(data[[column_name]], function(x) {
    split_string <- strsplit(as.character(x), " - ")[[1]]
    reversed_string <- paste(rev(split_string), collapse = " - ")
    return(reversed_string)
  })
  
  return(data)
}


# Function to write results to Excel
write_results_to_excel <- function(result_list, file_name) {
  sheet_names <- c()
  # Create a new workbook
  wb <- createWorkbook()
  
  # Loop through each result
  for (key in names(result_list)) {
    result <- result_list[[key]]
	cat(paste0("=========\n"))
    cat(paste0(key,"\n"))
    if (is.list(result)) {
      # If the result is a nested list, write each element into a separate sheet
      for (sub_key in names(result)) {
        if (grepl("N/A", sub_key)) {
            cat("Skipping sub_key with 'N/A':", sub_key, "\n")
            next # Skip this sub_key
        }
		 
        # Generate a unique sheet name
        sheet_name <- gsub("Q3_", "", paste0(key, "_", sub_key))
        sheet_name <- gsub("Neg_", "", sheet_name)
        sheet_name <- gsub("T_OS_00", "", sheet_name)
        sheet_name <- gsub("_High - Low", "", sheet_name)
        sheet_name <- gsub("Primary", "Pri", sheet_name)
        # sheet_name <- gsub("Seg", "FacSite", sheet_name) # Will delete this line later

        # Add worksheet and write data
        addWorksheet(wb, sheetName = sheet_name)
        writeData(wb, sheet = sheet_name, x = reverse_column(result[[sub_key]], "Contrast"))
		sheet_names <- c(sheet_names, sheet_name)
		cat(paste0(sheet_name,"\n"))
      }
    } else {
	
      # If the result is not a list, write it directly to one sheet
      sheet_name <- gsub("Q3_", "", key)
	  sheet_name <- gsub("Neg_", "", sheet_name)
	  sheet_name <- gsub("T_OS_00", "", sheet_name)
	  sheet_name <- gsub("_High - Low", "", sheet_name)
	  sheet_name <- gsub("Primary", "Pri", sheet_name)
	  # sheet_name <- gsub("Seg", "FacSite", sheet_name) # Will delete this line later
	  
      # Add worksheet and write data
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb, sheet = sheet_name, x =  reverse_column(result, "Contrast"))
	  sheet_names <- c(sheet_names, sheet_name)
	  cat(paste0(sheet_name,"\n"))
    }
  }
  
  # Save the workbook
  saveWorkbook(wb, file_name, overwrite = TRUE)
  
  # Write out contrast test list. 
  write.table(data.frame(Test = sheet_names), paste0(output_dir, "/contrast.list.tsv"), quote=F, row.names=F, col.names=T, sep="\t" )
}

# Write Neg results to an Excel file
write_results_to_excel(results_Neg, paste0(output_dir, "/Neg_Results.xlsx"))

# Write Q3 results to an Excel file
write_results_to_excel(results_Q3, paste0(output_dir, "/Q3_Results.xlsx"))


cat("Excel files 'Q3_Results.xlsx' and 'Neg_Results.xlsx' have been created.")


#################################################################
# Save image 
save.image(file = paste0(output_dir, "/GeoMx_data_test.RData"))



#################################################################
# -------------
# For all the comparisons, get the distribution in each class
# Factor Site
table(pData(target_geoMx_data)$Factor_Site)

seg="Vimentin"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg, ]$Factor_Site)

seg="CD163"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg, ]$Factor_Site)

# Factor Dog ID
table(pData(target_geoMx_data)$Factor_DogID)

seg="Vimentin"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg, ]$Factor_DogID)

seg="CD163"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg, ]$Factor_DogID)


site="Lung"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Factor_Site==site, ]$Factor_DogID)

site="Primary"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Factor_Site==site, ]$Factor_DogID)

seg="Vimentin"
site="Lung"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$Factor_DogID)

seg="Vimentin"
site="Primary"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$Factor_DogID)


# MonoEnrichment_HigherThreshold
table(pData(target_geoMx_data)$MonoEnrichment_HigherThreshold)

seg="Vimentin"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg, ]$MonoEnrichment_HigherThreshold)

seg="CD163"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg, ]$MonoEnrichment_HigherThreshold)


site="Lung"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Factor_Site==site, ]$MonoEnrichment_HigherThreshold)

site="Primary"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Factor_Site==site, ]$MonoEnrichment_HigherThreshold)


seg="Vimentin"
site="Lung"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$MonoEnrichment_HigherThreshold)

seg="Vimentin"
site="Primary"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$MonoEnrichment_HigherThreshold)


seg="CD163"
site="Lung"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$MonoEnrichment_HigherThreshold)

seg="CD163"
site="Primary"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$MonoEnrichment_HigherThreshold)

# MonoEnrichment_LowerThreshold
table(pData(target_geoMx_data)$MonoEnrichment_LowerThreshold)

seg="Vimentin"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg, ]$MonoEnrichment_LowerThreshold)

seg="CD163"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg, ]$MonoEnrichment_LowerThreshold)


site="Lung"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Factor_Site==site, ]$MonoEnrichment_LowerThreshold)

site="Primary"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Factor_Site==site, ]$MonoEnrichment_LowerThreshold)


seg="Vimentin"
site="Lung"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$MonoEnrichment_LowerThreshold)

seg="Vimentin"
site="Primary"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$MonoEnrichment_LowerThreshold)

seg="CD163"
site="Lung"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$MonoEnrichment_LowerThreshold)

seg="CD163"
site="Primary"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$MonoEnrichment_LowerThreshold)


# Tcell_Infiltration_HighT
table(pData(target_geoMx_data)$Tcell_Infiltration_HighT)

seg="Vimentin"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg, ]$Tcell_Infiltration_HighT)

seg="CD163"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg, ]$Tcell_Infiltration_HighT)


site="Lung"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Factor_Site==site, ]$Tcell_Infiltration_HighT)

site="Primary"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Factor_Site==site, ]$Tcell_Infiltration_HighT)


seg="Vimentin"
site="Lung"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$Tcell_Infiltration_HighT)

seg="Vimentin"
site="Primary"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$Tcell_Infiltration_HighT)

seg="CD163"
site="Lung"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$Tcell_Infiltration_HighT)

seg="CD163"
site="Primary"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$Tcell_Infiltration_HighT)



# Tcell_Infiltration_LowT
table(pData(target_geoMx_data)$Tcell_Infiltration_LowT)

seg="Vimentin"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg, ]$Tcell_Infiltration_LowT)

seg="CD163"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg, ]$Tcell_Infiltration_LowT)


site="Lung"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Factor_Site==site, ]$Tcell_Infiltration_LowT)

site="Primary"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Factor_Site==site, ]$Tcell_Infiltration_LowT)


seg="Vimentin"
site="Lung"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$Tcell_Infiltration_LowT)

seg="Vimentin"
site="Primary"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$Tcell_Infiltration_LowT)


seg="CD163"
site="Lung"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$Tcell_Infiltration_LowT)

seg="CD163"
site="Primary"
table(pData(target_geoMx_data)[pData(target_geoMx_data)$Segment__Name__Label_==seg&pData(target_geoMx_data)$Factor_Site==site, ]$Tcell_Infiltration_LowT)

