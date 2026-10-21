######################################################
###############################################################################
# GeoMx RNA-seq Analysis Pipeline (Part I)
#
# Description:
# This script performs initial processing and analysis of NanoString GeoMx
# RNA-seq data, including:
#   1. Library loading and data import
#   2. Quality control (QC) and preprocessing
#   3. Normalization
#   4. Unsupervised analysis (e.g., PCA, clustering)
#   5. Visualization and output of results
#
# Notes:
# - File paths in this script are provided as examples and should be updated
#   according to your local or HPC environment.
# - This script is intended as part of a multi-step workflow.
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
# relative paths (adjust to your environment)
library_path = "../software/Rlib"  # Path to R libraries (adjust as needed)
working_dir = "../geomx"  # Working directory
data_dir = "../K9_GeoMx/GEOMX_ANALYSIS_FILES/"  # Data directory
dcc_dir = "../K9_GeoMx/DCC Files Including Dummy Files/"  # DCC files directory
pkc_dir = "../K9_GeoMx/GEOMX_ANALYSIS_FILES/PKCFiles/"  # PKC files directory
annotation_file = "../geomx/annotation/K9_GeoMx_Annotations_V2_modified.xlsx"  # Annotation file


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
library(ggforce)
library(networkD3)
library(plotly)
library(htmlwidgets) # For saving interactive HTML
library(openxlsx)
library(knitr)
library(scales)

# Set working directory and define data paths
setwd(working_dir)  # Adjust path as needed
# Get today's date
output_dir <- paste0("out_", Sys.Date(), "/processed/")
dir.create(output_dir, showWarnings = FALSE)

# Load data
# Adjust these paths and file names to match your environment
DCCFiles <- dir(dcc_dir, pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
PKCFiles <- dir(pkc_dir, pattern = ".pkc$", full.names = TRUE, recursive = TRUE)

geoMx_data <- readNanoStringGeoMxSet(
  dccFiles = DCCFiles,
  pkcFiles = PKCFiles,
  phenoDataFile = annotation_file,
  phenoDataSheet = "Template",
  phenoDataDccColName = "Sample_ID",
  protocolDataColNames = c("aoi", "roi"),
  experimentDataColNames = c("panel")
)

# Sample Overview

# Define Sankey column names
sankeyCols <- c("source", "target", "value")

# Prepare links from your data
link1 <- count(pData(geoMx_data), `slide name`, Factor_DogID)
link2 <- count(pData(geoMx_data),  Factor_DogID, Factor_Site)
link3 <- count(pData(geoMx_data),  Factor_Site, Segment__Name__Label_)

# Rename columns to fit Sankey requirements
colnames(link1) <- sankeyCols
colnames(link2) <- sankeyCols
colnames(link3) <- sankeyCols

# Combine all links
links <- rbind(link1,link2,link3)

# Create unique nodes from the source and target columns
nodes <- unique(data.frame(name=c(links$source, links$target)))


# Sankey diagram uses 0-based indexing; update source and target indices
links$source <- as.integer(match(links$source,nodes$name)-1)
links$target <- as.integer(match(links$target,nodes$name)-1)

# Create the Sankey diagram
fig <- plot_ly(
  type = "sankey",
  arrangement = "fixed",
  node = list(
    label = nodes$name,
    pad = 15,
    thickness = 20
  ),
  link = list(
    source = links$source,
    target = links$target,
    value = links$value
  )
)


# Save as interactive HTML
htmlwidgets::saveWidget(fig, paste0(output_dir, "/sankey_plot_with_labels.html"))

###############################################################################
# 2. Quality Control (QC) and Preprocessing ----------------------------------

# Before we begin, we will shift any expression counts with a value of 0 to 1 to enable in downstream transformations.
geoMx_data <- shiftCountsOne(geoMx_data, useDALogic = TRUE)

#--------------------------# 
# 2.1 Segment QC           #
#--------------------------#
# We first assess sequencing quality and adequate tissue sampling for every ROI/AOI segment.
QC_params <- list(
  minSegmentReads = 1000,
  percentTrimmed = 80,
  percentStitched = 80,
  percentAligned = 80,
  percentSaturation = 50,
  minNegativeCount = 1,
  maxNTCCount = 9000,
  minArea = 1000
)
geoMx_data <- setSegmentQCFlags(geoMx_data, qcCutoffs = QC_params)

# collect QC results
QCResults <- protocolData(geoMx_data)[["QCFlags"]]
flag_columns <- colnames(QCResults)

QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
						 
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})

QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

print(table(QCResults$QCStatus))
#    PASS WARNING 
#    106      40 

# Visualize Segment QC
col_by <- "Segment__Name__Label_"
# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

# Now we explore each of the QC metrics for the segments.
# Save Trimmed (%) Plot
ggsave(file.path(output_dir, "Trimmed_Percentage.pdf"), 
       plot = QC_histogram(sData(geoMx_data), "Trimmed (%)", col_by, 80), 
       width = 8, height = 6)

# Save Stitched (%) Plot
ggsave(file.path(output_dir, "Stitched_Percentage.pdf"), 
       plot = QC_histogram(sData(geoMx_data), "Stitched (%)", col_by, 80), 
       width = 8, height = 6)

# Save Aligned (%) Plot
ggsave(file.path(output_dir, "Aligned_Percentage.pdf"), 
       plot = QC_histogram(sData(geoMx_data), "Aligned (%)", col_by, 80), 
       width = 8, height = 6)

# Save Saturated (%) Plot
ggsave(file.path(output_dir, "Saturated_Percentage.pdf"), 
       plot = QC_histogram(sData(geoMx_data), "Saturated (%)", col_by, 50) +
              labs(title = "Sequencing Saturation (%)", 
                   x = "Sequencing Saturation (%)"), 
       width = 8, height = 6)

# Save Area Plot (Log10 Scale)
ggsave(file.path(output_dir, "Area_Log10.pdf"), 
       plot = QC_histogram(sData(geoMx_data), "area", col_by, 1000, scale_trans = "log10"), 
       width = 8, height = 6)

cat("Plots have been saved to:", output_dir, "\n")

# calculate the negative geometric means for each module
# calculates and visualizes the geometric mean of negative control probes for each module in the dataset
negativeGeoMeans <- 
  esBy(negativeControlSubset(geoMx_data), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(geoMx_data)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
pkcs <- annotation(geoMx_data)
modules <- gsub(".pkc", "", pkcs)
negCols <- paste0("NegGeoMean_", modules)
pData(geoMx_data)[, negCols] <- sData(geoMx_data)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QC_histogram(pData(geoMx_data), ann, col_by, 2, scale_trans = "log10")
  ggsave(file.path(output_dir, paste0(ann, "_histogram.pdf")), plot = plt, width = 8, height = 6)
}


# detatch neg_geomean columns ahead of aggregateCounts call
pData(geoMx_data) <- pData(geoMx_data)[, !colnames(pData(geoMx_data)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
kable(table(NTC_Count = sData(geoMx_data)$NTC),
      col.names = c("NTC Count", "# of Segments"))
kable(QC_Summary, caption = "QC Summary Table for each Segment")

# Define the output Excel file name
output_file <- file.path(output_dir, "/QC_and_NTC_Tables.xlsx")

# Create NTC Count Table
ntc_table <- as.data.frame(table(NTC_Count = sData(geoMx_data)$NTC))
colnames(ntc_table) <- c("NTC Count", "# of Segments")

# Add row names to QC_Summary as a column
qc_summary_with_rownames <- cbind(Segment = rownames(QC_Summary), QC_Summary)

# Create a new Excel workbook
wb <- createWorkbook()

# Add NTC Count Table as the first sheet
addWorksheet(wb, "NTC Count Table")
writeData(wb, "NTC Count Table", ntc_table)

# Add QC Summary Table as the second sheet
addWorksheet(wb, "QC Summary Table")
writeData(wb, "QC Summary Table", qc_summary_with_rownames)

# Save the workbook to file
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("Tables saved to:", output_file, "\n")

# Remove flagged segments
# As the final step in Segment QC, we remove flagged segments that do not meet our QC cutoffs.
geoMx_data <- geoMx_data[, QCResults$QCStatus == "PASS"]
dim(geoMx_data)
# Features  Samples 
#    2010      106 

#-------------------------------------# 	
# 2.2 Probe QC                        #
#-------------------------------------#
# Before we summarize our data into gene-level count data, we will remove low-performing probes.
# After Probe QC, there will always remain at least one probe representing every gene target. In other words, Probe QC never removes genes from your data.
# Set Probe QC Flags
geoMx_data <- setBioProbeQCFlags(geoMx_data, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE) 

ProbeQCResults <- fData(geoMx_data)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
#   Passed Global Local
#    2004      0     6 

# Exclude Outlier Probes
# Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(geoMx_data, 
         fData(geoMx_data)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(geoMx_data)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)
# Features  Samples 
# 2010      106 
geoMx_data <- ProbeQCPassed 

#-------------------------------------# 
# 2.3 Create Gene-level Count Data.   #
#-------------------------------------# 
# With our Probe QC steps complete, we will generate a gene-level count matrix. 
# The count for any gene with multiple probes per segment is calculated as the geometric mean of those probes.
length(unique(featureData(geoMx_data)[["TargetName"]]))
#> [1] 1963

# collapse to targets
target_geoMx_data <- aggregateCounts(geoMx_data)
dim(target_geoMx_data)
#> Features  Samples 
#>    1963      106
exprs(target_geoMx_data)[1:5, 1:2]

#-------------------------------------# 
# 2.4 Limit of Quantification         #
#-------------------------------------# 
# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_geoMx_data))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_geoMx_data)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_geoMx_data)[, vars[1]] * 
             pData(target_geoMx_data)[, vars[2]] ^ cutoff)
  }
}
pData(target_geoMx_data)$LOQ <- LOQ

#-------------------------------------# 
# 2.5 Filtering                       #
#-------------------------------------# 
# Filtering out either segments and/or genes with abnormally low signal based on LOQ.
# We determine the number of genes detected in each segment across the dataset.

LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_geoMx_data)$Module == module
  Mat_i <- t(esApply(target_geoMx_data[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_geoMx_data)$TargetName, ]

# Segment Gene Detection
# We first filter out segments with exceptionally low signal. 
# These segments will have a small fraction of panel genes detected above the LOQ relative to the other segments in the study.
pData(target_geoMx_data)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_geoMx_data)$GeneDetectionRate <-
  pData(target_geoMx_data)$GenesDetected / nrow(target_geoMx_data)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_geoMx_data)$DetectionThreshold <- 
  cut(pData(target_geoMx_data)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))


# Define output filenames
figure_file <- file.path(output_dir, "Gene_Detection_Rate_Stacked_Barplot.pdf")
table_file <- file.path(output_dir, paste0("Gene_Detection_Rate_Table.xlsx"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
bar_plot <- ggplot(pData(target_geoMx_data),
       aes(x = DetectionThreshold)) +
    geom_bar(aes(fill = Factor_Site)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Gene Detection Rate",
         y = "Segments, #",
         fill = "Segment Type")
		 
detection_table <- table(pData(target_geoMx_data)$DetectionThreshold,
	pData(target_geoMx_data)$class)

# Save the plot to a PDF file
ggsave(figure_file, plot = bar_plot, width = 8, height = 6)

# Save the table to an Excel file
wb <- createWorkbook()
addWorksheet(wb, "Detection Table")
writeData(wb, "Detection Table", as.data.frame(detection_table))
saveWorkbook(wb, table_file, overwrite = TRUE)
cat("Table saved to:", table_file, "\n")


# 0.1 as cutoff
target_geoMx_data <-
    target_geoMx_data[, pData(target_geoMx_data)$GeneDetectionRate >= .05]

dim(target_geoMx_data)
# Features  Samples 
#    1963       79 

# Gene detection rate
# Next, we determine the detection rate for genes across the study. To illustrate this idea, we create a small gene list (goi) to review.
LOQ_Mat <- LOQ_Mat[, colnames(target_geoMx_data)]
fData(target_geoMx_data)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_geoMx_data)$DetectionRate <-
  fData(target_geoMx_data)$DetectedSegments / nrow(pData(target_geoMx_data))

# Gene of interest detection table
# I guess we need to change the genes of interest
goi <- c("FANCF", "SET", "NASP", "WNT4", "CD68", "TCF7L1",
         "AXIN1", "Negative Probe")
goi_df <- data.frame(
  Gene = goi,
  Number = fData(target_geoMx_data)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(target_geoMx_data)[goi, "DetectionRate"]))

# Gene filtering
# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
    unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                  function(x) {sum(fData(target_geoMx_data)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_geoMx_data))
rownames(plot_detect) <- plot_detect$Freq

gene_rate <- ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
              vjust = 1.6, color = "black", size = 4) +
    scale_fill_gradient2(low = "orange2", mid = "lightblue",
                         high = "dodgerblue3", midpoint = 0.65,
                         limits = c(0,1),
                         labels = scales::percent) +
    theme_bw() +
    scale_y_continuous(labels = scales::percent, limits = c(0,1),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "% of Segments",
         y = "Genes Detected, % of Panel > LOQ")

# Save the plot to a PDF file
ggsave(paste0(output_dir, "/Gene_detection_rate.pdf"), plot = gene_rate, width = 8, height = 6)

# Filtering
# We typically set a % Segment cutoff ranging from 5-20% based on the biological diversity of our dataset. For this study, we will select 10% as our cutoff. 
# Subset to target genes detected in at least 10% of the samples.
#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_geoMx_data), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_geoMx_data <- 
    target_geoMx_data[fData(target_geoMx_data)$DetectionRate >= 0.1 |
                        fData(target_geoMx_data)$TargetName %in% neg_probes, ]
dim(target_geoMx_data)
# Features  Samples 
#     549       79 

# retain only detected genes of interest
goi <- goi[goi %in% rownames(target_geoMx_data)]


###############################################################################
# 3. Normalization -----------------------------------------------------------
library(reshape2)  # for melt
library(cowplot)   # for plot_grid

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "Factor_Site"
Stat_data <- 
    data.frame(row.names = colnames(exprs(target_geoMx_data)),
               Segment = colnames(exprs(target_geoMx_data)),
               Annotation = pData(target_geoMx_data)[, ann_of_interest],
               Q3 = unlist(apply(exprs(target_geoMx_data), 2,
                                 quantile, 0.75, na.rm = TRUE)),
               NegProbe = exprs(target_geoMx_data)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")


plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
    geom_histogram(bins = 40) + theme_bw() +
    scale_x_continuous(trans = "log2") +
    facet_wrap(~Annotation, nrow = 1) + 
    scale_fill_brewer(palette = 3, type = "qual") +
    labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
    geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
    geom_point() + guides(color = "none") + theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
    geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
    geom_point() + theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))
combined_plot <- plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))

ggsave(paste0(output_dir, "/Q3_normalization.pdf"), plot = combined_plot, width = 10, height = 8)

# Q3 normalization
# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
target_geoMx_data <- normalize(target_geoMx_data,
                             norm_method = "quant", 
                             desiredQuantile = .75,
                             toElt = "q_norm")

# Background normalization for WTA/CTA without custom spike-in
target_geoMx_data <- normalize(target_geoMx_data ,
                             norm_method = "neg", 
                             fromElt = "exprs",
                             toElt = "neg_norm")
							 

# visualize the first 10 segments with Q3 normalization method
pdf(paste0(output_dir, "/Segment_Q3_Boxplots.pdf"), width = 8, height = 10)  # Specify dimensions

# Set up layout for two plots on one page
par(mfrow = c(2, 1))  # 2 rows, 1 column

boxplot(exprs(target_geoMx_data)[,1:10],
        col = "#9EDAE5", main = "Raw Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Raw")

boxplot(assayDataElement(target_geoMx_data[,1:10], elt = "q_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Q3 Normalized")

dev.off()

pdf(paste0(output_dir, "/Segment_Neg_Boxplots.pdf"), width = 8, height = 10)  # Specify dimensions

# Set up layout for two plots on one page
par(mfrow = c(2, 1))  # 2 rows, 1 column

boxplot(exprs(target_geoMx_data)[,1:10],
        col = "#9EDAE5", main = "Raw Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Raw")

boxplot(assayDataElement(target_geoMx_data[,1:10], elt = "neg_norm"),
        col = "#FF7F0E", main = "Neg Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Neg. Normalized")

dev.off()

###############################################################################
# 4. Unsupervised Analysis -----------------------------------------------------------

#-------------------------------------# 
# 4.1 UMAP&t-SNE.                     #
#-------------------------------------# 

library(umap)
library(Rtsne)

# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42

# Q3 Normalization
# run UMAP
umap_out <-
    umap(t(log2(assayDataElement(target_geoMx_data , elt = "q_norm"))),  
         config = custom_umap)

pData(target_geoMx_data)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
umap <- ggplot(pData(target_geoMx_data),
       aes(x = UMAP1, y = UMAP2, color = Factor_Site, shape = class)) +
    geom_point(size = 3) +
    theme_bw()


# run tSNE
set.seed(42) # set the seed for tSNE as well
tsne_out <-
    Rtsne(t(log2(assayDataElement(target_geoMx_data , elt = "q_norm"))),
          perplexity = ncol(target_geoMx_data)*.15)
pData(target_geoMx_data)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
tsne <- ggplot(pData(target_geoMx_data),
       aes(x = tSNE1, y = tSNE2, color = Factor_Site, shape = class)) +
    geom_point(size = 3) +
    theme_bw()
	
# Combine the UMAP and t-SNE plots
combined_plot <- plot_grid(umap, tsne, ncol = 1, labels = c("A", "B"))
	
# Save the combined plot as a PDF
ggsave(paste0(output_dir, "/unsupervise_clustering_Q3.pdf"), plot = combined_plot, width = 8, height = 10)

# Neg Normalization
# run UMAP
umap_out <-
    umap(t(log2(assayDataElement(target_geoMx_data , elt = "neg_norm"))),  
         config = custom_umap)

pData(target_geoMx_data)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
umap <- ggplot(pData(target_geoMx_data),
       aes(x = UMAP1, y = UMAP2, color = Factor_Site, shape = class)) +
    geom_point(size = 3) +
    theme_bw()


# run tSNE
set.seed(42) # set the seed for tSNE as well
tsne_out <-
    Rtsne(t(log2(assayDataElement(target_geoMx_data , elt = "neg_norm"))),
          perplexity = ncol(target_geoMx_data)*.15)
pData(target_geoMx_data)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
tsne <- ggplot(pData(target_geoMx_data),
       aes(x = tSNE1, y = tSNE2, color = Factor_Site, shape = class)) +
    geom_point(size = 3) +
    theme_bw()
	
# Combine the UMAP and t-SNE plots
combined_plot <- plot_grid(umap, tsne, ncol = 1, labels = c("A", "B"))
	
# Save the combined plot as a PDF
ggsave(paste0(output_dir, "/unsupervise_clustering_Neg.pdf"), plot = combined_plot, width = 8, height = 10)


#-------------------------------------# 
# 4.2 Clustering high CV Genes        #
#-------------------------------------# 
# We plot the results using unsupervised hierarchical clustering, displayed as a heatmap.
library(pheatmap)  # for pheatmap

# Q3 Normalization
# create a log2 transform of the data for analysis
assayDataElement(object = target_geoMx_data, elt = "log_q") <-
    assayDataApply(target_geoMx_data, 2, FUN = log, base = 2, elt = "q_norm")

# create CV function
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_geoMx_data,
                         elt = "log_q", MARGIN = 1, calc_CV)
# show the highest CD genes and their CV values
sort(CV_dat, decreasing = TRUE)[1:5]
#>   CAMK2N1    AKR1C1      AQP2     GDF15       REN 
#> 0.5886006 0.5114973 0.4607206 0.4196469 0.4193216

# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]
pdf(paste0(output_dir, "/High_CV_genes_clustering_Q3.pdf"),  width = 10, height = 8)
pheatmap(assayDataElement(target_geoMx_data[GOI, ], elt = "log_q"),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = 
             pData(target_geoMx_data)[, c("class", "Segment__Name__Label_", "Factor_Site")])

dev.off()


# Neg Normalization
# create a log2 transform of the data for analysis
assayDataElement(object = target_geoMx_data, elt = "log_neg") <-
    assayDataApply(target_geoMx_data, 2, FUN = log, base = 2, elt = "neg_norm")

# create CV function
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_geoMx_data,
                         elt = "log_neg", MARGIN = 1, calc_CV)
# show the highest CD genes and their CV values
sort(CV_dat, decreasing = TRUE)[1:5]
#>   CAMK2N1    AKR1C1      AQP2     GDF15       REN 
#> 0.5886006 0.5114973 0.4607206 0.4196469 0.4193216

# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]
pdf(paste0(output_dir, "/High_CV_genes_clustering_Neg.pdf"),  width = 10, height = 8)
pheatmap(assayDataElement(target_geoMx_data[GOI, ], elt = "log_neg"),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = 
             pData(target_geoMx_data)[, c("class", "Segment__Name__Label_", "Factor_Site")])

dev.off()
###############################################################################
# 5. Output expression data
expr.matrix.Neg <- assayDataElement(target_geoMx_data, elt = "neg_norm")
genes <- rownames(expr.matrix.Neg)
expr.matrix.Neg <- as.data.frame(expr.matrix.Neg)
expr.matrix.Neg$GENE <- genes
expr.matrix.Neg <- expr.matrix.Neg[,c(ncol(expr.matrix.Neg), 1:(ncol(expr.matrix.Neg)-1))]

expr.matrix.Q3 <- assayDataElement(target_geoMx_data, elt = "q_norm")
genes <- rownames(expr.matrix.Q3)
expr.matrix.Q3 <- as.data.frame(expr.matrix.Q3)
expr.matrix.Q3$GENE <- genes
expr.matrix.Q3 <- expr.matrix.Q3[,c(ncol(expr.matrix.Q3), 1:(ncol(expr.matrix.Q3)-1))]

write.table(expr.matrix.Neg, paste0(output_dir, "/Background_Norm_EXPR.csv"), quote=F, col.names=T, row.names=F, sep=",")
write.table(expr.matrix.Q3, paste0(output_dir, "/Q3_Norm_EXPR.csv"), quote=F, col.names=T, row.names=F, sep=",")


###############################################################################
# 5. Save RData --------------------------------------------------------------
# Save the R session
save.image(file = paste0(output_dir, "/GeoMx_data_processed.RData"))


