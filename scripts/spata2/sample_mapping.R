# ============================================================================
# Sample Mapping: Visium HD <-> IMC
# HCC Bone Metastasis Spatial Atlas
# ============================================================================
# 
# Patient ID: 27038840
# 
# Two tissue types:
#   1. Bone Metastasis (spine biopsy)
#   2. Primary HCC Tumor (surgical resection)
#
# ============================================================================

# --- Standardized Output Names ---
# Use these prefixes for ALL output files (plots, tables) to ensure consistency
OUTPUT_NAMES <- list(
  bone_metastasis = "BoneMet",
  primary_tumor = "PrimaryTumor"
)

# --- Sample Mapping Table ---
SAMPLE_MAPPING <- data.frame(
  tissue_type = c("bone_metastasis", "primary_tumor"),
  tissue_label = c("Bone Metastasis", "Primary HCC Tumor"),
  output_prefix = c("BoneMet", "PrimaryTumor"),
  
  # Visium HD samples (original names from SpaceRanger)
  visium_sample = c("D1_27038840_HCCboneMet_seg", 
                    "SRout_27038840_HCC_tumour"),
  visium_path = c("results/spaceranger/D1_27038840_HCCboneMet_seg_spaceranger_out",
                  "results/spaceranger/SRout_27038840_HCC_tumour"),
  
  # IMC samples (original names from MCD files)
  imc_sample = c("27038840_Biopsy_Spine", 
                 "27038840_Resection"),
  imc_mcd_file = c("27038840_Biopsy_Spine.mcd",
                   "27038840_Resection.mcd"),
  imc_n_rois = c(4, 5),
  
  stringsAsFactors = FALSE
)

# --- IMC ROI Details ---
IMC_ROIS <- data.frame(
  tissue_type = c(rep("bone_metastasis", 4), rep("primary_tumor", 5)),
  imc_sample = c(rep("27038840_Biopsy_Spine", 4), rep("27038840_Resection", 5)),
  roi_id = c(1:4, 1:5),
  roi_file = c(
    # Bone metastasis ROIs
    "27038840_Biopsy_Spine_1_1.txt",
    "27038840_Biopsy_Spine_2_2.txt",
    "27038840_Biopsy_Spine_3_3.txt",
    "27038840_Biopsy_Spine_4_4.txt",
    # Primary tumor ROIs
    "27038840_Resection_1_1.txt",
    "27038840_Resection_2_2.txt",
    "27038840_Resection_3_3.txt",
    "27038840_Resection_4_4.txt",
    "27038840_Resection_5_5.txt"
  ),
  stringsAsFactors = FALSE
)

# --- Helper Functions ---

#' Get Visium sample name from tissue type
get_visium_sample <- function(tissue_type) {
  SAMPLE_MAPPING$visium_sample[SAMPLE_MAPPING$tissue_type == tissue_type]
}

#' Get IMC sample name from tissue type
get_imc_sample <- function(tissue_type) {
  SAMPLE_MAPPING$imc_sample[SAMPLE_MAPPING$tissue_type == tissue_type]
}

#' Get tissue type from Visium sample name
get_tissue_from_visium <- function(visium_sample) {
  idx <- grep(visium_sample, SAMPLE_MAPPING$visium_sample, fixed = TRUE)
  if (length(idx) > 0) {
    return(SAMPLE_MAPPING$tissue_type[idx[1]])
  }
  return(NA)
}

#' Get tissue type from IMC sample name
get_tissue_from_imc <- function(imc_sample) {
  idx <- grep(imc_sample, SAMPLE_MAPPING$imc_sample, fixed = TRUE)
  if (length(idx) > 0) {
    return(SAMPLE_MAPPING$tissue_type[idx[1]])
  }
  return(NA)
}

#' Get matching IMC sample for a Visium sample
get_matched_imc <- function(visium_sample) {
  tissue <- get_tissue_from_visium(visium_sample)
  if (!is.na(tissue)) {
    return(get_imc_sample(tissue))
  }
  return(NA)
}

#' Get matching Visium sample for an IMC sample
get_matched_visium <- function(imc_sample) {
  tissue <- get_tissue_from_imc(imc_sample)
  if (!is.na(tissue)) {
    return(get_visium_sample(tissue))
  }
  return(NA)
}

#' Get all ROI files for an IMC sample
get_imc_rois <- function(imc_sample) {
  IMC_ROIS$roi_file[IMC_ROIS$imc_sample == imc_sample]
}

#' Get tissue label (pretty name)
get_tissue_label <- function(tissue_type) {
  SAMPLE_MAPPING$tissue_label[SAMPLE_MAPPING$tissue_type == tissue_type]
}

#' Get standardized output prefix for a sample (ST or IMC)
#' This ensures consistent naming across all outputs
#' @param sample_name Either a Visium or IMC sample name
#' @return Standardized prefix like "BoneMet" or "PrimaryTumor"
get_output_prefix <- function(sample_name) {
  # Try Visium first
  tissue <- get_tissue_from_visium(sample_name)
  if (is.na(tissue)) {
    # Try IMC
    tissue <- get_tissue_from_imc(sample_name)
  }
  if (!is.na(tissue)) {
    return(SAMPLE_MAPPING$output_prefix[SAMPLE_MAPPING$tissue_type == tissue])
  }
  # Fallback: clean up sample name
  return(gsub("[^A-Za-z0-9]", "_", sample_name))
}

#' Generate standardized output filename
#' @param sample_name Original sample name (Visium or IMC)
#' @param analysis_type e.g., "celltype", "cd8_spatial", "interaction"
#' @param suffix e.g., "heatmap", "umap", "stats"
#' @param ext File extension (default "png")
#' @return Standardized filename like "BoneMet_celltype_heatmap.png"
make_output_name <- function(sample_name, analysis_type, suffix = NULL, ext = "png") {
  prefix <- get_output_prefix(sample_name)
  if (!is.null(suffix) && suffix != "") {
    return(sprintf("%s_%s_%s.%s", prefix, analysis_type, suffix, ext))
  } else {
    return(sprintf("%s_%s.%s", prefix, analysis_type, ext))
  }
}

#' Get platform type (ST or IMC) from sample name
get_platform <- function(sample_name) {
  if (!is.na(get_tissue_from_visium(sample_name))) {
    return("ST")
  } else if (!is.na(get_tissue_from_imc(sample_name))) {
    return("IMC")
  }
  return("Unknown")
}

# --- Print mapping summary ---
print_sample_mapping <- function() {
  cat("\n=== HCC Bone Metastasis Sample Mapping ===\n")
  cat("Patient ID: 27038840\n\n")
  cat("Standardized Output Prefixes:\n")
  cat("  Bone Metastasis -> BoneMet_*\n")
  cat("  Primary Tumor   -> PrimaryTumor_*\n\n")
  
  for (i in 1:nrow(SAMPLE_MAPPING)) {
    cat(sprintf("--- %s (prefix: %s) ---\n", 
                SAMPLE_MAPPING$tissue_label[i],
                SAMPLE_MAPPING$output_prefix[i]))
    cat(sprintf("  Visium HD: %s\n", SAMPLE_MAPPING$visium_sample[i]))
    cat(sprintf("  IMC:       %s (%d ROIs)\n", 
                SAMPLE_MAPPING$imc_sample[i], 
                SAMPLE_MAPPING$imc_n_rois[i]))
    cat("\n")
  }
}

# Print on source
if (interactive()) {
  print_sample_mapping()
}
