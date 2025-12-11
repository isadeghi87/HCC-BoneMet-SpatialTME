#!/usr/bin/env Rscript
# 00_load_build_objects.R
# Load SpaceRanger data and build SPATA2 objects
# HIGH-QUALITY PLOTTING VERSION

suppressPackageStartupMessages({
  library(SPATA2)
  library(ggplot2)
  library(future)
})

source("scripts/spata2/plot_themes.R")

# Configuration
options(spata2r.verbose = TRUE)
plan("multisession", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2)

SAMPLE_INFO <- list(
  list(
    name = "D1_27038840_HCCboneMet_seg",
    directory = "results/spaceranger/D1_27038840_HCCboneMet_seg_spaceranger_out/outs",
    description = "HCC bone metastasis"
  ),
  list(
    name = "SRout_27038840_HCC_tumour",
    directory = "results/spaceranger/SRout_27038840_HCC_tumour/outs",
    description = "HCC primary tumour"
  )
)

OBJECTS_DIR <- "analysis/objects"
PLOTS_DIR <- "analysis/plots/qc"
LOGS_DIR <- "analysis/logs"

dir.create(OBJECTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOGS_DIR, paste0("00_load_build_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=== Load and Build SPATA2 Objects ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

for (sample in SAMPLE_INFO) {
  cat("\n--- Processing:", sample$name, "---\n")
  
  if (!dir.exists(sample$directory)) {
    cat("  Directory not found:", sample$directory, "\n")
    next
  }
  
  output_file <- file.path(OBJECTS_DIR, paste0(sample$name, "_raw.rds"))
  
  if (file.exists(output_file)) {
    cat("  Object already exists, skipping load.\n")
    next
  }
  
  tryCatch({
    cat("  Loading VisiumHD data (16um resolution)...\n")
    
    spata_obj <- initiateSpataObjectVisiumHD(
      sample_name = sample$name,
      directory_visium = sample$directory,
      square_res = "16um",
      unload = TRUE
    )
    
    cat("  ✓ Object created\n")
    cat("    Spots:", nObs(spata_obj), "\n")
    cat("    Genes:", nGenes(spata_obj), "\n")
    
    # Identify tissue outline
    cat("  Identifying tissue outline...\n")
    spata_obj <- identifyTissueOutline(spata_obj, method = "image")
    
    # Save raw object
    saveRDS(spata_obj, output_file)
    cat("  Saved:", output_file, "\n")
    
    rm(spata_obj)
    gc()
    
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
  })
}

cat("\n=== Object Loading Complete ===\n")
sink()
