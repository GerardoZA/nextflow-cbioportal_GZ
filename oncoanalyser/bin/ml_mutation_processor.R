#!/usr/bin/env Rscript

library(tidyverse)
library(httr)
library(jsonlite)

#' Fetch and cache cancer hotspots data
#' @param hotspots_file Optional path to a pre-downloaded copy of
#'   https://www.cancerhotspots.org/api/hotspots/single (for offline compute nodes)
#' @return DataFrame of hotspot mutations
#' @keywords internal
fetch_hotspots <- function(hotspots_file = NULL){

  if (!is.null(hotspots_file)) {
    # Read pre-staged single residue hotspots
    message("Reading hotspot data from ", hotspots_file, "...")
    single_hotspots <- fromJSON(hotspots_file, flatten = TRUE)
  } else {
    # Fetch data from cancerhotspots.org
    message("Fetching hotspot data from cancerhotspots.org...")

    # Single residue hotspots
    response <- GET("https://www.cancerhotspots.org/api/hotspots/single", config=config(ssl_verifypeer = FALSE))
    single_hotspots <- fromJSON(rawToChar(response$content), flatten = TRUE)
  }

  # Process single hotspots data
  hotspots <- single_hotspots %>%
    as_tibble() %>%
    # Select base columns
    select(
      hugoSymbol,
      residue,
      tumorTypeCount,
      tumorCount,
      qValue,
      starts_with("variantAminoAcid.")
    ) %>%
    # Convert wide amino acid format to long
    pivot_longer(
      cols = starts_with("variantAminoAcid."),
      names_to = "variant",
      values_to = "count",
      names_prefix = "variantAminoAcid."
    ) %>%
    # Remove variants with no counts
    filter(!is.na(count) & count > 0) %>%
    # Clean up variant names
    mutate(
      variant = str_remove(variant, "^\\.|del$|dup$|ins.*$"),
      type = "single"
    ) %>%
    # Filter significant hotspots
    filter(qValue < 0.05) %>%
    # Create unique identifier for each hotspot
    mutate(
      hotspot_id = paste(hugoSymbol, residue, variant, sep = "_"),
      residue = as.character(residue)
    )

  message(sprintf("Found %d significant hotspot mutations across %d genes",
                 nrow(hotspots),
                 length(unique(hotspots$hugoSymbol))))

  return(hotspots)
}

#' Helper function to strictly map MAF effect strings to biological weights
#' @keywords internal
get_mutation_weight <- function(effect) {
  case_when(
    effect %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site") ~ 4,
    effect %in% c("Splice_Site", "Splice_Region") ~ 3,
    effect %in% c("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins") ~ 2,
    effect %in% c("Silent") ~ 0,
    TRUE ~ 0 # Covers 3'Flank, 5'Flank, 3'UTR, 5'UTR, RNA, Targeted_Region, etc.
  )
}

#' Process mutation data for deep learning input
#' @param input_file Path to mutation result data
#' @param min_freq Minimum mutation frequency across samples to include gene (default: 0.01)
#' @param hotspots_file Optional path to pre-downloaded hotspot JSON (NULL = fetch from API)
#' @return DataFrame containing hybrid mutation encoding
#' @export
process_mutation_data <- function(input_file, min_freq = 0.01, hotspots_file = NULL) {

  # Fetch hotspot data
  hotspots <- fetch_hotspots(hotspots_file)

  # Read mutation data
  mutations <- read_tsv(input_file, show_col_types = FALSE)
  mutations <- mutations[mutations$effect != "RNA",]

  # Record initial dimensions
  initial_genes <- length(unique(mutations$gene))
  n_samples <- length(unique(mutations$sample))

  # Get unique samples and genes
  samples <- unique(mutations$sample)
  dataset_genes <- unique(mutations$gene)

  # 1. Filter genes by mutation frequency
  gene_freq <- table(mutations$gene) / n_samples
  frequent_genes <- names(gene_freq)[gene_freq >= min_freq]

  # 2. Hotspot Rescue: Keep known driver genes regardless of frequency
  hotspot_genes_list <- unique(hotspots$hugoSymbol)
  rescue_genes <- intersect(dataset_genes, hotspot_genes_list)

  # 3. Combine and deduplicate
  retained_genes <- unique(c(frequent_genes, rescue_genes))
  n_rescued <- length(setdiff(rescue_genes, frequent_genes))

  # Create different encodings
  binary_matrix <- encode_binary_mutations(mutations, samples, retained_genes)
  effect_matrix <- encode_effect_mutations(mutations, samples, retained_genes)
  vaf_matrix <- encode_vaf_mutations(mutations, samples, retained_genes)
  integrated_matrix <- encode_integrated_mutations(mutations, samples, retained_genes)

  # Create hybrid encoding with hotspots and other mutations
  hybrid_matrix <- create_hybrid_encoding(mutations, samples, retained_genes, hotspots)

  # Round numerical values to 3 decimal places
  matrices <- list(
    binary = binary_matrix,
    effect = effect_matrix,
    vaf = vaf_matrix,
    integrated = integrated_matrix,
    hybrid = hybrid_matrix
  )

  matrices <- lapply(matrices, function(mat) {
    mat %>% mutate(across(-sample, ~round(., 3)))
  })

  # Validate all matrices
  lapply(matrices, validate_mutation_data)

  # Save all versions
  file_suffixes <- c("binary", "effect", "vaf", "integrated", "hybrid")

  for (suffix in file_suffixes) {
    print(sprintf("Saving mutations_processed_%s.tsv...", suffix))
    write_tsv(
      matrices[[suffix]],
      sprintf("mutations_processed_%s.tsv", suffix)
    )
  }

  # Calculate summary statistics

  n_hotspot_features <- sum(grepl("_", colnames(hybrid_matrix)))
  n_other_features <- sum(endsWith(colnames(hybrid_matrix), "_other"))

  # Print processing summary
  message("\nMutation processing summary:")
  message(sprintf("- Initial number of genes: %d", initial_genes))
  message(sprintf("- Removed %d genes (frequency < %g%%)",
                 initial_genes - length(frequent_genes), min_freq * 100))
  message(sprintf("- Rescued %d rare hotspot genes below the frequency threshold", n_rescued))
  message(sprintf("- Final retained genes: %d", length(retained_genes)))
  message(sprintf("- Number of samples: %d", n_samples))
  message(sprintf("- Number of hotspot features (columns): %d", n_hotspot_features))
  message(sprintf("- Number of 'other' gene features (columns): %d", n_other_features))
  message("- Encodings generated:")
  message("  1. Binary (mutation presence/absence)")
  message("  2. Effect-based (functional impact weights)")
  message("  3. VAF-based (variant allele frequencies)")
  message("  4. Integrated (effect * VAF)")
  message("  5. Hybrid (hotspot-specific + other mutations)\n")

  # Return invisibly to prevent console flooding
  return(invisible(data.frame(matrices[["hybrid"]])))
}

#' Encode mutations using binary representation
#' @keywords internal
encode_binary_mutations <- function(mutations, samples, genes) {
  mutation_matrix <- matrix(0, nrow = length(samples), ncol = length(genes),
                          dimnames = list(samples, genes))

  for (i in seq_len(nrow(mutations))) {
    if (mutations$gene[i] %in% genes) {
      mutation_matrix[mutations$sample[i], mutations$gene[i]] <- 1
    }
  }

  as_tibble(mutation_matrix, rownames = "sample", .name_repair = "unique")
}

#' Encode mutations using effect categories
#' @keywords internal
encode_effect_mutations <- function(mutations, samples, genes) {
  mutation_matrix <- matrix(0, nrow = length(samples), ncol = length(genes),
                          dimnames = list(samples, genes))

  for (i in seq_len(nrow(mutations))) {
    if (mutations$gene[i] %in% genes) {
      weight <- get_mutation_weight(mutations$effect[i])
      current_weight <- mutation_matrix[mutations$sample[i], mutations$gene[i]]
      mutation_matrix[mutations$sample[i], mutations$gene[i]] <- max(current_weight, weight)
    }
  }

  as_tibble(mutation_matrix, rownames = "sample", .name_repair = "unique")
}

#' Encode mutations using variant allele frequency
#' @keywords internal
encode_vaf_mutations <- function(mutations, samples, genes) {
  mutation_matrix <- matrix(0, nrow = length(samples), ncol = length(genes),
                          dimnames = list(samples, genes))

  for (i in seq_len(nrow(mutations))) {
    if (mutations$gene[i] %in% genes) {
      current_vaf <- mutation_matrix[mutations$sample[i], mutations$gene[i]]
      mutation_matrix[mutations$sample[i], mutations$gene[i]] <-
        max(current_vaf, mutations$dna_vaf[i])
    }
  }

  as_tibble(mutation_matrix, rownames = "sample", .name_repair = "unique")
}

#' Encode mutations using integrated approach (effect * VAF)
#' @keywords internal
encode_integrated_mutations <- function(mutations, samples, genes) {
  mutation_matrix <- matrix(0, nrow = length(samples), ncol = length(genes),
                          dimnames = list(samples, genes))

  weight_matrix <- matrix(-1, nrow = length(samples), ncol = length(genes),
                          dimnames = list(samples, genes))

  for (i in seq_len(nrow(mutations))) {
    if (mutations$gene[i] %in% genes) {
      weight <- get_mutation_weight(mutations$effect[i])
      integrated_value <- weight * mutations$dna_vaf[i]

      sample_id <- mutations$sample[i]
      gene_id <- mutations$gene[i]
      stored_weight <- weight_matrix[sample_id, gene_id]

      if (weight > stored_weight) {
        weight_matrix[sample_id, gene_id] <- weight
        mutation_matrix[sample_id, gene_id] <- integrated_value
      } else if (weight == stored_weight) {
        mutation_matrix[sample_id, gene_id] <- max(
          mutation_matrix[sample_id, gene_id],
          integrated_value
        )
      }
    }
  }

  as_tibble(mutation_matrix, rownames = "sample", .name_repair = "unique")
}

#' Create hybrid encoding combining hotspots and other mutations with integrated encoding
#' @keywords internal
create_hybrid_encoding <- function(mutations, samples, genes, hotspots) {
  hotspot_genes <- unique(hotspots$hugoSymbol)
  hotspot_lookup <- unique(hotspots$hotspot_id)
  names(hotspot_lookup) <- hotspot_lookup

  features <- c(
    unique(hotspots$hotspot_id),
    paste0(hotspot_genes, "_other"),
    setdiff(genes, hotspot_genes)
  )

  mutation_matrix <- matrix(0, nrow = length(samples),
                          ncol = length(features),
                          dimnames = list(samples, features))

  weight_matrix <- matrix(-1, nrow = length(samples),
                          ncol = length(features),
                          dimnames = list(samples, features))

  for (i in seq_len(nrow(mutations))) {
    if (mutations$gene[i] %in% genes) {
      gene <- mutations$gene[i]

      weight <- get_mutation_weight(mutations$effect[i])
      integrated_value <- weight * mutations$dna_vaf[i]

      feature_col <- NA
      if (gene %in% hotspot_genes) {
        aa_change <- mutations$Amino_Acid_Change[i]

        if (!is.na(aa_change)) {
          residue <- sub("p\\.([A-Z]\\d+).*", "\\1", aa_change)
          variant <- sub("p\\.[A-Z]*\\d+([A-Za-z0-9_=*])", "\\1", aa_change)
          hotspot_id <- paste(gene, residue, variant, sep = "_")
        } else {
          hotspot_id <- "NA"
        }

        if (hotspot_id %in% names(hotspot_lookup)) {
          feature_col <- hotspot_id
        } else {
          feature_col <- paste0(gene, "_other")
        }
      } else {
        feature_col <- gene
      }

      stored_weight <- weight_matrix[mutations$sample[i], feature_col]

      if (weight > stored_weight) {
        weight_matrix[mutations$sample[i], feature_col] <- weight
        mutation_matrix[mutations$sample[i], feature_col] <- integrated_value
      } else if (weight == stored_weight) {
        mutation_matrix[mutations$sample[i], feature_col] <- max(
          mutation_matrix[mutations$sample[i], feature_col],
          integrated_value
        )
      }
    }
  }

  total_hotspots <- sum(mutation_matrix[, hotspots$hotspot_id] > 0)
  total_other <- sum(mutation_matrix[, paste0(hotspot_genes, "_other")] > 0)
  total_regular <- sum(mutation_matrix[, setdiff(genes, hotspot_genes)] > 0)

  message(sprintf("Hybrid encoding summary:"))
  message(sprintf("- Hotspot mutations: %d", total_hotspots))
  message(sprintf("- Other mutations in hotspot genes: %d", total_other))
  message(sprintf("- Regular gene mutations: %d", total_regular))

  as_tibble(mutation_matrix, rownames = "sample", .name_repair = "unique")
}

#' Validate mutation data structure and content
#' @param mutation_data Processed mutation data frame
#' @return Logical indicating if validation passed (invisible)
#' @keywords internal
validate_mutation_data <- function(mutation_data) {
  if (nrow(mutation_data) == 0) stop("No samples remained after processing")
  if (!"sample" %in% colnames(mutation_data)) stop("Sample column not found in mutation data")
  if (any(duplicated(mutation_data$sample))) stop("Duplicate samples found in mutation data")

  non_numeric_check <- mutation_data %>%
    select(-sample) %>%
    sapply(function(x) all(is.numeric(x)))

  if (!all(non_numeric_check)) stop("Non-numeric values found in mutation data")

  invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0){
  stop("Usage: Rscript ml_mutation_processor.R <path to mutation result data> [path to hotspots JSON]")
}

input_file <- args[1]
hotspots_file <- if (length(args) >= 2) args[2] else NULL

if (!file.exists(input_file)) {
  stop("Error : Input file does not exist: ", input_file)
}

if (!is.null(hotspots_file) && !file.exists(hotspots_file)) {
  stop("Error : Hotspots file does not exist: ", hotspots_file)
}

# The invisible wrap prevents it from printing the dataframe object back to console
invisible(process_mutation_data(input_file, hotspots_file = hotspots_file))
