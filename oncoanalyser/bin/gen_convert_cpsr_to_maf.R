#!/usr/bin/env Rscript

#' Convert CPSR TSV to MAF format and append to existing MAF file
#' Updated version that handles different CPSR column formats
#'
#' @param cpsr_file Path to the CPSR TSV file
#' @param maf_file Path to the existing MAF file
#' @param output_file Path to save the combined output MAF file
#'
#' @return None (writes to file)
convert_cpsr_to_maf <- function(cpsr_file, maf_file, output_file) {

  # Function to check if a string starts with a pattern
  starts_with <- function(x, pattern) {
    substr(x, 1, nchar(pattern)) == pattern
  }

  # Read the MAF file to get headers and extract sample information
  maf_lines <- readLines(maf_file)

  # Check if the first line is a version line
  if (starts_with(maf_lines[1], "#version")) {
    header_line <- maf_lines[2]
    data_start <- 3
  } else {
    header_line <- maf_lines[1]
    data_start <- 2
  }

  # Extract column names
  if (starts_with(header_line, "#")) {
    maf_columns <- strsplit(substring(header_line, 2), "\t")[[1]]
  } else {
    maf_columns <- strsplit(header_line, "\t")[[1]]
  }

  # Data rows of the somatic MAF (empty if the MAF has a header only)
  somatic_lines <- if (length(maf_lines) >= data_start) maf_lines[data_start:length(maf_lines)] else character(0)

  # vcf2maf leaves Mutation_Status blank: label the somatic MAF rows "Somatic" (CPSR rows below
  # are "Germline"). Non-blank values (e.g. RNA info from INTEGRATE_RNA_VARIANTS) are kept.
  split_maf_line <- function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    # strsplit drops trailing empty fields
    if (length(fields) < length(maf_columns)) {
      fields <- c(fields, rep("", length(maf_columns) - length(fields)))
    }
    fields
  }

  status_idx <- match("Mutation_Status", maf_columns)
  if (!is.na(status_idx)) {
    somatic_lines <- vapply(somatic_lines, function(line) {
      fields <- split_maf_line(line)
      if (fields[status_idx] == "") fields[status_idx] <- "Somatic"
      paste(fields, collapse = "\t")
    }, character(1), USE.NAMES = FALSE)
  }

  # Read the MAF data to extract sample barcodes
  maf_data <- read.delim(maf_file, comment.char="#", stringsAsFactors=FALSE)

  # Extract unique sample barcodes
  tumor_sample <- unique(maf_data$Tumor_Sample_Barcode)[1]
  normal_sample <- unique(maf_data$Matched_Norm_Sample_Barcode)[1]

  cat("Found tumor sample:", tumor_sample, "\n")
  cat("Found normal sample:", normal_sample, "\n")

  # Read the CPSR file
  cpsr_data <- read.delim(cpsr_file, stringsAsFactors=FALSE)

  # Check if CPSR file is empty (only header, no data)
  if (nrow(cpsr_data) == 0) {
    cat("WARNING: CPSR file contains only headers with no data rows. Creating empty output file.\n")

    # Write output with existing MAF data only
    output_con <- file(output_file, "w")
    if (starts_with(maf_lines[1], "#version")) {
      writeLines(maf_lines[1], output_con)
    }
    writeLines(header_line, output_con)
    writeLines(somatic_lines, output_con)
    close(output_con)

    cat("Output file created with existing MAF data only.\n")
    quit(status = 0)
  }

  # Check column names in the CPSR file to determine format
  cpsr_columns <- colnames(cpsr_data)
  cat("CPSR file columns detected:", length(cpsr_columns), "\n")

  # Define column mapping for possible variations in CPSR format
  column_mapping <- list(
    entrez_id = c("ENTREZ_ID", "ENTREZGENE"),
    dbsnp_rs = c("DBSNP", "DBSNP_RSID"),
    refseq_mrna = c("REFSEQ_MRNA", "REFSEQ_TRANSCRIPT_ID"),
    protein_change = c("PROTEIN_CHANGE", "ALTERATION")
  )

  # Function to get column from mapping
  get_column <- function(mapping_key) {
    possible_cols <- column_mapping[[mapping_key]]
    for (col in possible_cols) {
      if (col %in% cpsr_columns) {
        return(col)
      }
    }
    return(NULL)  # Return NULL if no matching column found
  }

  # Determine which columns to use
  entrez_col <- get_column("entrez_id")
  dbsnp_col <- get_column("dbsnp_rs")
  refseq_col <- get_column("refseq_mrna")
  protein_change_col <- get_column("protein_change")

  cat("Using column mappings:",
      "\nEntrez:", entrez_col,
      "\ndbSNP:", dbsnp_col,
      "\nRefSeq:", refseq_col,
      "\nProtein Change:", protein_change_col, "\n")

  # Create empty MAF data frame with all necessary columns
  maf_entries <- data.frame(matrix(ncol=length(maf_columns), nrow=0))
  colnames(maf_entries) <- maf_columns

  # Process each CPSR entry
  for (i in 1:nrow(cpsr_data)) {
    cpsr_entry <- cpsr_data[i,]

    # Create a new MAF entry with all columns initialized as empty
    maf_entry <- as.list(rep("", length(maf_columns)))
    names(maf_entry) <- maf_columns

    # Parse VAR_ID (e.g., "1_45332803_T_C" or "2_201232808_TAGTAAG_T")
    var_parts <- strsplit(cpsr_entry$VAR_ID, "_")[[1]]
    chrom <- var_parts[1]
    pos <- as.integer(var_parts[2])
    ref <- var_parts[3]
    alt <- var_parts[4]

    # Basic mappings that are common to all variants
    maf_entry$Hugo_Symbol <- cpsr_entry$SYMBOL
    maf_entry$Entrez_Gene_Id <- if(!is.null(entrez_col)) cpsr_entry[[entrez_col]] else ""
    maf_entry$Center <- "."
    maf_entry$NCBI_Build <- "GRCh38"
    maf_entry$Chromosome <- chrom  # CPSR VAR_ID has no "chr" prefix; add it to match the vcf2maf rows
    if (!grepl("^chr", maf_entry$Chromosome)) {
      maf_entry$Chromosome <- paste0("chr", maf_entry$Chromosome)
    }
    maf_entry$Strand <- "+"

    # Handle dbSNP RS ID
    if (!is.null(dbsnp_col) && !is.na(cpsr_entry[[dbsnp_col]]) && cpsr_entry[[dbsnp_col]] != "NA" && cpsr_entry[[dbsnp_col]] != "") {
      maf_entry$dbSNP_RS <- cpsr_entry[[dbsnp_col]]
    } else {
      maf_entry$dbSNP_RS <- "novel"
    }

    # Sample information from existing MAF
    maf_entry$Tumor_Sample_Barcode <- tumor_sample
    maf_entry$Matched_Norm_Sample_Barcode <- normal_sample

    # Set mutation status to Germline as requested
    maf_entry$Mutation_Status <- "Germline"

    # Map variant type and coordinates
    variant_class <- cpsr_entry$VARIANT_CLASS

    if (variant_class == "SNV") {
      # For SNVs (Single Nucleotide Variants)
      maf_entry$Start_Position <- as.character(pos)
      maf_entry$End_Position <- as.character(pos)
      maf_entry$Reference_Allele <- ref
      maf_entry$Tumor_Seq_Allele1 <- ref
      maf_entry$Tumor_Seq_Allele2 <- alt
      maf_entry$Variant_Type <- "SNP"

    } else if (variant_class == "insertion") {
      # For insertions
      common_prefix_len <- 0
      for (j in 1:min(nchar(ref), nchar(alt))) {
        if (substr(ref, j, j) == substr(alt, j, j)) {
          common_prefix_len <- common_prefix_len + 1
        } else {
          break
        }
      }

      inserted <- substr(alt, common_prefix_len + 1, nchar(alt))
      maf_entry$Start_Position <- as.character(pos + common_prefix_len - 1)
      maf_entry$End_Position <- as.character(pos + common_prefix_len - 1)
      maf_entry$Reference_Allele <- "-"
      maf_entry$Tumor_Seq_Allele1 <- "-"
      maf_entry$Tumor_Seq_Allele2 <- inserted
      maf_entry$Variant_Type <- "INS"

    } else if (variant_class == "deletion") {
      # For deletions
      common_prefix_len <- 0
      for (j in 1:min(nchar(ref), nchar(alt))) {
        if (substr(ref, j, j) == substr(alt, j, j)) {
          common_prefix_len <- common_prefix_len + 1
        } else {
          break
        }
      }

      deleted <- substr(ref, common_prefix_len + 1, nchar(ref))
      maf_entry$Start_Position <- as.character(pos + common_prefix_len)
      maf_entry$End_Position <- as.character(pos + nchar(ref) - 1)
      maf_entry$Reference_Allele <- deleted
      maf_entry$Tumor_Seq_Allele1 <- deleted
      maf_entry$Tumor_Seq_Allele2 <- "-"
      maf_entry$Variant_Type <- "DEL"
    }

    # Map matched normal alleles
    maf_entry$Match_Norm_Seq_Allele1 <- maf_entry$Reference_Allele
    maf_entry$Match_Norm_Seq_Allele2 <- maf_entry$Reference_Allele

    # Set RNA field placeholder
    if ("t_depth_rna" %in% maf_columns) maf_entry$t_depth_rna <- ""
    if ("t_ref_count_rna" %in% maf_columns) maf_entry$t_ref_count_rna <- ""
    if ("t_alt_count_rna" %in% maf_columns) maf_entry$t_alt_count_rna <- ""
    if ("t_vaf_rna" %in% maf_columns) maf_entry$t_vaf_rna <- ""
    if ("Flag_RNA_Expressed" %in% maf_columns) maf_entry$Flag_RNA_Expressed <- ""

    # Map HGVSc and HGVSp
    if ("HGVSc" %in% cpsr_columns && !is.na(cpsr_entry$HGVSc) && cpsr_entry$HGVSc != "NA" && cpsr_entry$HGVSc != "") {
      maf_entry$HGVSc <- cpsr_entry$HGVSc
    }

    if ("HGVSp" %in% cpsr_columns && !is.na(cpsr_entry$HGVSp) && cpsr_entry$HGVSp != "NA" && cpsr_entry$HGVSp != "") {
      maf_entry$HGVSp <- cpsr_entry$HGVSp
      # Extract short form from HGVSp
      if (grepl(":p.", cpsr_entry$HGVSp)) {
        parts <- strsplit(cpsr_entry$HGVSp, ":p.")[[1]]
        maf_entry$HGVSp_Short <- paste0("p.", parts[2])
      }
    }

    # Try to get protein change from other columns if HGVSp not available
    if ((!is.null(protein_change_col)) &&
        (is.na(maf_entry$HGVSp) || maf_entry$HGVSp == "") &&
        !is.na(cpsr_entry[[protein_change_col]]) &&
        cpsr_entry[[protein_change_col]] != "") {
      maf_entry$HGVSp_Short <- cpsr_entry[[protein_change_col]]
    }

    # Map Variant_Classification based on CONSEQUENCE
    csq <- cpsr_entry$CONSEQUENCE
    if (grepl("missense_variant", csq)) {
      maf_entry$Variant_Classification <- "Missense_Mutation"
    } else if (grepl("synonymous_variant", csq)) {
      maf_entry$Variant_Classification <- "Silent"
    } else if (grepl("splice_donor_variant", csq) || grepl("splice_acceptor_variant", csq)) {
      maf_entry$Variant_Classification <- "Splice_Site"
    } else if (grepl("3_prime_UTR_variant", csq)) {
      maf_entry$Variant_Classification <- "3'UTR"
    } else if (grepl("5_prime_UTR_variant", csq)) {
      maf_entry$Variant_Classification <- "5'UTR"
    } else if (grepl("upstream_gene_variant", csq)) {
      maf_entry$Variant_Classification <- "5'Flank"
    } else if (grepl("downstream_gene_variant", csq)) {
      maf_entry$Variant_Classification <- "3'Flank"
    } else if (grepl("intron_variant", csq)) {
      maf_entry$Variant_Classification <- "Intron"
    } else if (grepl("nonsense", csq) || grepl("stop_gained", csq)) {
      maf_entry$Variant_Classification <- "Nonsense_Mutation"
    } else if (grepl("frameshift", csq)) {
      maf_entry$Variant_Classification <- "Frame_Shift"
    } else if (csq == "NA" || csq == "") {
      maf_entry$Variant_Classification <- "IGR"  # Default to intergenic if no consequence
    } else {
      maf_entry$Variant_Classification <- "IGR"  # Default to intergenic
    }

    # Set transcript info
    maf_entry$Transcript_ID <- cpsr_entry$ENSEMBL_TRANSCRIPT_ID

    # Add VEP annotations
    maf_entry$Allele <- alt
    maf_entry$Gene <- cpsr_entry$ENSEMBL_GENE_ID
    maf_entry$Feature <- cpsr_entry$ENSEMBL_TRANSCRIPT_ID
    maf_entry$Feature_type <- "Transcript"

    # Extract first consequence if multiple are present
    if (grepl(",", csq)) {
      consequences <- strsplit(csq, ",")[[1]]
      maf_entry$Consequence <- consequences[1]
    } else {
      maf_entry$Consequence <- csq
    }

    maf_entry$IMPACT <- "MODIFIER"  # Default impact
    maf_entry$FILTER <- "PASS"

    if (!is.na(cpsr_entry$CCDS) && cpsr_entry$CCDS != "NA" && cpsr_entry$CCDS != "") {
      maf_entry$CCDS <- cpsr_entry$CCDS
    }

    if (!is.null(refseq_col) && !is.na(cpsr_entry[[refseq_col]]) &&
        cpsr_entry[[refseq_col]] != "NA" && cpsr_entry[[refseq_col]] != "") {
      maf_entry$RefSeq <- cpsr_entry[[refseq_col]]
    }

    # Add to results
    maf_entries <- rbind(maf_entries, as.data.frame(as.list(maf_entry), stringsAsFactors=FALSE))
  }

  # Germline takes priority: a somatic MAF row that is the same variant as a CPSR germline call
  # is dropped, so only the CPSR (Germline) row remains. Matched on MAF-style coordinates
  # (chr prefix ignored in case either source lacks it).
  variant_key <- function(chrom, start, ref, alt) {
    paste(sub("^chr", "", chrom), start, ref, alt, sep = "_")
  }
  germline_keys <- variant_key(maf_entries$Chromosome, maf_entries$Start_Position,
                               maf_entries$Reference_Allele, maf_entries$Tumor_Seq_Allele2)
  key_idx <- match(c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"), maf_columns)

  if (!anyNA(key_idx) && length(somatic_lines) > 0) {
    is_germline <- vapply(somatic_lines, function(line) {
      fields <- split_maf_line(line)
      variant_key(fields[key_idx[1]], fields[key_idx[2]], fields[key_idx[3]], fields[key_idx[4]]) %in% germline_keys
    }, logical(1), USE.NAMES = FALSE)
    somatic_lines <- somatic_lines[!is_germline]
    cat("Dropped", sum(is_germline), "somatic MAF rows that are the same variant as a CPSR germline call\n")
  }

  # Prepare to write the output
  # First, determine if there's a version line to preserve
  has_version <- FALSE
  version_line <- ""
  if (starts_with(maf_lines[1], "#version")) {
    has_version <- TRUE
    version_line <- maf_lines[1]
  }

  # Open output file for writing
  output_con <- file(output_file, "w")

  # Write version line if it exists
  if (has_version) {
    writeLines(version_line, output_con)
  }

  # Prepare column header line
  if (starts_with(header_line, "#")) {
    header_to_write <- header_line
  } else {
    header_to_write <- paste(maf_columns, collapse="\t")
    #header_to_write <- paste0("#", paste(maf_columns, collapse="\t"))
  }
  writeLines(header_to_write, output_con)

  # Write existing MAF data (skipping header)
  writeLines(somatic_lines, output_con)

  # Write new entries
  write.table(maf_entries, output_con, sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE, append=TRUE)

  close(output_con)

  cat("Conversion complete. Added", nrow(maf_entries), "germline variants to", output_file, "\n")
}

# Parse command line arguments if script is run directly
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 3) {
    cat("Usage: Rscript cpsr_to_maf.R <cpsr_file> <maf_file> <output_file>\n")
    quit(status = 1)
  }

  cpsr_file <- args[1]
  maf_file <- args[2]
  output_file <- args[3]

  convert_cpsr_to_maf(cpsr_file, maf_file, output_file)
}
