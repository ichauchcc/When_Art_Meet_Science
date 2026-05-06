#!/usr/bin/env Rscript

# ============================================================
# genome_art_clinvar.R
#
# Purpose:
#   Create artistic genome-based visualisation from ClinVar-style
#   human variant data.
#
# Expected input columns:
#   chr, pos, ref, alt, gene, consequence, classification
#
# Example usage:
#   Rscript Variant_art_clinvar.R clinvar_variants.tsv genome_art.png
#
# Author:
#   YuChen C / genomArtR draft
# ============================================================


# ------------------------------------------------------------
# 0. Load required packages
# ------------------------------------------------------------

required_packages <- c(
  "tidyverse",
  "ggforce",
  "scales"
)

install_missing_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]

  if (length(missing_pkgs) > 0) {
    message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
  }
}

install_missing_packages(required_packages)

library(tidyverse)
library(ggforce)
library(scales)


# ------------------------------------------------------------
# 1. Human chromosome lengths, GRCh38
# ------------------------------------------------------------

chr_lengths_grch38 <- tibble::tribble(
  ~chr,   ~length,
  "chr1",  248956422,
  "chr2",  242193529,
  "chr3",  198295559,
  "chr4",  190214555,
  "chr5",  181538259,
  "chr6",  170805979,
  "chr7",  159345973,
  "chr8",  145138636,
  "chr9",  138394717,
  "chr10", 133797422,
  "chr11", 135086622,
  "chr12", 133275309,
  "chr13", 114364328,
  "chr14", 107043718,
  "chr15", 101991189,
  "chr16", 90338345,
  "chr17", 83257441,
  "chr18", 80373285,
  "chr19", 58617616,
  "chr20", 64444167,
  "chr21", 46709983,
  "chr22", 50818468,
  "chrX",  156040895,
  "chrY",  57227415
)


# ------------------------------------------------------------
# 2. Read ClinVar-style variant file
# ------------------------------------------------------------

read_clinvar_variant_table <- function(file) {

  if (!file.exists(file)) {
    stop("Input file does not exist: ", file)
  }

  message("Reading input file: ", file)

  # Automatically read TSV or CSV
  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    variants <- readr::read_csv(file, show_col_types = FALSE)
  } else {
    variants <- readr::read_tsv(file, show_col_types = FALSE)
  }

  required_cols <- c(
    "chr",
    "pos",
    "ref",
    "alt",
    "gene",
    "consequence",
    "classification"
  )

  missing_cols <- setdiff(required_cols, colnames(variants))

  if (length(missing_cols) > 0) {
    stop(
      "Input file is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  variants <- variants %>%
    mutate(
      chr = as.character(chr),
      pos = as.numeric(pos),
      ref = as.character(ref),
      alt = as.character(alt),
      gene = as.character(gene),
      consequence = as.character(consequence),
      classification = as.character(classification)
    ) %>%
    filter(!is.na(chr), !is.na(pos))

  message("Loaded ", nrow(variants), " variants.")

  return(variants)
}


# ------------------------------------------------------------
# 3. Helper function: standardise chromosome names
# ------------------------------------------------------------

standardise_chromosomes <- function(variants) {

  variants %>%
    mutate(
      chr = as.character(chr),
      chr = str_replace(chr, "^chr", ""),
      chr = case_when(
        chr %in% as.character(1:22) ~ paste0("chr", chr),
        chr %in% c("X", "x") ~ "chrX",
        chr %in% c("Y", "y") ~ "chrY",
        chr %in% c("M", "MT", "m", "mt") ~ "chrM",
        TRUE ~ chr
      )
    )
}


# ------------------------------------------------------------
# 4. Helper function: classify variant impact for artwork
# ------------------------------------------------------------

classify_variant_impact <- function(consequence) {

  consequence <- tolower(consequence)

  case_when(
    str_detect(
      consequence,
      "frameshift|stop_gained|stop_lost|start_lost|splice_acceptor|splice_donor"
    ) ~ "High impact",

    str_detect(
      consequence,
      "missense|inframe_deletion|inframe_insertion|protein_altering"
    ) ~ "Protein altering",

    str_detect(
      consequence,
      "synonymous"
    ) ~ "Synonymous",

    str_detect(
      consequence,
      "intron|intergenic|upstream|downstream|utr"
    ) ~ "Non-coding / regulatory",

    TRUE ~ "Other"
  )
}


# ------------------------------------------------------------
# 5. Helper function: simplify ClinVar classification
# ------------------------------------------------------------

simplify_classification <- function(classification) {

  classification <- tolower(classification)

  case_when(
    str_detect(classification, "pathogenic") &
      !str_detect(classification, "conflicting|uncertain") ~ "Pathogenic / Likely pathogenic",

    str_detect(classification, "benign") &
      !str_detect(classification, "conflicting|uncertain") ~ "Benign / Likely benign",

    str_detect(classification, "uncertain|vus") ~ "Uncertain significance",

    str_detect(classification, "conflicting") ~ "Conflicting interpretations",

    TRUE ~ "Other / not provided"
  )
}


# ------------------------------------------------------------
# 6. Main function to create genome art
# ------------------------------------------------------------

create_genome_art <- function(
    variants,
    chr_lengths = chr_lengths_grch38,
    title = "ClinVar Genome Constellation",
    subtitle = "Variants mapped across the GRCh38 human genome",
    seed = 1571,
    highlight_genes = c(
      "TTN", "SCN5A", "RYR2", "MYH7", "MYBPC3",
      "LMNA", "DSP", "PKP2", "BAG3", "FLNC", "DES", "PLN"
    ),
    background = "black",
    show_chr_labels = TRUE
) {

  set.seed(seed)

  # Standardise and filter variants
  variants_clean <- variants %>%
    standardise_chromosomes() %>%
    filter(chr %in% chr_lengths$chr) %>%
    left_join(chr_lengths, by = "chr") %>%
    filter(pos >= 1, pos <= length)

  if (nrow(variants_clean) == 0) {
    stop("No valid variants remain after chromosome and position filtering.")
  }

  message("Using ", nrow(variants_clean), " variants for artwork.")

  # Create genome-wide coordinate layout
  genome_layout <- chr_lengths %>%
    mutate(
      chr = factor(chr, levels = chr_lengths$chr),
      cum_start = lag(cumsum(length), default = 0),
      cum_end = cum_start + length,
      midpoint = cum_start + length / 2
    )

  genome_size <- sum(chr_lengths$length)

  # Prepare chromosome ring data
  ring_df <- genome_layout %>%
    mutate(
      theta_start = 2 * pi * cum_start / genome_size,
      theta_end = 2 * pi * cum_end / genome_size,
      theta_mid = 2 * pi * midpoint / genome_size,
      label_x = 1.24 * cos(theta_mid),
      label_y = 1.24 * sin(theta_mid)
    )

  ring_segments <- purrr::map_dfr(seq_len(nrow(ring_df)), function(i) {
    theta_seq <- seq(
      ring_df$theta_start[i],
      ring_df$theta_end[i],
      length.out = 120
    )

    tibble(
      chr = ring_df$chr[i],
      theta = theta_seq,
      x_outer = cos(theta_seq),
      y_outer = sin(theta_seq),
      x_inner = 0.72 * cos(theta_seq),
      y_inner = 0.72 * sin(theta_seq)
    )
  })

  # Prepare variant coordinates
  plot_df <- variants_clean %>%
    left_join(
      genome_layout %>% select(chr, cum_start, cum_end, midpoint),
      by = "chr"
    ) %>%
    mutate(
      genome_pos = cum_start + pos,
      theta = 2 * pi * genome_pos / genome_size,

      impact_group = classify_variant_impact(consequence),
      clinvar_group = simplify_classification(classification),

      is_highlight_gene = gene %in% highlight_genes,

      base_radius = case_when(
        clinvar_group == "Pathogenic / Likely pathogenic" ~ 1.05,
        is_highlight_gene ~ 0.98,
        impact_group == "High impact" ~ 0.92,
        impact_group == "Protein altering" ~ 0.86,
        TRUE ~ 0.80
      ),

      radius = base_radius + rnorm(n(), mean = 0, sd = 0.025),

      x = radius * cos(theta),
      y = radius * sin(theta),

      visual_group = case_when(
        is_highlight_gene ~ gene,
        TRUE ~ "Other genes"
      ),

      point_size = case_when(
        clinvar_group == "Pathogenic / Likely pathogenic" ~ 3.8,
        clinvar_group == "Conflicting interpretations" ~ 2.8,
        clinvar_group == "Uncertain significance" ~ 2.1,
        impact_group == "High impact" ~ 2.4,
        impact_group == "Protein altering" ~ 1.6,
        TRUE ~ 0.8
      ),

      alpha_value = case_when(
        clinvar_group == "Pathogenic / Likely pathogenic" ~ 0.95,
        clinvar_group == "Conflicting interpretations" ~ 0.80,
        clinvar_group == "Uncertain significance" ~ 0.65,
        TRUE ~ 0.35
      )
    )

  # Create plot
  p <- ggplot() +

    # Inner orbit
    geom_path(
      data = ring_segments,
      aes(x = x_inner, y = y_inner, group = chr),
      linewidth = 0.25,
      alpha = 0.25,
      colour = "white"
    ) +

    # Outer orbit
    geom_path(
      data = ring_segments,
      aes(x = x_outer, y = y_outer, group = chr),
      linewidth = 0.35,
      alpha = 0.35,
      colour = "white"
    ) +

    # Soft chromosome radial dividers
    geom_segment(
      data = ring_df,
      aes(
        x = 0.70 * cos(theta_start),
        y = 0.70 * sin(theta_start),
        xend = 1.02 * cos(theta_start),
        yend = 1.02 * sin(theta_start)
      ),
      linewidth = 0.15,
      alpha = 0.18,
      colour = "white"
    ) +

    # Variant points
    geom_point(
      data = plot_df,
      aes(
        x = x,
        y = y,
        size = point_size,
        alpha = alpha_value,
        colour = visual_group,
        shape = clinvar_group
      )
    ) +

    # Optional chromosome labels
    {
      if (show_chr_labels) {
        geom_text(
          data = ring_df,
          aes(
            x = label_x,
            y = label_y,
            label = str_remove(as.character(chr), "chr")
          ),
          colour = "white",
          size = 3,
          alpha = 0.78
        )
      }
    } +

    coord_equal() +

    scale_size_identity() +
    scale_alpha_identity() +

    scale_shape_manual(
      values = c(
        "Pathogenic / Likely pathogenic" = 8,
        "Uncertain significance" = 16,
        "Conflicting interpretations" = 17,
        "Benign / Likely benign" = 1,
        "Other / not provided" = 3
      ),
      drop = FALSE
    ) +

    guides(
      colour = guide_legend(override.aes = list(size = 4, alpha = 1)),
      shape = guide_legend(override.aes = list(size = 4, alpha = 1))
    ) +

    labs(
      title = title,
      subtitle = subtitle,
      colour = "Gene highlight",
      shape = "ClinVar classification"
    ) +

    theme_void(base_size = 14) +

    theme(
      plot.background = element_rect(fill = background, colour = NA),
      panel.background = element_rect(fill = background, colour = NA),

      plot.title = element_text(
        colour = "white",
        hjust = 0.5,
        size = 24,
        face = "bold",
        margin = margin(b = 6)
      ),

      plot.subtitle = element_text(
        colour = "grey85",
        hjust = 0.5,
        size = 11,
        margin = margin(b = 16)
      ),

      legend.position = "right",
      legend.background = element_rect(fill = background, colour = NA),
      legend.key = element_rect(fill = background, colour = NA),
      legend.title = element_text(colour = "white", size = 10),
      legend.text = element_text(colour = "grey90", size = 8),

      plot.margin = margin(20, 20, 20, 20)
    )

  return(p)
}


# ------------------------------------------------------------
# 7. Export function
# ------------------------------------------------------------

export_genome_art <- function(
    plot,
    output_png = "clinvar_genome_art.png",
    output_pdf = "clinvar_genome_art.pdf",
    width = 12,
    height = 12,
    dpi = 600,
    background = "black"
) {

  message("Saving PNG: ", output_png)

  ggsave(
    filename = output_png,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = background
  )

  message("Saving PDF: ", output_pdf)

  ggsave(
    filename = output_pdf,
    plot = plot,
    width = width,
    height = height,
    bg = background
  )

  message("Export complete.")
}


# ------------------------------------------------------------
# 8. Run from command line or interactively
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {

  input_file <- args[1]

  output_png <- ifelse(
    length(args) >= 2,
    args[2],
    "clinvar_genome_art.png"
  )

  output_pdf <- str_replace(output_png, "\\.png$", ".pdf")

  variants <- read_clinvar_variant_table(input_file)

  p <- create_genome_art(
    variants = variants,
    title = "ClinVar Genome Constellation",
    subtitle = "Clinical variant classifications arranged across the GRCh38 human genome",
    seed = 1571
  )

  export_genome_art(
    plot = p,
    output_png = output_png,
    output_pdf = output_pdf,
    width = 12,
    height = 12,
    dpi = 600,
    background = "black"
  )

} else {

  message("No input file provided.")
  message("Creating an example mock artwork instead.")

  set.seed(1571)

  example_variants <- tibble(
    chr = sample(chr_lengths_grch38$chr, 2000, replace = TRUE),
    pos = purrr::map_dbl(
      chr,
      ~ sample(
        1:chr_lengths_grch38$length[chr_lengths_grch38$chr == .x],
        1
      )
    ),
    ref = sample(c("A", "C", "G", "T"), 2000, replace = TRUE),
    alt = sample(c("A", "C", "G", "T"), 2000, replace = TRUE),
    gene = sample(
      c(
        "TTN", "SCN5A", "RYR2", "MYH7", "MYBPC3",
        "LMNA", "DSP", "PKP2", "BAG3", "FLNC",
        "DES", "PLN", "Other"
      ),
      2000,
      replace = TRUE,
      prob = c(rep(0.025, 12), 0.70)
    ),
    consequence = sample(
      c(
        "missense_variant",
        "synonymous_variant",
        "frameshift_variant",
        "stop_gained",
        "splice_donor_variant",
        "intron_variant"
      ),
      2000,
      replace = TRUE
    ),
    classification = sample(
      c(
        "Pathogenic",
        "Likely pathogenic",
        "Uncertain significance",
        "Conflicting interpretations of pathogenicity",
        "Likely benign",
        "Benign"
      ),
      2000,
      replace = TRUE,
      prob = c(0.04, 0.06, 0.30, 0.10, 0.20, 0.30)
    )
  )

  p <- create_genome_art(
    variants = example_variants,
    title = "Mock ClinVar Genome Constellation",
    subtitle = "Example artwork generated from simulated ClinVar-style variant data",
    seed = 1571
  )

  print(p)

  export_genome_art(
    plot = p,
    output_png = "mock_clinvar_genome_art.png",
    output_pdf = "mock_clinvar_genome_art.pdf",
    width = 12,
    height = 12,
    dpi = 600,
    background = "black"
  )
}
