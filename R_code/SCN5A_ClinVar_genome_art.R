#!/usr/bin/env Rscript

# ============================================================
# SCN5A ClinVar Variant Art
#
# Input:
#   clinvar_SCN5A.tsv
#
# Required columns:
#   chr, pos, ref, alt, gene, consequence, classification,
#   conflicting_classification_detail
#
# Output:
#   SCN5A_ClinVar_art.png
#   SCN5A_ClinVar_art.pdf
#
# Concept:
#   Each ClinVar variant in SCN5A is represented as a point along
#   the SCN5A genomic locus. Classification controls colour;
#   consequence controls shape; conflicting classifications are
#   drawn as an outer "uncertainty halo".
# ============================================================


# ------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------

library(tidyverse)
library(ggforce)
library(scales)


# ------------------------------------------------------------
# 1. Input file
# ------------------------------------------------------------

input_file <- "clinvar_SCN5A.tsv"

if (!file.exists(input_file)) {
  stop("Cannot find input file: ", input_file)
}


# ------------------------------------------------------------
# 2. Read SCN5A ClinVar variants
# ------------------------------------------------------------

scn5a <- readr::read_tsv(
  input_file,
  show_col_types = FALSE
)

required_cols <- c(
  "chr",
  "pos",
  "ref",
  "alt",
  "gene",
  "consequence",
  "classification",
  "conflicting_classification_detail"
)

missing_cols <- setdiff(required_cols, colnames(scn5a))

if (length(missing_cols) > 0) {
  stop(
    "Input file is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

scn5a <- scn5a %>%
  mutate(
    chr = as.character(chr),
    chr = if_else(str_detect(chr, "^chr"), chr, paste0("chr", chr)),
    pos = as.numeric(pos),
    ref = as.character(ref),
    alt = as.character(alt),
    gene = as.character(gene),
    consequence = as.character(consequence),
    classification = as.character(classification),
    conflicting_classification_detail = as.character(conflicting_classification_detail)
  ) %>%
  filter(
    gene == "SCN5A",
    !is.na(pos)
  )

message("Loaded ", nrow(scn5a), " SCN5A variants.")


# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

classify_consequence <- function(x) {
  x <- tolower(x)

  case_when(
    str_detect(
      x,
      "frameshift|nonsense|stop_gained|stop_lost|start_lost|splice_acceptor|splice_donor"
    ) ~ "High-impact / LoF",

    str_detect(
      x,
      "missense|inframe_deletion|inframe_insertion|inframe_indel"
    ) ~ "Protein-altering",

    str_detect(x, "synonymous") ~ "Synonymous",

    str_detect(
      x,
      "intron|utr|non-coding|non_coding|3_prime|5_prime"
    ) ~ "Non-coding / transcript",

    TRUE ~ "Other"
  )
}


simplify_clinvar <- function(x) {
  x <- as.character(x)

  case_when(
    x %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic") ~
      "Pathogenic / likely pathogenic",

    x %in% c("Uncertain_significance") ~
      "Uncertain significance",

    x %in% c("Conflicting_classifications_of_pathogenicity") ~
      "Conflicting classifications",

    x %in% c("Benign", "Likely_benign", "Benign/Likely_benign") ~
      "Benign / likely benign",

    x %in% c("not_provided", ".") ~
      "Not provided",

    TRUE ~ "Other"
  )
}


variant_type_from_ref_alt <- function(ref, alt) {
  case_when(
    nchar(ref) == 1 & nchar(alt) == 1 ~ "SNV",
    nchar(ref) > nchar(alt) ~ "Deletion",
    nchar(ref) < nchar(alt) ~ "Insertion",
    TRUE ~ "Complex"
  )
}


# ------------------------------------------------------------
# 4. Prepare artistic coordinates
# ------------------------------------------------------------

scn5a_art <- scn5a %>%
  mutate(
    clinvar_group = simplify_clinvar(classification),
    consequence_group = classify_consequence(consequence),
    variant_type = variant_type_from_ref_alt(ref, alt),
    is_conflicting = clinvar_group == "Conflicting classifications",

    locus_start = min(pos, na.rm = TRUE),
    locus_end = max(pos, na.rm = TRUE),
    locus_fraction = (pos - locus_start) / (locus_end - locus_start),

    # Convert genomic position into an angle.
    theta = 2 * pi * locus_fraction,

    # Radial position is partly biological and partly aesthetic.
    base_radius = case_when(
      clinvar_group == "Pathogenic / likely pathogenic" ~ 1.08,
      clinvar_group == "Conflicting classifications" ~ 1.00,
      clinvar_group == "Uncertain significance" ~ 0.92,
      clinvar_group == "Benign / likely benign" ~ 0.82,
      TRUE ~ 0.74
    ),

    # Deterministic gentle jitter for visual separation.
    radius = base_radius + rnorm(n(), mean = 0, sd = 0.025),

    x = radius * cos(theta),
    y = radius * sin(theta),

    point_size = case_when(
      clinvar_group == "Pathogenic / likely pathogenic" ~ 3.5,
      clinvar_group == "Conflicting classifications" ~ 3.0,
      clinvar_group == "Uncertain significance" ~ 2.0,
      clinvar_group == "Benign / likely benign" ~ 1.2,
      TRUE ~ 0.9
    ),

    alpha_value = case_when(
      clinvar_group == "Pathogenic / likely pathogenic" ~ 0.95,
      clinvar_group == "Conflicting classifications" ~ 0.90,
      clinvar_group == "Uncertain significance" ~ 0.65,
      clinvar_group == "Benign / likely benign" ~ 0.35,
      TRUE ~ 0.25
    )
  )


# ------------------------------------------------------------
# 5. Create SCN5A locus rings
# ------------------------------------------------------------

ring_df <- tibble(
  theta = seq(0, 2 * pi, length.out = 1000),
  outer_x = cos(theta),
  outer_y = sin(theta),
  middle_x = 0.92 * cos(theta),
  middle_y = 0.92 * sin(theta),
  inner_x = 0.74 * cos(theta),
  inner_y = 0.74 * sin(theta)
)

# Add genomic tick labels around the circle
tick_df <- tibble(
  pos = pretty(scn5a_art$pos, n = 8)
) %>%
  filter(
    pos >= min(scn5a_art$pos),
    pos <= max(scn5a_art$pos)
  ) %>%
  mutate(
    locus_fraction = (pos - min(scn5a_art$pos)) /
      (max(scn5a_art$pos) - min(scn5a_art$pos)),
    theta = 2 * pi * locus_fraction,
    x_start = 1.13 * cos(theta),
    y_start = 1.13 * sin(theta),
    x_end = 1.18 * cos(theta),
    y_end = 1.18 * sin(theta),
    label_x = 1.29 * cos(theta),
    label_y = 1.29 * sin(theta),
    label = paste0(round(pos / 1e6, 2), " Mb")
  )


# ------------------------------------------------------------
# 6. Make artwork
# ------------------------------------------------------------

set.seed(1571)

p <- ggplot() +

  # Main SCN5A locus rings
  geom_path(
    data = ring_df,
    aes(x = outer_x, y = outer_y),
    colour = "white",
    linewidth = 0.4,
    alpha = 0.30
  ) +

  geom_path(
    data = ring_df,
    aes(x = middle_x, y = middle_y),
    colour = "white",
    linewidth = 0.25,
    alpha = 0.18
  ) +

  geom_path(
    data = ring_df,
    aes(x = inner_x, y = inner_y),
    colour = "white",
    linewidth = 0.25,
    alpha = 0.14
  ) +

  # Genomic coordinate ticks
  geom_segment(
    data = tick_df,
    aes(
      x = x_start,
      y = y_start,
      xend = x_end,
      yend = y_end
    ),
    colour = "white",
    linewidth = 0.25,
    alpha = 0.40
  ) +

  geom_text(
    data = tick_df,
    aes(
      x = label_x,
      y = label_y,
      label = label
    ),
    colour = "grey80",
    size = 2.6,
    alpha = 0.80
  ) +

  # Halo for conflicting variants
  geom_point(
    data = scn5a_art %>% filter(is_conflicting),
    aes(x = x, y = y),
    shape = 21,
    size = 5.2,
    stroke = 0.45,
    colour = "white",
    fill = NA,
    alpha = 0.55
  ) +

  # Variant points
  geom_point(
    data = scn5a_art,
    aes(
      x = x,
      y = y,
      colour = clinvar_group,
      shape = consequence_group,
      size = point_size,
      alpha = alpha_value
    )
  ) +

  # Centre title text
  annotate(
    "text",
    x = 0,
    y = 0.08,
    label = "SCN5A",
    colour = "white",
    size = 10,
    fontface = "bold"
  ) +

  annotate(
    "text",
    x = 0,
    y = -0.06,
    label = paste0(
      nrow(scn5a_art),
      " ClinVar variants\nchr3:",
      format(min(scn5a_art$pos), big.mark = ","),
      "–",
      format(max(scn5a_art$pos), big.mark = ",")
    ),
    colour = "grey85",
    size = 3.2,
    lineheight = 1.05
  ) +

  coord_equal(clip = "off") +

  scale_size_identity() +
  scale_alpha_identity() +

  scale_colour_manual(
    values = c(
      "Pathogenic / likely pathogenic" = "#ff4d6d",
      "Conflicting classifications" = "#ffd166",
      "Uncertain significance" = "#7bdff2",
      "Benign / likely benign" = "#80ed99",
      "Not provided" = "#adb5bd",
      "Other" = "#dee2e6"
    ),
    drop = FALSE
  ) +

  scale_shape_manual(
    values = c(
      "High-impact / LoF" = 8,
      "Protein-altering" = 16,
      "Synonymous" = 1,
      "Non-coding / transcript" = 3,
      "Other" = 4
    ),
    drop = FALSE
  ) +

  guides(
    colour = guide_legend(
      title = "ClinVar classification",
      override.aes = list(size = 4, alpha = 1)
    ),
    shape = guide_legend(
      title = "Variant consequence",
      override.aes = list(size = 4, alpha = 1)
    )
  ) +

  labs(
    title = "SCN5A ClinVar Variant Constellation",
    subtitle = "Clinical classifications and variant consequences mapped across the SCN5A locus",
    caption = "Data source: ClinVar-derived SCN5A variant table | Artwork generated in R"
  ) +

  theme_void(base_size = 14) +

  theme(
    plot.background = element_rect(fill = "black", colour = NA),
    panel.background = element_rect(fill = "black", colour = NA),

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
      margin = margin(b = 18)
    ),

    plot.caption = element_text(
      colour = "grey60",
      hjust = 0.5,
      size = 8,
      margin = margin(t = 16)
    ),

    legend.position = "right",
    legend.background = element_rect(fill = "black", colour = NA),
    legend.key = element_rect(fill = "black", colour = NA),
    legend.title = element_text(colour = "white", size = 10),
    legend.text = element_text(colour = "grey90", size = 8),

    plot.margin = margin(30, 40, 30, 30)
  )


# ------------------------------------------------------------
# 7. Display and export high-resolution artwork
# ------------------------------------------------------------

print(p)

ggsave(
  filename = "SCN5A_ClinVar_art.png",
  plot = p,
  width = 12,
  height = 12,
  dpi = 600,
  bg = "black"
)

ggsave(
  filename = "SCN5A_ClinVar_art.pdf",
  plot = p,
  width = 12,
  height = 12,
  bg = "black"
)

message("Saved: SCN5A_ClinVar_art.png")
message("Saved: SCN5A_ClinVar_art.pdf")


# ------------------------------------------------------------
# 8. Optional summary table
# ------------------------------------------------------------

summary_table <- scn5a_art %>%
  count(clinvar_group, consequence_group, sort = TRUE)

readr::write_tsv(
  summary_table,
  "SCN5A_ClinVar_art_summary.tsv"
)

message("Saved: SCN5A_ClinVar_art_summary.tsv")
