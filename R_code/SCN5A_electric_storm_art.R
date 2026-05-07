#!/usr/bin/env Rscript

# ============================================================
# SCN5A Electric Storm Art
#
# Generative artwork from ClinVar SCN5A variants.
#
# Input:
#   clinvar_SCN5A.tsv
#
# Required columns:
#   chr, pos, ref, alt, gene, consequence, classification,
#   conflicting_classification_detail
#
# Output:
#   SCN5A_electric_storm.png
#   SCN5A_electric_storm.pdf
#
# Artistic concept:
#   SCN5A encodes the cardiac sodium channel Nav1.5.
#   This artwork treats variants as charged particles moving
#   through an abstract electrical field.
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_packages <- c(
  "tidyverse",
  "ggforce",
  "ggfx",
  "ambient",
  "scales"
)

install_missing_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]

  if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
  }
}

install_missing_packages(required_packages)

library(tidyverse)
library(ggforce)
library(ggfx)
library(ambient)
library(scales)


# ------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------

input_file <- "clinvar_SCN5A.tsv"

output_png <- "SCN5A_electric_storm.png"
output_pdf <- "SCN5A_electric_storm.pdf"

seed <- 1571

set.seed(seed)


# ------------------------------------------------------------
# 2. Read data
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Cannot find input file: ", input_file)
}

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


variant_type_from_ref_alt <- function(ref, alt) {
  case_when(
    nchar(ref) == 1 & nchar(alt) == 1 ~ "SNV",
    nchar(ref) > nchar(alt) ~ "Deletion",
    nchar(ref) < nchar(alt) ~ "Insertion",
    TRUE ~ "Complex"
  )
}


# ------------------------------------------------------------
# 4. Prepare variant data for generative artwork
# ------------------------------------------------------------

pos_min <- min(scn5a$pos, na.rm = TRUE)
pos_max <- max(scn5a$pos, na.rm = TRUE)

scn5a_art <- scn5a %>%
  mutate(
    clinvar_group = simplify_clinvar(classification),
    consequence_group = classify_consequence(consequence),
    variant_type = variant_type_from_ref_alt(ref, alt),

    locus_fraction = (pos - pos_min) / (pos_max - pos_min),

    # Convert genomic coordinate to a base x-axis.
    x_base = rescale(locus_fraction, to = c(-5.8, 5.8)),

    # Random artistic displacement.
    # This keeps genomic ordering but makes the image organic.
    charge = case_when(
      clinvar_group == "Pathogenic / likely pathogenic" ~ 1.00,
      clinvar_group == "Conflicting classifications" ~ 0.75,
      clinvar_group == "Uncertain significance" ~ 0.55,
      clinvar_group == "Benign / likely benign" ~ 0.22,
      TRUE ~ 0.15
    ),

    y_wave = sin(x_base * 1.3) +
      0.5 * sin(x_base * 2.7) +
      0.25 * cos(x_base * 5.1),

    y_base = y_wave * charge,

    x = x_base + rnorm(n(), 0, 0.09 + charge * 0.10),
    y = y_base + rnorm(n(), 0, 0.25 + charge * 0.45),

    particle_size = case_when(
      clinvar_group == "Pathogenic / likely pathogenic" ~ runif(n(), 2.7, 6.2),
      clinvar_group == "Conflicting classifications" ~ runif(n(), 2.0, 4.8),
      clinvar_group == "Uncertain significance" ~ runif(n(), 1.1, 3.2),
      clinvar_group == "Benign / likely benign" ~ runif(n(), 0.4, 1.4),
      TRUE ~ runif(n(), 0.3, 1.1)
    ),

    alpha_value = case_when(
      clinvar_group == "Pathogenic / likely pathogenic" ~ runif(n(), 0.75, 1.00),
      clinvar_group == "Conflicting classifications" ~ runif(n(), 0.60, 0.90),
      clinvar_group == "Uncertain significance" ~ runif(n(), 0.35, 0.70),
      clinvar_group == "Benign / likely benign" ~ runif(n(), 0.12, 0.35),
      TRUE ~ runif(n(), 0.08, 0.22)
    )
  )


# ------------------------------------------------------------
# 5. Create an ambient noise background
# ------------------------------------------------------------

# Grid for background electric field
grid_n <- 420

field <- expand_grid(
  x = seq(-6.4, 6.4, length.out = grid_n),
  y = seq(-4.4, 4.4, length.out = grid_n)
) %>%
  mutate(
    noise1 = ambient::fracture(
      ambient::gen_simplex,
      x = x,
      y = y,
      frequency = 0.55,
      fractal = "fbm",
      octaves = 5
    ),
    noise2 = ambient::fracture(
      ambient::gen_simplex,
      x = x + 100,
      y = y - 100,
      frequency = 1.1,
      fractal = "billow",
      octaves = 4
    ),
    distance_to_channel = abs(
      y - (
        sin(x * 1.3) +
          0.5 * sin(x * 2.7) +
          0.25 * cos(x * 5.1)
      )
    ),
    current = exp(-distance_to_channel^2 / 0.65) *
      rescale(noise1 + 0.65 * noise2, to = c(0, 1))
  )


# ------------------------------------------------------------
# 6. Create flowing current lines
# ------------------------------------------------------------

n_lines <- 190

current_lines <- map_dfr(seq_len(n_lines), function(i) {
  x_seq <- seq(-6.2, 6.2, length.out = 300)

  phase <- runif(1, -pi, pi)
  offset <- rnorm(1, 0, 1.2)
  amplitude <- runif(1, 0.05, 0.35)

  tibble(
    line_id = i,
    x = x_seq,
    y = sin(x_seq * 1.3 + phase) +
      0.5 * sin(x_seq * 2.7 + phase / 2) +
      0.25 * cos(x_seq * 5.1 - phase) +
      offset +
      amplitude * sin(x_seq * runif(1, 3, 9)),
    alpha_line = runif(1, 0.03, 0.16),
    line_width = runif(1, 0.10, 0.55)
  )
})


# ------------------------------------------------------------
# 7. Create variant particle trails
# ------------------------------------------------------------

# Select stronger variants for trails.
trail_variants <- scn5a_art %>%
  filter(
    clinvar_group %in% c(
      "Pathogenic / likely pathogenic",
      "Conflicting classifications",
      "Uncertain significance"
    )
  ) %>%
  slice_sample(n = min(450, n()))

trails <- pmap_dfr(
  list(
    trail_variants$x,
    trail_variants$y,
    trail_variants$clinvar_group,
    seq_len(nrow(trail_variants))
  ),
  function(x0, y0, group, id) {

    n_steps <- sample(12:32, 1)

    tibble(
      trail_id = id,
      step = seq_len(n_steps),
      x = x0 + cumsum(rnorm(n_steps, 0.025, 0.050)),
      y = y0 + cumsum(rnorm(n_steps, 0.000, 0.070)),
      clinvar_group = group,
      alpha_trail = seq(0.02, 0.22, length.out = n_steps)
    )
  }
)


# ------------------------------------------------------------
# 8. Create decorative arcs around the channel
# ------------------------------------------------------------

arc_data <- tibble(
  arc_id = seq_len(80),
  x0 = runif(80, -6.0, 6.0),
  y0 = runif(80, -3.8, 3.8),
  r = runif(80, 0.15, 1.4),
  start = runif(80, 0, 2 * pi),
  end = start + runif(80, pi / 8, pi * 1.2),
  alpha_arc = runif(80, 0.025, 0.16),
  line_width = runif(80, 0.10, 0.50)
)


# ------------------------------------------------------------
# 9. Make the artwork
# ------------------------------------------------------------

p <- ggplot() +

  # Deep electric field background
  geom_raster(
    data = field,
    aes(x = x, y = y, fill = current),
    alpha = 0.95,
    interpolate = TRUE
  ) +

  scale_fill_gradientn(
    colours = c(
      "#030014",
      "#09002b",
      "#15104d",
      "#102a55",
      "#0a4f64",
      "#1f7a8c",
      "#7bdff2"
    ),
    guide = "none"
  ) +

  # Decorative arcs
  with_blur(
    geom_arc(
      data = arc_data,
      aes(
        x0 = x0,
        y0 = y0,
        r = r,
        start = start,
        end = end,
        alpha = alpha_arc,
        linewidth = line_width
      ),
      colour = "#7bdff2"
    ),
    sigma = 2
  ) +

  # Current lines
  with_blur(
    geom_path(
      data = current_lines,
      aes(
        x = x,
        y = y,
        group = line_id,
        alpha = alpha_line,
        linewidth = line_width
      ),
      colour = "#e0fbfc",
      lineend = "round"
    ),
    sigma = 1.2
  ) +

  # Trails behind important variants
  geom_path(
    data = trails,
    aes(
      x = x,
      y = y,
      group = trail_id,
      colour = clinvar_group,
      alpha = alpha_trail
    ),
    linewidth = 0.45,
    lineend = "round"
  ) +

  # Glow layer for all particles
  with_outer_glow(
    geom_point(
      data = scn5a_art,
      aes(
        x = x,
        y = y,
        colour = clinvar_group,
        size = particle_size,
        alpha = alpha_value
      )
    ),
    colour = "white",
    sigma = 7,
    expand = 3
  ) +

  # Main particle layer
  geom_point(
    data = scn5a_art,
    aes(
      x = x,
      y = y,
      colour = clinvar_group,
      shape = consequence_group,
      size = particle_size,
      alpha = alpha_value
    ),
    stroke = 0.25
  ) +

  # Additional white-hot core for pathogenic variants
  with_outer_glow(
    geom_point(
      data = scn5a_art %>%
        filter(clinvar_group == "Pathogenic / likely pathogenic"),
      aes(x = x, y = y),
      colour = "white",
      size = 1.2,
      alpha = 0.95
    ),
    colour = "#ff4d6d",
    sigma = 9,
    expand = 6
  ) +

  # Title, integrated as part of artwork
  annotate(
    "text",
    x = -5.95,
    y = 4.05,
    label = "SCN5A",
    colour = "white",
    size = 15,
    fontface = "bold",
    hjust = 0
  ) +

  annotate(
    "text",
    x = -5.90,
    y = 3.58,
    label = "ClinVar electric storm",
    colour = "#bde0fe",
    size = 5,
    hjust = 0
  ) +

  annotate(
    "text",
    x = -5.90,
    y = 3.25,
    label = paste0(
      nrow(scn5a_art),
      " variants  |  chr3:",
      format(pos_min, big.mark = ","),
      "–",
      format(pos_max, big.mark = ",")
    ),
    colour = "grey85",
    size = 3.2,
    hjust = 0
  ) +

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

  scale_size_identity() +
  scale_alpha_identity() +
  scale_linewidth_identity() +

  coord_equal(
    xlim = c(-6.4, 6.4),
    ylim = c(-4.4, 4.4),
    expand = FALSE
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
    caption = "Data-derived generative artwork from ClinVar SCN5A variants"
  ) +

  theme_void(base_size = 14) +

  theme(
    plot.background = element_rect(fill = "#030014", colour = NA),
    panel.background = element_rect(fill = "#030014", colour = NA),

    legend.position = c(0.82, 0.20),
    legend.background = element_rect(fill = alpha("#030014", 0.45), colour = NA),
    legend.key = element_rect(fill = alpha("#030014", 0.1), colour = NA),
    legend.title = element_text(colour = "white", size = 9),
    legend.text = element_text(colour = "grey90", size = 7),

    plot.caption = element_text(
      colour = "grey65",
      hjust = 0.5,
      size = 8,
      margin = margin(t = 8)
    ),

    plot.margin = margin(0, 0, 0, 0)
  )


# ------------------------------------------------------------
# 10. Export
# ------------------------------------------------------------

print(p)

ggsave(
  filename = output_png,
  plot = p,
  width = 14,
  height = 10,
  dpi = 600,
  bg = "#030014"
)

ggsave(
  filename = output_pdf,
  plot = p,
  width = 14,
  height = 10,
  bg = "#030014"
)

message("Saved: ", output_png)
message("Saved: ", output_pdf)


# ------------------------------------------------------------
# 11. Save summary table
# ------------------------------------------------------------

summary_table <- scn5a_art %>%
  count(clinvar_group, consequence_group, variant_type, sort = TRUE)

readr::write_tsv(
  summary_table,
  "SCN5A_electric_storm_summary.tsv"
)

message("Saved: SCN5A_electric_storm_summary.tsv")
