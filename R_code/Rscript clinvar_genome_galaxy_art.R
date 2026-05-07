############################################################
# A Few Known Stars in the Human Genome
# Generative scientific artwork using ClinVar P/LP variants
#
# Input:
#   clinvar_P_LP.tsv
#
# Expected columns:
#   chr, pos, ref, alt, gene, consequence, classification,
#   conflicting_classification_detail
#
# Output:
#   ClinVar_P_LP_genome_galaxy_GRCh38.png
#   ClinVar_P_LP_genome_galaxy_GRCh38.pdf
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

############################################################
# 1. Input and output
############################################################

input_file <- "clinvar_P_LP.tsv"
output_prefix <- "ClinVar_P_LP_genome_galaxy_GRCh38"

set.seed(20011012)

############################################################
# 2. GRCh38 chromosome lengths
############################################################

chr_lengths <- tibble(
  chr = as.character(c(1:22, "X", "Y")),
  chr_len = c(
    248956422, 242193529, 198295559, 190214555, 181538259,
    170805979, 159345973, 145138636, 138394717, 133797422,
    135086622, 133275309, 114364328, 107043718, 101991189,
    90338345, 83257441, 80373285, 58617616, 64444167,
    46709983, 50818468, 156040895, 57227415
  )
) %>%
  mutate(
    chr = factor(chr, levels = as.character(c(1:22, "X", "Y"))),
    chr_index = row_number(),
    chr_start = lag(cumsum(chr_len), default = 0),
    chr_end = chr_start + chr_len,
    chr_mid = chr_start + chr_len / 2,
    genome_fraction_start = chr_start / max(chr_end),
    genome_fraction_end = chr_end / max(chr_end),
    genome_fraction_mid = chr_mid / max(chr_end)
  )

genome_size <- max(chr_lengths$chr_end)

############################################################
# 3. Read ClinVar P/LP TSV
############################################################

clinvar <- read_tsv(
  input_file,
  show_col_types = FALSE,
  col_types = cols(.default = col_character())
)

# Clean column names
colnames(clinvar) <- colnames(clinvar) %>%
  str_replace_all("\\s+", "_") %>%
  str_replace_all("\\.", "_")

required_cols <- c(
  "chr", "pos", "ref", "alt", "gene", "consequence",
  "classification", "conflicting_classification_detail"
)

missing_cols <- setdiff(required_cols, colnames(clinvar))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

############################################################
# 4. Clean and classify variants
############################################################

clinvar_clean <- clinvar %>%
  mutate(
    chr = as.character(chr),
    chr = str_replace(chr, "^chr", ""),
    chr = case_when(
      chr %in% as.character(1:22) ~ chr,
      chr %in% c("23", "X") ~ "X",
      chr %in% c("24", "Y") ~ "Y",
      TRUE ~ chr
    ),
    pos = suppressWarnings(as.numeric(pos)),
    ref = replace_na(ref, ""),
    alt = replace_na(alt, ""),
    gene = replace_na(gene, ""),
    consequence = replace_na(consequence, ""),
    classification = replace_na(classification, ""),
    conflicting_classification_detail = replace_na(
      conflicting_classification_detail,
      ""
    )
  ) %>%
  filter(
    chr %in% as.character(c(1:22, "X", "Y")),
    !is.na(pos),
    chr != "chr"
  )

clinvar_plp <- clinvar_clean %>%
  filter(
    str_detect(classification, "Pathogenic|Likely_pathogenic") |
      str_detect(
        conflicting_classification_detail,
        "Pathogenic|Likely_pathogenic"
      )
  ) %>%
  mutate(
    classification_group = case_when(
      str_detect(classification, "Conflicting_classifications_of_pathogenicity") &
        str_detect(
          conflicting_classification_detail,
          "Pathogenic|Likely_pathogenic"
        ) ~ "Conflicting with P/LP evidence",

      str_detect(classification, "Pathogenic/Likely_pathogenic") ~
        "Pathogenic / Likely pathogenic",

      str_detect(classification, "Likely_pathogenic") &
        !str_detect(classification, "Pathogenic/Likely_pathogenic") ~
        "Likely pathogenic",

      str_detect(classification, "Pathogenic") ~
        "Pathogenic",

      TRUE ~ "Other P/LP-related"
    )
  )

############################################################
# 5. Map variants to traceable genome coordinates
############################################################

variants <- clinvar_plp %>%
  left_join(chr_lengths, by = "chr") %>%
  filter(
    !is.na(chr_len),
    pos <= chr_len
  ) %>%
  mutate(
    genome_pos = chr_start + pos,
    genome_fraction = genome_pos / genome_size,
    chr_fraction = pos / chr_len
  )

message("Total P/LP-related variants plotted: ", nrow(variants))

############################################################
# 6. Convert genomic coordinates into generative galaxy space
############################################################
# Scientific mapping:
#   chromosome determines the galaxy arm
#   position within chromosome determines distance along that arm
#
# Artistic transformation:
#   each chromosome becomes a spiral-like arm segment
#   controlled jitter creates a sky-like, non-track appearance

variants_art <- variants %>%
  mutate(
    chr_num = as.numeric(chr),

    # Each chromosome gets a slightly different angular sector
    arm_angle = 2 * pi * (chr_num - 1) / 24,

    # Genomic position along chromosome moves along the arm
    radial_base = 0.25 + 4.8 * chr_fraction,

    # Spiral twist within each chromosome arm
    theta = arm_angle +
      1.85 * chr_fraction +
      0.20 * sin(chr_fraction * 2 * pi * 3),

    # Controlled randomness for generative sky effect
    radial_jitter = rnorm(n(), mean = 0, sd = 0.055),
    theta_jitter = rnorm(n(), mean = 0, sd = 0.045),

    r = radial_base + radial_jitter,
    theta2 = theta + theta_jitter,

    x = r * cos(theta2),
    y = r * sin(theta2),

    # Star appearance
    star_size = case_when(
      classification_group == "Pathogenic" ~ runif(n(), 0.55, 0.40),
      classification_group == "Pathogenic / Likely pathogenic" ~ runif(n(), 0.45, 0.20),
      classification_group == "Likely pathogenic" ~ runif(n(), 0.35, 0.095),
      classification_group == "Conflicting with P/LP evidence" ~ runif(n(), 0.25, 0.070),
      TRUE ~ runif(n(), 0.20, 0.055)
    ),

    star_alpha = case_when(
      classification_group == "Pathogenic" ~ runif(n(), 0.55, 0.095),
      classification_group == "Pathogenic / Likely pathogenic" ~ runif(n(), 0.50, 0.085),
      classification_group == "Likely pathogenic" ~ runif(n(), 0.40, 0.075),
      classification_group == "Conflicting with P/LP evidence" ~ runif(n(), 0.25, 0.060),
      TRUE ~ runif(n(), 0.20, 0.045)
    )
  )

############################################################
# 7. Create subtle chromosome galaxy-arm guide curves
############################################################

arm_guides <- map_dfr(seq_len(nrow(chr_lengths)), function(i) {
  chr_i <- chr_lengths$chr[i]
  chr_num_i <- i
  arm_angle_i <- 2 * pi * (chr_num_i - 1) / 24

  tibble(
    chr = as.character(chr_i),
    chr_fraction = seq(0, 1, length.out = 250)
  ) %>%
    mutate(
      radial_base = 0.25 + 4.8 * chr_fraction,
      theta = arm_angle_i +
        1.85 * chr_fraction +
        0.20 * sin(chr_fraction * 2 * pi * 3),
      x = radial_base * cos(theta),
      y = radial_base * sin(theta),
      chr_index = chr_num_i
    )
})

############################################################
# 8. Background distant stars and genome dust
############################################################

# Distant unlabelled stars represent the vast uninterpreted genome.
# These are decorative only, not ClinVar variants.

n_dust <- 9000

background_dust <- tibble(
  theta = runif(n_dust, 0, 2 * pi),
  r = sqrt(runif(n_dust, 0, 1)) * 5.7,
  x = r * cos(theta),
  y = r * sin(theta),
  size = runif(n_dust, 0.03, 0.18),
  alpha = runif(n_dust, 0.015, 0.12)
)

# A faint central glow
n_nebula <- 3500

nebula <- tibble(
  theta = runif(n_nebula, 0, 2 * pi),
  r = rgamma(n_nebula, shape = 2.2, rate = 1.2),
  x = r * cos(theta) + rnorm(n_nebula, 0, 0.15),
  y = r * sin(theta) + rnorm(n_nebula, 0, 0.15),
  size = runif(n_nebula, 0.05, 0.40),
  alpha = runif(n_nebula, 0.015, 0.08)
)

############################################################
# 9. Chromosome labels
############################################################

chr_labels <- arm_guides %>%
  group_by(chr) %>%
  slice_max(chr_fraction, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    label = paste0("chr", chr),
    label_x = x * 1.08,
    label_y = y * 1.08
  )

############################################################
# 10. Optional gene labels for highly represented genes
############################################################

top_genes <- variants_art %>%
  filter(gene != "", gene != ".") %>%
  count(gene, sort = TRUE) %>%
  slice_head(n = 20)

gene_labels <- variants_art %>%
  semi_join(top_genes, by = "gene") %>%
  group_by(gene) %>%
  summarise(
    x = median(x, na.rm = TRUE),
    y = median(y, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(n)) %>%
  mutate(
    label = paste0(gene)
  )

############################################################
# 11. Plot
############################################################

p <- ggplot() +

  # Deep-space background dust
  geom_point(
    data = background_dust,
    aes(x = x, y = y, size = size, alpha = alpha),
    colour = "white"
  ) +

  # Nebula-like central glow
  geom_point(
    data = nebula,
    aes(x = x, y = y, size = size, alpha = alpha),
    colour = "#6baed6"
  ) +

  # Subtle chromosome galaxy arms
  geom_path(
    data = arm_guides,
    aes(x = x, y = y, group = chr),
    colour = "#506080",
    linewidth = 0.22,
    alpha = 0.25
  ) +

  # Star glow layer: larger and transparent
  geom_point(
    data = variants_art,
    aes(
      x = x,
      y = y,
      colour = classification_group,
      size = star_size
    ),
    alpha = 0.13,
    shape = 16
  ) +

  # Main star layer
  geom_point(
    data = variants_art,
    aes(
      x = x,
      y = y,
      colour = classification_group,
      size = star_size,
      alpha = star_alpha
    ),
    shape = 8,
    stroke = 0.22
  ) +

  # Tiny bright cores
  geom_point(
    data = variants_art %>% slice_sample(n = min(8000, nrow(variants_art))),
    aes(
      x = x,
      y = y,
      colour = classification_group
    ),
    size = 0.12,
    alpha = 0.85
  ) +

  # Chromosome labels
  geom_text(
    data = chr_labels,
    aes(x = label_x, y = label_y, label = label),
    colour = "grey72",
    size = 2.7,
    alpha = 0.75
  ) +

  # Gene labels, deliberately subtle
  geom_text(
    data = gene_labels,
    aes(x = x, y = y, label = label),
    colour = "grey88",
    size = 2.4,
    alpha = 0.50
  ) +

  coord_equal(xlim = c(-6.3, 6.3), ylim = c(-6.3, 6.3), expand = FALSE) +

  scale_colour_manual(
    values = c(
      "Pathogenic" = "#fff7bc",
      "Likely pathogenic" = "#fec44f",
      "Pathogenic / Likely pathogenic" = "#f03b20",
      "Conflicting with P/LP evidence" = "#9ecae1",
      "Other P/LP-related" = "#c994c7"
    )
  ) +

  scale_size_identity() +
  scale_alpha_identity() +

  labs(
    title = "A Few Known Stars in the Human Genome",
    subtitle = paste0(
      comma(nrow(variants_art)),
      " ClinVar variants with Pathogenic or Likely pathogenic evidence mapped across GRCh38"
    ),
    colour = "ClinVar evidence",
    caption = paste0(
      "Each star represents a ClinVar variant with Pathogenic or Likely pathogenic evidence. ",
      "Chromosome and genomic position determine each star's traceable location before artistic transformation. ",
      "The surrounding darkness represents the vast genome still limited by ascertainment, interpretation, and discovery."
    )
  ) +

  theme_void(base_size = 13) +
  theme(
    plot.background = element_rect(fill = "#02030a", colour = NA),
    panel.background = element_rect(fill = "#02030a", colour = NA),

    plot.title = element_text(
      colour = "white",
      size = 25,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 8)
    ),

    plot.subtitle = element_text(
      colour = "grey78",
      size = 11,
      hjust = 0.5,
      margin = margin(b = 18)
    ),

    plot.caption = element_text(
      colour = "grey55",
      size = 8.2,
      hjust = 0.5,
      margin = margin(t = 18, b = 4)
    ),

    legend.position = "bottom",
    legend.background = element_rect(fill = "#02030a", colour = NA),
    legend.key = element_rect(fill = "#02030a", colour = NA),
    legend.title = element_text(colour = "grey90", size = 10),
    legend.text = element_text(colour = "grey75", size = 9)
  )

############################################################
# 12. Save high-resolution files
############################################################

ggsave(
  filename = paste0(output_prefix, ".png"),
  plot = p,
  width = 14,
  height = 14,
  dpi = 600,
  bg = "#02030a"
)

ggsave(
  filename = paste0(output_prefix, ".pdf"),
  plot = p,
  width = 14,
  height = 14,
  bg = "#02030a"
)

print(p)

message("Saved PNG: ", paste0(output_prefix, ".png"))
message("Saved PDF: ", paste0(output_prefix, ".pdf"))

############################################################
# 13. Basic summary table for traceability
############################################################

summary_by_chr <- variants_art %>%
  count(chr, classification_group, name = "n") %>%
  arrange(chr, classification_group)

write_tsv(
  summary_by_chr,
  paste0(output_prefix, "_summary_by_chr.tsv")
)

message("Saved summary: ", paste0(output_prefix, "_summary_by_chr.tsv"))
