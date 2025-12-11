# Mar-2025 ==== author: Brooke Genovese
##============================================================================##
# * EXPLORATORY SEX DIFFERENCES ANALYSIS: rCV CALCULATIONS *
# * ERB specimens *
# * 01_msdap.R *
##============================================================================##

# --- load libraries --- #
library(dplyr)
library(tidyverse)
library(readxl)
library(ggplot2)
library(cowplot)
library(scales)  
library(ggrepel)

##============================================================================##
## publication theme across plots
##============================================================================##

my_theme <- function(base_size = 10, base_family = "Arial") {
  theme_cowplot(font_size = base_size, font_family = base_family) +
    theme(
      axis.text = element_text(size = base_size, family = base_family),
      axis.title = element_text(size = base_size + 1, family = base_family),
      plot.title = element_text(size = base_size + 2, face = "bold", family = base_family),
      legend.title = element_text(size = base_size, family = base_family),
      legend.text = element_text(size = base_size - 1, family = base_family),
      legend.position = "right",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
} 

##============================================================================##
## 
##============================================================================##

# --- import data --- #

dat0 <- read_tsv("./protein_abundance__filter by group independently.tsv")
dat0 <- as.data.frame(dat0)
dat0$leading_accession <- sapply(strsplit(dat0$protein_id, ";"), `[`, 1)
dat0 <- dat0[, c("leading_accession", setdiff(names(dat0), "leading_accession"))]

anno <- read_csv("") # read in annotation file

sex_samples <- tibble(
  sample = c("FL030524_BGbatPlasm-Dia_90m_c1", 
             "FL030524_BGbatPlasm-Dia_90m_c2", 
             "FL030524_BGbatPlasm-Dia_90m_c3",
             "FL030524_BGbatPlasm-Dia_90m_c4", 
             "FL030524_BGbatPlasm-Dia_90m_c5", 
             "FL030524_BGbatPlasm-Dia_90m_c6"),
  sex = c("Female", "Female", "Female", "Male", "Male", "Male")  
)

##============================================================================##
## calculating MAD
##============================================================================## 

intensity_cols <- grep("^FL", names(dat0), value = TRUE)
dat_long <- dat0 %>%
  select(leading_accession, all_of(intensity_cols)) %>%
  pivot_longer(
    cols = all_of(intensity_cols),
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  left_join(sex_samples, by = "sample")

# compute MAD within each sex for each protein
mad_by_sex <- dat_long %>%
  group_by(leading_accession, sex) %>%
  summarise(mad = mad(abundance, na.rm = TRUE), .groups = "drop")

# optional formats
mad_wide <- mad_by_sex %>%
  pivot_wider(names_from = sex, values_from = mad, names_prefix = "MAD_")

# average across sexes for overall within-sex variability
mad_summary <- mad_wide %>%
  mutate(MAD_mean = rowMeans(across(starts_with("MAD_")), na.rm = TRUE))

# what are the top variable proteins?
mad_summary %>%
  arrange(desc(MAD_mean)) %>%
  head(20)

ggplot(mad_wide, aes(x = MAD_Female, y = MAD_Male)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "MAD (Female)",
    y = "MAD (Male)",
    title = "Protein variability within sexes (MAD)"
  ) +
  theme_minimal() 


# now compute MAD + median abundance for each protein/sex
mad_median_by_sex <- dat_long %>%
  group_by(leading_accession, sex) %>%
  summarise(
    MAD = mad(abundance, na.rm = TRUE),
    Median = median(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(rCV = MAD / Median)

# per sex
var_wide <- mad_median_by_sex %>%
  select(leading_accession, sex, MAD, Median, rCV) %>%
  pivot_wider(
    names_from = sex,
    values_from = c(MAD, Median, rCV),
    names_sep = "_"
  )

# means across sexes
var_summary <- var_wide %>%
  mutate(
    MAD_mean = rowMeans(across(starts_with("MAD_")), na.rm = TRUE),
    rCV_mean = rowMeans(across(starts_with("rCV_")), na.rm = TRUE)
  )


### 

# tag proteins by detection pattern
tag_sex_detection <- function(df,
                              median_female_col = "Median_Female",
                              median_male_col   = "Median_Male",
                              tag_col_name = "sex_detection_tag") {
  if (!all(c(median_female_col, median_male_col) %in% names(df))) {
    stop("Median column names not found in dataframe.")
  }
  
  df %>%
    mutate(
      detected_female = !is.na(.data[[median_female_col]]),
      detected_male   = !is.na(.data[[median_male_col]]),
      !!tag_col_name := case_when(
        detected_female & detected_male ~ "both_sexes",
        detected_female & !detected_male ~ "female_only",
        !detected_female & detected_male ~ "male_only",
        TRUE ~ "undetected"
      )
    ) %>%
    select(-detected_female, -detected_male)
}

# apply the function to var_summary table
# (var_summary should have columns Median_Female and Median_Male)
var_summary_tagged <- tag_sex_detection(var_summary)

# quick check
var_summary_tagged %>%
  count(sex_detection_tag) %>%
  arrange(desc(n)) -> detection_counts
print(detection_counts)

# example: get top 20 proteins by rCV_mean that are detected in both sexes
top20_both_by_rCV <- var_summary_tagged %>%
  filter(sex_detection_tag == "both_sexes") %>%
  arrange(desc(rCV_mean)) %>%
  slice_head(n = 20)

# example: top 20 proteins by MAD_mean that are male-only 
top20_male_only_by_MAD <- var_summary_tagged %>%
  filter(sex_detection_tag == "male_only") %>%
  arrange(desc(MAD_mean)) %>%
  slice_head(n = 20)


# label extreme outliers 
label_n <- 10
top_labels <- var_summary_tagged %>%
  filter(sex_detection_tag == "both_sexes") %>%
  arrange(desc(rCV_mean)) %>%
  slice_head(n = label_n)

# restrict to proteins quantified in both sexes
both_sex_only <- var_summary_tagged %>%
  filter(sex_detection_tag == "both_sexes")

# --- Most variable proteins ---
top10_female_var <- both_sex_only_anno %>%
  arrange(desc(rCV_Female)) %>%
  slice_head(n = 10) %>%
  mutate(tag = "Most variable (Female)")

top10_male_var <- both_sex_only_anno %>%
  arrange(desc(rCV_Male)) %>%
  slice_head(n = 10) %>%
  mutate(tag = "Most variable (Male)")

top_labels_var <- bind_rows(top10_female_var, top10_male_var)

p_var <- ggplot(both_sex_only_anno, aes(x = rCV_Female, y = rCV_Male)) +
  geom_point(aes(color = (rCV_Male - rCV_Female)), alpha = 0.8, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_color_gradient2(
    low = "#F8766D", mid = "grey80", high = "#00afbb",
    midpoint = 0, name = "rCV(M - F)"
  ) +
  geom_text_repel(
    data = top_labels_var,
    aes(label = human_gene_name, color = NULL),
    size = 3, max.overlaps = 25
  ) +
  labs(
    x = "rCV (Female)",
    y = "rCV (Male)",
    title = "Coefficient of Variation comparison between sexes",
    subtitle = "Top 10 most variable proteins per sex (both quantified)"
  ) +
  theme_cowplot()

p_var 

#output_file2 <- file.path("./results/ctrl", "05_top_variable_sex.xlsx")
#write_xlsx(list("top10vari_female" = top10_female_var, "top10vari_male" = top10_male_var, "top10stable_female" = top10_female_stable, "top10stable_male" = top10_male_stable), path = output_file2)
# --- Most stable proteins ---
top10_female_stable <- both_sex_only_anno %>%
  arrange(rCV_Female) %>%      # ascending
  slice_head(n = 10) %>%
  mutate(tag = "Most stable (Female)")

top10_male_stable <- both_sex_only_anno %>%
  arrange(rCV_Male) %>%
  slice_head(n = 10) %>%
  mutate(tag = "Most stable (Male)")

top_labels_stable <- bind_rows(top10_female_stable, top10_male_stable)

p_stable <- ggplot(both_sex_only_anno, aes(x = rCV_Female, y = rCV_Male)) +
  geom_point(aes(color = (rCV_Male - rCV_Female)), alpha = 0.8, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_color_gradient2(
    low = "#F8766D", mid = "grey80", high = "#00afbb",
    midpoint = 0, name = "rCV(M - F)"
  ) +
  geom_text_repel(
    data = top_labels_stable,
    aes(label = human_gene_name, color = NULL),
    size = 3, max.overlaps = 25
  ) +
  labs(
    x = "rCV (Female)",
    y = "rCV (Male)",
    title = "Coefficient of Variation comparison between sexes",
    subtitle = "Top 10 most stable proteins per sex (both quantified)"
  ) +
  theme_cowplot()


p_stable_zoom <- p_stable +
  coord_cartesian(xlim = c(0, 0.05), ylim = c(0, 0.05)) +
  labs(subtitle = "Top 10 most stable proteins per sex (zoomed-in to very low rCV)")

p_stable_zoom 
summary(both_sex_only$rCV_Female)



