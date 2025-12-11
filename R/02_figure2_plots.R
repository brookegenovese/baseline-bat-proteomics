# Mar-2025 ==== author: Brooke Genovese
##============================================================================##
# * GENERATING PLOTS AND FIGURES *
# * ERB specimens + batch covariate *
# * 01_msdap.R *
##============================================================================##

# --- load libraries --- #
library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl)
library(writexl)
library(ggplot2)
library(proteomicsCV)
library(cowplot)
library(viridis) 
library(gt)
library(boot)
library(webshot2)


##============================================================================##
## publication theme across plots
##============================================================================##

my_theme <- function(base_size = 12, base_family = "Helvetica") {
  theme_cowplot(font_size = base_size, font_family = base_family) +
    theme(
      axis.text  = element_text(size = base_size, family = base_family, face = "bold"),
      axis.title = element_text(
        size = base_size + 1,
        family = base_family,
        face = "bold"     
      ),
      plot.title = element_text(size = base_size + 2, face = "bold", family = base_family),
      legend.title = element_text(size = base_size, family = base_family),
      legend.text = element_text(size = base_size - 1, family = base_family),
      legend.position = "right",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

##============================================================================##
## sample correlation matrix
##============================================================================##

# --- import data --- #

indep <- read_tsv("./protein_abundance__filter by group independently.tsv")
indep <- indep %>%
  filter(!str_detect(protein_id, "Cont_")) %>%
  separate(protein_id, into = c("leading_accession", "other_ids"), sep = ";", extra = "drop") %>%
  dplyr::select(-other_ids) 

fl_columns <- grep("^FL", names(indep), value = TRUE)  
# column names that start with "FL" are samples
datInd <- data.frame(
  row.names = indep$leading_accession,
  indep[, fl_columns]
) 

colnames(datInd)[1] <- "C1"
colnames(datInd)[2] <- "C2"
colnames(datInd)[3] <- "C3"
colnames(datInd)[4] <- "C4"
colnames(datInd)[5] <- "C5"
colnames(datInd)[6] <- "C6"

correlation_matrix <- cor(datInd, use = "pairwise.complete.obs", method = "pearson")
correlation_matrix_df <- as.data.frame(as.table(correlation_matrix))

# --- plot --- #

P1 <- ggplot(correlation_matrix_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0.5, 1), option = "inferno") +
  coord_fixed() +
  my_theme() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 8))
  ) +
  labs(x = NULL, y = "Samples", fill = "Pearson corr. coeff.")
P1 

##============================================================================##
## peptide / protein counts (barplot)
##============================================================================##

# --- import data --- #
pdata <- read_xlsx("./samples.xlsx") # derived from MS-DAP output

# --- calc median peptides and proteins with CI --- #
# func to compute median for bootstrapping
median_boot <- function(data, indices, colname) {
  sample_data <- data[indices, colname, drop = TRUE]
  return(median(sample_data, na.rm = TRUE))
}

# func to get median, CI, and half-width of CI
get_median_ci <- function(data, colname, R = 1000) {
  boot_out <- boot(data = data, statistic = median_boot, R = R, colname = colname)
  ci <- boot.ci(boot_out, type = "perc")$percent[4:5]
  plus_minus <- (ci[2] - ci[1]) / 2
  
  list(
    median = median(data[[colname]], na.rm = TRUE),
    lower_CI = ci[1],
    upper_CI = ci[2],
    plus_minus = plus_minus
  )
}

peptides_stats <- get_median_ci(pdata, "detected_peptides")
proteins_stats <- get_median_ci(pdata, "detected_proteins")
peptides_stats
proteins_stats

# --- define max values ---
max_peptides <- max(pdata %>% pull(detected_peptides))
max_proteins <- 400  # explicitly set to 400 for scaling

# --- process pdata for plot ---
pdata_long <- pdata %>%
  dplyr::select(sample_id, shortname, group, detected_proteins, detected_peptides) %>%
  pivot_longer(cols = c("detected_proteins", "detected_peptides"), 
               names_to = "Type", 
               values_to = "Count") %>%
  mutate(Type = ifelse(Type == "detected_proteins", "Protein", "Peptide")) %>%
  mutate(Count_Scaled = Count / ifelse(Type == "Protein", max_proteins, max_peptides) * ifelse(Type == "Protein", -1, 1))

# --- set explicit ticks ---
ticks_peptides <- pretty(1:max_peptides, n = 5, min.n = 4)
ticks_proteins <- c(0, 100, 200, 300, 400)  

# --- compute mirrored coords ---
ticks_coord <- c(-1 * rev(ticks_proteins / max_proteins), ticks_peptides / max_peptides)
ticks_label <- c(rev(ticks_proteins), ticks_peptides)

# --- plot ---
P2 <- ggplot(pdata_long, aes(x = shortname, y = Count_Scaled, fill = Type)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.6) +
  geom_text(aes(label = Count, y = Count_Scaled / 2), colour = "black", size = 4.5) +
  scale_y_continuous(
    breaks = ticks_coord,
    labels = ticks_label,
    limits = c(-1, 1)  # explicit axis limits 
  ) +
  scale_fill_manual(
    values = c("Peptide" = "#d8b365", "Protein" = "#5ab4ac"),
    breaks = c("Protein", "Peptide")) +
  coord_flip() +
  labs(y = "Protein / Peptide Counts", x = "Samples", fill = "Type", title = NULL) +
  my_theme() +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 12),
    legend.position = "top",
    legend.justification = "center",
    legend.box.just = "center",
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.text = element_text(size = 12, face = "bold"), 
    legend.title = element_blank()
  )

P2
##============================================================================##
## rank protein abundance + uniprot keyword roundup
##============================================================================##
# --- import data --- #
c <- read.csv("") # .csv spectronaut output ranking

# --- keyword roundup --- #
process_keywords <- function(df, keyword_col, protein_col) {
  total_proteins <- df %>%
    dplyr::select({{ protein_col }}) %>%
    dplyr::distinct() %>%
    nrow()
  
  df %>%
    dplyr::select({{ keyword_col }}, {{ protein_col }}) %>%
    dplyr::filter(!is.na({{ keyword_col }})) %>%
    tidyr::separate_rows({{ keyword_col }}, sep = ";") %>%
    dplyr::group_by({{ keyword_col }}) %>%
    dplyr::summarise(
      total_occurrences = dplyr::n(),
      unique_protein_count = dplyr::n_distinct({{ protein_col }}),
      proteins = paste(unique({{ protein_col }}), collapse = ", "),
      percent_proteins_represented = round((unique_protein_count / total_proteins) * 100, 2)
    ) %>%
    dplyr::arrange(dplyr::desc(total_occurrences))
}

# bat keywords
bat_keywords_summary <- process_keywords(c, bat_keywords, bat_protein_description)
# human keywords
human_keywords_summary <- process_keywords(c, human_keywords, human_gene_name)

head(bat_keywords_summary)
head(human_keywords_summary)

# --- data wrangle for plot --- #
#function
get_proteins_by_keyword <- function(df, keyword_col, protein_col, keyword_of_interest, save_as_character = FALSE) {
  result <- df %>%
    dplyr::select({{ keyword_col }}, {{ protein_col }}) %>%
    dplyr::filter(!is.na({{ keyword_col }})) %>%
    tidyr::separate_rows({{ keyword_col }}, sep = ";") %>%
    dplyr::filter(trimws({{ keyword_col }}) == keyword_of_interest) %>%
    dplyr::distinct({{ protein_col }}) %>%
    dplyr::pull({{ protein_col }})
  
  if (save_as_character) {
    result <- paste(result, collapse = ", ")
  }
  
  return(result)
}

comp_path <- get_proteins_by_keyword(c, human_keywords, leading_accession, "Complement pathway", save_as_character = FALSE)
proteasome <- get_proteins_by_keyword(c, human_keywords, leading_accession, "Proteasome", save_as_character = FALSE)
oxidoreductase <- get_proteins_by_keyword(c, human_keywords, leading_accession, "Oxidoreductase", save_as_character = FALSE)

# to specifically gather all apolipoproteins
apo_accessions <- c %>%
  filter(grepl("^APO", human_gene_name)) %>%
  pull(leading_accession) %>%
  unique()

c <- c %>%
  mutate(highlight = case_when(
    leading_accession %in% apo_accessions ~ "apolipoprotein",  # prioritized first
    leading_accession %in% proteasome ~ "proteasome",
    leading_accession %in% comp_path ~ "complement pathway",
    leading_accession %in% oxidoreductase ~ "oxidoreductase",
    TRUE ~ "other"
  ))

apo_labels <- c %>% filter(highlight == "apolipoprotein")

# --- plot --- #

P3 <- ggplot(c, aes(x = rank, y = log10_medianquantity)) +
  geom_line(color = "darkgrey", linewidth = 0.50) +
  geom_point(aes(fill = "other"), shape = 21, size = 2, stroke = 0.5, color = "grey") +
  geom_point(data = subset(c, highlight != "other"), aes(fill = highlight),
             shape = 21, size = 3.0, stroke = 0.5, color = "black") +
  scale_fill_manual(
    values = c(
      "complement pathway" = "#F8766D",
      "proteasome" = "#ffd92f",
      "oxidoreductase" = "#7CAE00",
      "apolipoprotein" = "#C77CFF",
      "other" = "darkgrey"
    ),
    name = "UniProt Keywords",
    breaks = c("complement pathway", "proteasome", "oxidoreductase", "apolipoprotein")
  ) +
  labs(title = "", x = "Protein Rank", y = "Log10 Median Quantity") +
  my_theme() +
  theme(
    legend.position = c(.6, .75),
    legend.text = element_text(size = 14, face = "bold"),   
    legend.title = element_text(size = 15, face = "bold")   
  )

P3 

###

# filter and split by group
filtered_tbl <- c %>%
  filter(highlight %in% c("complement pathway", "proteasome", "oxidoreductase", "apolipoprotein")) %>%
  dplyr::select(
    Rank = rank,
    `Gene Symbol` = human_gene_name,
    `Log10 Median Quantity` = log10_medianquantity,
    Group = highlight
  ) %>%
  arrange(Group, Rank) %>%
  group_split(Group)

# create a list of gt tables 
group_tables <- lapply(filtered_tbl, function(df) {
  group_name <- unique(df$Group)
  
  df %>%
    dplyr::select(-Group) %>%
    gt() %>%
    tab_spanner(
      label = toupper(group_name),
      columns = everything()
    ) %>%
    cols_label(
      Rank = "Rank",
      `Gene Symbol` = "Gene Symbol",
      `Log10 Median Quantity` = "Log10 Median Quantity"
    ) %>%
    fmt_number(
      columns = `Log10 Median Quantity`,
      decimals = 2
    ) %>%
    tab_options(
      table.font.size = px(11),
      table.font.names = "Helvetica",
      heading.align = "left",
      column_labels.font.weight = "bold",
      table.border.top.width = px(0),
      table.border.bottom.width = px(0),
      heading.border.bottom.width = px(0),
      table_body.border.bottom.width = px(0),
      table_body.border.top.width = px(0),
      table_body.hlines.width = px(0.5),
      table_body.hlines.color = "#D3D3D3",
      data_row.padding = px(3)
    )
})

group_tables[[1]]
group_tables[[2]]
group_tables[[3]]
group_tables[[4]] 

if (!dir.exists("./figs/temp/tables")) {
  dir.create("./figs/temp/tables", recursive = TRUE)
}
for (i in seq_along(group_tables)) {
  file_name <- paste0("./figs/temp/tables/table_", i, ".png")
  gtsave(data = group_tables[[i]], filename = file_name)
}

##============================================================================##
## distribution of top 20 proteins
##============================================================================##

# extract top 20 proteins based on rank
top_20_accessions <- c %>%
  arrange(rank) %>%  
  slice_head(n = 20) %>%  
  pull(leading_accession)  

by_contr <- read_tsv("./protein_abundance__filter by contrast; female vs male # condition_variable_ group # additional_variables_ batch.tsv")
by_contr <- by_contr %>%
  mutate(leading_accession = sapply(strsplit(protein_id, ";"), `[`, 1))
top_20 <- by_contr %>% # subset for top 20
  filter(leading_accession %in% top_20_accessions)
top_20 <- top_20 %>%
  left_join(c %>% dplyr::select(leading_accession, human_gene_name), by = "leading_accession")
top_20_long <- top_20 %>%
  pivot_longer(cols = 4:9, names_to = "sample", values_to = "abundance") %>%
  mutate(log10_abundance = abundance / log2(10))  # convert

top_20_long <- top_20_long %>%
  left_join(c %>% dplyr::select(leading_accession, log10_medianquantity), by = "leading_accession")

median_log10 <- top_20_long %>% # compute median log10 abundance for each protein
  group_by(leading_accession) %>%
  summarize(median_log10_abundance = median(log10_abundance, na.rm = TRUE))

# compute scaling factors
scaling_factors <- median_log10 %>%
  left_join(c %>% dplyr::select(leading_accession, log10_medianquantity), by = "leading_accession") %>%
  mutate(scaling_factor = log10_medianquantity / median_log10_abundance)
top_20_long <- top_20_long %>%
  left_join(scaling_factors %>% dplyr::select(leading_accession, scaling_factor), by = "leading_accession") %>%
  mutate(scaled_log10_abundance = log10_abundance * scaling_factor)

# --- plot --- #

P4 <- ggplot(top_20_long, aes(x = reorder(human_gene_name, -scaled_log10_abundance, median), 
                              y = scaled_log10_abundance)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = NA, alpha = 0.7) + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 3.0, color = "#5ab4ac") +
  scale_y_continuous(breaks = seq(9, 5, by = -0.5), limits = c(6, 9)) +
  labs(title = "", x = "Protein (Gene Name)", y = "Abundance (Log10)") +
  my_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
P4

# --- CoV of top 20 --- #

cv <- top_20 %>% 
  dplyr::select(leading_accession, dplyr::contains("FL030524"))
cv <- as.data.frame(cv)

rownames(cv) <- cv$leading_accession
cv <- cv[, -1]  

# calculate CVs
cvs <- protLogCV(cv, 2)
print(cvs)
cv_df <- data.frame(Protein = names(cvs), CV = cvs)

# "bin" CV values 
cv_df <- cv_df %>%
  mutate(CV_category = factor(case_when(
    CV < 20 ~ "< 20%",
    CV >= 20 & CV < 30 ~ "20-30%",
    CV >= 30 ~ "≥ 30%"
  ), levels = c("< 20%", "20-30%", "≥ 30%")))  

cv_summary <- cv_df %>%
  group_by(CV_category) %>%
  summarise(Protein_Count = n(), .groups = "drop")

# --- plot --- #
Pcv <- ggplot(cv_summary, aes(x = Protein_Count, y = CV_category)) +
  geom_bar(stat = "identity", fill = "grey60", width = 0.6) +
  labs(x = "# Proteins", y = "CV %") +
  my_theme() +
  theme(axis.title.y = element_text(margin = margin(r = 8)),
        axis.text.y = element_text(size = 10),
        axis.text = element_text(face = "plain"))

P4_inset <- ggdraw() +
  draw_plot(P4) +  # main plot
  draw_plot(Pcv, x = 0.65, y = 0.65, width = 0.3, height = 0.2)

P4_inset

#####

pan_AB <- plot_grid(NULL, P2, NULL, P1,
                    rel_heights = c(0.05, 1, 0.05, 1),  
                    ncol = 1,
                    label_y = 1.07,
                    labels = c("","A","","B"),
                    label_size = 18,
                    label_fontfamily = "sans")
pan_AB

pan_C <- plot_grid(NULL, P3,
                   rel_heights = c(0.1, 1.85),
                   ncol = 1,
                   label_y = 1 + 0.07 /1.85,
                   labels = c("","C"),
                   label_size = 18,
                   label_fontfamily = "sans")
pan_C

pan_ABC <- plot_grid(
  pan_AB,
  NULL,             
  pan_C,
  nrow = 1,
  rel_widths = c(2, 0.15, 3)   
)
pan_ABC

fig_P4inset <- P4_inset + theme(plot.margin = unit(c(5.5, 5.5, 0, 5.5), units = "pt"))
fig_P4inset
pan_D <- plot_grid(fig_P4inset,
                   rel_widths = c(0.09, 0.93),
                   nrow = 1,
                   labels = c("D",""),
                   label_size = 18,
                   label_fontfamily = "sans")
pan_D

pan_full <- plot_grid(pan_ABC, pan_D,
                      ncol = 1,
                      rel_heights = c(2, 2),  
                      label_size = 18,
                      label_fontfamily = "sans")
pan_full 

output_path <- "" # define

ggsave(filename = output_path,
       plot = pan_full,
       width = 11,
       height = 10,
       units = "in",
       dpi = 300,
       bg = "white")  



