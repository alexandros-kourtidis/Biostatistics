# ------------------
# Environment setup
# ------------------
library(tidyverse)
library(janitor)
library(lubridate)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(broom.mixed)
library(stringr)
library(dplyr)
library(ggplot2)
library(forcats)
library(ggsignif)

rm(list = ls())

# ------------------------------------
# 1. Read files
# ------------------------------------

long_b123 <- readRDS("data_processed/batch123_long.rds")
long_b123MK <- readRDS("data_processed/batch123MK_long.rds")

# treatment and ponds : fixed effect variables
# clonal line, date, observer : random effect variables


# ------------------------------------
# 2. Model fitting
# ------------------------------------
# 1. Alex + MK: WITH observer as random
long_b123MK_olre <- long_b123MK %>% mutate(obs_id = factor(row_number()))
m_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment * pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b123MK_olre
)

# 2. Alex: WITHOUT observer as random
long_b123_olre <- long_b123 %>% mutate(obs_id = factor(row_number()))
m_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment * pond_code +
    (1 | clonal_line) + (1 | date) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b123_olre
)

# ------------------------------------
# 3. Estimated Marginal Means per pond and population comparisons
# ------------------------------------
# 1. Alex + MK
emm_b123MK <- emmeans(m_b123MK, ~ pond_code | treatment, type = "response")
summary(emm_b123MK)
pairs_C5_b123MK <- contrast(emm_b123MK, 
                            method = "trt.vs.ctrl", 
                            ref = "C5", 
                            by = "treatment", 
                            adjust = "none")   # or "tukey"
summary(pairs_C5_b123MK)

# 2. Alex
emm_b123 <- emmeans(m_b123, ~ pond_code | treatment , type = "response")
summary(emm_b123)
pairs_C5_b123 <- contrast(emm_b123, 
                            method = "trt.vs.ctrl", 
                            ref = "C5", 
                            adjust = "none")   # or "tukey"
summary(pairs_C5_b123)


# ------------------------------------
# 4. Barplots of pond probabilities (colored by type) 
# ------------------------------------
pond_map <- long_b123 %>% distinct(pond_code, pond_type)

# Function to convert emmeans output to a tidy data frame
emm_to_df <- function(emm_obj) {
  df <- as_tibble(emm_obj)
  
  # Identify the response (estimate) column
  resp_candidates <- c("response", ".response", "prob", ".emmean", "rate")
  resp_col <- resp_candidates[resp_candidates %in% names(df)][1]
  if (is.na(resp_col)) {
    stop("Could not find the response column in emmeans output. ",
         "Columns were: ", paste(names(df), collapse = ", "))
  }
  
  # Identify lower/upper confidence limit cols (several possible namings)
  lcl_candidates <- c("asymp.LCL", "lower.CL", "LCL")
  ucl_candidates <- c("asymp.UCL", "upper.CL", "UCL")
  lcl_col <- lcl_candidates[lcl_candidates %in% names(df)][1]
  ucl_col <- ucl_candidates[ucl_candidates %in% names(df)][1]
  
  # Build a standardized tibble
  out <- df %>%
    mutate(prob = .data[[resp_col]])
  
  if (!is.na(lcl_col) && !is.na(ucl_col)) {
    out <- out %>%
      mutate(
        lower = .data[[lcl_col]],
        upper = .data[[ucl_col]]
      )
  } 
  out
}
# Function to convert pairwise comparisons to a tidy data frame
pairs_to_df <- function(pairs_obj) {
  df <- as_tibble(pairs_obj)
  
  df %>%
    mutate(
      pond_code = str_replace(contrast, " / C5", ""),
      p.value = p.value,
      sig = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE ~ ""
      )
    ) %>%
    select(treatment, pond_code, p.value, sig)
}
# Function to make a barplot with all three levels
plot_emm_bar <- function(df, pond_map, title = "") {
  df %>%
    left_join(pond_map, by = "pond_code") %>%
    mutate(
      pond_code = fct_reorder(pond_code, prob),
      treatment = factor(treatment, levels = c("control","0.5","0.7"))
    ) %>%
    ggplot(aes(x = pond_code, y = prob, fill = pond_type)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.8) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.9),
                  width = 0.2) +
    
    # ADD SIGNIFICANCE ASTERISKS
    geom_text(
      aes(
        y = upper + 0.05,   # shifts the asterisk above CI bar
        label = sig
      ),
      size = 6,
      vjust = 0,
      fontface = "bold"
    ) +
    
    facet_wrap(~ treatment) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = title,
      x = "Pond code",
      y = "Predicted survival probability"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Convert emm to dataframes
df_b123MK <- emm_to_df(emm_b123MK)
df_b123   <- emm_to_df(emm_b123)
# Convert pw comparisons to df
sig_b123MK <- pairs_to_df(pairs_C5_b123MK)
sig_b123   <- pairs_to_df(pairs_C5_b123)
# Merge significance annotations into EMM dataframe
df_b123MK_plot <- df_b123MK %>%
  left_join(sig_b123MK, by = c("pond_code", "treatment"))
df_b123_plot <- df_b123 %>%
  left_join(sig_b123, by = c("pond_code", "treatment"))
# Plot 
p_b123MK <- plot_emm_bar(df_b123MK_plot, pond_map, "Batch123MK — Survival Probabilities - Significance vs C5")
p_b123 <- plot_emm_bar(df_b123_plot, pond_map, "Batch123 — Survival Probabilities - Significance vs C5")
# Save as png
ggsave("plots/p_b123MK.png", p_b123MK, width = 10, height = 7, dpi = 300)
ggsave("plots/p_b123.png",   p_b123,   width = 10, height = 7, dpi = 300)
# Save as pdf
ggsave("plots/p_b123MK.pdf", p_b123MK, width = 10, height = 7)
ggsave("plots/p_b123.pdf",   p_b123,   width = 10, height = 7)
