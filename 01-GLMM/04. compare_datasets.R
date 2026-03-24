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
library(patchwork)

rm(list = ls())

# ------------------
# 1. Read files
# ------------------
long_b123MK <- read_csv("data_processed/batch123MK_long.csv")

long_b123MK <- long_b123MK %>%
  mutate(
    treatment  = factor(treatment, levels = c("control","0.5","0.7")),
    pond_code  = factor(pond_code, levels = c("C5", setdiff(unique(pond_code), "C5"))),
    pond_type  = factor(pond_type, levels = c("natu","agri","city")),
    observer   = factor(observer),
    clonal_line = factor(clonal_line)
  )
summary(long_b123MK)
saveRDS(long_b123MK, "batch123MK_long.rds")

long_b123 <- readRDS("data_processed/batch123_long.rds")
long_b123MK <- readRDS("data_processed/batch123MK_long.rds")
long_MK <- readRDS("data_processed/MK_long.rds")
long_b3 <- readRDS("data_processed/batch3_long.rds")
long_b12 <- readRDS("data_processed/batch12_long.rds")


# ------------------------------------
# 2. Model fitting
# ------------------------------------
# 1. batch123 + MK: WITH observer as random
long_b123MK_olre <- long_b123MK %>% mutate(obs_id = factor(row_number()))
m_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment * pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b123MK_olre
)
# 2. batch123: WITHOUT observer as random
long_b123_olre <- long_b123 %>% mutate(obs_id = factor(row_number()))
m_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment * pond_code +
    (1 | clonal_line) + (1 | date) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b123_olre
)
# 3. batch12: WITHOUT observer as random
long_b12_olre <- long_b12 %>% mutate(obs_id = factor(row_number()))
m_b12 <- glmmTMB(
  cbind(alive, dead) ~ treatment * pond_code +
    (1 | clonal_line) + (1 | date) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b12_olre
)
# 4. batch3: WITHOUT observer as random
long_b3_olre <- long_b3 %>% mutate(obs_id = factor(row_number()))
m_b3 <- glmmTMB(
  cbind(alive, dead) ~ treatment * pond_code +
    (1 | clonal_line) + (1 | date) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b3_olre
)
# 5. MK: WITHOUT POND x TREATMENT interaction
long_MK_olre <- long_MK %>% mutate(obs_id = factor(row_number()))
m_MK <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_MK_olre
)

# ------------------------------------
# 3. Estimated Marginal Means per pond and population comparisons
# ------------------------------------
# 1. batch123 + MK
emm_b123MK <- emmeans(m_b123MK, ~ pond_code | treatment, type = "response")
# 2. batch123
emm_b123 <- emmeans(m_b123, ~ pond_code | treatment , type = "response")
# 3. batch12
emm_b12 <- emmeans(m_b12, ~ pond_code | treatment, type = "response")
# 4. batch3
emm_b3 <- emmeans(m_b3, ~ pond_code | treatment , type = "response")
# 5. MK
emm_MK <- emmeans(m_MK, ~ pond_code | treatment , type = "response")

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
# ----
# Convert emm to dataframes
# ----
df_b123MK <- emm_to_df(emm_b123MK)
df_b123   <- emm_to_df(emm_b123)
df_b12 <- emm_to_df(emm_b12)
df_b3   <- emm_to_df(emm_b3)
df_MK <- emm_to_df(emm_MK)
# ----
# Plot 
# ----
p_b123MK <- plot_emm_bar(df_b123MK, pond_map, "Batch 123 + MK — Survival Probabilities ")
p_b123 <- plot_emm_bar(df_b123, pond_map, "Batch 123 — Survival Probabilities ")
p_b12 <- plot_emm_bar(df_b12, pond_map, "Batch 12 — Survival Probabilities ")
p_b3 <- plot_emm_bar(df_b3, pond_map, "Batch 3 — Survival Probabilities ")
p_MK <- plot_emm_bar(df_MK, pond_map, "Batch MK — Survival Probabilities ")

combined_plot <- 
  p_b3     /
  p_b12    /
  p_b123   /
  p_MK /
  p_b123MK +
  plot_layout(ncol = 1)

# ----
# Save as png
# ----
ggsave(
  filename = "plots/compare_datasets_survival.png",
  plot = combined_plot,
  width = 10,
  height = 30,   # increase height to avoid compression
  dpi = 300
)

# ----
# Save as pdf
# ----
ggsave(
  filename = "plots/compare_datasets_survival.png",
  plot = combined_plot,
  width = 10,
  height = 30,   # increase height to avoid compression
)



