
# ------------------
# Environment setup
# ------------------
# Install once:
#install.packages(c("tidyverse","janitor","lubridate","glmmTMB","DHARMa","emmeans","broom.mixed"))

library(tidyverse)
library(janitor)
library(lubridate)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(broom.mixed)
library(stringr)
library(dplyr)

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

# treatment and ponds : fixed effect variables
# clonal line, date, observer : random effect variables

# just batch 3
long_b3 <- read_csv("data_processed/batch3_long.csv")

long_b3 <- long_b3 %>%
  mutate(
    treatment  = factor(treatment, levels = c("control","0.5","0.7")),
    pond_code  = factor(pond_code, levels = c("C5", setdiff(unique(pond_code), "C5"))),
    pond_type  = factor(pond_type, levels = c("natu","agri","city")),
    observer   = factor(observer),
    clonal_line = factor(clonal_line)
  )
saveRDS(long_b3, "batch3_long.rds")

# just batch 1 and 2
long_b12 <- read_csv("data_processed/batch12_long.csv")

long_b12 <- long_b12 %>%
  mutate(
    treatment  = factor(treatment, levels = c("control","0.5","0.7")),
    pond_code  = factor(pond_code, levels = c("C5", setdiff(unique(pond_code), "C5"))),
    pond_type  = factor(pond_type, levels = c("natu","agri","city")),
    observer   = factor(observer),
    clonal_line = factor(clonal_line)
  )
saveRDS(long_b12, "batch12_long.rds")

# just the MK 
long_MK <- read_csv("data_processed/matkris_processed.csv")

long_MK <- long_MK %>%
  mutate(
    treatment  = factor(treatment, levels = c("control","0.5","0.7")),
    pond_code  = factor(pond_code, levels = c("C5", setdiff(unique(pond_code), "C5"))),
    pond_type  = factor(pond_type, levels = c("natu","agri","city")),
    observer   = factor(observer),
    clonal_line = factor(clonal_line)
  )
saveRDS(long_MK, "MK_long.rds")
  
# =====================================
#  OVERALL WINNERS: Binomial GLMM + OLRE
# =====================================
# W/ Interaction between pond and treatment
# -------------------------------------
# Alex + MK: WITH observer as random
long_b123MK_olre <- long_b123MK %>% mutate(obs_id = factor(row_number()))

m_olre_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment * pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b123MK_olre
)

# Alex: WITHOUT observer as random
long_b123_olre <- long_b123 %>% mutate(obs_id = factor(row_number()))

m_olre_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment * pond_code +
    (1 | clonal_line) + (1 | date) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b123_olre
)
AIC(m_beta_b123MK, m_olre_b123MK, mnoOb_beta_b123, m_olre_b123)


# =====================================
# WINNERS FOR THE BETA-BINOMIAL:
# =====================================
# MK + Alex: WITH observer as random
m_beta_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer),
  family = betabinomial(link = "logit"),
  data = long_b123MK
)

# Alex: WITHOUT observer as random
mnoOb_beta_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date),
  family = betabinomial(link = "logit"),
  data = long_b123
)


# =====================================
#  Beta-binomial GLMM
# =====================================
# -------------------------
# observer effect
# -------------------------
# Alex + MK
mnoOb_beta_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date),
  family = betabinomial(link = "logit"),
  data = long_b123MK
)

m_beta_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer),
  family = betabinomial(link = "logit"),
  data = long_b123MK
)

AIC(mnoOb_beta_b123MK, m_beta_b123MK)
# Alex
mnoOb_beta_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date),
  family = betabinomial(link = "logit"),
  data = long_b123
)
m_beta_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer),
  family = betabinomial(link = "logit"),
  data = long_b123
)
AIC(mnoOb_beta_b123, m_beta_b123)

# Result: including MK (Matthias and Kristina's) results, observer becomes significant!!!
# ALSO: a lot more variation unexplained (higher AICs)


# -------------------------
# pond as a fixed effect
# -------------------------
# Alex + MK
mnoPo_beta_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment + 
    (1 | clonal_line) + (1 | date) + (1 | observer),
  family = betabinomial(link = "logit"),
  data = long_b123MK
)

m_beta_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer),
  family = betabinomial(link = "logit"),
  data = long_b123MK
)

AIC(mnoPo_beta_b123MK, m_beta_b123MK)
# Alex
mnoPo_beta_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment +
    (1 | clonal_line) + (1 | date) + (1 | observer),
  family = betabinomial(link = "logit"),
  data = long_b123
)
m_beta_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer),
  family = betabinomial(link = "logit"),
  data = long_b123
)
AIC(mnoPo_beta_b123, m_beta_b123)

# results: Pond code important for both datasets!

# -------------------------
# Date
# -------------------------
# Alex + MK
mnoDa_beta_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | observer),
  family = betabinomial(link = "logit"),
  data = long_b123MK
)

m_beta_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer),
  family = betabinomial(link = "logit"),
  data = long_b123MK
)

# Alex
mnoDa_beta_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | observer),
  family = betabinomial(link = "logit"),
  data = long_b123
)
m_beta_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer),
  family = betabinomial(link = "logit"),
  data = long_b123
)

AIC(mnoDa_beta_b123MK, m_beta_b123MK, mnoDa_beta_b123, m_beta_b123)

# result: date important for both!

# -------------------------
# Pond type
# -------------------------
# Alex + MK
mwPT_beta_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_type +
    (1 | clonal_line) + (1 | date) + (1 | observer)  + (1 | pond_code),
  family = betabinomial(link = "logit"),
  data = long_b123MK
)

m_beta_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code
    (1 | clonal_line) + (1 | date) + (1 | observer),
  family = betabinomial(link = "logit"),
  data = long_b123MK
)
mwPTr_beta_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer)  + (1 | pond_type),
  family = betabinomial(link = "logit"),
  data = long_b123MK
)
# Alex
mwPT_beta_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_type +
    (1 | clonal_line) + (1 | date) + (1 | observer)  + (1 | pond_code),
  family = betabinomial(link = "logit"),
  data = long_b123
)
m_beta_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer),
  family = betabinomial(link = "logit"),
  data = long_b123
)
mwPTr_beta_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer)  + (1 | pond_type),
  family = betabinomial(link = "logit"),
  data = long_b123
)
AIC(mwPT_beta_b123MK, m_beta_b123MK, mwPTr_beta_b123MK, mwPT_beta_b123, m_beta_b123, mwPTr_beta_b123)

# Result: in both datasets,  no pond_type does NOT help either as fixed or as random variable!!
# note: you can not have both pond_type and pond_code as fixed effects because they are nested


# =====================================
#  Binomial GLMM + OLRE
# =====================================
# Alex + MK
long_b123MK_olre <- long_b123MK %>% mutate(obs_id = factor(row_number()))

m_olre_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b123MK_olre
)

# Alex
long_b123_olre <- long_b123 %>% mutate(obs_id = factor(row_number()))

m_olre_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b123_olre
)
AIC(m_beta_b123MK, m_olre_b123MK, mnoOb_beta_b123, m_olre_b123)


# -------------------------
# Observer
# -------------------------
# Alex + MK
mnoOb_olre_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment * pond_code +
    (1 | clonal_line) + (1 | date) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b123MK_olre
)
AIC(m_olre_b123MK, mnoOb_olre_b123MK)

# Result: the observer is important for the combined dataset!!

# -------------------------
# INTERACTION pond*treatment 
# -------------------------
# (allows for separate statistical testings between 
# ponds for each treatment)
# -------------------------
# Alex + MK
mPTInter_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment * pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b123MK_olre
)
mPT_b123MK <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b123MK_olre
)
AIC(mPT_b123MK, mPTInter_b123MK)
# Alex: 
mTP_olre_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b123_olre
)
mTPInter_olre_b123 <- glmmTMB(
  cbind(alive, dead) ~ treatment * pond_code +
    (1 | clonal_line) + (1 | date) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b123_olre
)
AIC(mTP_olre_b123, mTPInter_olre_b123)

# Result: It is in both cases BETTER WITH INTERACTION BETWEEN POND*TREATMENT !!!
# !!!

# =====================================
#  Check for only MK or b3
# =====================================
summary(long_MK)
long_MK_olre <- long_MK %>% mutate(obs_id = factor(row_number()))
summary(long_MK_olre)

# MK
m_MK <- glmmTMB(
  cbind(alive, dead) ~ treatment * pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_MK_olre
)
mnoOb_MK <- glmmTMB(
  cbind(alive, dead) ~ treatment * pond_code +
    (1 | clonal_line) + (1 | date)  + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_MK_olre
)
mnoInt_MK <- glmmTMB(
  cbind(alive, dead) ~ treatment + pond_code +
    (1 | clonal_line) + (1 | date)  + (1 | observer)  + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_MK_olre
)
AIC(m_MK, mnoOb_MK, mnoInt_MK)
# RESULT: the MK alone needs simpler model (NO interaction pond x treatment)

# b3
long_b3
long_b3_olre <- long_b3 %>% mutate(obs_id = factor(row_number()))
summary(long_b3)
m_b3 <- glmmTMB(
  cbind(alive, dead) ~ treatment * pond_code +
    (1 | clonal_line) + (1 | date) + (1 | observer) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b3_olre
)
mnoOb_b3 <- glmmTMB(
  cbind(alive, dead) ~ treatment * pond_code +
    (1 | clonal_line) + (1 | date)  + (1 | obs_id),
  family = binomial(link = "logit"),
  data = long_b3_olre
)
AIC(m_b3, mnoOb_b3)
#RESULT: the b3 alone can handle the complex model (with interaction)

# =====================================
# DHARMA TESTS
# =====================================

# Alex + MK:
# Simulate residuals
res_beta_b123MK <- simulateResiduals(m_beta_b123MK, plot = TRUE)  # opens quick plots
# Formal tests
testDispersion(res_beta_b123MK)            # over/under-dispersion
testZeroInflation(res_beta_b123MK)         # should be NA / not applicable for (beta)binomial
#testTemporalAutocorrelation(res_beta_b123MK, time = long_b123MK_olre$date)  # if date is ordered

# Alex:
# Simulate residuals
res_beta_b123 <- simulateResiduals(m_beta_b123, plot = TRUE)  # opens quick plots
# Formal tests
testDispersion(res_beta_b123)            # over/under-dispersion
testZeroInflation(res_beta_b123)         # should be NA / not applicable for (beta)binomial

