
# ------------------
# Environment setup
# ------------------
# Install once:
#install.packages(c("tidyverse","janitor","lubridate","glmmTMB","DHARMa","emmeans","broom.mixed"))

library(tidyverse)
library(janitor)
library(lubridate)
library(broom.mixed)
library(stringr)

rm(list = ls())
# ------------------
# 1. Read file
# ------------------
raw <- read_csv("batch123_cleaned.csv") %>% clean_names()

# Harmonize names
raw <- raw %>%
  mutate(
    date = as.character(date),
    date = if_else(is.na(date) | date %in% c("NA",""), NA_character_, date),
    observer = factor(observer),
    pond_type = factor(pond_type, levels = c("natu","agri","city")),
    # basic sanity coercions
    ctrl      = as.numeric(ctrl),
    `0_5ug_L` = as.numeric(`x0_5ug`),
    `0_7ug_L` = as.numeric(`x0_7ug`)
  )

# Date parsing with year info
long_months <- c("Oct","Nov","Dec","Jan","Feb","Mar")

dat <- raw %>%
  mutate(
    # 1) Clean raw date text
    date_chr = str_trim(as.character(date)),
    date_chr = if_else(date_chr %in% c("", "NA"), NA_character_, date_chr),
    
    # strip any trailing punctuation (e.g., "06-Mar," -> "06-Mar")
    date_chr = str_replace_all(date_chr, "[,;]+$", ""),
    
    # 2) Extract day ("07" from "07-Nov")
    day_txt = if_else(!is.na(date_chr), str_sub(date_chr, 1, 2), NA_character_),
    day     = suppressWarnings(as.integer(day_txt)),
    
    # 3) Extract month abbreviation ("Nov" from "07-Nov")
    month_abbr = if_else(!is.na(date_chr), str_sub(date_chr, 4, 6), NA_character_),
    
    # 4) Map month to numeric (10..12, 1..3)
    month = case_when(
      month_abbr == "Oct" ~ 10L,
      month_abbr == "Nov" ~ 11L,
      month_abbr == "Dec" ~ 12L,
      month_abbr == "Jan" ~  1L,
      month_abbr == "Feb" ~  2L,
      month_abbr == "Mar" ~  3L,
      TRUE ~ NA_integer_
    ),
    
    # 5) Assign year based on month
    year = case_when(
      month %in% c(10L, 11L, 12L) ~ 2025L,  # Oct/Nov/Dec 2025
      month %in% c( 1L,  2L,  3L) ~ 2026L,  # Jan/Feb/Mar 2026
      TRUE ~ NA_integer_
    ),
    
    # 6) Build a proper Date only when all components are present
    date_full = if_else(!is.na(day) & !is.na(month) & !is.na(year),
                        make_date(year = year, month = month, day = day),
                        as_date(NA))
  ) %>%
  # 7) Keep necessary columns (use the new date_full)
  select(
    clonal_line, pond_code, pond_name, line_number, pond_type, clone_code,
    ctrl, `0_5ug_L`, `0_7ug_L`, observer, date_full, notes
  ) %>%
  rename(date = date_full)  # from now on, use `date` as the proper Date

# save file as a csv and an rds for future use
write_csv(dat, "batch123_cleaned_processed.csv")
saveRDS(dat, "batch123_cleaned_processed.rds")

# ------------------
# 2. QC check
# ------------------

# Values outside [0,10]
bad_vals <- dat %>%
  pivot_longer(cols = c(ctrl, `0_5ug_L`, `0_7ug_L`),
               names_to = "treatment", values_to = "alive") %>%
  filter(!is.na(alive) & (alive < 0 | alive > 10))

# Missing dates
n_date_na <- sum(is.na(dat$date))

# Notes that indicate irregularities
notes_nonempty <- dat %>%
  filter(!is.na(notes) & str_trim(notes) != "")

# Possible duplicates: same clonal_line + date + observer with identical counts
dups <- dat %>%
  group_by(clonal_line, date, observer) %>%
  filter(n() > 1) %>%
  summarise(n = n(),
            ctrl_vals = paste(ctrl, collapse = ","),
            d05_vals  = paste(`0_5ug_L`, collapse = ","),
            d07_vals  = paste(`0_7ug_L`, collapse = ","),
            .groups = "drop") %>%
  mutate(identical_counts = ctrl_vals == unique(ctrl_vals)[1] &
           d05_vals  == unique(d05_vals)[1]  &
           d07_vals  == unique(d07_vals)[1])

qc_summary <- list(
  n_rows = nrow(dat),
  n_date_na = n_date_na,
  any_out_of_bounds = nrow(bad_vals) > 0,
  n_notes = nrow(notes_nonempty),
  duplicate_keys = dups %>% filter(identical_counts)
)

qc_summary
dups

# ------------------
# 3. Long format
# ------------------
long <- dat %>%
  pivot_longer(
    cols = c(ctrl, `0_5ug_L`, `0_7ug_L`),
    names_to = "treatment",
    values_to = "alive"
  ) %>%
  mutate(
    treatment = recode(treatment,
                       ctrl = "control",
                       `0_5ug_L` = "0.5",
                       `0_7ug_L` = "0.7"),
    treatment = factor(treatment, levels = c("control","0.5","0.7")),
    dose = as.numeric(recode(as.character(treatment),
                             control = "0", `0.5` = "0.5", `0.7` = "0.7")),
    n = 10L,
    dead = n - alive,
    # optional: flags
    note_flag = !is.na(notes) & str_trim(notes) != "",
    copied_flag = note_flag & str_detect(notes, regex("copied", ignore_case = TRUE)),
    over10_flag = note_flag & str_detect(notes, regex("11 alive|12 alive|>10|more than 10", ignore_case = TRUE)),
    control_fail_flag = (treatment == "control" & alive < 9)  # OECD-202 control validity <10% immob.
  )

# Optional: think of dropping suspect rows for primary model
long_clean <- long %>%
  filter(!copied_flag) %>%                     # drop "copied" lines
  group_by(clonal_line, date, observer, treatment) %>% 
  slice_head(n = 1) %>%                        # in case of exact duplicates
  ungroup()


# Set C5 as reference pond
long_clean <- long_clean %>%
  mutate(
    pond_code = relevel(factor(pond_code), ref = "C5")
  )

# save file as a csv and an rds for future use
write_csv(long_clean, "batch123_long.csv")
saveRDS(long_clean, "batch123_long.rds")
