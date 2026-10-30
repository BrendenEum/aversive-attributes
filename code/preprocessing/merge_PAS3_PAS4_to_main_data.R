# =============================================================================
# merge_PAS3_PAS4_to_main_data.R
#
# Combine the two mouse-tracking experiments (PAS3 = Study 1, PAS4 = Study 2)
# into one trial-level file and one fixation-level file:
#
#   data/processed/main_data.csv       one row per trial
#   data/processed/main_data_long.csv  one row per fixation (mouse hover)
#
# Both files get a `Study` column (1 or 2) and a study-unique `Subject` id.
# Everything downstream (notebooks/Main_Analysis.Rmd) reads these two files.
# =============================================================================

library(tidyverse)

paths    <- yaml::yaml.load_file(here::here("config/paths.yaml"))
data_dir <- file.path(paths$root, paths$data$processed)


# -----------------------------------------------------------------------------
# 1. Read the four input files
# -----------------------------------------------------------------------------

pas3      <- read_csv(file.path(data_dir, "PAS3.csv"),     show_col_types = FALSE)
pas3_long <- read_csv(file.path(data_dir, "PAS3long.csv"), show_col_types = FALSE)
pas4      <- read_csv(file.path(data_dir, "PAS4.csv"),     show_col_types = FALSE)
pas4_long <- read_csv(file.path(data_dir, "PAS4long.csv"), show_col_types = FALSE)


# -----------------------------------------------------------------------------
# 2. Study 2 exclusions
#
# Six PAS4 subjects have a subset of trials where the P/N attribute labels are
# swapped relative to the sign of the values (e.g. an "N" row carrying positive
# values). We drop them entirely, from both the trial and the fixation file.
# -----------------------------------------------------------------------------

swapped_label_subjects <- c(420, 561, 679, 714, 967, 1015)

pas4      <- pas4      %>% filter(!(SubjectNumber %in% swapped_label_subjects))
pas4_long <- pas4_long %>% filter(!(SubjectNumber %in% swapped_label_subjects))


# -----------------------------------------------------------------------------
# 3. Make the two studies' columns line up
#
# The only naming difference is a within-block trial counter called Trial10 in
# PAS3 and Trial20 in PAS4. It isn't used anywhere, so drop it from both.
# -----------------------------------------------------------------------------

pas3      <- pas3      %>% select(-Trial10)
pas3_long <- pas3_long %>% select(-Trial10)
pas4      <- pas4      %>% select(-Trial20)
pas4_long <- pas4_long %>% select(-Trial20)

stopifnot(identical(names(pas3),      names(pas4)))
stopifnot(identical(names(pas3_long), names(pas4_long)))

# PAS4's long file stores Trial as text; make it numeric like everything else.
pas4_long <- pas4_long %>% mutate(Trial = as.numeric(Trial))


# -----------------------------------------------------------------------------
# 4. Stack the two studies
#
# Study = 1 for PAS3, 2 for PAS4. SubjectNumber is NOT unique across studies
# (hundreds of IDs repeat), so we also build a study-unique Subject id, e.g.
# "S1_5" and "S2_5" are two different people.
# -----------------------------------------------------------------------------

main_data <- bind_rows(
  pas3 %>% mutate(Study = 1),
  pas4 %>% mutate(Study = 2)
) %>%
  mutate(Subject = paste0("S", Study, "_", SubjectNumber)) %>%
  relocate(Study, Subject)

main_data_long <- bind_rows(
  pas3_long %>% mutate(Study = 1),
  pas4_long %>% mutate(Study = 2)
) %>%
  mutate(Subject = paste0("S", Study, "_", SubjectNumber)) %>%
  relocate(Study, Subject)


# -----------------------------------------------------------------------------
# 5. Fixation-level fixes
#
# (a) Units. FixDur was recorded in 10 ms ticks (100 Hz mouse sampling): within
#     every trial, sum(FixDur) == 100 * (DwellL + DwellR), and the Dwell columns
#     are in seconds. Convert to milliseconds so the units are obvious.
# (b) Fixation order. Rows are already in temporal order within a trial, but we
#     add an explicit counter so first/last fixations never depend on row order.
# -----------------------------------------------------------------------------

main_data_long <- main_data_long %>%
  mutate(FixDur = FixDur * 10) %>%
  group_by(Study, Subject, Trial) %>%
  mutate(FixNum = row_number()) %>%
  ungroup()

# Sanity check on the unit conversion: summed fixation time per trial should
# now equal total dwell time (both in seconds) up to rounding.
unit_check <- main_data_long %>%
  group_by(Study, Subject, Trial) %>%
  summarise(
    fix_s   = sum(FixDur) / 1000,
    dwell_s = first(DwellL + DwellR),
    .groups = "drop"
  )
stopifnot(all(abs(unit_check$fix_s - unit_check$dwell_s) < 0.01))


# -----------------------------------------------------------------------------
# 6. Save, and print a summary so the console output documents the merge
# -----------------------------------------------------------------------------

write_csv(main_data,      file.path(data_dir, "main_data.csv"))
write_csv(main_data_long, file.path(data_dir, "main_data_long.csv"))

summary_tbl <- main_data %>%
  group_by(Study) %>%
  summarise(subjects = n_distinct(Subject), trials = n(), .groups = "drop") %>%
  left_join(
    main_data_long %>% group_by(Study) %>% summarise(fixations = n(), .groups = "drop"),
    by = "Study"
  )

print(summary_tbl)
